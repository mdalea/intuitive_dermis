# Texture 
#CUDA_VISIBLE_DEVICES=3 python intuitivespikes_snn_sg.py --datasetPath ../spikegen/global_outfile/spikegen/ALL_N_lcadc-0-uV-v2-0-uV-0-uV-0-pct-4-mr_cnt-192-multitrial-Kylberg_filt_6_scan_actualdimscale_chip_lchpf--1to40-1to40_sorted/ --epochs 150 --hidden_size 256 --num_classes 6 --batch_size 32 --tmax_ms 400 --Ts 10 --kfolds 6 --cache yes --save_bestmodel yes


"""
Train 1-Hidden Layer LIF SNN using Surrogate Gradient
Trained model is used for SVM classification
Author: Ali Safa - IMEC- KU Leuven, Federico Corradi - IMEC-NL
Modified by: Mark Alea - KU Leuven
04/02/2023
"""

import sys, os
CURRENT_TEST_DIR = os.getcwd()
sys.path.append(CURRENT_TEST_DIR + "/../../src")

import numpy as np
from eNetworks import mini_eMLP
import torch.nn as nn
import torch
from torch.utils.data import TensorDataset, Dataset, DataLoader, ConcatDataset
import matplotlib.pyplot as plt
from sklearn.utils import shuffle
from sklearn.model_selection import KFold

import tonic
from tonic.io import read_mnist_file
from tonic import CachedDataset

import argparse

path = os.getcwd()
parser = argparse.ArgumentParser("Ali Safa's Surrogate Gradient Training with cross-entropy loss (For REAL chip spikes)")
parser.add_argument("--datasetPath",type=str,help="Dataset path", required=True)
parser.add_argument("--epochs",type=int,help="No. of epochs of training",required=True)
parser.add_argument("--hidden_size",type=int,help="MLP hidden size",required=True)
parser.add_argument("--ch_size",type=int,help="No. of input channels",required=False)
parser.add_argument("--x_size",type=int,help="Input x-dim",required=False)
parser.add_argument("--y_size",type=int,help="Input y-dim",required=False)
parser.add_argument("--num_classes",type=int,help="No. of output classes", required=True)
parser.add_argument("--batch_size",type=int,help="Training batch size",required=True)
parser.add_argument("--lr",type=float,help="Learning rate",required=False)
parser.add_argument("--tmax_ms",type=float,help="Maximum sensor recording considered",required=True)
parser.add_argument("--Ts",type=float,help="Sampling time",required=True)
parser.add_argument("--show_inputspikes",type=str,help="Visualize input spikes",required=False)
parser.add_argument("--kfolds",type=int,help="No. of k-folds",required=True)
parser.add_argument("--save_bestmodel",type=str,help="Save best model and run k-fold with it",required=False)
parser.add_argument("--cache",type=str,help="Cache dataset for faster access",required=True)
parser.add_argument("--retrainPath",type=str,help="Re-train existing network. Enter model path .pt",required=False)
parser.add_argument("--retrainPath_id",type=str,help="Custom ID to differentiate which model was loaded",required=False)

args = parser.parse_args()

# Hyper-parameters - match orig example
if args.ch_size:
    ch_size = args.ch_size
else:
    ch_size = 2

if args.x_size:
    x_size = args.x_size
else:
    x_size = 16

if args.y_size:
    y_size = args.y_size
else:
    y_size = 12

#use custom ID to name retrain path
if args.retrainPath_id:
    retrainPath = '-' + args.retrainPath_id

else:
    if args.retrainPath:
        tokens_slash = args.retrainPath.split('/')
        print(tokens_slash[-2])
        print(tokens_slash[-1])
        tokens_dot = tokens_slash[-1].split('.')
        tokens = tokens_dot[0].split('-')
        #print(modelname[0])
        retrainPath = '-retrainPath-' + tokens_slash[-2] + '-' + tokens[-2] + '-' + tokens[-1]
        print(retrainPath)
    else:
        retrainPath = ''


input_size = ch_size*x_size*y_size 
hidden_size = args.hidden_size 
num_classes = args.num_classes
num_epochs = args.epochs
batch_size = args.batch_size
print('hidden_size')
print(hidden_size)
print('num_classes')
print(num_classes)
print('batch_size')
print(batch_size)
if args.lr:
    learning_rate = args.lr
else:
    learning_rate = 1e-3

if args.datasetPath:
    tokens_slash = args.datasetPath.split('/')
    datasetPath = '-dsetPath-' + tokens_slash[-2]
else:
    datasetPath = ''

if args.retrainPath:
    extraPath =  'ep-' + str(args.epochs) +  retrainPath 
else:
    extraPath =  'ep-' + str(args.epochs) + '-isize-' + str(ch_size) + '-' + str(x_size) + '-' + str(y_size) + '-hsize-' + str(args.hidden_size) + '-bsize-' + str(args.batch_size) + '-lr-' + str(learning_rate) + '-tmax_ms-' + str(args.tmax_ms) + '-Ts-' + str(args.Ts) + "-k-" + str(args.kfolds) + datasetPath 

sampleFile_train = args.datasetPath + "train1K.txt"
sampleFile_test = args.datasetPath + "test100.txt"
savePath = path + "/SNN_" + extraPath

print(args.datasetPath)
print(savePath)
mode=0o777
if os.path.exists(savePath) == False:
    os.mkdir(savePath)

# Dataset loading
# Load events into frames for GPU training
class SpkDataset_toFrameTimeBins(Dataset):
    def __init__(self, datasetPath, sampleFile, samplingTime):
        self.path = datasetPath 
        self.samples = np.loadtxt(sampleFile).astype('int')
        self.samplingTime = float(samplingTime)

    def __getitem__(self, index):
        inputIndex  = self.samples[index, 0]
        classLabel  = self.samples[index, 1]


        dtype = np.dtype([('x', '<i8'), ('y', '<i8'), ('t', 'int32'), ('p', '<i8')]) 
        ordering = dtype.names

        events = read_mnist_file(self.path+str(inputIndex.item())+".bs2", dtype=dtype)


        frame_transform = tonic.transforms.ToFrame(sensor_size=(x_size,y_size,ch_size),n_time_bins=int(args.tmax_ms/args.Ts))
        frames = frame_transform(events).astype(np.float32)

        if args.show_inputspikes:
            return frames,classLabel
        else:
            return frames.reshape(-1, ch_size*x_size*y_size),classLabel

    def __len__(self):
        return self.samples.shape[0]

            
criterion = nn.NLLLoss()


def train(model, train_loader, test_loader, optimizer, epochs, batch_size, acc_hist_test, fold):

    for e in range(epochs):
        train_acc = 0
        train_loss_sum = 0
        total = 0
        for i, (images, labels) in enumerate(train_loader):

            optimizer.zero_grad()
            hidden_mem =  torch.zeros(model.hidden_size_1) # random init
            output_mem =  torch.zeros(model.output_size) #random init
            hidden_spike = torch.zeros(model.hidden_size_1)
            output_spike = torch.zeros(model.output_size)
            
            predictions, train_loss, nbr_events_avrg = model(images, labels, hidden_mem, hidden_spike, output_mem, output_spike)
            predicted = predictions
            
                
            train_loss.sum().backward(retain_graph=True)
            train_loss_sum += train_loss.data.cpu().numpy()
            
                
            optimizer.step()
            predicted = predicted.t()
            train_acc += (predicted.T[0] == labels).sum() 
            total += labels.size(0)
        
        train_acc = train_acc.data.cpu().numpy()
        acc_hist_train.append(train_acc/total)
        loss_hist_train.append(train_loss_sum.item()/total)
        print("Fold: " + str(fold) + " Epoch: " + str(e) + " Loss: " + str(train_loss_sum.item()/total) + " Accuracy: " + str(train_acc/total))

        # Test per epoch -> comment if too slow
        acc_hist_test = test(model, test_loader, batch_size, acc_hist_test, fold)

        # SAVE best model so far (based on train acc)
        if args.save_bestmodel:
            if train_acc/total >= np.max(acc_hist_train):
                print(" -->  Saving best trained model : epoch " + str(e) + " Accuracy: " + str(train_acc/total))
                torch.save(model.state_dict(), savePath + '/chip-classifier-' + str(fold) + '-best.pt') 

    return acc_hist_test
        
def test(model, testloader, batch_size, acc_hist_test, fold):
    correct = 0
    total = 0
    correct_best = 0
    test_loss_sum = 0
    loss_hist_test = []

    y_pred = []
    y_true = []

    # Evaluationfor this fold
    for i, (images, labels) in enumerate(testloader):

      hidden_mem =  torch.zeros(model.hidden_size_1) # random init
      output_mem =  torch.zeros(model.output_size) #random init
      hidden_spike = torch.zeros(model.hidden_size_1)
      output_spike = torch.zeros(model.output_size)
            
      predictions, test_loss, nbr_events_avrg = model(images, labels, hidden_mem, hidden_spike, output_mem, output_spike)
      predicted = predictions
            
      test_loss_sum += test_loss.data.cpu().numpy()
            
      predicted = predicted.t()
      correct += (predicted.T[0] == labels).sum() 
      total += labels.size(0) 
      print(predicted.T[0])
      print(labels)

      y_pred.extend(predicted.T[0].data.cpu().numpy())
      y_true.extend(labels.long().cpu().numpy())

    test_acc = correct.data.cpu().numpy()
    acc_hist_test.append(correct.data.cpu().numpy()/total)
    print("--> Test:  Loss: " + str(test_loss_sum.item()/total) + " Accuracy: " + str(test_acc/total))

    # SAVE best model so far (based on test acc)
    if args.save_bestmodel:
        if test_acc/total >= np.max(acc_hist_test):
            print(" -->  Saving best trained model (based on test): Accuracy: " + str(test_acc/total))
            torch.save(model.state_dict(), savePath + '/chip-classifier-test-' + str(fold) + '-best.pt') 


    return acc_hist_test

if __name__ == '__main__':      

    # Define the cuda device to run the code on.
    torch.cuda.empty_cache()
    device = torch.device('cuda')

    trainingSet = SpkDataset_toFrameTimeBins(datasetPath =args.datasetPath, 
                                sampleFile  =sampleFile_train,
                                samplingTime=args.Ts)

    testingSet = SpkDataset_toFrameTimeBins(datasetPath  =args.datasetPath, 
                                sampleFile  =sampleFile_test,
                                samplingTime=args.Ts)


    trainingSet_cached = CachedDataset(trainingSet, cache_path='./cache-' + extraPath + '/spkdataset/train',reset_cache=True)
    trainLoader = DataLoader(dataset=trainingSet_cached, batch_size=batch_size, shuffle=True, num_workers=4)

    testingSet_cached = CachedDataset(testingSet, cache_path='./cache-' + extraPath + '/spkdataset/test',reset_cache=True)
    testLoader = DataLoader(dataset=testingSet_cached, batch_size=batch_size, shuffle=True, num_workers=4)


    if args.show_inputspikes:
        for i in range(1):
            input, classLabel = trainingSet[i]
            print(input.shape)
            print(classLabel)
            plt.subplot(5,5,i+1)
            #plt.ylabel("Class: " + str(classLabel))
            #print ON channel spikes
            im = plt.imshow(input[0][0])
        plt.colorbar(im, ax=plt.gcf().axes)
        plt.show()

        for i in range(40):
            input, classLabel = trainingSet[0]
            plt.subplot(5,8,i+1)
            #plt.ylabel("Class: " + str(classLabel))
            #print ON channel spikes
            im = plt.imshow(input[i][0])
        plt.colorbar(im, ax=plt.gcf().axes)
        plt.show()

    # Configuration options
    k_folds = args.kfolds
  
    # For fold results
    results = {}
    results_best = {}
    results_max = {}
  
    # Set fixed random number seed
    torch.manual_seed(42)
  
    # Prepare dataset by concatenating Train/Test part; we split later.
    dataset = ConcatDataset([trainingSet_cached, testingSet_cached])

    # Define the K-fold Cross Validator
    kfold = KFold(n_splits=k_folds, shuffle=True)
    
    # Start print
    print('--------------------------------')

    # store accuracy per class for all folds
    accs_all = [];

    # Store mean accuracy for all folds for top K
    top_KK_acc_folds = [];

    # Store mean accuracy for all folds for top K (chosen K)
    top_KK_acc_chosenK_folds = [];

    # K-fold Cross Validation model evaluation
    for fold, (train_ids, test_ids) in enumerate(kfold.split(dataset)):
    
      # Print
      print(f'FOLD {fold}')
      print('--------------------------------')
    
      # Sample elements randomly from a given list of ids, no replacement.
      train_subsampler = torch.utils.data.SubsetRandomSampler(train_ids)
      test_subsampler = torch.utils.data.SubsetRandomSampler(test_ids)
    
      # Define data loaders for training and testing data in this fold
      trainloader = torch.utils.data.DataLoader(
                            dataset, 
                            batch_size=batch_size, sampler=train_subsampler)
      testloader = torch.utils.data.DataLoader(
                            dataset,
                            batch_size=batch_size, sampler=test_subsampler)
    
      # Init the neural network
      model = mini_eMLP(ch_size, x_size, y_size, hidden_size, num_classes, criterion=criterion, batch_size=batch_size)
      model_best = mini_eMLP(ch_size, x_size, y_size, hidden_size, num_classes, criterion=criterion, batch_size=batch_size)

      if args.retrainPath:
          model.load_state_dict(torch.load(args.retrainPath))     

      #base_params = [model.i2h.weight, model.i2h.bias, model.h2o.weight, model.h2o.bias]
      base_params = [model.i2h.weight, model.i2h.bias, model.h2o.weight, model.h2o.bias, model.h2h.weight, model.h2h.bias]
      optimizer = torch.optim.Adam([{'params': base_params},], lr=learning_rate)



      loss_hist_train = []
      loss_hist_test = []
      acc_hist_train = []
      acc_hist_test = []


      # Training
      acc_hist_test_ = train(model, trainloader, testloader, optimizer, num_epochs, batch_size, acc_hist_test, fold) 
      acc_hist_test_max = np.max(acc_hist_test_)
      epoch_max = np.argmax(acc_hist_test_)


      ## Save model
      torch.save(model.state_dict(), savePath + '/chip-classifier-' + str(fold) + '.pt')               
        

      # Testing
      # Load best model
      if args.save_bestmodel:
          model_best.load_state_dict(torch.load( savePath + '/chip-classifier-' + str(fold) + '-best.pt'))


      correct = 0
      total = 0
      correct_best = 0
      #test_acc = 0
      test_loss_sum = 0
      acc_hist_test = []
      loss_hist_test = []

      #if epoch == num_epochs-1:
      y_pred = []
      y_true = []

      # Last Evaluationfor this fold 
      for i, (images, labels) in enumerate(testloader):

        hidden_mem =  torch.zeros(model.hidden_size_1) # random init
        output_mem =  torch.zeros(model.output_size) #random init
        hidden_spike = torch.zeros(model.hidden_size_1)
        output_spike = torch.zeros(model.output_size)
            
        predictions, test_loss, nbr_events_avrg = model(images, labels, hidden_mem, hidden_spike, output_mem, output_spike)
        predicted = predictions
            
        test_loss_sum += test_loss.data.cpu().numpy()
            
        predicted = predicted.t()
        correct += (predicted.T[0] == labels).sum() 
        total += labels.size(0) 
        print(predicted.T[0])
        print(labels)

        #if epoch == num_epochs-1:
        y_pred.extend(predicted.T[0].data.cpu().numpy())
        y_true.extend(labels.long().cpu().numpy())

        # test best model
        if args.save_bestmodel:
          hidden_mem =  torch.zeros(model.hidden_size_1) # random init
          output_mem =  torch.zeros(model.output_size) #random init
          hidden_spike = torch.zeros(model.hidden_size_1)
          output_spike = torch.zeros(model.output_size)
            
          predictions, _, _ = model_best(images, labels, hidden_mem, hidden_spike, output_mem, output_spike)
          predicted = predictions
            
            
          predicted = predicted.t()
          correct_best += (predicted.T[0] == labels).sum() 
          print(predicted.T[0])
          print(labels)


      test_acc_ = correct / total
      acc_hist_test.append(test_acc_)
      loss_hist_test.append(test_loss_sum.item()/total)


      print('Accuracy for fold %d: %f %%' % (fold, 100.0 * test_acc_))
      print('--------------------------------')
      print('--> At epoch = %d, Max: %f %%' % (epoch_max, 100.0 * acc_hist_test_max))       
      results[fold] = 100.0 * (test_acc_)
      results_max[fold] =  100.0 * (acc_hist_test_max)

      #----------------------------------------------
      # Testing (BEST trained model)
      if args.save_bestmodel:
        test_acc_best = correct_best / total
        print('[BEST model] Accuracy for fold %d: %f %%' % (fold, 100.0 * test_acc_best))
        print('--------------------------------')
        results_best[fold] = 100.0 * (test_acc_best) 

      with open(savePath + '/loss_hist_train-' + str(fold) + '.txt', 'w') as f:
          for item in loss_hist_train:
            f.write("%f\n" % item)

      with open(savePath + '/loss_hist_test-' + str(fold) + '.txt', 'w') as f:
          for item in loss_hist_test:
              f.write("%f\n" % item)

      with open(savePath + '/acc_hist_train-' + str(fold) + '.txt', 'w') as f:
          for item in acc_hist_train:
              f.write("%f\n" % item)

      with open(savePath + '/acc_hist_test-' + str(fold) + '.txt', 'w') as f:
          for item in acc_hist_test:
              f.write("%f\n" % item)

    with open(savePath + '/results.txt', 'w') as f:
      sum = 0.0
      for key, value in results.items():
        f.write("Fold %d: %f \n" % (key, value))
        sum += value
      f.write("Average: %f %%\n" % (sum/len(results.items())))

    with open(savePath + '/results_max.txt', 'w') as f:
      sum = 0.0
      for key, value in results_max.items():
        f.write("Fold %d: %f \n" % (key, value))
        sum += value
      f.write("Average: %f %%\n" % (sum/len(results_max.items())))


    # Print fold results
    print(f'K-FOLD CROSS VALIDATION RESULTS FOR {k_folds} FOLDS')
    print('--------------------------------')
    sum = 0.0
    for key, value in results.items():
      print(f'Fold {key}: {value} %')
      sum += value
    print(f'Average: {sum/len(results.items())} %')

    # Print fold results
    print(f'MAX:   K-FOLD CROSS VALIDATION RESULTS FOR {k_folds} FOLDS')
    print('--------------------------------')
    sum = 0.0
    for key, value in results_max.items():
      print(f'Fold {key}: {value} %')
      sum += value
    print(f'MAX Average: {sum/len(results.items())} %')


    # Print fold results [BEST model]
    if args.save_bestmodel:
      print(f'[BEST model] K-FOLD CROSS VALIDATION RESULTS FOR {k_folds} FOLDS')
      print('--------------------------------')
      sum = 0.0
      for key, value in results_best.items():
        print(f'Fold {key}: {value} %')
        sum += value
      print(f'Average: {sum/len(results_best.items())} %')

      with open(savePath + '/results_best.txt', 'w') as f:
        sum = 0.0
        for key, value in results_best.items():
          f.write("Fold %d: %f \n" % (key, value))
          sum += value
        f.write("Average: %f %%\n" % (sum/len(results_best.items())))


    # Plot the results.
    # Learning loss
    plt.figure(3)
    plt.semilogy(loss_hist_train, label='Training (last fold)')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title(savePath + '   chip-classifier')
    plt.legend()
    plt.savefig(savePath + '/loss-lastfold.png')

    # Learning accuracy
    plt.figure(4)
    plt.plot(acc_hist_train, label='Training (last fold)')
    #plt.yscale('log')
    plt.xlabel('Epoch')
    plt.ylabel('Accuracy')
    plt.title(savePath + '   chip-classifier')
    plt.legend()
    plt.savefig(savePath + '/acc-lastfold.png')

    #plt.show()

    cmd = 'rm -rf cache-' + extraPath + '/'
    os.system(cmd)





