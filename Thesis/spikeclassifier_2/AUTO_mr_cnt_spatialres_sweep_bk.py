import os
import numpy as np
import argparse

path = os.getcwd()
parser = argparse.ArgumentParser("Automate SNN and SVM")
parser.add_argument("--cuda",type=int,help="CUDA core used", required=False)

args = parser.parse_args()
'''
parser.add_argument("--epochs",type=int,help="No. of epochs of training",required=False)
parser.add_argument("--batch_size",type=int,help="Training batch size",required=False)
parser.add_argument("--hidden_size",type=int,help="MLP hidden size",required=False)
parser.add_argument("--tmax_ms",type=float,help="Maximum sensor recording considered",required=False)
parser.add_argument("--Ts",type=float,help="Sampling time",required=False)
parser.add_argument("--kfolds",type=int,help="No. of k-folds",required=False)
parser.add_argument("--svm_window_size",type=int,help="SVM window size",required=False)



if args.epochs:
    epochs = args.epochs
else:
    epochs = 100

if args.svm_window_size:
    window_size = args.svm_window_size
else:
    window_size = 100
'''

if args.cuda:
    cuda = args.cuda
else:
    cuda = 0

kfolds=6
runID = 'sensor_sweep'


## CLEAN
datasetID = 'lcadc-0-uV2-0-uV2-0-uV2-0-mult-0-uV2-3'
mr_cnt = '9'
spatialres = '0p025'

os.system(f"CUDA_VISIBLE_DEVICES={cuda} python intuitivespikes_snn_sg.py --datasetPath ../spikegen/global_outfile/spikegen/ALL_N_{datasetID}-multitrial-Kylberg_filt_6_scan_oversampled20x_actualdimscale--1to40-1to40_sorted-mr_cnt-{mr_cnt}-spatialres-{spatialres}/ --epochs 150 --hidden_size 256 --num_classes 6 --batch_size 32 --tmax_ms 400 --Ts 10 --kfolds 6 --cache yes --save_bestmodel yes")

for fold in range(kfolds):
    os.system(f"CUDA_VISIBLE_DEVICES={cuda} python intuitivespikes_svmclassifier.py --datasetPath ../spikegen/global_outfile/spikegen/ALL_N_{datasetID}-multitrial-Kylberg_filt_6_scan_oversampled20x_actualdimscale--1to40-1to40_sorted-mr_cnt-{mr_cnt}-spatialres-{spatialres}/ --epochs 5 --hidden_size 256 --num_classes 6 --batch_size 32 --tmax_ms 400 --Ts 10 --kfolds 10 --cache yes --svm_window_size 40 --retrainPath SNN_ep-150-isize-2-16-12-hsize-256-bsize-32-lr-0.001-tmax_ms-400.0-Ts-10.0-k-6-dsetPath-ALL_N_{datasetID}-multitrial-Kylberg_filt_6_scan_oversampled20x_actualdimscale--1to40-1to40_sorted-mr_cnt-{mr_cnt}-spatialres-{spatialres}/chip-classifier-test-{fold}-best.pt --retrainPath_id {runID}-{datasetID}-mr_cnt-{mr_cnt}-spatialres-{spatialres}-best-{fold}")


#-----------------------------

