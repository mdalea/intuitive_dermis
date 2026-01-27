"""
Train SVM Classifier on the Filtered Output Spikes of the Hidden Layer
Texture and Flutter Frequency Classification
Author: Mark Alea - KU Leuven
Further modified: GPU-accelerated feature extraction (hidden spikes), faster dataloaders, less sync.
04/02/2023
"""

import sys, os
CURRENT_TEST_DIR = os.getcwd()
sys.path.append(CURRENT_TEST_DIR + "/../../src")

import numpy as np
from eNetworks_svmonly import mini_eMLP
import torch
from torch.utils.data import Dataset, DataLoader, ConcatDataset
import matplotlib.pyplot as plt
from sklearn.model_selection import KFold

import tonic
from tonic.io import read_mnist_file
from tonic import CachedDataset

from sklearn import svm
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV

import argparse


# ---- compatibility (older torch) ----
_INFERENCE_CM = getattr(torch, "inference_mode", torch.no_grad)
_HAS_AMP = hasattr(torch.cuda, "amp") and hasattr(torch.cuda.amp, "autocast")


def str2bool(v):
    if v is None:
        return False
    if isinstance(v, bool):
        return v
    v = str(v).strip().lower()
    if v in ("1", "true", "t", "yes", "y", "on"):
        return True
    if v in ("0", "false", "f", "no", "n", "off"):
        return False
    raise argparse.ArgumentTypeError(f"Expected bool-like value, got: {v}")


def choose_device(device_arg: str) -> torch.device:
    device_arg = (device_arg or "auto").lower()
    if device_arg == "cpu":
        return torch.device("cpu")
    if device_arg == "cuda":
        return torch.device("cuda")
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def choose_amp(amp_arg: str, device: torch.device) -> bool:
    amp_arg = (amp_arg or "no").lower()
    if device.type != "cuda":
        return False
    if not _HAS_AMP:
        return False
    return amp_arg in ("yes", "auto")


path = os.getcwd()
parser = argparse.ArgumentParser("SVM Classifier on Filtered Hidden Layer Spikes from Trained SNN")
parser.add_argument("--datasetPath", type=str, required=True)
parser.add_argument("--epochs", type=int, required=True)
parser.add_argument("--hidden_size", type=int, required=True)
parser.add_argument("--ch_size", type=int, required=False)
parser.add_argument("--x_size", type=int, required=False)
parser.add_argument("--y_size", type=int, required=False)
parser.add_argument("--num_classes", type=int, required=True)
parser.add_argument("--batch_size", type=int, required=True)
parser.add_argument("--lr", type=float, required=False)
parser.add_argument("--tmax_ms", type=float, required=True)
parser.add_argument("--Ts", type=float, required=True)
parser.add_argument("--show_inputspikes", type=str, required=False)
parser.add_argument("--retrainPath", type=str, required=False)
parser.add_argument("--kfolds", type=int, required=True)
parser.add_argument("--save_bestmodel", type=str, required=False)
parser.add_argument("--cache", type=str, required=True)
parser.add_argument("--svm_window_size", type=int, required=False)
parser.add_argument("--scaler", type=str, required=False)
parser.add_argument("--retrainPath_id", type=str, required=False)

# GPU knobs (CUDA core selection is via CUDA_VISIBLE_DEVICES outside)
parser.add_argument("--device", type=str, default="auto", choices=["auto", "cpu", "cuda"])
parser.add_argument("--amp", type=str, default="no", choices=["no", "yes", "auto"])
parser.add_argument("--num_workers", type=int, default=4)
parser.add_argument("--pin_memory", type=str, default="auto", choices=["auto", "yes", "no"])
parser.add_argument("--persistent_workers", type=str, default="yes")
parser.add_argument("--reset_cache", type=str, default="yes")
parser.add_argument("--keep_cache", type=str, default="no")

# Debug prints
parser.add_argument("--print_preds", type=str, default="yes", help="Print labels + predicted (SLOW)")

args = parser.parse_args()

# Hyper-parameters
ch_size = args.ch_size if args.ch_size else 2
x_size = args.x_size if args.x_size else 16
y_size = args.y_size if args.y_size else 12
window_size = args.svm_window_size if args.svm_window_size else 100

use_scaler = bool(args.scaler)
print_preds = str2bool(args.print_preds)

# run-name tags
if args.retrainPath_id:
    retrainPath_tag = '-' + args.retrainPath_id
else:
    if args.retrainPath:
        tokens_slash = args.retrainPath.split('/')
        tokens_dot = tokens_slash[-1].split('.')
        tokens = tokens_dot[0].split('-')
        retrainPath_tag = '-retrainPath-' + tokens_slash[-2] + '-' + tokens[-2] + '-' + tokens[-1]
    else:
        retrainPath_tag = ''

learning_rate = args.lr if args.lr else 1e-3

if args.datasetPath:
    tokens_slash = args.datasetPath.split('/')
    datasetPath_tag = '-dsetPath-' + tokens_slash[-2]
else:
    datasetPath_tag = ''

if args.retrainPath:
    extraPath = 'svmep-' + str(args.epochs) + '-svm_wsize-' + str(window_size) + ('-scaler-' if use_scaler else '') + retrainPath_tag
else:
    extraPath = (
        'svmep-' + str(args.epochs)
        + '-isize-' + str(ch_size) + '-' + str(x_size) + '-' + str(y_size)
        + '-hsize-' + str(args.hidden_size)
        + '-bsize-' + str(args.batch_size)
        + '-swsize-' + str(window_size)
        + '-lr-' + str(learning_rate)
        + '-tmax_ms-' + str(args.tmax_ms)
        + '-Ts-' + str(args.Ts)
        + "-k-" + str(args.kfolds)
        + ('-scaler-' if use_scaler else '')
        + retrainPath_tag
        + datasetPath_tag
    )

sampleFile_train = args.datasetPath + "train1K.txt"
sampleFile_test = args.datasetPath + "test100.txt"
savePath = path + "/SVM_" + extraPath

print(args.datasetPath)
print(savePath)
if not os.path.exists(savePath):
    os.mkdir(savePath)


# Dataset loading
class SpkDataset_toFrameTimeBins(Dataset):
    def __init__(self, datasetPath, sampleFile, samplingTime):
        self.path = datasetPath
        self.samples = np.loadtxt(sampleFile).astype('int')
        self.samplingTime = float(samplingTime)

    def __getitem__(self, index):
        inputIndex = self.samples[index, 0]
        classLabel = self.samples[index, 1]

        dtype = np.dtype([('x', '<i8'), ('y', '<i8'), ('t', 'int32'), ('p', '<i8')])
        events = read_mnist_file(self.path + str(int(inputIndex)) + ".bs2", dtype=dtype)

        frame_transform = tonic.transforms.ToFrame(
            sensor_size=(x_size, y_size, ch_size),
            n_time_bins=int(args.tmax_ms / args.Ts)
        )
        frames = frame_transform(events).astype(np.float32)

        return frames.reshape(-1, ch_size * x_size * y_size), classLabel

    def __len__(self):
        return self.samples.shape[0]


def extract_features(model, loader, device, use_amp):
    """
    Run the SNN forward on GPU to get hidden_spike_filt features.
    Return:
      X: np.ndarray [N, hidden_size * timewindows]
      y: np.ndarray [N]
    """
    model.eval()
    feats = []
    ys = []

    with _INFERENCE_CM():
        for images, labels in loader:
            # images: [B,T,input]
            if device.type == "cuda":
                images = images.to(device, non_blocking=True)
                labels_t = labels.to(device, non_blocking=True)
            else:
                labels_t = labels

            # states on same device
            hidden_mem = torch.zeros(model.hidden_size_1, device=device)
            output_mem = torch.zeros(model.output_size, device=device)
            hidden_spike = torch.zeros(model.hidden_size_1, device=device)
            output_spike = torch.zeros(model.output_size, device=device)

            if use_amp and device.type == "cuda" and _HAS_AMP:
                with torch.cuda.amp.autocast():
                    hidden_spike_filt = model(images, labels_t, hidden_mem, hidden_spike, output_mem, output_spike)
            else:
                hidden_spike_filt = model(images, labels_t, hidden_mem, hidden_spike, output_mem, output_spike)

            # [B, H, Tw] -> [B, H*Tw]
            Xb = hidden_spike_filt.detach().float().cpu().numpy()
            Xb = Xb.reshape(Xb.shape[0], -1)

            feats.append(Xb)
            ys.append(labels.detach().cpu().numpy())

    X = np.concatenate(feats, axis=0) if len(feats) else np.zeros((0, 1), dtype=np.float32)
    y = np.concatenate(ys, axis=0) if len(ys) else np.zeros((0,), dtype=np.int64)
    return X, y


def train_test(model, train_loader, test_loader, epochs, grid_search_linear, device, use_amp):
    acc_hist_train = []
    acc_hist_test = []

    scaler = StandardScaler()

    for e in range(epochs):
        # --- Train features ---
        X_train, y_train = extract_features(model, train_loader, device, use_amp)
        if use_scaler:
            X_train = scaler.fit_transform(X_train)

        # Fit SVM (CPU)
        grid_search_linear.fit(X_train, y_train)
        pred_train = grid_search_linear.predict(X_train)

        if print_preds:
            print(torch.from_numpy(y_train))
            print(pred_train)

        train_acc = (pred_train == y_train).mean()
        print(f"Epoch {e} best_params={grid_search_linear.best_params_} train_acc={train_acc:.4f} best_cv={grid_search_linear.best_score_:.4f}")
        acc_hist_train.append(train_acc)

        # --- Test features ---
        X_test, y_test = extract_features(model, test_loader, device, use_amp)
        if use_scaler:
            X_test = scaler.transform(X_test)

        pred_test = grid_search_linear.predict(X_test)
        if print_preds:
            print(torch.from_numpy(y_test))
            print(pred_test)

        test_acc = (pred_test == y_test).mean()
        print(f"Epoch {e} test_acc={test_acc:.4f}")
        acc_hist_test.append(test_acc)

    return acc_hist_train, acc_hist_test


if __name__ == '__main__':

    # Device / AMP
    device = choose_device(args.device)
    use_amp = choose_amp(args.amp, device)

    reset_cache = str2bool(args.reset_cache)
    keep_cache = str2bool(args.keep_cache)

    pin_memory = (device.type == "cuda") if args.pin_memory == "auto" else str2bool(args.pin_memory)
    persistent_workers = str2bool(args.persistent_workers) and (args.num_workers > 0)

    print(f"[INFO] device={device} amp={use_amp} num_workers={args.num_workers} pin_memory={pin_memory} persistent_workers={persistent_workers}")

    trainingSet = SpkDataset_toFrameTimeBins(args.datasetPath, sampleFile_train, args.Ts)
    testingSet = SpkDataset_toFrameTimeBins(args.datasetPath, sampleFile_test, args.Ts)

    trainingSet_cached = CachedDataset(trainingSet, cache_path='./cache-' + extraPath + '/spkdataset/train', reset_cache=reset_cache)
    testingSet_cached  = CachedDataset(testingSet,  cache_path='./cache-' + extraPath + '/spkdataset/test',  reset_cache=reset_cache)

    # Optional visualize
    if args.show_inputspikes:
        inp, classLabel = trainingSet[0]
        inp = inp.reshape(-1, ch_size, y_size, x_size)
        print(inp.shape, classLabel)
        plt.imshow(inp[0][0])
        plt.show()

    k_folds = args.kfolds
    torch.manual_seed(42)
    np.random.seed(42)

    dataset = ConcatDataset([trainingSet_cached, testingSet_cached])
    kfold = KFold(n_splits=k_folds, shuffle=True, random_state=42)

    # SVM grid search
    parameters_linear = {'C': [0.1, 1, 10, 100]}
    linear_svm = svm.SVC(kernel='linear')
    grid_search_linear = GridSearchCV(linear_svm, parameters_linear, cv=kfold)

    results = {}
    results_max = {}

    print('--------------------------------')

    for fold, (train_ids, test_ids) in enumerate(kfold.split(dataset)):
        print(f'FOLD {fold}')
        print('--------------------------------')

        train_subsampler = torch.utils.data.SubsetRandomSampler(train_ids)
        test_subsampler  = torch.utils.data.SubsetRandomSampler(test_ids)

        trainloader = DataLoader(
            dataset,
            batch_size=args.batch_size,
            sampler=train_subsampler,
            num_workers=args.num_workers,
            pin_memory=pin_memory,
            persistent_workers=persistent_workers,
        )
        testloader = DataLoader(
            dataset,
            batch_size=args.batch_size,
            sampler=test_subsampler,
            num_workers=args.num_workers,
            pin_memory=pin_memory,
            persistent_workers=persistent_workers,
        )

        # Init SNN (used only for feature extraction)
        model = mini_eMLP(ch_size, x_size, y_size, args.hidden_size, args.num_classes, window_size=window_size).to(device)

        if args.retrainPath:
            model.load_state_dict(torch.load(args.retrainPath, map_location=device))

        acc_hist_train, acc_hist_test = train_test(model, trainloader, testloader, args.epochs, grid_search_linear, device, use_amp)

        acc_max = float(np.max(acc_hist_test)) if len(acc_hist_test) else 0.0
        epoch_max = int(np.argmax(acc_hist_test)) if len(acc_hist_test) else 0

        results[fold] = 100.0 * float(acc_hist_test[-1])
        results_max[fold] = 100.0 * acc_max

        print('Linear SVM Accuracy for fold %d: %f %%' % (fold, results[fold]))
        print('--> At epoch = %d, Linear SVM Max: %f %%' % (epoch_max, results_max[fold]))
        print('--------------------------------')

    print(f'K-FOLD CROSS VALIDATION RESULTS FOR {k_folds} FOLDS - SVM Linear')
    print('--------------------------------')
    avg = sum(results.values()) / len(results) if len(results) else 0.0
    for k, v in results.items():
        print(f'Fold {k}: {v} %')
    print(f'Average: {avg} %')

    print(f'MAX:   K-FOLD CROSS VALIDATION RESULTS FOR {k_folds} FOLDS - SVM Linear')
    print('--------------------------------')
    avgm = sum(results_max.values()) / len(results_max) if len(results_max) else 0.0
    for k, v in results_max.items():
        print(f'Fold {k}: {v} %')
    print(f'MAX Average: {avgm} %')

    with open(savePath + '/results_svm_linear.txt', 'w') as f:
        for k, v in results.items():
            f.write("Fold %d: %f \n" % (k, v))
        f.write("Average: %f %%\n" % avg)

    with open(savePath + '/results_svm_linear_max.txt', 'w') as f:
        for k, v in results_max.items():
            f.write("Fold %d: %f \n" % (k, v))
        f.write("Average: %f %%\n" % avgm)

    # Plots (last fold)
    plt.figure(1)
    plt.plot(acc_hist_train, label='Training (last fold) - SVM Linear')
    plt.xlabel('Epoch'); plt.ylabel('Accuracy')
    plt.title(savePath + '   svm-linear-train')
    plt.legend()
    plt.savefig(savePath + '/acc-svm-linear-train-lastfold.png')

    plt.figure(3)
    plt.plot(acc_hist_test, label='Test (last fold) - SVM Linear')
    plt.xlabel('Epoch'); plt.ylabel('Accuracy')
    plt.title(savePath + '   svm-linear-test')
    plt.legend()
    plt.savefig(savePath + '/acc-svm-linear-test-lastfold.png')

    if not keep_cache:
        cmd = 'rm -rf cache-' + extraPath + '/'
        os.system(cmd)

