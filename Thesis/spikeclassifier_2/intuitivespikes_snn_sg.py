# ============================================
# intuitivespikes_snn_sg.py  (FULL, MODIFIED)
# ============================================
# Texture
# Example:
# CUDA_VISIBLE_DEVICES=3 python intuitivespikes_snn_sg.py \
#   --datasetPath ../spikegen/global_outfile/spikegen/ALL_N_lcadc-0-uV-v2-0-uV-0-uV-0-uV-0-pct-4-mr_cnt-192-multitrial-Kylberg_filt_6_scan_actualdimscale_chip_lchpf--1to40-1to40_sorted/ \
#   --epochs 150 --hidden_size 256 --num_classes 6 --batch_size 32 --tmax_ms 400 --Ts 10 --kfolds 6 \
#   --cache yes --save_bestmodel yes \
#   --device cuda --amp no --num_workers 4
#
# NEW:
#   --outdir_prefix "SNN"              # default "SNN" (keeps old behavior)
#   --outdir_suffix ""                # default "" (keeps old behavior)
#
# To get your requested shorter output folder name:
#   --outdir_prefix "SNN_dsetPath" --outdir_suffix ""   (and it will create:)
#     <cwd>/SNN_dsetPath-<datasetFolderBasename>
#
# NEW (THIS REQUEST): if output directory is in short mode (outdir_prefix != "SNN"),
# the cache path will ALSO be short and include ONLY the dataset basename:
#   ./cache-dsetPath-<datasetFolderBasename>/spkdataset/{train,test}
#
# (Dataset folder name on disk is NOT changed.)

"""
Train 1-Hidden Layer LIF SNN using Surrogate Gradient
Trained model is used for SVM classification
Author: Ali Safa - IMEC- KU Leuven, Federico Corradi - IMEC-NL
Modified by: Mark Alea - KU Leuven
Further modified to actually use CUDA + speed up loaders + reduce sync.
Backward-compatible with older PyTorch (no inference_mode requirement).
04/02/2023
"""

import sys, os
CURRENT_TEST_DIR = os.getcwd()
sys.path.append(CURRENT_TEST_DIR + "/../../src")

import numpy as np
from eNetworks import mini_eMLP
import torch.nn as nn
import torch
from torch.utils.data import Dataset, DataLoader, ConcatDataset
import matplotlib.pyplot as plt
from sklearn.model_selection import KFold

import tonic
from tonic.io import read_mnist_file
from tonic import CachedDataset

import argparse

# -----------------------------
# Compatibility helpers (older torch)
# -----------------------------
_INFERENCE_CM = getattr(torch, "inference_mode", torch.no_grad)

_HAS_AMP = False
try:
    _HAS_AMP = hasattr(torch.cuda, "amp") and hasattr(torch.cuda.amp, "autocast") and hasattr(torch.cuda.amp, "GradScaler")
except Exception:
    _HAS_AMP = False

# -----------------------------
# Helpers
# -----------------------------
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
    # Default is "no" in this script; AMP only if supported + cuda
    amp_arg = (amp_arg or "no").lower()
    if device.type != "cuda":
        return False
    if not _HAS_AMP:
        return False
    if amp_arg == "yes":
        return True
    if amp_arg == "auto":
        return True
    return False

# -----------------------------
# Args
# -----------------------------
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

# CUDA speed / reproducibility knobs (CUDA core selection is via CUDA_VISIBLE_DEVICES outside)
parser.add_argument("--device", type=str, default="auto", choices=["auto", "cpu", "cuda"], help="Training device")
parser.add_argument("--amp", type=str, default="no", choices=["no", "yes", "auto"], help="Mixed precision (CUDA only, if supported)")
parser.add_argument("--deterministic", type=str, default="no", help="Try to make CUDA ops deterministic (slower)")
parser.add_argument("--num_workers", type=int, default=4, help="DataLoader workers")
parser.add_argument("--pin_memory", type=str, default="auto", choices=["auto", "yes", "no"], help="Pin memory for faster H2D copies")
parser.add_argument("--print_preds", type=str, default="no", help="Print batch preds/labels during test (SLOW)")
parser.add_argument("--reset_cache", type=str, default="yes", help="Reset CachedDataset (yes/no)")
parser.add_argument("--keep_cache", type=str, default="no", help="Keep cache folder at end (yes/no)")

# -----------------------------
# NEW: output folder naming controls
# -----------------------------
parser.add_argument(
    "--outdir_prefix",
    type=str,
    default="SNN",
    help="Output folder prefix. Default keeps old behavior: 'SNN' -> savePath = <cwd>/SNN_<extraPath>."
)
parser.add_argument(
    "--outdir_suffix",
    type=str,
    default="",
    help="Optional suffix appended to output folder name (after dataset tag). Default ''."
)

args = parser.parse_args()

# -----------------------------
# Hyper-parameters
# -----------------------------
ch_size = args.ch_size if args.ch_size else 2
x_size  = args.x_size if args.x_size else 16
y_size  = args.y_size if args.y_size else 12

# use custom ID to name retrain path
if args.retrainPath_id:
    retrainPath_tag = '-' + args.retrainPath_id
else:
    if args.retrainPath:
        tokens_slash = args.retrainPath.split('/')
        tokens_dot = tokens_slash[-1].split('.')
        tokens = tokens_dot[0].split('-')
        retrainPath_tag = '-retrainPath-' + tokens_slash[-2] + '-' + tokens[-2] + '-' + tokens[-1]
        print(retrainPath_tag)
    else:
        retrainPath_tag = ''

hidden_size = args.hidden_size
num_classes = args.num_classes
num_epochs = args.epochs
batch_size = args.batch_size

print('hidden_size'); print(hidden_size)
print('num_classes'); print(num_classes)
print('batch_size'); print(batch_size)

learning_rate = args.lr if args.lr else 1e-3

# dataset folder basename (used for output tag)
dset_basename = ""
if args.datasetPath:
    tokens_slash = args.datasetPath.split('/')
    dset_basename = tokens_slash[-2]
    datasetPath_tag = '-dsetPath-' + dset_basename
else:
    datasetPath_tag = ''

# ----------------------------------------
# Output folder naming (SHORTER option)
# ----------------------------------------
# OLD behavior:
#   savePath = <cwd> + "/SNN_" + extraPath
#   where extraPath includes ep/isize/.../k/... + datasetPath_tag
#
# NEW behavior requested by you:
#   savePath = <cwd> + "/SNN_dsetPath-" + <datasetFolderBasename>
#
# Implemented via --outdir_prefix:
#   --outdir_prefix "SNN_dsetPath"
# and ignore extraPath in that mode by building:
#   savePath = <cwd> + "/" + outdir_prefix + "-" + <datasetFolderBasename> + outdir_suffix
#
# Default keeps full naming (backward-compatible).
# ----------------------------------------

# Build extraPath as before (kept for backwards compatibility)
if args.retrainPath:
    extraPath = 'ep-' + str(args.epochs) + retrainPath_tag
else:
    extraPath = (
        'ep-' + str(args.epochs)
        + '-isize-' + str(ch_size) + '-' + str(x_size) + '-' + str(y_size)
        + '-hsize-' + str(args.hidden_size)
        + '-bsize-' + str(args.batch_size)
        + '-lr-' + str(learning_rate)
        + '-tmax_ms-' + str(args.tmax_ms)
        + '-Ts-' + str(args.Ts)
        + "-k-" + str(args.kfolds)
        + datasetPath_tag
    )

sampleFile_train = args.datasetPath + "train1K.txt"
sampleFile_test  = args.datasetPath + "test100.txt"

# NEW savePath logic
short_mode = (args.outdir_prefix != "SNN")

if not short_mode:
    # legacy mode (unchanged)
    savePath = path + "/SNN_" + extraPath
else:
    # requested short mode: "<cwd>/SNN_dsetPath-<datasetFolderBasename>"
    if args.datasetPath:
        savePath = path + "/" + args.outdir_prefix + "-" + dset_basename + str(args.outdir_suffix or "")
    else:
        savePath = path + "/" + args.outdir_prefix + str(args.outdir_suffix or "")

print(args.datasetPath)
print(savePath)
if not os.path.exists(savePath):
    os.mkdir(savePath)

# -----------------------------
# Dataset loading
# -----------------------------
class SpkDataset_toFrameTimeBins(Dataset):
    def __init__(self, datasetPath, sampleFile, samplingTime):
        self.path = datasetPath
        self.samples = np.loadtxt(sampleFile).astype('int')
        self.samplingTime = float(samplingTime)

    def __getitem__(self, index):
        inputIndex  = self.samples[index, 0]
        classLabel  = self.samples[index, 1]

        dtype = np.dtype([('x', '<i8'), ('y', '<i8'), ('t', 'int32'), ('p', '<i8')])
        events = read_mnist_file(self.path + str(int(inputIndex)) + ".bs2", dtype=dtype)

        frame_transform = tonic.transforms.ToFrame(
            sensor_size=(x_size, y_size, ch_size),
            n_time_bins=int(args.tmax_ms / args.Ts)
        )
        frames = frame_transform(events).astype(np.float32)

        if args.show_inputspikes:
            return frames, classLabel
        else:
            return frames.reshape(-1, ch_size * x_size * y_size), classLabel

    def __len__(self):
        return self.samples.shape[0]

criterion = nn.NLLLoss()

# -----------------------------
# Train / Test
# -----------------------------
def train(model, train_loader, test_loader, optimizer, epochs, batch_size,
          acc_hist_train, loss_hist_train, acc_hist_test, fold,
          device, use_amp, scaler):

    for e in range(epochs):
        model.train()

        train_correct = torch.zeros((), device=device, dtype=torch.long)
        train_loss_sum = torch.zeros((), device=device, dtype=torch.float32)
        total = 0

        for images, labels in train_loader:
            images = images.to(device, non_blocking=True)
            labels = labels.to(device, non_blocking=True)

            optimizer.zero_grad()

            # IMPORTANT: keep same shapes as your original code (1D state)
            hidden_mem   = torch.zeros(model.hidden_size_1, device=device)
            output_mem   = torch.zeros(model.output_size, device=device)
            hidden_spike = torch.zeros(model.hidden_size_1, device=device)
            output_spike = torch.zeros(model.output_size, device=device)

            if use_amp and _HAS_AMP:
                with torch.cuda.amp.autocast():
                    predictions, train_loss, _ = model(images, labels, hidden_mem, hidden_spike, output_mem, output_spike)
                    loss = train_loss.sum()
                scaler.scale(loss).backward()
                scaler.step(optimizer)
                scaler.update()
            else:
                predictions, train_loss, _ = model(images, labels, hidden_mem, hidden_spike, output_mem, output_spike)
                loss = train_loss.sum()
                loss.backward()
                optimizer.step()

            train_loss_sum = train_loss_sum + train_loss.detach().sum()

            predicted = predictions.t()
            train_correct = train_correct + (predicted.T[0] == labels).sum()
            total += labels.size(0)

        train_acc = (train_correct.float() / float(total)).item()
        train_loss_epoch = (train_loss_sum / float(total)).item()

        acc_hist_train.append(train_acc)
        loss_hist_train.append(train_loss_epoch)

        print(f"Fold: {fold} Epoch: {e} Loss: {train_loss_epoch} Accuracy: {train_acc}")

        # Test per epoch -> keep behavior, but use compatible inference context
        acc_hist_test = test(model, test_loader, batch_size, acc_hist_test, fold, device, use_amp)

        # SAVE best model so far (based on train acc)
        if args.save_bestmodel:
            if train_acc >= float(np.max(acc_hist_train)):
                print(f" -->  Saving best trained model : epoch {e} Accuracy: {train_acc}")
                torch.save(model.state_dict(), savePath + f'/chip-classifier-{fold}-best.pt')

    return acc_hist_test

def test(model, testloader, batch_size, acc_hist_test, fold, device, use_amp):
    model.eval()

    correct = torch.zeros((), device=device, dtype=torch.long)
    total = 0
    test_loss_sum = torch.zeros((), device=device, dtype=torch.float32)

    do_print = str2bool(args.print_preds)

    with _INFERENCE_CM():
        for images, labels in testloader:
            images = images.to(device, non_blocking=True)
            labels = labels.to(device, non_blocking=True)

            hidden_mem   = torch.zeros(model.hidden_size_1, device=device)
            output_mem   = torch.zeros(model.output_size, device=device)
            hidden_spike = torch.zeros(model.hidden_size_1, device=device)
            output_spike = torch.zeros(model.output_size, device=device)

            if use_amp and _HAS_AMP:
                with torch.cuda.amp.autocast():
                    predictions, test_loss, _ = model(images, labels, hidden_mem, hidden_spike, output_mem, output_spike)
            else:
                predictions, test_loss, _ = model(images, labels, hidden_mem, hidden_spike, output_mem, output_spike)

            test_loss_sum = test_loss_sum + test_loss.detach().sum()

            predicted = predictions.t()
            correct = correct + (predicted.T[0] == labels).sum()
            total += labels.size(0)

            if do_print:
                print(predicted.T[0])
                print(labels)

    test_acc = (correct.float() / float(total)).item()
    test_loss_epoch = (test_loss_sum / float(total)).item()

    acc_hist_test.append(test_acc)
    print(f"--> Test:  Loss: {test_loss_epoch} Accuracy: {test_acc}")

    # SAVE best model so far (based on test acc)
    if args.save_bestmodel:
        if test_acc >= float(np.max(acc_hist_test)):
            print(f" -->  Saving best trained model (based on test): Accuracy: {test_acc}")
            torch.save(model.state_dict(), savePath + f'/chip-classifier-test-{fold}-best.pt')

    return acc_hist_test

# -----------------------------
# Main
# -----------------------------
if __name__ == '__main__':

    # Device / AMP / determinism
    device = choose_device(args.device)
    use_amp = choose_amp(args.amp, device)
    deterministic = str2bool(args.deterministic)

    # Repro seed
    torch.manual_seed(42)
    np.random.seed(42)
    if device.type == "cuda":
        try:
            torch.cuda.manual_seed_all(42)
        except Exception:
            pass

    if device.type == "cuda":
        try:
            torch.cuda.empty_cache()
        except Exception:
            pass
        try:
            torch.backends.cudnn.benchmark = not deterministic
            torch.backends.cudnn.deterministic = deterministic
        except Exception:
            pass
        if deterministic:
            try:
                torch.use_deterministic_algorithms(True)
            except Exception:
                pass
        try:
            torch.set_float32_matmul_precision("high")
        except Exception:
            pass

    reset_cache = str2bool(args.reset_cache)
    keep_cache = str2bool(args.keep_cache)

    pin_memory = (device.type == "cuda") if args.pin_memory == "auto" else str2bool(args.pin_memory)

    print(f"[INFO] device={device} amp={use_amp} (amp_supported={_HAS_AMP}) deterministic={deterministic}")
    print(f"[INFO] num_workers={args.num_workers} pin_memory={pin_memory}")
    print(f"[INFO] short_mode={short_mode} outdir_prefix={args.outdir_prefix}")

    trainingSet = SpkDataset_toFrameTimeBins(
        datasetPath=args.datasetPath,
        sampleFile=sampleFile_train,
        samplingTime=args.Ts
    )
    testingSet = SpkDataset_toFrameTimeBins(
        datasetPath=args.datasetPath,
        sampleFile=sampleFile_test,
        samplingTime=args.Ts
    )

    # -----------------------------
    # Cache path naming (UPDATED)
    # -----------------------------
    # Legacy mode: cache tag uses full extraPath (old behavior)
    # Short mode: cache tag uses ONLY dataset basename (requested)
    if not short_mode:
        cache_tag = extraPath
    else:
        # Exactly what you asked: only include datasetPath basename
        cache_tag = "dsetPath-" + dset_basename if dset_basename else "dsetPath-UNKNOWN"

    cache_train_path = "./cache-" + cache_tag + "/spkdataset/train"
    cache_test_path  = "./cache-" + cache_tag + "/spkdataset/test"

    trainingSet_cached = CachedDataset(trainingSet, cache_path=cache_train_path, reset_cache=reset_cache)
    testingSet_cached  = CachedDataset(testingSet,  cache_path=cache_test_path,  reset_cache=reset_cache)

    # Optional visualization
    if args.show_inputspikes:
        for i in range(1):
            _input, classLabel = trainingSet[i]
            print(_input.shape)
            print(classLabel)
            plt.subplot(5,5,i+1)
            im = plt.imshow(_input[0][0])
        plt.colorbar(im, ax=plt.gcf().axes)
        plt.show()

        for i in range(40):
            _input, classLabel = trainingSet[0]
            plt.subplot(5,8,i+1)
            im = plt.imshow(_input[i][0])
        plt.colorbar(im, ax=plt.gcf().axes)
        plt.show()

    # K-fold
    k_folds = args.kfolds
    results = {}
    results_best = {}
    results_max = {}

    dataset = ConcatDataset([trainingSet_cached, testingSet_cached])
    kfold = KFold(n_splits=k_folds, shuffle=True)

    print('--------------------------------')

    for fold, (train_ids, test_ids) in enumerate(kfold.split(dataset)):

        print(f'FOLD {fold}')
        print('--------------------------------')

        train_subsampler = torch.utils.data.SubsetRandomSampler(train_ids)
        test_subsampler  = torch.utils.data.SubsetRandomSampler(test_ids)

        trainloader = DataLoader(
            dataset,
            batch_size=batch_size,
            sampler=train_subsampler,
            num_workers=args.num_workers,
            pin_memory=pin_memory,
        )
        testloader = DataLoader(
            dataset,
            batch_size=batch_size,
            sampler=test_subsampler,
            num_workers=args.num_workers,
            pin_memory=pin_memory,
        )

        model = mini_eMLP(ch_size, x_size, y_size, hidden_size, num_classes, criterion=criterion, batch_size=batch_size)
        model_best = mini_eMLP(ch_size, x_size, y_size, hidden_size, num_classes, criterion=criterion, batch_size=batch_size)

        model = model.to(device)
        model_best = model_best.to(device)

        if args.retrainPath:
            model.load_state_dict(torch.load(args.retrainPath, map_location=device))

        base_params = [
            model.i2h.weight, model.i2h.bias,
            model.h2o.weight, model.h2o.bias,
            model.h2h.weight, model.h2h.bias
        ]
        optimizer = torch.optim.Adam([{'params': base_params}], lr=learning_rate)

        scaler = None
        if use_amp and _HAS_AMP and device.type == "cuda":
            scaler = torch.cuda.amp.GradScaler()

        loss_hist_train = []
        loss_hist_test = []
        acc_hist_train = []
        acc_hist_test = []

        # Training
        acc_hist_test_ = train(
            model, trainloader, testloader, optimizer, num_epochs, batch_size,
            acc_hist_train, loss_hist_train, acc_hist_test, fold,
            device, use_amp, scaler
        )

        acc_hist_test_max = float(np.max(acc_hist_test_)) if len(acc_hist_test_) else 0.0
        epoch_max = int(np.argmax(acc_hist_test_)) if len(acc_hist_test_) else 0

        # Save model
        torch.save(model.state_dict(), savePath + f'/chip-classifier-{fold}.pt')

        # Load best model
        if args.save_bestmodel:
            model_best.load_state_dict(torch.load(savePath + f'/chip-classifier-{fold}-best.pt', map_location=device))

        # Final evaluation for this fold
        model.eval()
        model_best.eval()

        correct = torch.zeros((), device=device, dtype=torch.long)
        correct_best = torch.zeros((), device=device, dtype=torch.long)
        total = 0

        with _INFERENCE_CM():
            for images, labels in testloader:
                images = images.to(device, non_blocking=True)
                labels = labels.to(device, non_blocking=True)

                hidden_mem   = torch.zeros(model.hidden_size_1, device=device)
                output_mem   = torch.zeros(model.output_size, device=device)
                hidden_spike = torch.zeros(model.hidden_size_1, device=device)
                output_spike = torch.zeros(model.output_size, device=device)

                if use_amp and _HAS_AMP:
                    with torch.cuda.amp.autocast():
                        predictions, _, _ = model(images, labels, hidden_mem, hidden_spike, output_mem, output_spike)
                else:
                    predictions, _, _ = model(images, labels, hidden_mem, hidden_spike, output_mem, output_spike)

                predicted = predictions.t()
                correct = correct + (predicted.T[0] == labels).sum()
                total += labels.size(0)

                if args.save_bestmodel:
                    hidden_mem   = torch.zeros(model_best.hidden_size_1, device=device)
                    output_mem   = torch.zeros(model_best.output_size, device=device)
                    hidden_spike = torch.zeros(model_best.hidden_size_1, device=device)
                    output_spike = torch.zeros(model_best.output_size, device=device)

                    if use_amp and _HAS_AMP:
                        with torch.cuda.amp.autocast():
                            predictions_best, _, _ = model_best(images, labels, hidden_mem, hidden_spike, output_mem, output_spike)
                    else:
                        predictions_best, _, _ = model_best(images, labels, hidden_mem, hidden_spike, output_mem, output_spike)

                    predicted_best = predictions_best.t()
                    correct_best = correct_best + (predicted_best.T[0] == labels).sum()

        test_acc_ = (correct.float() / float(total)).item()
        results[fold] = 100.0 * test_acc_
        results_max[fold] = 100.0 * acc_hist_test_max

        print('Accuracy for fold %d: %f %%' % (fold, 100.0 * test_acc_))
        print('--------------------------------')
        print('--> At epoch = %d, Max: %f %%' % (epoch_max, 100.0 * acc_hist_test_max))

        if args.save_bestmodel:
            test_acc_best = (correct_best.float() / float(total)).item()
            results_best[fold] = 100.0 * test_acc_best
            print('[BEST model] Accuracy for fold %d: %f %%' % (fold, 100.0 * test_acc_best))
            print('--------------------------------')

        # Save histories
        with open(savePath + f'/loss_hist_train-{fold}.txt', 'w') as f:
            for item in loss_hist_train:
                f.write("%f\n" % item)

        with open(savePath + f'/acc_hist_train-{fold}.txt', 'w') as f:
            for item in acc_hist_train:
                f.write("%f\n" % item)

        with open(savePath + f'/acc_hist_test-{fold}.txt', 'w') as f:
            for item in acc_hist_test:
                f.write("%f\n" % item)

    with open(savePath + '/results.txt', 'w') as f:
        s = 0.0
        for key, value in results.items():
            f.write("Fold %d: %f \n" % (key, value))
            s += value
        f.write("Average: %f %%\n" % (s / len(results.items())))

    with open(savePath + '/results_max.txt', 'w') as f:
        s = 0.0
        for key, value in results_max.items():
            f.write("Fold %d: %f \n" % (key, value))
            s += value
        f.write("Average: %f %%\n" % (s / len(results_max.items())))

    print(f'K-FOLD CROSS VALIDATION RESULTS FOR {k_folds} FOLDS')
    print('--------------------------------')
    s = 0.0
    for key, value in results.items():
        print(f'Fold {key}: {value} %')
        s += value
    print(f'Average: {s / len(results.items())} %')

    print(f'MAX:   K-FOLD CROSS VALIDATION RESULTS FOR {k_folds} FOLDS')
    print('--------------------------------')
    s = 0.0
    for key, value in results_max.items():
        print(f'Fold {key}: {value} %')
        s += value
    print(f'MAX Average: {s / len(results_max.items())} %')

    if args.save_bestmodel:
        print(f'[BEST model] K-FOLD CROSS VALIDATION RESULTS FOR {k_folds} FOLDS')
        print('--------------------------------')
        s = 0.0
        for key, value in results_best.items():
            print(f'Fold {key}: {value} %')
            s += value
        print(f'Average: {s / len(results_best.items())} %')

        with open(savePath + '/results_best.txt', 'w') as f:
            s = 0.0
            for key, value in results_best.items():
                f.write("Fold %d: %f \n" % (key, value))
                s += value
            f.write("Average: %f %%\n" % (s / len(results_best.items())))

    # Plots (last fold histories)
    plt.figure(3)
    plt.semilogy(loss_hist_train, label='Training (last fold)')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title(savePath + '   chip-classifier')
    plt.legend()
    plt.savefig(savePath + '/loss-lastfold.png')

    plt.figure(4)
    plt.plot(acc_hist_train, label='Training (last fold)')
    plt.xlabel('Epoch')
    plt.ylabel('Accuracy')
    plt.title(savePath + '   chip-classifier')
    plt.legend()
    plt.savefig(savePath + '/acc-lastfold.png')

    # Remove cache folder unless requested
    if not keep_cache:
        cmd = 'rm -rf cache-' + cache_tag + '/'
        os.system(cmd)

