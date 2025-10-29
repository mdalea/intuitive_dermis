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
runID = 'fnoise_sweep_trainonclean'

datasetID_clean = 'lcadc-0-uV2-0-uV2-0-uV2-0-mult-0-uV2-4'
## FLICKER NOISE SWEEP --------------------------
#-----------------------------
datasetID = 'lcadc-0-uV2-0-uV2-0-uV2-0.10857-mult-0-uV2-4'

#os.system(f"CUDA_VISIBLE_DEVICES={cuda} python intuitivespikes_snn_sg.py --datasetPath ../spikegen/global_outfile_tbcas24/spikegen/ALL_N_{datasetID}-mr_cnt-192-multitrial-Kylberg_filt_6_scan_actualdimscale_chip_lchpf--1to40-1to40_sorted/ --epochs 150 --hidden_size 256 --num_classes 6 --batch_size 32 --tmax_ms 400 --Ts 10 --kfolds 6 --cache yes --save_bestmodel yes")

for fold in range(kfolds):
    os.system(f"CUDA_VISIBLE_DEVICES={cuda} python intuitivespikes_svmclassifier.py --datasetPath ../spikegen/global_outfile_tbcas24/spikegen/ALL_N_{datasetID}-mr_cnt-192-multitrial-Kylberg_filt_6_scan_actualdimscale_chip_lchpf--1to40-1to40_sorted/ --epochs 5 --hidden_size 256 --num_classes 6 --batch_size 32 --tmax_ms 400 --Ts 10 --kfolds 10 --cache yes --svm_window_size 40 --retrainPath SNN_ep-150-isize-2-16-12-hsize-256-bsize-32-lr-0.001-tmax_ms-400.0-Ts-10.0-k-6-dsetPath-ALL_N_{datasetID_clean}-mr_cnt-192-multitrial-Kylberg_filt_6_scan_actualdimscale_chip_lchpf--1to40-1to40_sorted/chip-classifier-test-{fold}-best.pt --retrainPath_id {runID}-{datasetID}-best-{fold}")
#-----------------------------
datasetID = 'lcadc-0-uV2-0-uV2-0-uV2-10.8574-mult-0-uV2-4'

#os.system(f"CUDA_VISIBLE_DEVICES={cuda} python intuitivespikes_snn_sg.py --datasetPath ../spikegen/global_outfile_tbcas24/spikegen/ALL_N_{datasetID}-mr_cnt-192-multitrial-Kylberg_filt_6_scan_actualdimscale_chip_lchpf--1to40-1to40_sorted/ --epochs 150 --hidden_size 256 --num_classes 6 --batch_size 32 --tmax_ms 400 --Ts 10 --kfolds 6 --cache yes --save_bestmodel yes")

for fold in range(kfolds):
    os.system(f"CUDA_VISIBLE_DEVICES={cuda} python intuitivespikes_svmclassifier.py --datasetPath ../spikegen/global_outfile_tbcas24/spikegen/ALL_N_{datasetID}-mr_cnt-192-multitrial-Kylberg_filt_6_scan_actualdimscale_chip_lchpf--1to40-1to40_sorted/ --epochs 5 --hidden_size 256 --num_classes 6 --batch_size 32 --tmax_ms 400 --Ts 10 --kfolds 10 --cache yes --svm_window_size 40 --retrainPath SNN_ep-150-isize-2-16-12-hsize-256-bsize-32-lr-0.001-tmax_ms-400.0-Ts-10.0-k-6-dsetPath-ALL_N_{datasetID_clean}-mr_cnt-192-multitrial-Kylberg_filt_6_scan_actualdimscale_chip_lchpf--1to40-1to40_sorted/chip-classifier-test-{fold}-best.pt --retrainPath_id {runID}-{datasetID}-best-{fold}")

#-----------------------------
datasetID = 'lcadc-0-uV2-0-uV2-0-uV2-1085.7362-mult-0-uV2-4'

#os.system(f"CUDA_VISIBLE_DEVICES={cuda} python intuitivespikes_snn_sg.py --datasetPath ../spikegen/global_outfile_tbcas24/spikegen/ALL_N_{datasetID}-mr_cnt-192-multitrial-Kylberg_filt_6_scan_actualdimscale_chip_lchpf--1to40-1to40_sorted/ --epochs 150 --hidden_size 256 --num_classes 6 --batch_size 32 --tmax_ms 400 --Ts 10 --kfolds 6 --cache yes --save_bestmodel yes")

for fold in range(kfolds):
    os.system(f"CUDA_VISIBLE_DEVICES={cuda} python intuitivespikes_svmclassifier.py --datasetPath ../spikegen/global_outfile_tbcas24/spikegen/ALL_N_{datasetID}-mr_cnt-192-multitrial-Kylberg_filt_6_scan_actualdimscale_chip_lchpf--1to40-1to40_sorted/ --epochs 5 --hidden_size 256 --num_classes 6 --batch_size 32 --tmax_ms 400 --Ts 10 --kfolds 10 --cache yes --svm_window_size 40 --retrainPath SNN_ep-150-isize-2-16-12-hsize-256-bsize-32-lr-0.001-tmax_ms-400.0-Ts-10.0-k-6-dsetPath-ALL_N_{datasetID_clean}-mr_cnt-192-multitrial-Kylberg_filt_6_scan_actualdimscale_chip_lchpf--1to40-1to40_sorted/chip-classifier-test-{fold}-best.pt --retrainPath_id {runID}-{datasetID}-best-{fold}")

## FLICKER NOISE SWEEP --------------------------
#-----------------------------


