# spatialres_array = [0.025 0.05 0.1 0.2 0.4 1 2 5];   % example spatial resolutions
# %spatialres_array = [2 5];   % example spatial resolutions
# mr_cnt_array = [3*3 8*8 14*14 16*16 18*18];

# Python equivalents (keeping the same sweep values)
spatialres_array = [0.025, 0.05, 0.1, 0.2, 0.4, 1, 2, 5]  # example spatial resolutions
#spatialres_array = [1]  # example spatial resolutions

#mr_cnt_array = [3*3, 8*8, 14*14, 16*16, 18*18]
mr_cnt_array = [18*18]
#mr_cnt_array = [3*3, 8*8, 14*14]  # quick test

import os
import sys
import math
import argparse
import subprocess

# -----------------------------
# Args
# -----------------------------
parser = argparse.ArgumentParser("Automate SNN and SVM sweeps")
parser.add_argument("--cuda", type=int, help="CUDA device index to use", required=False)
args = parser.parse_args()

cuda = args.cuda if args.cuda is not None else 0

# -----------------------------
# Fixed experiment params
# -----------------------------
epochs = 150
hidden_size = 256
batch_size = 32
lr = 0.001
tmax_ms = 400.0
Ts = 10.0
kfolds = 6

runID = "sensor_sweep"
datasetID = "lcadc-0-uV2-0-uV2-0-uV2-0-mult-0-uV2-3"

def fmt_float_1(x: float) -> str:
    # matches your example: 400.0, 10.0
    return f"{x:.1f}"

def fmt_spatialres(val: float) -> str:
    """
    Folder naming uses DOTS: 0.025 -> '0.025', 0.05 -> '0.05', 1 -> '1'
    Strip trailing zeros so 1.0 becomes '1' and 2.0 becomes '2'.
    """
    return f"{val:.12f}".rstrip("0").rstrip(".")

def run_cmd_or_die(cmd: str) -> None:
    """
    Run a shell command and terminate the WHOLE sweep immediately on any failure.
    - Prints the command before running.
    - Exits with the same non-zero return code if it fails.
    """
    print("\n[RUN]", cmd, flush=True)
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"\n[ERROR] Command failed with return code {e.returncode}:\n{cmd}\n", file=sys.stderr, flush=True)
        sys.exit(e.returncode if e.returncode is not None else 1)
    except Exception as e:
        print(f"\n[ERROR] Unexpected exception while running:\n{cmd}\nException: {e}\n", file=sys.stderr, flush=True)
        sys.exit(1)

for mr_cnt_val in mr_cnt_array:
    mr_cnt_int = int(mr_cnt_val)
    mr_cnt = str(mr_cnt_int)

    # -------------------------------------------------------------------------
    # !!! ADDRESSING MODEL FIX (THIS IS THE PART YOU MESSED UP BEFORE) !!!
    #
    # You HARD-FIXED the X address range to:
    #   x starts at 0 and x_max is ALWAYS 15  -> 16 X positions total.
    #
    # That means spikes are packed like this (row-major in X):
    #   addr = y * 16 + x
    #   x    = addr % 16
    #   y    = addr // 16
    #
    # BIG CONSEQUENCE:
    #   You only move to the next Y after X has gone 0..15.
    #   So each Y-row can hold 16 addresses.
    #
    # Therefore the CORRECT sizing is:
    #   x_size = 16  (ALWAYS)
    #   y_size = ceil(mr_cnt / 16)  (MINIMUM rows needed to fit mr_cnt addresses)
    #
    # This fixes the specific case you pointed out:
    #   mr_cnt = 9  -> ceil(9/16)=1  -> y_size=1 ✅ (fits entirely in first row)
    #
    # NOTE ABOUT YOUR EARLIER 14*14 EXAMPLE:
    #   mr_cnt = 196 -> ceil(196/16)=13 -> y_size=13 -> y_max=12 (0..12) ✅
    # If you previously wrote "y_max=13", that would be *padding* (14 rows) and wastes space.
    # If you ever WANT padding again, do it explicitly; don't bake it into y_size.
    # -------------------------------------------------------------------------

    X_MAX = 15
    x_size = X_MAX + 1  # ALWAYS 16 columns because x is fixed to [0..15]

    # Minimal number of Y rows required so total capacity >= mr_cnt
    y_size = math.ceil(mr_cnt_int / x_size)

    # Defensive sanity check so you never silently under-allocate addresses again.
    if x_size * y_size < mr_cnt_int:
        raise RuntimeError(
            f"Address capacity too small: x_size*y_size={x_size*y_size} < mr_cnt={mr_cnt_int}. "
            f"(This should never happen with ceil(mr_cnt/16).)"
        )

    # Optional: keep this just as a warning because your sweeps are usually squares.
    side = int(round(math.sqrt(mr_cnt_int)))
    if side * side != mr_cnt_int:
        print(f"[WARN] mr_cnt={mr_cnt_int} is not a perfect square (fine with fixed X=16 addressing).", flush=True)

    # Build the SAME run-name your SNN training creates (now uses computed x_size/y_size)
    snn_run_name = (
        f"SNN_ep-{epochs}-isize-2-{x_size}-{y_size}"
        f"-hsize-{hidden_size}-bsize-{batch_size}-lr-{lr}"
        f"-tmax_ms-{fmt_float_1(tmax_ms)}-Ts-{fmt_float_1(Ts)}-k-{kfolds}-"
    )

    for spatialres_val in spatialres_array:
        spatialres = fmt_spatialres(spatialres_val)

        # Your dataset path pattern (note the "_sorted/" suffix)
        datasetPath = (
            f"../spikegen/global_outfile/spikegen/ALL_N_{datasetID}"
            f"-multitrial-Kylberg_filt_6_scan_oversampled20x_actualdimscale--1to40-1to40"
            f"-mr_cnt-{mr_cnt}-spatialres-{spatialres}_sorted/"
        )

        # ---------------- SNN ----------------
        run_cmd_or_die(
            f"CUDA_VISIBLE_DEVICES={cuda} python intuitivespikes_snn_sg.py "
            f"--datasetPath {datasetPath} "
            f"--epochs {epochs} --hidden_size {hidden_size} --num_classes 6 --batch_size {batch_size} "
            f"--tmax_ms {tmax_ms} --Ts {Ts} --kfolds {kfolds} --cache yes --save_bestmodel yes "
            f"--x_size {x_size} --y_size {y_size}"
        )

        # ---------------- SVM (per fold) ----------------
        for fold in range(kfolds):
            retrainPath = (
                f"{snn_run_name}"
                f"dsetPath-ALL_N_{datasetID}-multitrial-Kylberg_filt_6_scan_oversampled20x_actualdimscale--1to40-1to40"
                f"-mr_cnt-{mr_cnt}-spatialres-{spatialres}_sorted/"
                f"chip-classifier-test-{fold}-best.pt"
            )

            run_cmd_or_die(
                f"CUDA_VISIBLE_DEVICES={cuda} python intuitivespikes_svmclassifier.py "
                f"--datasetPath {datasetPath} "
                f"--epochs 5 --hidden_size {hidden_size} --num_classes 6 --batch_size {batch_size} "
                f"--tmax_ms {tmax_ms} --Ts {Ts} --kfolds {kfolds} --cache yes --svm_window_size 40 "
                f"--x_size {x_size} --y_size {y_size} "
                f"--retrainPath {retrainPath} "
                f"--retrainPath_id {runID}-{datasetID}-mr_cnt-{mr_cnt}-spatialres-{spatialres}-best-{fold}"
            )

#-----------------------------

