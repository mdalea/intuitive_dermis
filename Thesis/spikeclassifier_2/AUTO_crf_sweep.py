# -----------------------------
# New dataset sweep dimensions
# -----------------------------
crf_array = [4, 8, 16, 64]

# (crf_scope, dist_tag) pairs matching your folder names exactly
dataset_variants = [
    ("global_crf", "skewdist"),
    ("global_crf", "unidist-max8"),
    ("local_crf",  "skewdist"),
    ("local_crf",  "unidist-max8"),
]

runID = "sensor_sweep_crf"

# base part (everything up to "...-3"), same as before
datasetID_base = "lcadc-0-uV2-0-uV2-0-uV2-0-mult-0-uV2-3"

# NOTE: your new datasets use --1to10-1to10 (not --1to40-1to40)
SCAN_TAG = "multitrial-Kylberg_filt_6_scan_oversampled20x_actualdimscale--1to10-1to10"

# ... keep your addressing code above unchanged ...

for mr_cnt_val in mr_cnt_array:
    mr_cnt_int = int(mr_cnt_val)
    mr_cnt = str(mr_cnt_int)

    X_MAX = 15
    x_size = X_MAX + 1
    y_size = math.ceil(mr_cnt_int / x_size)

    if x_size * y_size < mr_cnt_int:
        raise RuntimeError(
            f"Address capacity too small: x_size*y_size={x_size*y_size} < mr_cnt={mr_cnt_int}."
        )

    snn_run_name = (
        f"SNN_ep-{epochs}-isize-2-{x_size}-{y_size}"
        f"-hsize-{hidden_size}-bsize-{batch_size}-lr-{lr}"
        f"-tmax_ms-{fmt_float_1(tmax_ms)}-Ts-{fmt_float_1(Ts)}-k-{kfolds}-"
    )

    # -----------------------------
    # NEW: sweep CRF + dataset variants
    # -----------------------------
    for crf_x in crf_array:
        for crf_scope, dist_tag in dataset_variants:

            # this becomes: lcadc-...-3-crf-X-global_crf-skewdist  (etc.)
            datasetID = f"{datasetID_base}-crf-{crf_x}-{crf_scope}-{dist_tag}"

            for spatialres_val in spatialres_array:
                spatialres = fmt_spatialres(spatialres_val)

                # EXACT folder name you provided (plus variable mr_cnt/spatialres)
                datasetFolderName = (
                    f"ALL_N_{datasetID}-{SCAN_TAG}"
                    f"-mr_cnt-{mr_cnt}-spatialres-{spatialres}_sorted"
                )

                datasetPath = f"../spikegen/global_outfile/spikegen/{datasetFolderName}/"

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
                    # IMPORTANT: match how your SNN saved folder encodes datasetPath in the checkpoint path:
                    # "dsetPath-" + <datasetFolderName> + "/chip-classifier-test-...-best.pt"
                    retrainPath = (
                        f"{snn_run_name}"
                        f"dsetPath-{datasetFolderName}/"
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

