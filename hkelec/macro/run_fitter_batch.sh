#!/bin/bash
#
# id: run_fitter_batch.sh
# Place: ~/hkelec/DiscreteSoftware/Analysis/macro/
# Created: 2025-11-05
#
# 概要: gausfit / meanfinder / peakfinder のワークフローを統一して実行する
#       引数で method を指定 (gaus|mean|peak|all)。実行ログを出力する。
#

BASE_DIR="$(dirname "$0")"

# executables
GAUS_FITTER="${BASE_DIR}/fit_results/gausfit"
MEAN_FITTER="${BASE_DIR}/fit_results/meanfinder"
PEAK_FITTER="${BASE_DIR}/fit_results/peakfinder"
GAIN_SELECTOR="${BASE_DIR}/fit_results/select_gain"
GAIN_SELECTOR_MEAN="${BASE_DIR}/fit_results/select_gain_mean"
CT_PLOT_CREATOR="${BASE_DIR}/fit_results/create_ct_plot"

usage(){
    cat <<EOF
Usage: $0 <target_dir> <method> [pedestal_dir] [options]
  method: gaus | mean | peak | all
  options:
    --fit-charge | --fit-time | --fit-all
    --no-pdf
    --make-ct-plot
EOF
}

if [ "$#" -lt 2 ]; then
    usage; exit 1
fi

# TARGET_DIR=$1
# METHOD=$2
# PEDESTAL_DIR=$3 ##########$TARGET_DIR
# shift 2
# if [ "$#" -ge 1 ] && ! [[ "$1" =~ ^-- ]]; then
#     PEDESTAL_DIR=$3; shift #### $1
# fi

# FIT_OPTION="--fit-charge"

# 1. TARGET_DIR と METHOD を取得
TARGET_DIR=$1
METHOD=$2

# 2. PEDESTAL_DIR のデフォルト値を TARGET_DIR に設定
PEDESTAL_DIR="$TARGET_DIR" 
shift 2 # $1 と $2 を消費

# 3. 残りの第1引数($1)が存在し、かつオプション(--で始まらない)場合、
#    それを PEDESTAL_DIR として採用し、引数リストから消費($1)する
if [ "$#" -ge 1 ] && ! [[ "$1" =~ ^-- ]]; then
    PEDESTAL_DIR=$1; shift 
fi

# 4. FIT_OPTION 以降のオプション解析は変更なし
FIT_OPTION="--fit-charge"
PDF_OPTION=""
MAKE_CT_PLOT="no"
for arg in "$@"; do
    case $arg in
        --fit-charge|--fit-time|--fit-all) FIT_OPTION=$arg ;; 
        --no-pdf) PDF_OPTION=$arg ;;
        --make-ct-plot) MAKE_CT_PLOT="yes" ;;
        *) echo "Unknown option: $arg"; usage; exit 1 ;;
    esac
done

# prepare logfile
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOGFILE="$TARGET_DIR/run_fitter_${METHOD}_${TIMESTAMP}.log"
mkdir -p "$(dirname "$LOGFILE")"
exec > >(tee -a "$LOGFILE") 2>&1

echo "=== run_fitter_batch start: $(date) ==="
echo "target: $TARGET_DIR  method: $METHOD  fit: $FIT_OPTION  make_ct: $MAKE_CT_PLOT"

# helper: run per-file command if method matches
run_for_method(){
    local m=$1; shift
    case "$METHOD" in
        all) $@ ;;
        $m) $@ ;;
        *) return 0 ;;
    esac
}

echo "\n--- Pedestal step ---"
PEDESTAL_ROOT_FILE="hkelec_pedestal_hithist.root"
PEDESTAL_ROOT_PATH="$PEDESTAL_DIR/$PEDESTAL_ROOT_FILE"
PEDESTAL_TXT_PATH="$PEDESTAL_DIR/${PEDESTAL_ROOT_FILE/.root/_fits.txt}"
if [ -f "$PEDESTAL_ROOT_PATH" ]; then
    if [ -x "${BASE_DIR}/fit_results/fit_pedestal" ]; then
        echo "running pedestal fitter..."
        "${BASE_DIR}/fit_results/fit_pedestal" "$PEDESTAL_ROOT_PATH" "$PDF_OPTION"
    else
        echo "pedestal fitter not found, skipping"
        touch "$PEDESTAL_TXT_PATH"
    fi
else
    echo "pedestal root not found: $PEDESTAL_ROOT_PATH"
    touch "$PEDESTAL_TXT_PATH"
fi

# Process eventhist files
echo "\n--- Processing eventhist files ---"
for evf in "$TARGET_DIR"/*_eventhist.root; do
    [ -f "$evf" ] || continue
    echo "Processing: $evf"

    # run fit-charge / fit-time with the correct executable depending on METHOD
    if [[ "$FIT_OPTION" == "--fit-charge" || "$FIT_OPTION" == "--fit-all" ]]; then
        case "$METHOD" in
            gaus)
                "$GAUS_FITTER" "$evf" "--fit-charge" "$PDF_OPTION" ;;
            mean)
                "$MEAN_FITTER" "$evf" "--fit-charge" "$PDF_OPTION" ;;
            peak)
                "$PEAK_FITTER" "$evf" "--fit-charge" "$PDF_OPTION" ;;
            all)
                "$GAUS_FITTER" "$evf" "--fit-charge" "$PDF_OPTION"
                "$MEAN_FITTER" "$evf" "--fit-charge" "$PDF_OPTION"
                "$PEAK_FITTER" "$evf" "--fit-charge" "$PDF_OPTION" ;;
        esac
    fi

    if [[ "$FIT_OPTION" == "--fit-time" || "$FIT_OPTION" == "--fit-all" ]]; then
        case "$METHOD" in
            gaus)
                "$GAUS_FITTER" "$evf" "--fit-time" "$PDF_OPTION" ;;
            mean)
                "$MEAN_FITTER" "$evf" "--fit-time" "$PDF_OPTION" ;;
            peak)
                "$PEAK_FITTER" "$evf" "--fit-time" "$PDF_OPTION" ;;
            all)
                "$GAUS_FITTER" "$evf" "--fit-time" "$PDF_OPTION"
                "$MEAN_FITTER" "$evf" "--fit-time" "$PDF_OPTION"
                "$PEAK_FITTER" "$evf" "--fit-time" "$PDF_OPTION" ;;
        esac
    fi
done

echo "\n--- Creating summaries and HV vs Charge files ---"
if [[ "$METHOD" == "gaus" || "$METHOD" == "all" ]]; then
    SUMMARY_FILE_CHARGE="$TARGET_DIR/summary_gausfit_all.txt"
    rm -f "$SUMMARY_FILE_CHARGE"
    echo "# ch,type,voltage,peak,peak_err,sigma,sigma_err,chi2_ndf,rough_sigma" > "$SUMMARY_FILE_CHARGE"
    cat "$TARGET_DIR"/*_gausfit.txt | grep -v '^#' >> "$SUMMARY_FILE_CHARGE"
    echo "Running gain selector (gaus)"
    if [ -x "$GAIN_SELECTOR" ]; then
        "$GAIN_SELECTOR" "$SUMMARY_FILE_CHARGE" "$PEDESTAL_TXT_PATH" "$TARGET_DIR" "gaus"
    else
        echo "GAIN_SELECTOR not found: $GAIN_SELECTOR"
    fi
fi

if [[ "$METHOD" == "mean" || "$METHOD" == "all" ]]; then
    SUMMARY_FILE_CHARGE_MEAN="$TARGET_DIR/summary_HV_vs_Charge_mean.txt"
    # meanfinder outputs *_mean.txt; select_gain_mean produces selected files and summary
    # create aggregated mean summary so select_gain_mean can read a single CSV-like file
    rm -f "$SUMMARY_FILE_CHARGE_MEAN"
    echo "# ch,type,voltage,mean,mean_err,rms" > "$SUMMARY_FILE_CHARGE_MEAN"
    # Some *_mean.txt files may contain wrapped lines; pick only lines that start with channel number
    grep -h '^[0-9]' "$TARGET_DIR"/*_mean.txt 2>/dev/null >> "$SUMMARY_FILE_CHARGE_MEAN" || true
    if [ -x "$GAIN_SELECTOR_MEAN" ]; then
        echo "Running gain selector (mean)"
        # pass: <summary_mean_all.txt> <pedestal_fits.txt> <output_dir>
        "$GAIN_SELECTOR_MEAN" "$SUMMARY_FILE_CHARGE_MEAN" "$PEDESTAL_TXT_PATH" "$TARGET_DIR"
    else
        echo "GAIN_SELECTOR_MEAN not found: $GAIN_SELECTOR_MEAN"
    fi
fi

if [[ "$METHOD" == "peak" || "$METHOD" == "all" ]]; then
    # peak workflow may create HV_vs_ChargeSubtracted files; keep existing behaviour
    echo "Note: peak workflow HV file creation should be handled by existing batch scripts before calling this unified script if needed."
fi

echo "\n--- Time summary aggregation ---"
SUMMARY_FILE_TIME="$TARGET_DIR/summary_timefit_all.txt"
rm -f "$SUMMARY_FILE_TIME"
echo "# ch,type,voltage,tts,sigma,fwhm,peak,tau,chi2_ndf" > "$SUMMARY_FILE_TIME"
cat "$TARGET_DIR"/*_timefit.txt 2>/dev/null | grep -v '^#' >> "$SUMMARY_FILE_TIME" || true

# Create Charge vs Time files if requested
if [ "$MAKE_CT_PLOT" == "yes" ]; then
    echo "\n--- Creating Charge vs Time data and plots ---"
    # call CT plot creator if summaries exist
    if [ -f "$SUMMARY_FILE_TIME" ]; then
        # For CT plot creation, we need:
        # 1. charge summary (from gausfit/mean/peak)
        # 2. time summary (already in SUMMARY_FILE_TIME)
        # 3. pedestal file
        # 4. output directory
        CHARGE_SUMMARY=""
        case "$METHOD" in
            gaus|all)
                if [ -f "$SUMMARY_FILE_CHARGE" ]; then
                    CHARGE_SUMMARY="$SUMMARY_FILE_CHARGE"
                fi ;;
            mean)
                if [ -f "$SUMMARY_FILE_CHARGE_MEAN" ]; then
                    CHARGE_SUMMARY="$SUMMARY_FILE_CHARGE_MEAN"
                fi ;;
            peak)
                # Try to find first *_peak.txt in target dir
                PEAK_SUMMARY=$(find "$TARGET_DIR" -name "*_peak.txt" | head -n 1)
                if [ -n "$PEAK_SUMMARY" ]; then
                    CHARGE_SUMMARY="$PEAK_SUMMARY"
                fi ;;
        esac
        
        if [ -n "$CHARGE_SUMMARY" ] && [ -x "$CT_PLOT_CREATOR" ]; then
            # decide which method the charge summary corresponds to so CT filenames reflect it
            CT_METHOD="$METHOD"
            case "$CHARGE_SUMMARY" in
                *gaus* ) CT_METHOD="gaus" ;;
                *mean* ) CT_METHOD="mean" ;;
                *peak* ) CT_METHOD="peak" ;;
                * ) CT_METHOD="$METHOD" ;;
            esac
            # pass method_for_ct as last arg
            "$CT_PLOT_CREATOR" "$CHARGE_SUMMARY" "$SUMMARY_FILE_TIME" "$PEDESTAL_TXT_PATH" "$TARGET_DIR" "$CT_METHOD"
        else
            echo "CT_PLOT_CREATOR not found or charge summary missing"
        fi
    else
        echo "summary_timefit_all.txt not found; cannot create CT plots"
    fi
fi

echo "\n=== run_fitter_batch end: $(date) ==="
echo "Log saved to: $LOGFILE"

exit 0
