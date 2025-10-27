#!/bin/bash
#
# id: run_gausfit_batch.sh
# Place: ~/hkelec/DiscreteSoftware/Analysis/macro/
# Last Edit: 2025-10-18 Gemini
#
# 概要: gausfit系列の解析フロー全体を自動で実行するスクリプト。
#       オプションでCharge vs Timeグラフ用データの作成も可能。
# 実行可能
#

# --- 1. 設定 ---
BASE_DIR="$(dirname "$0")"
PEDESTAL_FITTER="${BASE_DIR}/fit_results/fit_pedestal"
GAUS_FITTER="${BASE_DIR}/fit_results/gausfit"
GAIN_SELECTOR="${BASE_DIR}/fit_results/select_gain"
CT_PLOT_CREATOR="${BASE_DIR}/fit_results/create_ct_plot" # ★★★ 新規追加 ★★★
PEDESTAL_ROOT_FILE="hkelec_pedestal_hithist.root"

# --- 2. 引数の解釈 ---
if [ "$#" -lt 1 ]; then
    echo "使い方: $0 <対象のデータディレクトリ> [ペデスタルディレクトリ] [オプション]"
    echo "オプション: --fit-charge | --fit-time | --fit-all | --no-pdf | --make-ct-plot"
    exit 1
fi
TARGET_DIR=$1
PEDESTAL_DIR=$TARGET_DIR
if [ "$#" -ge 2 ] && ! [[ "$2" =~ ^-- ]]; then
    PEDESTAL_DIR=$2
fi

PDF_OPTION=""
FIT_OPTION="--fit-charge"
MAKE_CT_PLOT="no" # ★★★ 新規追加 ★★★
for arg in "$@"; do
    case $arg in
        --fit-charge|--fit-time|--fit-all) FIT_OPTION=$arg ;;
        --no-pdf) PDF_OPTION=$arg ;;
        --make-ct-plot) MAKE_CT_PLOT="yes" ;; # ★★★ 新規追加 ★★★
    esac
done

# ... (ステップ1-3は変更なし) ...
echo "--- フル解析バッチ処理を開始します ---"
echo "対象ディレクトリ: $TARGET_DIR"
echo "フィットオプション: $FIT_OPTION"

echo ""
echo "--- ステップ1: ペデスタルファイルをフィットしています... ---"
PEDESTAL_ROOT_PATH="$PEDESTAL_DIR/$PEDESTAL_ROOT_FILE"
PEDESTAL_TXT_PATH="$PEDESTAL_DIR/${PEDESTAL_ROOT_FILE/.root/_fits.txt}"
if [ -f "$PEDESTAL_ROOT_PATH" ]; then
    "$PEDESTAL_FITTER" "$PEDESTAL_ROOT_PATH" "$PDF_OPTION"
else
    echo "警告: ペデスタルファイル $PEDESTAL_ROOT_PATH が見つかりません。"
    touch "$PEDESTAL_TXT_PATH"
fi

# ... 信号フィット、サマリー作成 ...
if [[ "$FIT_OPTION" == "--fit-charge" || "$FIT_OPTION" == "--fit-all" ]]; then
    echo ""
    echo "--- ステップA-1: 各信号ROOTファイルに電荷フィットを実行しています... ---"
    for file in "$TARGET_DIR"/*_eventhist.root; do
        if [ -f "$file" ]; then
            "$GAUS_FITTER" "$file" "--fit-charge" "$PDF_OPTION"
        fi
    done
    echo ""
    echo "--- ステップA-2: 電荷フィットのサマリーを作成しています... ---"
    SUMMARY_FILE_CHARGE="$TARGET_DIR/summary_gausfit_all.txt"
    rm -f "$SUMMARY_FILE_CHARGE"
    echo "# ch,type,voltage,peak,peak_err,sigma,sigma_err,chi2_ndf,rough_sigma" > "$SUMMARY_FILE_CHARGE"
    cat "$TARGET_DIR"/*_gausfit.txt | grep -v '^#' >> "$SUMMARY_FILE_CHARGE"
    
    echo ""
    echo "--- ステップA-3: 最適なADCを選択し、HV vs Chargeグラフファイルを作成しています... ---"
    "$GAIN_SELECTOR" "$SUMMARY_FILE_CHARGE" "$PEDESTAL_TXT_PATH" "$TARGET_DIR"
fi

if [[ "$FIT_OPTION" == "--fit-time" || "$FIT_OPTION" == "--fit-all" ]]; then
    echo ""
    echo "--- ステップB-1: 各信号ROOTファイルに時間フィットを実行しています... ---"
    for file in "$TARGET_DIR"/*_eventhist.root; do
        if [ -f "$file" ]; then
            "$GAUS_FITTER" "$file" "--fit-time" "$PDF_OPTION"
        fi
    done
    echo ""
    echo "--- ステップB-2: 時間フィットのサマリーを作成しています... ---"
    SUMMARY_FILE_TIME="$TARGET_DIR/summary_timefit_all.txt"
    rm -f "$SUMMARY_FILE_TIME"
    echo "# ch,type,voltage,tts,sigma,fwhm,peak,tau,chi2_ndf" > "$SUMMARY_FILE_TIME"
    cat "$TARGET_DIR"/*_timefit.txt | grep -v '^#' >> "$SUMMARY_FILE_TIME"
fi


# --- ★★★ 新規ステップ: Charge vs Time グラフ用データの作成 ★★★ ---
if [ "$MAKE_CT_PLOT" == "yes" ]; then
    echo ""
    echo "--- 最終ステップ: Charge vs Time グラフ用データを作成しています... ---"
    # 電荷と時間の両方のサマリーファイルが存在する場合のみ実行
    if [ -f "$SUMMARY_FILE_CHARGE" ] && [ -f "$SUMMARY_FILE_TIME" ]; then
        "$CT_PLOT_CREATOR" "$SUMMARY_FILE_CHARGE" "$SUMMARY_FILE_TIME" "$PEDESTAL_TXT_PATH" "$TARGET_DIR"
    else
        echo "警告: Charge vs Time グラフ作成に必要なサマリーファイルが見つかりません。"
        echo "(--fit-all オプションでスクリプトを再実行してください)"
    fi
fi


echo ""
echo "✅ 全ての処理が完了しました。"
