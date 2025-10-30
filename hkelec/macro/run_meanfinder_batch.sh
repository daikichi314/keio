#!/bin/bash
#
# id: run_meanfinder_batch.sh
# Place: ~/hkelec/DiscreteSoftware/Analysis/macro/
# Last Edit: 2025-10-29 Gemini
#
# 概要: meanfinder (電荷平均値 + 時間EMGフィット) の解析フロー全体を自動実行する。
#       Charge vs Time グラフ用データの作成機能も含む。
#

# --- 1. 設定 (★ gausfit から変更) ---
BASE_DIR="$(dirname "$0")"
PEDESTAL_FITTER="${BASE_DIR}/fit_results/fit_pedestal"
# 1. GAUS_FITTER -> MEAN_FINDER
MEAN_FINDER="${BASE_DIR}/fit_results/meanfinder"
# 2. GAIN_SELECTOR -> GAIN_SELECTOR_MEAN
GAIN_SELECTOR_MEAN="${BASE_DIR}/fit_results/select_gain_mean"
CT_PLOT_CREATOR="${BASE_DIR}/fit_results/create_ct_plot"
PEDESTAL_ROOT_FILE="hkelec_pedestal_hithist.root"

# --- 2. 引数の解釈 (変更なし) ---
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
# 3. デフォルトを --fit-charge に変更
FIT_OPTION="--fit-charge" 
MAKE_CT_PLOT="no"
for arg in "$@"; do
    case $arg in
        --fit-charge|--fit-time|--fit-all) FIT_OPTION=$arg ;;
        --no-pdf) PDF_OPTION=$arg ;;
        --make-ct-plot) MAKE_CT_PLOT="yes" ;;
    esac
done

echo "--- 電荷平均値・時間フィット バッチ処理を開始します ---"
echo "対象ディレクトリ: $TARGET_DIR"
echo "フィットオプション: $FIT_OPTION"

# --- 3. ステップ1: ペデスタルフィット (変更なし) ---
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

# --- 4. ステップA: 電荷平均値の計算 (★ gausfit から変更) ---
if [[ "$FIT_OPTION" == "--fit-charge" || "$FIT_OPTION" == "--fit-all" ]]; then
    echo ""
    echo "--- ステップA-1: 各信号ROOTファイルの電荷平均値を計算しています... ---"
    for file in "$TARGET_DIR"/*_eventhist.root; do
        if [ -f "$file" ]; then
            # 4. MEAN_FINDER を --fit-charge で呼び出す
            # (電荷計算はPDFオプション不要)
            "$MEAN_FINDER" "$file" "--fit-charge"
        fi
    done
    echo ""
    echo "--- ステップA-2: 電荷平均値のサマリーを作成しています... ---"
    # 5. 出力ファイル名を summary_mean_all.txt に変更
    SUMMARY_FILE_CHARGE="$TARGET_DIR/summary_mean_all.txt"
    rm -f "$SUMMARY_FILE_CHARGE"
    # 6. ヘッダーを mean 用に変更
    echo "# ch,type,voltage,mean,mean_err,rms" > "$SUMMARY_FILE_CHARGE"
    # 7. *_mean.txt を収集
    cat "$TARGET_DIR"/*_mean.txt | grep -v '^#' >> "$SUMMARY_FILE_CHARGE"
    
    echo ""
    echo "--- ステップA-3: HV vs Charge (Mean) グラフファイルを作成しています... ---"
    # 8. select_gain_mean を呼び出す
    "$GAIN_SELECTOR_MEAN" "$SUMMARY_FILE_CHARGE" "$PEDESTAL_TXT_PATH" "$TARGET_DIR"
fi

# --- 5. ステップB: 時間フィット (★ gausfit から変更) ---
if [[ "$FIT_OPTION" == "--fit-time" || "$FIT_OPTION" == "--fit-all" ]]; then
    echo ""
    echo "--- ステップB-1: 各信号ROOTファイルに時間フィットを実行しています... ---"
    for file in "$TARGET_DIR"/*_eventhist.root; do
        if [ -f "$file" ]; then
            # 9. MEAN_FINDER を --fit-time で呼び出す
            "$MEAN_FINDER" "$file" "--fit-time" "$PDF_OPTION"
        fi
    done
    echo ""
    echo "--- ステップB-2: 時間フィットのサマリーを作成しています... ---"
    SUMMARY_FILE_TIME="$TARGET_DIR/summary_timefit_all.txt"
    rm -f "$SUMMARY_FILE_TIME"
    # 10. ヘッダーを (gausfit/ttshistofit と) 同一に
    echo "# ch,type,voltage,tts(sigma),sigma,fwhm(calc),peak(calc),tau(1/lambda),chi2_ndf" > "$SUMMARY_FILE_TIME"
    cat "$TARGET_DIR"/*_timefit.txt | grep -v '^#' >> "$SUMMARY_FILE_TIME"
fi


# --- 6. 最終ステップ: Charge vs Time グラフ用データの作成 (★ gausfit から変更) ---
if [ "$MAKE_CT_PLOT" == "yes" ]; then
    echo ""
    echo "--- 最終ステップ: Charge vs Time グラフ用データを作成しています... ---"
    
    # 11. 電荷サマリーファイル名を _mean_all.txt に変更
    SUMMARY_FILE_CHARGE_FOR_CT="$TARGET_DIR/summary_HV_vs_Charge_mean.txt"
    
    # 12. 時間サマリーファイル名 (timefit_all) は変更なし
    SUMMARY_FILE_TIME_FOR_CT="$TARGET_DIR/summary_timefit_all.txt" 
    
    # 13. C-Tグラフ作成 (create_ct_plot を呼び出す)
    # (入力ファイル名が _mean になっている点に注意)
    if [ -f "$SUMMARY_FILE_CHARGE_FOR_CT" ] && [ -f "$SUMMARY_FILE_TIME_FOR_CT" ]; then
        "$CT_PLOT_CREATOR" "$SUMMARY_FILE_CHARGE_FOR_CT" "$SUMMARY_FILE_TIME_FOR_CT" "$PEDESTAL_TXT_PATH" "$TARGET_DIR"
    else
        echo "警告: Charge vs Time グラフ作成に必要なサマリーファイルが見つかりません。"
        echo "($SUMMARY_FILE_CHARGE_FOR_CT または $SUMMARY_FILE_TIME_FOR_CT)"
        echo "(--fit-all オプションでスクリプトを再実行してください)"
    fi
fi


echo ""
echo "✅ 全ての処理が完了しました。"