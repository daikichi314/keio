#!/bin/bash
#
# id: run_hv_fitter.sh
# Place: ~/hkelec/DiscreteSoftware/Analysis/macro/
# Last Edit: 2025-10-16 Gemini
#
# 概要: HV vs Charge のグラフ用テキストファイルを一括でフィットし、
#       結果をPDFとサマリーテキストファイルに出力する。
# 実行可能
#

# 1. 設定
BASE_DIR="$(dirname "$0")"
FITTER_EXECUTABLE="${BASE_DIR}/fit_results/fit_hv_gain"

# 2. 引数の検証
if [ "$#" -ne 2 ]; then
    echo "使い方: $0 <対象のデータディレクトリ> <データタイプ>"
    echo "データタイプ: "
    echo "  'gaus' -> gausfit系列 (HV_vs_Charge_chXX.txt)"
    echo "  'peak' -> peakfinder系列 (HV_vs_ChargeSubtracted_hgain_chXX.txt)"
    echo "  'mean' -> meanfinder系列 (HV_vs_Charge_mean_chXX.txt)"
    exit 1
fi

TARGET_DIR=$1
DATA_TYPE=$2

# 3. データタイプに応じてファイル名パターンを決定
FILE_PATTERN=""
SUMMARY_FILE=""
if [ "$DATA_TYPE" == "gaus" ]; then
    FILE_PATTERN="HV_vs_Charge_ch*.txt"
    SUMMARY_FILE="$TARGET_DIR/summary_HV_vs_Charge_gaus_fit_results.txt"
elif [ "$DATA_TYPE" == "peak" ]; then
    FILE_PATTERN="HV_vs_Charge_peak_hgain_ch*.txt"
    SUMMARY_FILE="$TARGET_DIR/summary_HV_vs_Charge_peak_fit_results.txt"
elif [ "$DATA_TYPE" == "mean" ]; then
    # select_gain_mean.C の出力ファイルを使用
    FILE_PATTERN="HV_vs_ChargeSelected_mean_ch*.txt"
    SUMMARY_FILE="$TARGET_DIR/summary_HV_vs_Charge_mean_fit_results.txt"
else
    echo "エラー: 不正なデータタイプです。"
    exit 1
fi

# 4. 最終的なサマリーファイルの準備
rm -f "$SUMMARY_FILE"
echo "# ch,param_b,param_b_err,param_a,param_a_err,HV for Gain = 10e7 (if 1 p.e.)" > "$SUMMARY_FILE"

# 5. 対象となるテキストファイルをループ処理
for file in "$TARGET_DIR"/$FILE_PATTERN; do
    if [ -f "$file" ]; then
        echo "処理中: $file"
        # C++プログラムを実行し、結果をサマリーファイルに追記
        "$FITTER_EXECUTABLE" "$file" >> "$SUMMARY_FILE"
    fi
done

echo ""
echo "✅ フィット処理が完了しました。"

