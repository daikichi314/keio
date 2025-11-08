#!/bin/bash
#
# id: run_ct_plotter.sh
# Place: ~/hkelec/DiscreteSoftware/Analysis/macro/
# Last Edit: 2025-11-04
#
# 概要: Charge vs Time プロット（エラーバー付き）を
#       チャンネルごとに一括作成する
# 実行可能
#

# --- 1. 設定 ---
BASE_DIR="$(dirname "$0")"
PLOTTER="${BASE_DIR}/fit_results/plot_ct"

# --- 2. 引数の検証 ---
if [ "$#" -ne 2 ]; then
    echo "使い方: $0 <対象のデータディレクトリ> <データタイプ>"
    echo "データタイプ: "
    echo "  'gaus'  -> gausfit系列 (Charge_vs_Time_gaus_ch*.txt)"
    echo "  'peak'  -> peakfinder系列 (Charge_vs_Time_peak_ch*.txt)"
    echo "  'mean'  -> meanfinder系列 (Charge_vs_Time_mean_ch*.txt)"
    exit 1
fi

TARGET_DIR=$1
DATA_TYPE=$2

# --- 3. データタイプに応じてファイル名パターンを決定 ---
FILE_PATTERN=""
if [ "$DATA_TYPE" == "gaus" ]; then
    FILE_PATTERN="Charge_vs_Time_gaus_ch*.txt"
elif [ "$DATA_TYPE" == "peak" ]; then
    FILE_PATTERN="Charge_vs_Time_peak_ch*.txt"
elif [ "$DATA_TYPE" == "mean" ]; then
    FILE_PATTERN="Charge_vs_Time_mean_ch*.txt"
else
    echo "エラー: 不正なデータタイプです。"
    exit 1
fi

# --- 4. 全チャンネルのプロットを作成 ---
for file in "$TARGET_DIR"/$FILE_PATTERN; do
    if [ -f "$file" ]; then
        echo "処理中: $file"
        "$PLOTTER" "$file" "$TARGET_DIR"
    fi
done

echo ""
echo "✅ 全てのプロットが作成されました。"