#!/bin/bash
#
# id: run_peakfinder_batch.sh
# Place: ~/hkelec/DiscreteSoftware/Analysis/macro/
# Last Edit: 2025-10-17 Gemini
#
# 概要: peakfinder(電荷)とttsfit(時間)を一括実行するスクリプト。
# 実行可能
#

# --- 1. 設定 ---
BASE_DIR="$(dirname "$0")"
PEDESTAL_FITTER="${BASE_DIR}/fit_results/fit_pedestal"
PEAK_FINDER="${BASE_DIR}/fit_results/peakfinder" # peakfinder.Cは内部でttsfitも呼び出す
PEDESTAL_ROOT_FILE="hkelec_pedestal_hithist.root"
K_HGAIN=0.073 # pC/ADC for high gain

# --- 2. 引数の解釈 ---
if [ "$#" -lt 1 ]; then
    echo "使い方: $0 <対象のデータディレクトリ> [ペデスタルディレクトリ] [オプション]"
    echo "オプション: --fit-charge | --fit-time | --fit-all | --no-pdf"
    exit 1
fi
TARGET_DIR=$1
PEDESTAL_DIR=$TARGET_DIR
if [ "$#" -ge 2 ] && ! [[ "$2" =~ ^-- ]]; then
    PEDESTAL_DIR=$2
fi

FIT_OPTION="--fit-charge"
PDF_OPTION=""
for arg in "$@"; do
    case $arg in
        --fit-charge|--fit-time|--fit-all) FIT_OPTION=$arg ;;
        --no-pdf) PDF_OPTION=$arg ;;
    esac
done

echo "--- Peak Finder / TTS Fit バッチ処理を開始します ---"

# --- 3. (電荷処理時のみ) ペデスタルファイルをフィット ---
if [[ "$FIT_OPTION" == "--fit-charge" || "$FIT_OPTION" == "--fit-all" ]]; then
    PEDESTAL_ROOT_PATH="$PEDESTAL_DIR/$PEDESTAL_ROOT_FILE"
    PEDESTAL_TXT_PATH="$PEDESTAL_DIR/${PEDESTAL_ROOT_FILE/.root/_fits.txt}"
    if [ -f "$PEDESTAL_ROOT_PATH" ]; then
        "$PEDESTAL_FITTER" "$PEDESTAL_ROOT_PATH" --no-pdf
    else
        echo "警告: ペデスタルファイルが見つかりません。"
        touch "$PEDESTAL_TXT_PATH"
    fi
fi

# --- 4. 全ての eventhist.root ファイルに対して peakfinder を実行 ---
for file in "$TARGET_DIR"/*_eventhist.root; do
    if [ -f "$file" ]; then
        "$PEAK_FINDER" "$file" "$FIT_OPTION" "$PDF_OPTION"
    fi
done

# --- 5. (電荷処理時のみ) サマリー作成とペデスタル減算 ---
if [[ "$FIT_OPTION" == "--fit-charge" || "$FIT_OPTION" == "--fit-all" ]]; then
    SUMMARY_FILE="$TARGET_DIR/summary_peak_all.txt"
    rm -f "$SUMMARY_FILE"
    echo "# ch,type,voltage,peak_pos" > "$SUMMARY_FILE"
    cat "$TARGET_DIR"/*_peak.txt | grep -v '^#' >> "$SUMMARY_FILE"
    
    TYPE="hgain" # hgainのみ対象
    for ch in {0..11}; do
        GRAPH_FILE="$TARGET_DIR/HV_vs_ChargeSubtracted_${TYPE}_ch${ch}.txt"
        rm -f "$GRAPH_FILE"
        echo "# HV(V), Charge(pC)" > "$GRAPH_FILE"
        awk -F, -v ch="$ch" -v type="$TYPE" -v k_hgain="$K_HGAIN" '
            NR==FNR { if ($1==ch && $2==type) ped[$1]=$3; next }
            { if ($1==ch && $2==type) print $3, ($4 - ped[$1]) * k_hgain }
        ' "$PEDESTAL_TXT_PATH" "$SUMMARY_FILE" >> "$GRAPH_FILE"
        if [ -s "$GRAPH_FILE" ] && [ $(wc -l < "$GRAPH_FILE") -le 1 ]; then
            rm "$GRAPH_FILE"
        fi
    done
fi

# --- 6. (時間処理時のみ) サマリー作成 ---
if [[ "$FIT_OPTION" == "--fit-time" || "$FIT_OPTION" == "--fit-all" ]]; then
    SUMMARY_FILE_TIME="$TARGET_DIR/summary_timefit_all.txt"
    rm -f "$SUMMARY_FILE_TIME"
    echo "# ch,type,voltage,tts,sigma,fwhm,peak,tau,chi2_ndf" > "$SUMMARY_FILE_TIME"
    # peakfinderは時間フィットの結果を _timefit.txt に出力する
    cat "$TARGET_DIR"/*_timefit.txt | grep -v '^#' >> "$SUMMARY_FILE_TIME"
    echo "時間フィットのサマリーファイルを作成しました: $SUMMARY_FILE_TIME"
fi

echo ""
echo "✅ 全ての処理が完了しました。"
