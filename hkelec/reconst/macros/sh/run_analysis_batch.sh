#!/bin/bash
#
# id: run_analysis_batch.sh
# Place: /home/daiki/keio/hkelec/reconst/macros/sh/
# Created: 2025-11-21
#
# 概要: 指定されたディレクトリ群に対して、以下の解析フローを一括実行する
#   1. fit_pedestal (ペデスタル解析)
#   2. meanfinder (イベント電荷・時間解析)
#   3. plot_summary (集計・グラフ作成・フィッティング)
#

# --- 設定: 実行ファイルの場所 ---
# このスクリプトからの相対パス、または絶対パスで指定
BASE_DIR="$(cd "$(dirname "$0")" && pwd)"
CPP_DIR="${BASE_DIR}/../cpp"

EXE_PEDESTAL="${CPP_DIR}/fit_pedestal"
EXE_MEANFINDER="${CPP_DIR}/meanfinder"
EXE_PLOT="${CPP_DIR}/plot_summary"

# --- ヘルプ表示 ---
usage() {
    cat <<EOF
===============================================================================
  Batch Analysis Runner
===============================================================================
[概要]
  指定されたディレクトリごとに、Pedestal Fit -> Charge/Time Fit -> Summary Plot
  の一連の解析を実行します。

[使い方]
  $0 <target_dirs...> [options]

[例]
  $0 /data/no1 /data/no2 --fit-all

[オプション]
  --fit-charge : meanfinderで電荷のみ計算
  --fit-time   : meanfinderで時間のみ計算
  --fit-all    : 両方計算 (デフォルト)
  --no-pdf     : PDFを出力しない

[入出力ファイルの仕様]
  -----------------------------------------------------------------------------
  | 入力 | hkelec_pedestal_hithist.root | 必須 | ペデスタルデータ           |
  | 入力 | *_eventhist.root             | 必須 | イベントデータ             |
  -----------------------------------------------------------------------------
  | 出力 | _fits.txt, _mean.txt, ...    | 自動 | 各ステップの解析結果       |
  | 出力 | Charge_vs_Time_chXX.pdf      | 自動 | 最終的なフィット結果グラフ |
  -----------------------------------------------------------------------------
===============================================================================
EOF
}

# 引数がない場合はヘルプを表示
if [ "$#" -lt 1 ]; then
    usage
    exit 1
fi

# --- オプション解析 ---
DIRS=()
FIT_OPTION="--fit-all"
PDF_OPTION=""

while (( "$#" )); do
    case "$1" in
        --help|-h)
            usage
            exit 0
            ;;
        --fit-charge|--fit-time|--fit-all)
            FIT_OPTION="$1"
            shift
            ;;
        --no-pdf)
            PDF_OPTION="$1"
            shift
            ;;
        -*)
            echo "エラー: 不明なオプション $1"
            usage
            exit 1
            ;;
        *)
            # オプションでないものはディレクトリとしてリストに追加
            DIRS+=("$1")
            shift
            ;;
    esac
done

if [ ${#DIRS[@]} -eq 0 ]; then
    echo "エラー: 対象ディレクトリが指定されていません。"
    exit 1
fi

# --- 実行ファイルの存在確認 ---
for cmd in "$EXE_PEDESTAL" "$EXE_MEANFINDER" "$EXE_PLOT"; do
    if [ ! -x "$cmd" ]; then
        echo "エラー: 実行ファイルが見つかりません: $cmd"
        echo "cppディレクトリで 'make' を実行してください。"
        exit 1
    fi
done

# --- メインループ ---
for DIR in "${DIRS[@]}"; do
    echo "----------------------------------------------------------------"
    echo "Processing Directory: $DIR"
    echo "----------------------------------------------------------------"

    if [ ! -d "$DIR" ]; then
        echo "警告: ディレクトリが存在しません。スキップします: $DIR"
        continue
    fi

    # 1. Pedestal Fit
    PED_FILE="$DIR/hkelec_pedestal_hithist.root"
    if [ -f "$PED_FILE" ]; then
        echo "[Step 1] Running Pedestal Fit..."
        "$EXE_PEDESTAL" "$PED_FILE" "$PDF_OPTION"
    else
        echo "[Step 1] Skip: ペデスタルファイルが見つかりません ($PED_FILE)"
        # ペデスタルがないと後続のCharge計算が正しく行えない可能性があるが、
        # meanfinder側で処理(欠損時の挙動)に任せて続行する
    fi

    # 2. Meanfinder (Event Fit)
    echo "[Step 2] Running Meanfinder ($FIT_OPTION)..."
    # eventhistファイルを検索してループ
    count=0
    find "$DIR" -maxdepth 1 -name "*_eventhist.root" | sort | while read evf; do
        echo "  - $evf"
        "$EXE_MEANFINDER" "$evf" "$FIT_OPTION" "$PDF_OPTION"
        ((count++))
    done
    
    # パイプ内でのカウントは反映されないため、再度確認
    num_files=$(find "$DIR" -maxdepth 1 -name "*_eventhist.root" | wc -l)
    if [ "$num_files" -eq 0 ]; then
        echo "  警告: eventhistファイルが見つかりません。スキップします。"
        continue
    fi

    # 3. Plot Summary
    echo "[Step 3] Creating Summary Plots..."
    "$EXE_PLOT" "$DIR" "$PDF_OPTION"

    echo "Done: $DIR"
    echo ""
done

echo "全工程が完了しました。"
exit 0