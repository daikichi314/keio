#!/bin/bash

# id: run_reconstruction_batch.sh
# Place: /home/daiki/keio/hkelec/reconst/reco/
# Author: Gemini 3 Pro
# Last Edit: 2025-12-06
#
# 概要:
# 指定ディレクトリ内の *eventhist.root ファイルに対して
# reconstructor を一括実行するバッチスクリプト

# --- スクリプトの使い方を表示する関数 ---
usage() {
    echo "=========================================================================="
    echo " [概要]"
    echo " 指定されたディレクトリ内のすべての *eventhist.root ファイルに対して"
    echo " 再構成プログラム (reconstructor) を実行します。"
    echo ""
    echo " [使い方]"
    echo " $0 <TargetDirectory> [Suffix]"
    echo ""
    echo " [引数]"
    echo " <TargetDirectory> : 処理対象のディレクトリパス"
    echo "                     この中にROOTファイルとペデスタルファイルが必要です。"
    echo " [Suffix]          : 出力ファイル名のサフィックス（省略可）"
    echo "                     指定時: *_eventhist.root → *_reconst_<Suffix>.root"
    echo "                     省略時: *_eventhist.root → *_reconst.root"
    echo ""
    echo " [前提条件]"
    echo " 1. ディレクトリ内に 'hkelec_pedestal_hithist_means.txt' が存在すること。"
    echo " 2. このスクリプトと同じディレクトリに実行ファイル 'reconstructor' が存在すること。"
    echo "    スクリプト位置から絶対パスを解決して実行します。"
    echo "=========================================================================="
}

# --- 引数チェック ---
if [ "$#" -lt 1 ]; then
    usage
    exit 1
fi

TARGET_DIR="$1"
SUFFIX="${2:-}"

# --- ディレクトリの存在チェック ---
if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: Directory '$TARGET_DIR' does not exist."
    exit 1
fi

# --- ペデスタルファイルの存在チェック ---
PEDESTAL_FILE="${TARGET_DIR}/hkelec_pedestal_hithist_means.txt"
if [ ! -f "$PEDESTAL_FILE" ]; then
    echo "Error: Pedestal file not found at: $PEDESTAL_FILE"
    exit 1
fi

# --- 実行ファイルの場所 (スクリプト位置から絶対パスを解決) ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RECONSTRUCTOR="${SCRIPT_DIR}/reconstructor"
if [ ! -f "$RECONSTRUCTOR" ]; then
    echo "Error: Executable not found: $RECONSTRUCTOR"
    echo "Please compile the code using 'make' first."
    exit 1
fi

echo "Start processing in: $TARGET_DIR"

# --- ループ処理 ---
# findコマンドで *eventhist.root を検索して処理
find "$TARGET_DIR" -name "*eventhist.root" | while read input_file; do
    
    # 出力ファイル名の生成 (_eventhist.root を除去して _reconst[_SUFFIX] を付与)
    # 例: run001_eventhist.root -> run001_reconst または run001_reconst_XXX
    base_name=$(basename "$input_file" _eventhist.root)
    if [ -n "$SUFFIX" ]; then
        output_base="${TARGET_DIR}/${base_name}_reconst_${SUFFIX}"
    else
        output_base="${TARGET_DIR}/${base_name}_reconst"
    fi
    
    echo "--------------------------------------------------------"
    echo "Processing: $input_file"
    echo "Output:     $output_base.root / .csv"
    
    # reconstructor の実行
    $RECONSTRUCTOR "$input_file" "$output_base"
    
    if [ $? -eq 0 ]; then
        echo "Done."
    else
        echo "Error occurred during processing $input_file"
    fi

done

echo "--------------------------------------------------------"
echo "All tasks finished."