#!/bin/bash

# --- スクリプトの使い方を表示する関数 ---
usage() {
    echo "=========================================================================="
    echo " [概要]"
    echo " 指定されたディレクトリ内のすべての *eventhist.root ファイルに対して"
    echo " 再構成プログラム (reconstructor) を実行します。"
    echo ""
    echo " [使い方]"
    echo " $0 <TargetDirectory>"
    echo ""
    echo " [引数]"
    echo " <TargetDirectory> : 処理対象のディレクトリパス"
    echo "                     この中にROOTファイルとペデスタルファイルが必要です。"
    echo ""
    echo " [前提条件]"
    echo " 1. ディレクトリ内に 'hkelec_pedestal_hithist_means.txt' が存在すること。"
    echo " 2. 同じディレクトリに実行ファイル 'reconstructor' が存在するか、"
    echo "    パスが通っていること (このスクリプトではカレントの ./reconstructor を想定)。"
    echo "=========================================================================="
}

# --- 引数チェック ---
if [ "$#" -lt 1 ]; then
    usage
    exit 1
fi

TARGET_DIR="$1"

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

# --- 実行ファイルの場所 (カレントディレクトリを想定) ---
RECONSTRUCTOR="./reconstructor"
if [ ! -f "$RECONSTRUCTOR" ]; then
    echo "Error: Executable '$RECONSTRUCTOR' not found in current directory."
    echo "Please compile the code using 'make' first."
    exit 1
fi

echo "Start processing in: $TARGET_DIR"

# --- ループ処理 ---
# findコマンドで *eventhist.root を検索して処理
find "$TARGET_DIR" -name "*eventhist.root" | while read input_file; do
    
    # 出力ファイル名の生成 (拡張子を除去して _reconst を付与)
    # 例: run001_eventhist.root -> run001_eventhist_reconst
    base_name=$(basename "$input_file" .root)
    output_base="${TARGET_DIR}/${base_name}_reconst"
    
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