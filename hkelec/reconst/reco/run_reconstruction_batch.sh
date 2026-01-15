#!/bin/bash

# id: run_reconstruction_batch.sh
# Place: /home/daiki/keio/hkelec/reconst/reco/
# Author: Gemini (Modified based on user request)
# Last Edit: 2025-01-08
#
# 概要:
# 指定ディレクトリ内の *eventhist.root ファイルに対して
# reconstructor を一括実行するバッチスクリプト
#
# [重要]
# 出力ファイル名はC++プログラム側で、指定されたオプション(-m, -q, -tなど)
# に基づいて自動的に決定されます（例: run01_reconst_3hits_bc_func_f_goodness.root）。
# これにより、条件違いでの結果の上書きを防ぎます。
#
# 使い方:
# $ ./run_reconstruction_batch.sh <TargetDirectory> [Options...]
# 例:
# $ ./run_reconstruction_batch.sh ./data -u 1 -m func_g -q bc -t goodness
#
# 引数:
# <TargetDirectory> : 処理対象のディレクトリパス (必須)
#                     この中にROOTファイルとペデスタルファイルが必要です。
# [Options...]      : reconstructor にそのまま渡されるオプション群
#                     -u <0/1>        : 0=4本必須(default), 1=3本許容(Unhit補完)
#                     -m <model>      : func_f=半径28.5cmモデル(default)
#                                      func_g=半径23.5cmモデル(1/r^2)
#                     -q <model>      : gaus=ガウス(default), bc=Baker-Cousins
#                                      none=電荷情報を使用しない (時間のみでフィット)
#                     -t <model>      : gaus=ガウス(default), emg=EMG, goodness=SK Goodness
#                                      none=時間情報を使用しない (電荷のみでフィット)
#
# 前提条件:
# 1. 指定ディレクトリに 'hkelec_pedestal_hithist_means.txt' が存在すること
# 2. このスクリプトと同じディレクトリに 'reconstructor' 実行ファイルが存在すること


# --- スクリプトの使い方を表示する関数 ---
usage() {
    echo "=========================================================================="
    echo " [概要]"
    echo " 指定されたディレクトリ内のすべての *eventhist.root ファイルに対して"
    echo " 再構成プログラム (reconstructor) を実行します。"
    echo ""
    echo " [使い方]"
    echo " $0 <TargetDirectory> [Options...]"
    echo ""
    echo " [引数]"
    echo " <TargetDirectory> : 処理対象のディレクトリパス (必須)"
    echo "                     この中にROOTファイルとペデスタルファイルが必要です。"
    echo " [Options...]      : reconstructor にそのまま渡されるオプション群"
    echo ""
    echo "   -u <0/1>        : 0=4本必須(default), 1=3本許容(Unhit補完)"
    echo "   -m <model>      : func_f=半径28.5cmモデル(default)"
    echo "                     func_g=半径23.5cmモデル(1/r^2)"
    echo "                     all=func_f と func_g の両方で解析"
    echo "   -q <model>      : gaus=ガウス(default), bc=Baker-Cousins"
    echo "                     none=電荷情報を使用しない (時間のみでフィット)"
    echo "   -t <model>      : gaus=ガウス(default), emg=EMG, goodness=SK Goodness"
    echo "                     none=時間情報を使用しない (電荷のみでフィット)"
    echo ""
    echo " [実行例]"
    echo " 1. デフォルト設定 (4本, FuncF, Gaussian):"
    echo "    $0 ./data"
    echo ""
    echo " 2. 3本ヒット許容、モデルG、Baker-Cousins、Goodness:"
    echo "    $0 ./data -u 1 -m func_g -q bc -t goodness"
    echo ""
    echo " 3. 両方のモデル(FuncF と FuncG)で解析を実行:"
    echo "    $0 ./data -m all"
    echo ""
    echo " 4. 時間情報のみを使ってフィット (電荷不使用):"
    echo "    $0 ./data -q none"
    echo ""
    echo " [前提条件]"
    echo " 1. ディレクトリ内に 'hkelec_pedestal_hithist_means.txt' が存在すること。"
    echo " 2. このスクリプトと同じディレクトリに実行ファイル 'reconstructor' が存在すること。"
    echo "=========================================================================="
}

# --- 引数チェック ---
if [ "$#" -lt 1 ]; then
    usage
    exit 1
fi

TARGET_DIR="$1"
shift # ディレクトリ引数をずらす

# 残りの引数をすべてオプションとして格納
RECO_OPTIONS="$@"

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

echo "========================================================"
echo " Start Batch Reconstruction"
echo " Directory : $TARGET_DIR"
echo " Options   : ${RECO_OPTIONS:-(Default)}"
echo "========================================================"

# --- ループ処理 ---
# findで検索し、ファイル名順にソートして処理
find "$TARGET_DIR" -name "*eventhist.root" | sort | while read input_file; do
    
    echo "--------------------------------------------------------"
    echo "Processing: $(basename "$input_file")"
    
    # reconstructor の実行
    # 引数: 入力ファイル オプション...
    # ※出力ファイル名はプログラム内部で自動生成されるため指定不要
    $RECONSTRUCTOR "$input_file" $RECO_OPTIONS
    
    if [ $? -eq 0 ]; then
        echo " -> Done."
    else
        echo " -> Error occurred!"
    fi

done

echo "--------------------------------------------------------"
echo "All tasks finished."