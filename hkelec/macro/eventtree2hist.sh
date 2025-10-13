#!/bin/bash

# --- 1. 変数の定義と引数の検証 ---

# カスタムライブラリのパス
CUSTOM_LIB_PATH="/home/hkpd/hkelec/DiscreteSoftware/build/lib"

# 実行ファイル名
EXECUTABLE_NAME="./eventtree2hist"

# 引数のチェック: 入力ファイルと出力ファイル、計2つ必要
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file.root> <output_file.root>"
    echo "例: $0 data/input.root output/hists.root"
    exit 1
fi

INPUT_FILE=$1
OUTPUT_FILE=$2

# --- 2. 環境変数の設定 ---

echo "--- 1. LD_LIBRARY_PATH を設定しています ---"
# LD_LIBRARY_PATHにカスタムライブラリのパスを追加
export LD_LIBRARY_PATH=${CUSTOM_LIB_PATH}:$LD_LIBRARY_PATH
echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

# --- 3. クリーンとコンパイル ---

echo ""
echo "--- 2. make clean で既存のビルドを削除します ---"
make clean

echo ""
echo "--- 3. make で再コンパイルします ---"
if make; then
    echo ""
    echo "✅ コンパイル成功: ${EXECUTABLE_NAME} が生成されました。"
else
    echo "❌ コンパイル失敗。スクリプトを終了します。"
    exit 1
fi

# --- 4. 実行 ---

echo ""
echo "--- 4. 解析を実行します ---"
echo "入力ファイル: $INPUT_FILE"
echo "出力ファイル: $OUTPUT_FILE"

# 実行ファイルを実行し、成功/失敗メッセージを表示
${EXECUTABLE_NAME} "$INPUT_FILE" "$OUTPUT_FILE"

if [ $? -eq 0 ]; then
    echo ""
    echo "🎉 解析実行が正常に完了しました。"
    echo "結果はファイル ${OUTPUT_FILE} に保存されました。"
else
    echo "⚠️ 解析実行中にエラーが発生しました。"
    exit 1
fi
