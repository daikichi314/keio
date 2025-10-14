#!/bin/bash

# --- 1. 変数の定義と引数の検証 ---

# カスタムライブラリのパス
CUSTOM_LIB_PATH="/home/hkpd/hkelec/DiscreteSoftware/build/lib"

# 実行ファイル名
EXECUTABLE_NAME="./eventtree2hist"

# 引数のチェック: 対象ディレクトリ、計1つ必要
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <target_directory>"
    echo "例: $0 /home/hkpd/hkelec/DiscreteSoftware/data/20251009/100pe"
    exit 1
fi

TARGET_DIR=$1

# 指定された引数がディレクトリか確認
if [ ! -d "$TARGET_DIR" ]; then
    echo "エラー: 指定されたパス '$TARGET_DIR' はディレクトリではありません。"
    exit 1
fi

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

# --- 4. 複数ファイルの実行 ---

echo ""
echo "--- 4. ディレクトリ内の複数ファイルの解析を実行します ---"
echo "対象ディレクトリ: $TARGET_DIR"

# 指定されたディレクトリ内の *eventtree.root ファイルをループ処理
for INPUT_FILE in "$TARGET_DIR"/*eventtree.root; do
    
    # マッチするファイルが一つもなかった場合の処理
    if [ ! -f "$INPUT_FILE" ]; then
        echo "警告: 対象ディレクトリ内に *eventtree.root ファイルが見つかりませんでした。"
        break
    fi

    # 出力ファイル名を生成 (eventtree.root -> eventhist.root)
    OUTPUT_FILE="${INPUT_FILE/eventtree.root/eventhist.root}"

    echo ""
    echo "------------------------------------------------------------"
    echo "-> 処理中: $INPUT_FILE"
    echo "   出力先: $OUTPUT_FILE"
    echo "------------------------------------------------------------"

    # 実行ファイルを実行し、成功/失敗メッセージを表示
    ${EXECUTABLE_NAME} "$INPUT_FILE" "$OUTPUT_FILE"

    if [ $? -eq 0 ]; then
        echo "🎉 正常に完了しました。"
    else
        echo "⚠️ このファイルの処理中にエラーが発生しました。"
    fi
done

echo ""
echo "✅ 全てのファイルの処理が完了しました。"
