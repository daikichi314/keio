#!/bin/bash

# --- 1. 変数の定義 ---

# このスクリプトファイル自体の絶対パスを取得し、それを基準に各パスを定義
# これにより、どこからスクリプトを実行しても正しく動作するようになります
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
echo "スクリプトの場所: ${SCRIPT_DIR}"

# カスタムライブラリのパス（これは絶対パスなので変更不要）
CUSTOM_LIB_PATH="/home/hkpd/hkelec/DiscreteSoftware/build/lib"

# 実行ファイルのフルパスを定義
EXECUTABLE_PATH="${SCRIPT_DIR}/eventtree2hist"

# --- 2. 引数の検証 ---

# 引数のチェック: 対象ディレクトリ、計1つ必要
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <target_directory>"
    echo "例: $0 /home/hkpd/hkelec/DiscreteSoftware/data/20251009/100pe"
    echo "例 (相対パス): $0 ../data/20251009/100pe"
    exit 1
fi

TARGET_DIR=$1

# 指定された引数がディレクトリか確認
if [ ! -d "$TARGET_DIR" ]; then
    echo "エラー: 指定されたパス '$TARGET_DIR' はディレクトリではありません。"
    exit 1
fi

# --- 3. 環境変数の設定 ---

echo "--- 1. LD_LIBRARY_PATH を設定しています ---"
# LD_LIBRARY_PATHにカスタムライブラリのパスを追加
export LD_LIBRARY_PATH=${CUSTOM_LIB_PATH}:$LD_LIBRARY_PATH
echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

# --- 4. クリーンとコンパイル ---

echo ""
echo "--- 2. make clean を実行します (in ${SCRIPT_DIR}) ---"
# -C オプションで、Makefileのあるディレクトリを指定してmakeを実行
if ! make -C "${SCRIPT_DIR}" clean; then
    echo "警告: 'make clean' に失敗しましたが、処理を続行します。"
fi


echo ""
echo "--- 3. make で再コンパイルします (in ${SCRIPT_DIR}) ---"
# 同様に、-C オプションでディレクトリを指定
if make -C "${SCRIPT_DIR}"; then
    echo ""
    echo "✅ コンパイル成功: ${EXECUTABLE_PATH} が生成されました。"
else
    echo "❌ コンパイル失敗。スクリプトを終了します。"
    exit 1
fi

# --- 5. 複数ファイルの実行 ---

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

    # 実行ファイルをフルパスで指定して実行
    "${EXECUTABLE_PATH}" "$INPUT_FILE" "$OUTPUT_FILE"

    if [ $? -eq 0 ]; then
        echo "🎉 正常に完了しました。"
    else
        echo "⚠️ このファイルの処理中にエラーが発生しました。"
    fi
done

echo ""
echo "✅ 全てのファイルの処理が完了しました。"
