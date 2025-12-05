# 1. コンパイラの設定
# 通常は g++ を使用します
CXX = g++

# 2. コンパイルオプションの設定
# -O2: 最適化レベル2
# -Wall: 警告を表示
# root-config --cflags: ROOTのヘッダファイルのパスなどを自動設定
CXXFLAGS = -O2 -Wall $(shell root-config --cflags)

# 3. リンクオプションの設定
# root-config --libs: ROOTのライブラリを自動設定
# -lMinuit: TMinuitクラスを使用するために明示的にリンク
LDFLAGS = $(shell root-config --libs) -lMinuit

# 4. 生成する実行ファイル名
# 任意の名前で構いませんが、ここでは reconstructor としています
TARGET = reconstructor

# 5. ソースファイル
# main.cc 内で他の .cc ファイルを include しているため、
# コンパイル対象は main.cc のみとなります
SRCS = main.cc

# 6. ビルドターゲット (make コマンドで実行される部分)
all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# 7. クリーンアップ (make clean で生成物を削除)
clean:
	rm -f $(TARGET)