# コンパイラの指定 (g++)
CXX = g++

# コンパイルオプション
# -O2: 最適化レベル2
# -Wall: 警告を全て表示
# root-config --cflags: ROOTライブラリのインクルードパスなどを展開
CXXFLAGS = -O2 -Wall $(shell root-config --cflags)

# リンクオプション
# root-config --libs: ROOTライブラリのリンク設定を展開
# -lMinuit: Minuitライブラリを明示的にリンク
LDFLAGS = $(shell root-config --libs) -lMinuit

# 生成する実行ファイルの名前
TARGET = reconstructor

# ソースファイルのリスト
SRCS = main.cc readData.cc onemPMTfit.cc

# オブジェクトファイル名 (.cc を .o に置換)
OBJS = $(SRCS:.cc=.o)

# デフォルトターゲット (make と打つとここが実行される)
all: $(TARGET)

# 実行ファイルの生成ルール
# $@ はターゲット名(reconstructor), $^ は依存ファイルリスト(OBJS)
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# 各ソースファイルのコンパイルルール
# $< は最初の依存ファイル(.cc), $@ はターゲット(.o)
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 生成ファイルを削除するターゲット
clean:
	rm -f $(TARGET) $(OBJS)