# 実行ファイル名
TARGET1 := gausfit
TARGET2 := peakfinder
# ソースファイル名
SRC1 := gausfit.C
SRC2 := peakfinder.C
# ROOTのコンパイルフラグとリンクフラグを取得
ROOTCFLAGS := $(shell root-config --cflags)
ROOTGLIBS := $(shell root-config --glibs)
# コンパイラとオプション
CXX := g++
CXXFLAGS := -O2 -Wall -fPIC $(ROOTCFLAGS)
# リンカフラグ
LDFLAGS := $(ROOTGLIBS)
.PHONY: all clean
# 'make all' または 'make' で両方の実行ファイルを作成
all: $(TARGET1) $(TARGET2)
# 高度なフィッターのためのルール
$(TARGET1): $(SRC1)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)
# シンプルなピークファインダーのためのルール
$(TARGET2): $(SRC2)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)
clean:
	rm -f $(TARGET1) $(TARGET2) *.o