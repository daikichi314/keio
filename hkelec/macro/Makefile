# 実行ファイル名
TARGET := eventtree2hist
# ソースファイル名（ユーザー指定のファイル名に変更）
SRC := eventtree2hist.C

# 1. ヘッダファイルのあるパス (Hit.h, MetaData.hの親ディレクトリ)
CUSTOM_INCLUDE_DIR := /home/hkpd/hkelec/DiscreteSoftware/RootInterface
# 2. 共有ライブラリ (.soファイル) のあるパス
CUSTOM_LIB_DIR := /home/hkpd/hkelec/DiscreteSoftware/build/lib
# 3. リンクするライブラリ名
CUSTOM_LIBS := -lRootInterfaceLib

# ROOTのコンパイルフラグとリンクフラグを取得
ROOTCFLAGS := $(shell root-config --cflags)
ROOTGLIBS := $(shell root-config --glibs)

# コンパイラとオプション
CXX := g++
# 💡 CXXFLAGS: ROOTのフラグに加え、カスタムヘッダパスを追加
CXXFLAGS := -O2 -Wall -fPIC $(ROOTCFLAGS) -I$(CUSTOM_INCLUDE_DIR)

# 💡 LDFLAGS: ROOTのライブラリに加え、カスタムライブラリのパスと名前を追加
LDFLAGS := $(ROOTGLIBS) -L$(CUSTOM_LIB_DIR) $(CUSTOM_LIBS)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGET) *.o