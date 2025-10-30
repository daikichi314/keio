# 実行ファイル名
TARGET := eventtree2hist
# ソースファイル名
SRC := eventtree2hist.C

# ---------------------------------------------------------------
# 1. 現在のユーザー名を取得
# ---------------------------------------------------------------
USER_NAME = $(shell whoami)

# ---------------------------------------------------------------
# 2. ユーザー名に基づいてカスタムパスを自動で切り替える
# ---------------------------------------------------------------
ifeq ($(USER_NAME), hkpd)
    # --- hkpd の場合のパス設定 ---
    CUSTOM_INCLUDE_DIR := /home/hkpd/hkelec/DiscreteSoftware/RootInterface
    CUSTOM_LIB_DIR := /home/hkpd/hkelec/DiscreteSoftware/build/lib

else ifeq ($(USER_NAME), daiki)
    # --- daiki の場合のパス設定 ---
    # (daiki のコピー先: /home/daiki/lab/hkelec_deps)
    CUSTOM_INCLUDE_DIR := /home/daiki/lab/hkelec_deps/RootInterface
    CUSTOM_LIB_DIR := /home/daiki/lab/hkelec_deps/build/lib

else
    # どちらでもない場合 (エラーメッセージを出して停止)
    $(error Unknown user $(USER_NAME). Please update Makefile.)
endif

# 3. リンクするライブラリ名
CUSTOM_LIBS := -lRootInterfaceLib

# ROOTのコンパイルフラグとリンクフラグを取得
# (注意) daiki, hkpd 両環境で `root-config` コマンドが正しく設定されている必要があります
ROOTCFLAGS := $(shell root-config --cflags)
ROOTGLIBS := $(shell root-config --glibs)

# コンパイラとオプション
CXX := g++
# CXXFLAGS: ROOTのフラグに加え、上記で設定したカスタムヘッダパスを追加
CXXFLAGS := -O2 -Wall -fPIC $(ROOTCFLAGS) -I$(CUSTOM_INCLUDE_DIR)

# LDFLAGS: ROOTのライブラリに加え、上記で設定したカスタムライブラリのパスと名前を追加
LDFLAGS := $(ROOTGLIBS) -L$(CUSTOM_LIB_DIR) $(CUSTOM_LIBS)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGET) *.o

