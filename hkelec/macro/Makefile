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

# LDFLAGS: リンク順序とオプションを再修正しました
# [修正点]
# 1. -Wl,--no-as-needed: 不要とみなされたライブラリの除外を防ぐ
# 2. -Wl,--allow-shlib-undefined: 自作ライブラリ内の未解決シンボル(TROOT::SetBatch等)によるエラーを抑制し、実行時解決に委ねる
LDFLAGS := -Wl,--no-as-needed -Wl,--allow-shlib-undefined -L$(CUSTOM_LIB_DIR) $(CUSTOM_LIBS) $(ROOTGLIBS)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGET) *.o