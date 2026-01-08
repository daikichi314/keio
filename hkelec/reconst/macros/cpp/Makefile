#
# id: Makefile
# Place: /home/daiki/keio/hkelec/reconst/macros/cpp/
# Last Edit: 2025-11-21 Gemini
#
# 概要: cppディレクトリ内のC++ソースコードをコンパイルするためのMakefile。
#       fit_pedestal, meanfinder, plot_summary を対象とする。
#
# 実行可能 (makeコマンド)
#

# --- 設定 ---
CXX := g++

# ROOTのコンパイルフラグとリンクフラグを取得
ROOTCFLAGS := $(shell root-config --cflags)
ROOTGLIBS := $(shell root-config --glibs)

CXXFLAGS := -O2 -Wall -fPIC $(ROOTCFLAGS)
LDFLAGS := $(ROOTGLIBS)

# --- ターゲット定義 ---
TARGET_PEDESTAL := fit_pedestal
TARGET_MEANFINDER := meanfinder
TARGET_PLOT := plot_summary
TARGET_PLOT_CSV := plot_from_csv

# ソースファイル名
SRC_PEDESTAL := fit_pedestal.C
SRC_MEANFINDER := meanfinder.C
SRC_PLOT := plot_summary.C
SRC_PLOT_CSV := plot_from_csv.C

.PHONY: all clean

# --- ルール ---

# 'make all' または 'make' で全ての実行ファイルを作成
all: $(TARGET_PEDESTAL) $(TARGET_MEANFINDER) $(TARGET_PLOT) $(TARGET_PLOT_CSV)

# fit_pedestal
$(TARGET_PEDESTAL): $(SRC_PEDESTAL)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# meanfinder
$(TARGET_MEANFINDER): $(SRC_MEANFINDER)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# plot_summary
$(TARGET_PLOT): $(SRC_PLOT)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# plot_from_csv
$(TARGET_PLOT_CSV): $(SRC_PLOT_CSV)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# クリーンアップ (生成された実行ファイルとオブジェクトファイルを削除)
clean:
	rm -f $(TARGET_PEDESTAL) $(TARGET_MEANFINDER) $(TARGET_PLOT) $(TARGET_PLOT_CSV) *.o