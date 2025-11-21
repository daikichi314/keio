#
# id: Makefile
# Place: ~/hkelec/DiscreteSoftware/Analysis/macro/fit_results/
# Last Edit: 2025-10-18 Gemini
#
# 概要: fit_results/ ディレクトリ内の全てのC++ソースコードをコンパイルし、
#       実行ファイルを作成するための設計図。
# 実行可能 (makeコマンド)
#
TARGET1 := gausfit
TARGET2 := fit_pedestal
TARGET3 := select_gain
TARGET4 := peakfinder
TARGET5 := fit_hv_gain
TARGET6 := create_ct_plot
TARGET7 := meanfinder
TARGET8 := select_gain_mean
TARGET9 := plot_ct

TARGET10 := fit_hv_gain_2

# ソースファイル名
SRC1 := gausfit.C
SRC2 := fit_pedestal.C
SRC3 := select_gain.C
SRC4 := peakfinder.C
SRC5 := fit_hv_gain.C
SRC6 := create_ct_plot.C
SRC7 := meanfinder.C
SRC8 := select_gain_mean.C
SRC9 := plot_ct.C

SRC10 := fit_hv_gain_2.C
# ヘッダーファイル
HEADER1 := tts_fitter.h

# ROOTのコンパイルフラグとリンクフラグを取得
ROOTCFLAGS := $(shell root-config --cflags)
ROOTGLIBS := $(shell root-config --glibs)
CXX := g++
CXXFLAGS := -O2 -Wall -fPIC $(ROOTCFLAGS)
LDFLAGS := $(ROOTGLIBS)
.PHONY: all clean

# 'make all' または 'make' で全ての実行ファイルを作成
all: $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7) $(TARGET8) $(TARGET9)

# gausfit: tts_fitter.hにも依存<-no
$(TARGET1): $(SRC1) #$(HEADER1)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# fit_pedestal
$(TARGET2): $(SRC2)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# select_gain
$(TARGET3): $(SRC3)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# peakfinder: tts_fitter.hにも依存<-no
$(TARGET4): $(SRC4) #$(HEADER1)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# fit_hv_gain
$(TARGET5): $(SRC5)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# create_ct_plot
$(TARGET6): $(SRC6)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# meanfinder
$(TARGET7): $(SRC7)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# select_gain_mean
$(TARGET8): $(SRC8)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# plot_ct
$(TARGET9): $(SRC9)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# fit_hv_gain_2
$(TARGET10): $(SRC10)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7) $(TARGET8) $(TARGET9) $(TARGET10) *.o
