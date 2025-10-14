#$Id$

# C++コンパイラをg++に設定
CXX = g++
# Cコンパイラも同じものを使う
CC = $(CXX)

OBJ = simplefit

# ROOTのコンパイルとリンクフラグを取得
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

# 共通のコンパイラフラグ
CFLAGS = -O0 -g -Wall -D rcomp
CPPFLAGS = $(ROOTCFLAGS) $(CFLAGS)

# 共通のライブラリ（lSpectrumはroot-config --libsに含まれない場合があるため明示的に追加）
LIBS = $(ROOTLIBS) -lSpectrum

# デバッグ用フラグ
DEBUGFLAGS = 

all: $(OBJ)

simplefit: simplefit.o
	$(CXX) $(CPPFLAGS) $< -o $@ $(LIBS)

# .c (小文字) ファイル用のルール
.c.o:
	$(CC) $(CPPFLAGS) -c -o $@ $<

# .C (大文字) ファイル用のルール
.C.o:
	$(CXX) $(CPPFLAGS) -c -o $@ $<

# .cc (小文字) ファイル用のルール
.cc.o:
	$(CXX) $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJ) *.o *~
