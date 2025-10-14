#$Id$
## Example to run
##    make
##    ./ttsfit 1000
## Note: Upper sigma is not Gaussian sigma, that shows better values because it's defined as sigma in the exponential modified Gaussian. 
CC = gcc
OBJ = ttsfit

CFLAGS= -O0 -g -Wall
CPPFLAGS += `root-config --cflags` -O0 -g -Wall
ROOTLIB += `root-config --libs` 
#DEBUGFLAGS = -DDEBUG
DEBUGFLAGS = 

### C++
ROOTCFLAGS = $(shell root-config --cflags)
ROOTMAIN = $(shell root-config --libs)
ROOTMAIN += -lSpectrum 
ROOTLIBS = -L$(ROOTSYS)/lib $(ROOTMAIN) 
LIBS = $(ROOTLIBS) $(CLHEPLIBS)
XXDEBUGFLAGS = -Wall -g -O0

CXXOPTIONS = $(CXXDEBUGFLAGS) $(ROOTCFLAGS) $(CLHEPFLAGS) $(CXXDEFINES)

#STDINC     = -I../../include


INCS = $(STDINC)

all	: $(OBJ)

ttsfit:	ttsfit.o
		$(CC) $(CFLAGS) $(CXXOPTIONS) $(LIBS) -lstdc++ $(DEBUGFLAGS) $(ROOTLIB) $< $(DAQLIB) $(STDINC) -lm -o $@

.c.o : 
	$(CC) $(INCS) $(ROOTLIB) $(CFLAGS) -c -o $@ $<

.cc.o : 
	$(CXX) $(INCS) $(ROOTLIB) $(CPPFLAGS) -c -o $@ $<

clean		:
		rm -f $(OBJ) *.o *~
