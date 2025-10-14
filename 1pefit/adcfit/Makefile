#$Id$
CC = gcc
CXX = g++ #####0829追加
OBJ = performance1pe

CFLAGS= -O0 -g -Wall -D rcomp
CPPFLAGS += `root-config --cflags` -O0 -g -Wall -D rcomp
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



INCS = $(STDINC)

all	: $(OBJ)

# performance1pe:	performance1pe.o
# 		$(CC) $(CFLAGS) $(CXXOPTIONS) $(LIBS) -lstdc++ $(DEBUGFLAGS) $(ROOTLIB) $< $(DAQLIB) $(STDINC) -lm -o $@ ##0829変更

performance1pe: performance1pe.o
	$(CXX) $(CXXOPTIONS) $^ $(LIBS) -o $@

.c.o : 
	$(CC) $(INCS) $(ROOTLIB) $(CFLAGS) -c -o $@ $<

.cc.o : 
	$(CXX) $(INCS) $(ROOTLIB) $(CPPFLAGS) -c -o $@ $<

clean		:
		rm -f $(OBJ) *.o *~

# #$Id$
# CC = gcc
# OBJ = simplefit

# CFLAGS= -O0 -g -Wall -D rcomp
# CPPFLAGS += `root-config --cflags` -O0 -g -Wall -D rcomp
# ROOTLIB += `root-config --libs` 
# #DEBUGFLAGS = -DDEBUG
# DEBUGFLAGS = 

# ### C++
# ROOTCFLAGS = $(shell root-config --cflags)
# ROOTMAIN = $(shell root-config --libs)
# ROOTMAIN += -lSpectrum 
# ROOTLIBS = -L$(ROOTSYS)/lib $(ROOTMAIN) 
# LIBS = $(ROOTLIBS) $(CLHEPLIBS)
# XXDEBUGFLAGS = -Wall -g -O0

# CXXOPTIONS = $(CXXDEBUGFLAGS) $(ROOTCFLAGS) $(CLHEPFLAGS) $(CXXDEFINES)



# INCS = $(STDINC)

# all	: $(OBJ)

# performance1pe:	performance1pe.o
# 		$(CC) $(CFLAGS) $(CXXOPTIONS) $(LIBS) -lstdc++ $(DEBUGFLAGS) $(ROOTLIB) $< $(DAQLIB) $(STDINC) -lm -o $@

# simplefit:	simplefit.o
# 		$(CC) $< $(CFLAGS) $(CXXOPTIONS) $(LIBS) -lstdc++ $(DEBUGFLAGS) $(ROOTLIB) $(DAQLIB) $(STDINC) -lm -o $@

# .c.o : 
# 	$(CC) $(INCS) $(ROOTLIB) $(CFLAGS) -c -o $@ $<

# .cc.o : 
# 	$(CXX) $(INCS) $(ROOTLIB) $(CPPFLAGS) -c -o $@ $<

# clean		:
# 		rm -f $(OBJ) *.o *~

