#CC=icpc
#OPTS = -openmp -O3 -lgsl -lgslcblas -lm $(EXT)
CC=g++
OPTS = -O3 -lgsl -lgslcblas -lsqlite3 $(EXT)
#OPTS = -fopenmp -O3 -lgsl -lgslcblas -lm $(EXT)
#ARCH=$(shell ARC=`uname -a |grep "x86_64"`;\
#	     if [ -n "$ARC" ]; then \
#		ARC="em64t"; \
#	     else \
#		ARC="ia32"; \
#	     fi \
#	     echo $ARC)
#OPTS=-lmkl_$(ARCH) -openmp -O3
#OPTS = -lmkl_em64t -lmkl_ia32 -openmp -O3 $(EXT)
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
HEADERS = $(wildcard *.h)
EXE = hamf.x

default: hamf.x

hamf.x: $(OBJECTS) $(HEADERS)
	$(CC) $(OPTS) $(OBJECTS) -o $(EXE) 

debug:
	make clean
	make "EXT=-g"
	-rm *.o #next make is clean

profile:
	make clean
	make "EXT=-g -p"
	-rm *.o #next make is clean

%.o: %.cpp
	$(CC) $(OPTS) -c $(@:.o=.cpp) 

clean:
	-rm *.o	
	-rm *.x
	-rm *.dat


