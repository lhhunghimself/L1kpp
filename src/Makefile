CXX=g++
CC=gcc
FF=gfortran
L_PATH = /usr/local/lib
LD=g++ 
LOPENMP= -lgomp 
OPT_FLAGS= -static -Ofast -fopenmp -ftree-vectorize -funroll-all-loops -m64 -msse3 -mavx
ifdef DEBUG
 OPT_FLAGS= -ggdb -m64 -msse3 -mavx
endif
FFLAGS=-fimplicit-none
CFLAGS= -I ./ -I ./FCSTools -I ./tclap-1.2.1/include -I $(OPT_FLAGS) -std=c++11 -Wno-write-strings -fexpensive-optimizations -ffast-math -finline-functions -frerun-loop-opt -static $(HDRS) $(DEFINES)
LIBS= -lm -lc -lgcc -lboost_system -lboost_filesystem -lboost_math_tr1 -lopenblas -lgfortran -lrt $(LOPENMP) -L$(L_PATH)

OBJ_L1k++ = L1k++.o
	
L1k++ : $(OBJ_L1k++)
	$(LD) -o L1k++ $(LFLAGS) $(OBJ_L1k++) $(LIBS)
L1k++.o: L1k++.cpp
	$(CXX) $(CFLAGS) -c L1k++.cpp

clean:
	rm *.o
	rm core.*

# all together now...
default: L1k++
