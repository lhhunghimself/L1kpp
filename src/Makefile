CXX=g++
CC=gcc
FF=gfortran
L_PATH = /usr/local/lib
LD=g++ 
LOPENMP= -lgomp 
OPT_FLAGS= -static -Ofast -fopenmp -ftree-vectorize -funroll-all-loops -m64 -msse3 -mavx
ifdef DEBUG
 OPT_FLAGS= -g -ggdb -m64 -msse3 -mavx
endif
FFLAGS=-fimplicit-none
CFLAGS= -I ./  -I ./tclap-1.2.1/include -I $(OPT_FLAGS) -std=c++11 -Wno-write-strings -fexpensive-optimizations -ffast-math -finline-functions -frerun-loop-opt -static $(HDRS) $(DEFINES)
LIBS= -lm -lc -lgcc -lboost_system -lboost_filesystem -lboost_math_tr1 -lopenblas -lgfortran -lrt $(LOPENMP) -L$(L_PATH)

OBJ_L1k++ = L1k++.o
	
L1k++ : $(OBJ_L1k++)
	$(LD) -o L1k++ $(LFLAGS) $(OBJ_L1k++) $(LIBS)
L1k++.o: L1k++.cpp gmm.hpp qNorm.hpp common.hpp
	$(CXX) $(CFLAGS) -c L1k++.cpp

clean:
	-@rm -rf *.o 2>/dev/null || true
	-@rm -rf core.* 2>/dev/null || true
	-@rm -rf 'L1k++' 2>/dev/null || true

# all together now...
default: L1k++
