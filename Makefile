CXX = g++ -g -std=c++11
CX = gcc
CXXFLAGS= -Wall -O0 #-O3
C_FLAGS= -std=c99 -Wall -Werror -pedantic -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls
C_DEFINES = "-DSFF=1"
LINKFLAGS = -lpthread -lz 
DEBUG=
OBJECTS = BOBHash.o cmlsketch.o ErrorCorrection.o GetKmers.o

# For Windows pthreads library: http://www.sourceware.org/pthreads-win32/
ifneq (,$(findstring MINGW,$(shell uname)))
	LINKFLAGS = -L. -lpthreadGC2
endif

all: torch

torch: main.o $(OBJECTS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJECTS) main.o $(LINKFLAGS)

main.o: main.cpp utils.h Reads.h Store.h StoreCML.h File.h KmerCode.h cmlsketch.h bloom_filter.hpp
ErrorCorrection.o: ErrorCorrection.cpp ErrorCorrection.h utils.h
GetKmers.o: GetKmers.cpp GetKmers.h utils.h
cmlsketch.o: cmlsketch.cpp cmlsketch.h bench_config.h bench_common.h define.h BOBHash.h
BOBHash.o: BOBHash.cpp BOBHash.h

clean:
	rm -f *.o *.gch lighter torch
