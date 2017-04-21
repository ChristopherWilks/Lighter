CXX = g++ -g -I./
CX = gcc -I./
CXXFLAGS= -std=c++11 -Wall -O3 #-O0 #-O3
C_FLAGS= -std=c99 -Wall -Werror -pedantic -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls
C_DEFINES = "-DSFF=1"
LINKFLAGS = -L./ -lpthread -lz -lbf
DEBUG=
OBJECTS = ErrorCorrection.o GetKmers.o

# For Windows pthreads library: http://www.sourceware.org/pthreads-win32/
ifneq (,$(findstring MINGW,$(shell uname)))
	LINKFLAGS = -L. -lpthreadGC2
endif

all: torch

torch: main.o $(OBJECTS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJECTS) main.o $(LINKFLAGS)

main.o: main.cpp utils.h Reads.h Store.h StoreBF.h File.h KmerCode.h bloom_filter.hpp bf.h
ErrorCorrection.o: ErrorCorrection.cpp ErrorCorrection.h utils.h
GetKmers.o: GetKmers.cpp GetKmers.h utils.h

clean:
	rm -f *.o *.gch lighter torch
