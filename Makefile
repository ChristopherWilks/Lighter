CXX = g++ -g
CX = gcc
CXXFLAGS= -Wall -O0 #-O3
C_FLAGS= -std=c99 -Wall -Werror -pedantic -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls
C_DEFINES = "-DSFF=1"
LINKFLAGS = -lpthread -lz 
DEBUG=
OBJECTS = BOBHash.o ErrorCorrection.o GetKmers.o

# For Windows pthreads library: http://www.sourceware.org/pthreads-win32/
ifneq (,$(findstring MINGW,$(shell uname)))
	LINKFLAGS = -L. -lpthreadGC2
endif

all: torch

torch: main.o $(OBJECTS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJECTS) main.o $(LINKFLAGS)

main.o: main.cpp utils.h Reads.h Store.h File.h KmerCode.h bloom_filter.hpp
ErrorCorrection.o: ErrorCorrection.cpp ErrorCorrection.h utils.h
GetKmers.o: GetKmers.cpp GetKmers.h utils.h
BOBHash.o: BOBHash.c BOBHash.h
	$(CX) $(C_DEFINES) $(C_FLAGS) -o BOBHash.o -c BOBHash.c

clean:
	rm -f *.o *.gch lighter torch
