CXX = g++ -g -std=c++11
CX = gcc -g
CXXFLAGS= -Wall -O3
#CXXFLAGS= -Wall #-O3
C_FLAGS= -std=c99 -Wall -Werror -pedantic -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls
C_DEFINES = "-DSFF=1"
LINKFLAGS = -lpthread -lz -lrt
DEBUG=
ARCH=
OBJECTS = gqf.o hashutil.o ErrorCorrection.o GetKmers.o

# For Windows pthreads library: http://www.sourceware.org/pthreads-win32/
ifneq (,$(findstring MINGW,$(shell uname)))
	LINKFLAGS = -L. -lpthreadGC2
endif

all: torch torchq

torch: main.o $(OBJECTS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJECTS) main.o $(LINKFLAGS)
torchq: mainq.o $(OBJECTS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJECTS) mainq.o $(LINKFLAGS)

mainq.o: mainq.cpp utils.h Reads.h Store.h hashutil.h StoreCQF.h File.h KmerCode.h sketch.h bloom_filter.hpp
main.o: main.cpp utils.h Reads.h Store.h hashutil.h StoreCQF.h File.h KmerCode.h sketch.h bloom_filter.hpp
ErrorCorrection.o: ErrorCorrection.cpp ErrorCorrection.h utils.h
GetKmers.o: GetKmers.cpp GetKmers.h utils.h
gqf.o: gqf.c gqf.h
	$(CXX)  -o gqf.o $(CXXFLAGS) -c gqf.c
hashutil.o: hashutil.cc hashutil.h
	$(CXX)  -o hashutil.o $(CXXFLAGS) -c hashutil.cc

clean:
	rm -f *.o *.gch lighter torch
