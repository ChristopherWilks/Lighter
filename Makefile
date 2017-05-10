CXX = g++ -g
CX = gcc
CXXFLAGS= -Wall -O3
C_FLAGS= -std=c99 -Wall -Werror -pedantic -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls
C_DEFINES = "-DSFF=1"
LINKFLAGS = -lpthread -lz -lrt
DEBUG=
OBJECTS = BOBHash.o sketch_config.o sffsketch.o ErrorCorrection.o GetKmers.o

# For Windows pthreads library: http://www.sourceware.org/pthreads-win32/
ifneq (,$(findstring MINGW,$(shell uname)))
	LINKFLAGS = -L. -lpthreadGC2
endif

all: torch torchq

torch: main.o $(OBJECTS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJECTS) main.o $(LINKFLAGS)

torchq: mainq.o $(OBJECTS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJECTS) mainq.o $(LINKFLAGS)

mainq.o: mainq.cpp utils.h Reads.h Store.h StoreSF.h File.h KmerCode.h sketch.h bloom_filter.hpp
main.o: main.cpp utils.h Reads.h Store.h StoreSF.h File.h KmerCode.h sketch.h bloom_filter.hpp
ErrorCorrection.o: ErrorCorrection.cpp ErrorCorrection.h utils.h
GetKmers.o: GetKmers.cpp GetKmers.h utils.h
sketch_config.o: sketch_config.c sketch_config.h sketch.h
	$(CX) $(C_DEFINES) $(C_FLAGS) -o sketch_config.o -c sketch_config.c
sffsketch.o: sffsketch.c sffsketch.h sketch_config.h sketch.h
	$(CX) $(C_DEFINES) $(C_FLAGS) -o sffsketch.o -c sffsketch.c
BOBHash.o: BOBHash.c BOBHash.h
	$(CX) $(C_DEFINES) $(C_FLAGS) -o BOBHash.o -c BOBHash.c

clean:
	rm -f *.o *.gch lighter torch torchq
