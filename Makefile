CXX = g++ -g -I./
CX = gcc -I./
CXXFLAGS= -std=c++11 -Wall -O3 #-O0 #-O3
C_FLAGS= -std=c99 -Wall -c -Werror -pedantic -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -fno-strict-aliasing -I./ -I/usr/include/
C_DEFINES = "-DSFF=1"
LINKFLAGS = -Wall -L./ -lpthread -lssl -lcrypto -lz
DEBUG=
OBJECTS = ErrorCorrection.o GetKmers.o hashutil.o #printutil.o 

HEADERS = $(wildcard ./*.h)


all: lighter

lighter: main.o $(OBJECTS) ${HEADERS}
	$(CXX) -o $@ $(CXXFLAGS) $(OBJECTS) main.o $(LINKFLAGS)

%.o: %.cpp ${HEADERS}

#Makefile
#$(CXX) $(CXXFLAGS) $< -o $@
#main.o: main.cpp utils.h Reads.h Store.h StoreCUCKOO.h File.h KmerCode.h bloom_filter.hpp
#ErrorCorrection.o: ErrorCorrection.cpp ErrorCorrection.h utils.h
#GetKmers.o: GetKmers.cpp GetKmers.h utils.h

clean:
	rm -f *.o *.gch lighter
