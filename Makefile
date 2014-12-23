CXX=g++
CXXFLAGS=-fPIC -std=c++11 -Wall -Wextra -O2 -s
CPPFLAGS=-I. -DHAVE_INLINE
LDFLAGS=-shared
LDLIBS=-lgsl -lgslcblas
DIR=build
EXAMPLES=examples

.PHONY: examples install clean

libLinRealRecursiveSeq.so: LinRealRecursiveSeq.o
	$(CXX) $(LDFLAGS) -o libLinRealRecursiveSeq.so $(LDLIBS) LinRealRecursiveSeq.o

LinRealRecursiveSeq.o: LinRealRecursiveSeq.cc
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o LinRealRecursiveSeq.o LinRealRecursiveSeq.cc

examples:
	mkdir -p $(EXAMPLES)/bin
	g++ -std=c++11 -Wall -Wextra -O2 -s -I. $(shell pwd)/libLinRealRecursiveSeq.so -o $(EXAMPLES)/bin/fibonacci $(EXAMPLES)/fibonacci.cc
	g++ -std=c++11 -Wall -Wextra -O2 -s -I. $(shell pwd)/libLinRealRecursiveSeq.so -o $(EXAMPLES)/bin/perrin    $(EXAMPLES)/perrin.cc

install:
	mkdir -p $(DIR)
	cp -f LinRealRecursiveSeq.h $(DIR)
	cp -f libLinRealRecursiveSeq.so $(DIR)

clean:
	rm -rf *.o *.so $(EXAMPLES)/bin
