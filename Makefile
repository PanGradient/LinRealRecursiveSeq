CXX=g++
CXXFLAGS=-std=c++11 -Wall -Wextra -O2 -s
CPPFLAGS=-I. -DHAVE_INLINE
LDFLAGS=-shared
LDLIBS=-lgsl -lgslcblas
DIR=build
EXAMPLES_DIR=examples

libLinRealRecursiveSeq.so: LinRealRecursiveSeq.o
	$(CXX) $(LDFLAGS) -o libLinRealRecursiveSeq.so $(LDLIBS) LinRealRecursiveSeq.o

LinRealRecursiveSeq.o: LinRealRecursiveSeq.cc
	$(CXX) -c -fPIC $(CXXFLAGS) $(CPPFLAGS) -o LinRealRecursiveSeq.o LinRealRecursiveSeq.cc

examples:
	mkdir -p $(EXAMPLES_DIR)/bin
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDLIBS) -o $(EXAMPLES_DIR)/bin/fibonacci $(EXAMPLES_DIR)/fibonacci.cc LinRealRecursiveSeq.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDLIBS) -o $(EXAMPLES_DIR)/bin/perrin    $(EXAMPLES_DIR)/perrin.cc    LinRealRecursiveSeq.cc

clean:
	rm -rf *.o *.so $(EXAMPLES_DIR)/bin

.PHONY: clean examples
