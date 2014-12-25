LinRealRecursiveSeq
===================

A simple C++ class allowing for relatively fast calculation of a recursive sequence elements (with constant coefficients). Elements and recurrence relation coefficients are assumed to be `double` variables.

This package provides a header file, defining the class, a source file containing its implementation and a Makefile for building a shared library.

# Prerequisites

* [GCC 4.7](http://gcc.gnu.org/) (or later) - provides support for constructor delegation (part of C++11 standard)
* [GSL 1.9](http://www.gnu.org/software/gsl/) (or later) - provides nonsymmetric matrix diagonalization method, GSL CBLAS Library and GSL BLAS Interface
* (optional) other BLAS library ([OpenBLAS](http://www.openblas.net/), [ATLAS](http://math-atlas.sourceforge.net/), ...) - alternative for GSL CBLAS Library

# Installation

Simplify type `make` in the main directory. This will locally build a shared library `libLinRealREcursiveSeq.so`. Both library and header file must be put in appropriate directories (i.e. those being part of `CPLUS_INCLUDE_PATH` and `LD_LIBRARY_PATH` variables).

If you want to compile example programs, type `make examples`. This will build executable files in the `examples/bin` directory.

If you want to clear all compiled files, type `make clear`.

This class uses GSL CBLAS Library. If user wants to use another BLAS library, please change `-lgslcblas` in `LDLIBS` variable in Makefile to appropriate one.

# Usage

`LinRealRecursiveSeq` class requires two `std::vector` objects or two `double` arrays (with their sizes provided by user), both of equal length. The first has to be composed of the first N elements of a sequence. The second has to contain constant coefficients of the recurrence relation. For instance, for Fibonnaci numbers `fib(n) = fib(n-1) + fib(n-2)`, the first array would be `{0.0, 1.0}` and the second `{1.0, 1.0}`. For details please read `doc/reference.pdf`.

## Compiling with the shared library

If user wants use the shared library, he/she has to compile his/her program in one of possible ways:
* `g++ <compilation/linker flags> -lLinRealRecursiveSeq <source file>` - if the header file and the shared library are in directories provided by `CPLUS_INCLUDE_PATH` and `LD_LIBRARY_PATH` variables
* `g++ <compilation/linker flags> -I<header directory> -L<shared library directory> -lLinRealRecursiveSeq <source file>` - if the header file and the shared library aren't in those directories (unrecommended, since `LD_LIBRARY_PATH` would have to be exported to contain shared library directory each time the program is executed)
* `g++ <compilation/linker flags> -I<header directory> <shared library directory>/libLinRealRecursiveSeq.so <source file>` - if the header file and the shared library aren't in those directories and user wants to specify the path (absolute or relative) to the shared library (unrecommended)

## Compiling with the source file

If user wants compile his/her program with provided source code as part of his/her own project, it can be done if appropriate g++ flags are used:

```
g++ -c <some compilation flags> -std=c++11 -DHAVE_INLINE -I<header directory> -o LinRealRecursiveSeq.o LinRealRecursiveSeq.cc
g++ <some linker flags> -lgsl -lgslcblas -o <output name> LinRealRecursiveSeq.o <other object files>
```
Optionally, `-lgslcblas` can be changed to other BLAS library.

# Other

For more information please read `doc/reference.pdf`.
