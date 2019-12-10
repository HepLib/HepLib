include ../opt.mk

all: SD.o SecDecG.o Helpers.o ErrMin.o cseParser.o CppFormat.o HCubature.o CUBA.o MinUit.o HookeJeeves.o Lib3/HCubature.o

SD.o : SD.cpp
	g++ $(opt) -c -o $@ $^

SecDecG.o : SecDecG.cpp
	g++ $(opt) -c -o $@ $^

Helpers.o : Helpers.cpp
	g++  $(opt) -fPIC -c -o $@ $^

ErrMin.o : ErrMin.cpp
	g++  $(opt) -fPIC -c -o $@ $^

cseParser.o : cseParser.cpp
	g++  $(opt) -fPIC -c -o $@ $^


CppFormat.o : CppFormat.cpp
	g++ $(opt) -c -o $@ $^

HCubature.o : HCubature.cpp
	g++ $(opt) -I Lib3 -fopenmp -c -o $@ $^

CUBA.o : CUBA.cpp
	g++ $(opt) -I Lib3 -fopenmp -c -o $@ $^

MinUit.o : MinUit.cpp
	g++ $(opt) -fopenmp -c -o $@ $^

HookeJeeves.o : HookeJeeves.cpp
	gcc $(opt) -c -o $@ $^

Lib3/HCubature.o : Lib3/HCubature.c
	gcc $(opt) -c -o $@ $^

clean:
	rm -f *.o
	rm -f Lib3/*.o




