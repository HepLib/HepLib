include ../opt.mk

all: SD.o SecDecG.o Helpers.o ErrMin.o cseParser.o CppFormat.o HCubature.o CUBA.o MinUit.o HookeJeeves.o Lib3/HCubature.o KillPowers.o ChengWu.o MatDet.o

%.o : %.cpp
	g++ $(opt) -c -o $@ $^

Lib3/HCubature.o : Lib3/HCubature.c
	gcc $(opt) -c -o $@ $^

HCubature.o : HCubature.cpp
	g++ $(opt) -I Lib3 -fopenmp -c -o $@ $^

CUBA.o : CUBA.cpp
	g++ $(opt) -I Lib3 -fopenmp -c -o $@ $^

MinUit.o : MinUit.cpp
	g++ $(opt) -fopenmp -c -o $@ $^

clean:
	rm -f *.o
	rm -f Lib3/*.o




