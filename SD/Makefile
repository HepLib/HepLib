include ../opt.mk

all: SD.o SecDecG.o Helpers.o ErrMin.o cseParser.o CppFormat.o HCubature.o CUBA.o MinUit.o HookeJeeves.o Lib3_HCubature.o KillPowers.o ChengWu.o MatDet.o CIPrepares.o Contours.o Integrates.o AsyMB.o InitEval.o

%.o : %.cpp
	g++ $(opt) -c -o $@ $^

Lib3_HCubature.o : Lib3_HCubature.c
	gcc $(opt) -c -o $@ $^

HCubature.o : HCubature.cpp
	g++ $(opt) -fopenmp -c -o $@ $^

clean:
	rm -f *.o
	rm -f Lib3/*.o




