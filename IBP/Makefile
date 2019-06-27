include ../opt.mk

all: Fermat.o IBP.o

Fermat.o : Fermat.cpp
	g++ $(opt) -c -o $@ $^

IBP.o : IBP.cpp
	g++ $(opt) -c -o $@ $^ 

clean:
	rm -f *.o




