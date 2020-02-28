include ../opt.mk

ifdef WSTP_PATH
all: ExGiNaC.o WSKernel.o
else
all: ExGiNaC.o
endif

ExGiNaC.o : ExGiNaC.cpp
	g++ $(opt) -c -o $@ $^

WSKernel.o : WSKernel.cpp
	g++ $(opt) -c -o $@ $^
	
clean:
	rm -f *.o




