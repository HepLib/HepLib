include opt.mk

default: install

dep: 
	make -j -C ExGiNaC
	make -j -C SD
	make -j -C MB

lib: libHepLib.so
	
libHepLib.so: ExGiNaC/*.o SD/*.o
	g++ $(opt) -shared -lgomp -lquadmath -ldl -lqhullstatic -lMinuit2 -lginac -lcln -lcubaq -lmpfr -lgmp -o $@ $^

install: libHepLib.so
	make -j -C bin
	cp -r include/* $(prefix)/include/
	cp libHepLib.so $(prefix)/lib/
	cp bin/garview $(prefix)/bin/
	cp bin/UF $(prefix)/bin/
	cp bin/gFF $(prefix)/bin/
	cp bin/gFFz $(prefix)/bin/
	cp bin/qFF $(prefix)/bin/
	cp bin/qFFz $(prefix)/bin/

clean:
	make -j -C ExGiNaC clean
	make -j -C SD clean
	make -j -C MB clean
	make -j -C IBP clean
	rm -f libHepLib.so
	make -j -C bin clean



