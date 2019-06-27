include opt.mk

default: install

dep: 
	make -C ExGiNaC
	make -C SD
	make -C IBP

libHepLib.so: ExGiNaC/*.o SD/*.o SD/Lib3/*.o IBP/*.o
	g++ $(opt) -shared -lgomp -lquadmath -ldl -lqhullstatic -lMinuit2 -lginac -lcln -o $@ $^

install: libHepLib.so
	cp -r include/* $(prefix)/include/
	cp libHepLib.so $(prefix)/lib/

clean:
	make -C ExGiNaC clean
	make -C SD clean
	make -C IBP clean
	rm -f libHepLib.so



