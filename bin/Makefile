include ../opt.mk

all: garview UF

% : %.cpp
	g++ $(opt) -o $@ $^ -lHepLib -lginac -lcln -ldl -lqhullstatic -lcubaq -lquadmath

clean:
	rm -f garview FF FFz UF




