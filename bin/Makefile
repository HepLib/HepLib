include ../opt.mk

all: garview

garview : garview.cpp
	g++ $(opt) -o $@ $^ -L .. -lHepLib -lginac -lcln -ldl -lqhullstatic


clean:
	rm -f garview




