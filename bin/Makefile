include ../opt.mk

all: garview FF FFz

% : %.cpp
	g++ $(opt) -o $@ $^ -L .. -lHepLib -lginac -lcln -ldl -lqhullstatic

clean:
	rm -f garview FF FFz




