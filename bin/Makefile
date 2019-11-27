include ../opt.mk

all: garview UF gFF gFFz qFF qFFz

% : %.cpp
	g++ $(opt) -o $@ $^ -lHepLib -lginac -lcln -ldl -lqhullstatic -lcubaq -lquadmath

clean:
	rm -f garview UF gFF gFFz qFF qFFz




