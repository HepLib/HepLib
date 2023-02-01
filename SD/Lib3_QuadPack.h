// modified version from https://github.com/drjerry/quadpackpp
#pragma once
#include "Lib3_GaussKronrod.h"

namespace {

enum {
	GSL_SUCCESS  = 0,
	GSL_EMAXITER = 11,  /* exceeded max number of iterations */
	GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
	GSL_EROUND   = 18,  /* failed because of roundoff error */
};

class Workspace : public GaussKronrod {
public:
    typedef void (*PrintHookerType) (Real*, Real*, size_t *, void *);
	Workspace();
	Workspace(size_t limit);
	Workspace(size_t limit, size_t m);
	~Workspace();
	int qag(FtnBase & f, Real_t a, Real_t b, Real_t epsabs, Real_t epsrel, Real& result, Real& abserr, PrintHookerType PrintHooker=NULL);

private:
	void	allocate(size_t limit);

	/* data from gsl_integration_workspace */
	size_t limit;
	size_t size;
	size_t nrmax;
	size_t i_work;
	size_t maximum_level;
	Real   *alist;
	Real   *blist;
	Real   *rlist;
	Real   *elist;
	size_t *order;
	size_t *level;

	/* auxillary functions for adaptive quadrature */
	void  append_interval(Real_t a1, Real_t b1, Real_t area1, Real_t error1);
	void  initialise(Real_t a, Real_t b);
	void  set_initial_result(Real_t result, Real_t error);
	void  qpsrt();
	void  sort_results();
	void  retrieve(Real& a, Real& b, Real& r, Real& e);
	Real  sum_results();
	void  update(Real_t a1, Real_t b1, Real_t area1, Real_t error1, Real_t a2, Real_t b2, Real_t area2, Real_t error2);

};

void Workspace::allocate(size_t n) {
	alist = new Real[n];
	blist = new Real[n];
	rlist = new Real[n];
	elist = new Real[n];
	order = new size_t[n];
	level = new size_t[n];
	size = 0;
	limit = n;
	maximum_level = 0;
}

Workspace::Workspace() : GaussKronrod() { allocate(1); }

Workspace::Workspace(size_t n) : GaussKronrod() {
	if (n == 0) n = 1;	// ensure that workspace has positive size
	allocate(n);
}

Workspace::Workspace(size_t n, size_t m) : GaussKronrod(m) {
	if (n == 0) n = 1;	// ensure that workspace has positive size
	allocate(n);
}

Workspace::~Workspace() {
	delete[] alist;
	delete[] blist;
	delete[] rlist;
	delete[] elist;
	delete[] order;
	delete[] level;
}

void Workspace::append_interval (Real_t a1, Real_t b1, Real_t area1, Real_t error1) {
	alist[size] = a1;
	blist[size] = b1;
	rlist[size] = area1;
	elist[size] = error1;
	order[size] = size;
	level[size] = 0;
	size++;
}

void Workspace::initialise (Real_t a, Real_t b) {
	size = 0;
	nrmax = 0;
	i_work = 0;
	alist[0] = a;
	blist[0] = b;
	rlist[0] = Real(0);
	elist[0] = Real(0);
	order[0] = 0;
	level[0] = 0;
	maximum_level = 0;
}

void Workspace::sort_results () {
	size_t i;

	for (i = 0; i < size; i++) {
      size_t i1 = order[i];
      Real e1 = elist[i1];
      size_t i_max = i1;
      size_t j;

      for (j = i + 1; j < size; j++) {
			size_t i2 = order[j];
			Real e2 = elist[i2];

			if (e2 >= e1) {
				i_max = i2;
				e1 = e2;
			}
       }

      if (i_max != i1) {
			order[i] = order[i_max];
			order[i_max] = i1;
		}
	}

	i_work = order[0] ;
}

void Workspace::qpsrt () {
	const size_t last = size - 1;

	Real errmax ;
	Real errmin ;
	int i, k, top;

	size_t i_nrmax = nrmax;
	size_t i_maxerr = order[i_nrmax] ;

	/* Check whether the list contains more than two error estimates */

	if (last < 2) {
      order[0] = 0 ;
      order[1] = 1 ;
      i_work = i_maxerr ;
      return ;
	}

	errmax = elist[i_maxerr] ;

	/* This part of the routine is only executed if, due to a difficult
	 integrand, subdivision increased the error estimate. In the normal
	 case the insert procedure should start after the nrmax-th largest
	 error estimate. */

	while (i_nrmax > 0 && errmax > elist[order[i_nrmax - 1]]) {
      order[i_nrmax] = order[i_nrmax - 1] ;
      i_nrmax-- ;
	}

	/* Compute the number of elements in the list to be maintained in
	 descending order. This number depends on the number of
	 subdivisions still allowed. */

	if(last < (limit/2 + 2)) top = last ;
	else top = limit - last + 1;

	/* Insert errmax by traversing the list top-down, starting
	 comparison from the element elist(order(i_nrmax+1)). */

	i = i_nrmax + 1 ;

	/* The order of the tests in the following line is important to
	 prevent a segmentation fault */

	while (i < top && errmax < elist[order[i]]) {
      order[i-1] = order[i] ;
      i++ ;
	}

	order[i-1] = i_maxerr ;

	/* Insert errmin by traversing the list bottom-up */

	errmin = elist[last] ;

	k = top - 1 ;

	while (k > i - 2 && errmin >= elist[order[k]]) {
      order[k+1] = order[k] ;
      k-- ;
	}

	order[k+1] = last ;

	/* Set i_max and e_max */

	i_maxerr = order[i_nrmax] ;

	i_work = i_maxerr ;
	nrmax = i_nrmax ;
}

void Workspace::set_initial_result (Real_t result, Real_t error) {
	size = 1;
	rlist[0] = result;
	elist[0] = error;
}

void Workspace::retrieve (Real& a, Real& b, Real& r, Real& e) {
	a = alist[i_work] ;
	b = blist[i_work] ;
	r = rlist[i_work] ;
	e = elist[i_work] ;
}

Real Workspace::sum_results () {
	size_t k;
	Real result_sum = Real(0);
	for (k = 0; k < size; k++) result_sum += rlist[k];
	return result_sum;
}

void Workspace::update (Real_t a1, Real_t b1, Real_t area1, Real_t error1, Real_t a2, Real_t b2, Real_t area2, Real_t error2) {
	const size_t i_max = i_work ;
	const size_t i_new = size ;

	const size_t new_level = level[i_max] + 1;

	/* append the newly-created intervals to the list */

	if (error2 > error1) {
      alist[i_max] = a2;        /* blist[maxerr] is already == b2 */
      rlist[i_max] = area2;
      elist[i_max] = error2;
      level[i_max] = new_level;

      alist[i_new] = a1;
      blist[i_new] = b1;
      rlist[i_new] = area1;
      elist[i_new] = error1;
      level[i_new] = new_level;
	} else {
      blist[i_max] = b1;        /* alist[maxerr] is already == a1 */
      rlist[i_max] = area1;
      elist[i_max] = error1;
      level[i_max] = new_level;

      alist[i_new] = a2;
      blist[i_new] = b2;
      rlist[i_new] = area2;
      elist[i_new] = error2;
      level[i_new] = new_level;
	}

	size++;

	if (new_level > maximum_level) maximum_level = new_level;

	qpsrt () ;
}

/* -------------------------------------------------------------------------- */
int Workspace::qag(FtnBase& f, Real_t a, Real_t b, Real_t epsabs, Real_t epsrel, Real& result, Real& abserr, PrintHookerType PrintHooker) {
	Real area, errsum;
	Real result0, abserr0, resabs0, resasc0;
	Real tolerance;
	size_t iteration = 0;

	initialise (a, b);
	result = Real(0);
	abserr = Real(0);

	this->qk(f, a, b, result0, abserr0, resabs0, resasc0);
	set_initial_result (result0, abserr0);

	/* Test on accuracy */
	tolerance = max (epsabs, epsrel * abs(result0));

	if ((abserr0 <= tolerance && abserr0 != resasc0) || abserr0 == Real(0)) {
      result = result0;
      abserr = abserr0;

if(PrintHooker) {
    std::cout << result << std::endl;
    std::cout << abserr << std::endl << std::endl;
}

      return GSL_SUCCESS;
	}

	area = result0;
	errsum = abserr0;
	iteration = 1;

	do {
      Real a1, b1, a2, b2;
      Real a_i, b_i, r_i, e_i;
      Real area1 = Real(0), area2 = Real(0), area12 = Real(0);
      Real error1 = Real(0), error2 = Real(0), error12 = Real(0);
      Real resasc1, resasc2;
      Real resabs1, resabs2;

      /* Bisect the subinterval with the largest error estimate */
      
      retrieve(a_i, b_i, r_i, e_i);

      a1 = a_i;
      b1 = (a_i + b_i) / Real(2);
      a2 = b1;
      b2 = b_i;

      this->qk(f, a1, b1, area1, error1, resabs1, resasc1);
      this->qk(f, a2, b2, area2, error2, resabs2, resasc2);

      area12 = area1 + area2;
      error12 = error1 + error2;

      errsum += (error12 - e_i);
      area += area12 - r_i;

      tolerance = max (epsabs, epsrel * abs(area));

      update (a1, b1, area1, error1, a2, b2, area2, error2);
      retrieve(a_i, b_i, r_i, e_i);
      iteration++;
    
if(PrintHooker) {
    result = sum_results();
    abserr = errsum;
    std::cout << result << std::endl;
    std::cout << abserr << std::endl << std::endl;
}

	} while (iteration < limit && errsum > tolerance);

	result = sum_results();
	abserr = errsum;

	if (errsum <= tolerance) return GSL_SUCCESS;
	else return GSL_EMAXITER;
}

}
