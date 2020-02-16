/* Adaptive multidimensional integration of a vector of integrands.
 *
 * Copyright (c) 2005-2013 Steven G. Johnson
 *
 * Portions (see comments) based on HIntLib (also distributed under
 * the GNU GPL, v2 or later), copyright (c) 2002-2005 Rudolf Schuerer.
 *     (http://www.cosy.sbg.ac.at/~rschuer/hintlib/)
 *
 * Portions (see comments) based on GNU GSL (also distributed under
 * the GNU GPL, v2 or later), copyright (c) 1996-2000 Brian Gough.
 *     (http://www.gnu.org/software/gsl/)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#pragma once

//typedef long double cubareal;
typedef __float128 REAL;

#include <stdlib.h> /* for size_t */
#include <quadmath.h>

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/* USAGE: Call hcubature or pcubature with your function as described
          in the README file. */

/* a vector integrand - evaluates the function at the given point x
   (an array of length ndim) and returns the result in fval (an array
   of length fdim).   The void* parameter is there in case you have
   to pass any additional data through to your function (it corresponds
   to the fdata parameter you pass to cubature).  Return 0 on
   success or nonzero to terminate the integration. */
typedef int (*integrand) (unsigned ndim, const REAL *x, void *, unsigned fdim, REAL *fval);

/* a vector integrand of a vector of npt points: x[i*ndim + j] is the
   j-th coordinate of the i-th point, and the k-th function evaluation
   for the i-th point is returned in fval[i*fdim + k].  Return 0 on success
   or nonzero to terminate the integration. */
typedef int (*integrand_v) (unsigned ndim, long long npt, const REAL *x, void *, unsigned fdim, REAL *fval);

/* Integrate the function f from xmin[dim] to xmax[dim], with at most
   maxEval function evaluations (0 for no limit), until the given
   absolute or relative error is achieved.  val returns the integral,
   and err returns the estimate for the absolute error in val; both
   of these are arrays of length fdim, the dimension of the vector
   integrand f(x). The return value of the function is 0 on success
   and non-zero if there  was an error. */

/* adapative integration by partitioning the integration domain ("h-adaptive")
   and using the same fixed-degree quadrature in each subdomain, recursively,
   until convergence is achieved. */
int hcubature(unsigned fdim, integrand f, void *fdata, unsigned dim, const REAL *xmin, const REAL *xmax, long long minEval, long long runEval, long long maxEval, REAL reqAbsError, REAL reqRelError, REAL *val, REAL *err);

typedef void (* PrintHookerType) (REAL*, REAL*, long long int*, void *);

/* as hcubature, but vectorized integrand */
int hcubature_v(unsigned fdim, integrand_v f, void *fdata, unsigned dim, const REAL *xmin, const REAL *xmax, long long minEval, long long runEval, long long maxEval, REAL reqAbsError, REAL reqRelError, REAL *val, REAL *err, PrintHookerType PrintHooker);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

