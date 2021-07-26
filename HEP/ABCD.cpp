/**
 * @file
 * @brief Functions to perform Tensor Index Reduction
 */
 
#include "HEP.h"
#include <cmath>

namespace HepLib {

    /**
     * @brief scalar integral A0, devided by (2pi)^(4-2ep)
     * @param m2 mass squared
     * @param n exponent for the propagator l^2-m^2
     * @return A0
     */
    ex A0(const ex m2, int n, ex d) {
        return I*pow(ex(-1),n)*pow(2,-d)*pow(Pi, -d/2)*pow(m2, -1 + d/2)*pow(tgamma(ex(n)), -1)*tgamma(-d/2 + n);
    }

}

