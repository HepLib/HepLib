// modified version from https://github.com/drjerry/quadpackpp
// GaussKronrod Rules
#pragma once
#include "mpreal.h"
#include <functional>

/* error return codes */
#define SUCCESS 0
#define FAILURE 1

namespace {
    typedef mpfr::mpreal Real;
    typedef const mpfr::mpreal & Real_t;
    

    class GaussKronrod {
    protected:
        size_t m_;  // Gauss-Legendre degree
        size_t n_;  // size of Gauss-Kronrod arrays
        Real*  xgk_;  // Gauss-Kronrod abscissae
        Real*  wg_;   // Gauss-Legendre weights
        Real*  wgk_;  // Gauss-Kronrod weights

    private:
        Real *coefs;  // Chebyshev coefficients
        Real *zeros;  // zeros of Legendre polynomial

        void legendre_zeros();
        void chebyshev_coefs();
        void gauss_kronrod_abscissae();
        void gauss_kronrod_weights();

        Real legendre_err(int deg, Real_t x, Real& err);
        Real legendre_deriv(int deg, Real_t x);
        Real chebyshev_series(Real_t x, Real& err);
        Real chebyshev_series_deriv(Real_t x);
        
    public:
        GaussKronrod(size_t m = 10);
        ~GaussKronrod();
        size_t size() { return n_; };
        Real xgk(int k) { return (0 <= k && k < n_) ? xgk_[k] : Real(0); }
        Real wgk(int k) { return (0 <= k && k < n_) ? wgk_[k] : Real(0); }
        Real wg(int k) { return (0 <= k && k < n_/2) ? wg_[k] : Real(0); }
    };

    GaussKronrod::GaussKronrod(size_t m) {
        m_ = m;
        n_ = m_ + 1;
        xgk_ = new Real[n_];
        wg_  = new Real[n_ / 2];
        wgk_ = new Real[n_];
        coefs = new Real[n_ + 1];
        zeros = new Real[m_ + 2];
        
        legendre_zeros();
        chebyshev_coefs();
        gauss_kronrod_abscissae();
        gauss_kronrod_weights();
    }

    GaussKronrod::~GaussKronrod() {
        delete[] xgk_;
        delete[] wgk_;
        delete[] wg_;
        delete[] coefs;
        delete[] zeros;
    }

    void GaussKronrod::legendre_zeros() {
        Real* temp = new Real[m_+1];
        zeros[0] = Real(-1);
        zeros[1] = Real(1);
        Real delta, epsilon;

        for (int k = 1; k <= m_; ++k) {
            // loop to locate zeros of P_k interlacing z_0,...,z_k
            for (int j = 0; j < k; ++j) {
                // Newton's method for P_k :
                // initialize solver at midpoint of (z_j, z_{j+1})
                delta = 1;
                Real x_j = (zeros[j] + zeros[j+1]) / 2;
                Real P_k = legendre_err(k, x_j, epsilon);
                while (abs(P_k) > epsilon && abs(delta) > mpfr::machine_epsilon()) {
                    delta = P_k / legendre_deriv(k, x_j);
                    x_j -= delta;
                    P_k = legendre_err(k, x_j, epsilon);
                }
                temp[j] = x_j;
            }

            // copy roots tmp_0,...,tmp_{k-1} to z_1,...z_k:
            zeros[k+1] = zeros[k];
            for (int j = 0; j < k; ++j) zeros[j+1] = temp[j];
        }
        delete[] temp;
    }

    void GaussKronrod::chebyshev_coefs() {
        size_t ell = (m_ + 1)/2;
        Real* alpha = new Real[ell+1];
        Real* f = new Real[ell+1];

        /* Care must be exercised in initalizing the constants in the definitions.
         * Compilers interpret expressions like "(2*k + 1.0)/(k + 2.0)" as floating
         * point precision, before casting to Real.
         */
        f[1] = Real(m_+1)/Real(2*m_ + 3);
        alpha[0] = Real(1); // coefficient of T_{m+1}
        alpha[1] = -f[1];

        for (int k = 2; k <= ell; ++k) {
            f[k] = f[k-1] * (2*k - 1) * (m_ + k) / (k * (2*m_ + 2*k + 1));
            alpha[k] = -f[k];
            for (int i = 1; i < k; ++i) alpha[k] -= f[i] * alpha[k-i];
        }

        for (int k = 0; k <= ell; ++k) {
            coefs[m_ + 1 - 2*k] = alpha[k];
            if (m_  >= 2*k) coefs[m_ - 2*k] = Real(0);
        }

        delete[] alpha;
        delete[] f;
    }

    void GaussKronrod::gauss_kronrod_weights() {
        Real err;
        // Gauss-Legendre weights:
        for (int k = 0; k < n_ / 2; ++k) {
            Real x = xgk_[2*k + 1];
            wg_[k] = (Real(-2) / ((m_ + 1) * legendre_deriv(m_, x) * legendre_err(m_+1, x, err)));
        }

        /* The ratio of leading coefficients of P_n and T_{n+1} is computed
         * from the recursive formulae for the respective polynomials.
         */
        Real F_m = Real(2) / Real(2*m_ + 1);
        for (int k = 1; k <= m_; ++k) F_m *= (Real(2*k) / Real(2*k - 1));

        // Gauss-Kronrod weights:
        for (size_t k = 0; k < n_; ++k) {
            Real x = xgk_[k];
            if (k % 2 == 0) {
                wgk_[k] = F_m / (legendre_err(m_, x, err) * chebyshev_series_deriv(x));
            } else {
                wgk_[k] = (wg_[k/2] + F_m / (legendre_deriv(m_, x) * chebyshev_series(x, err)));
            }
        }
    }

    void GaussKronrod::gauss_kronrod_abscissae() {
        Real delta, epsilon;

        for (int k = 0; 2*k < n_; ++k) {
            delta = 1;
            // Newton's method for E_{n+1} :
            Real x_k = (zeros[m_-k] + zeros[m_+1-k])/Real(2);
            Real E = chebyshev_series(x_k, epsilon);
            while (abs(E) > epsilon && abs(delta) > mpfr::machine_epsilon()) {
                delta = E / chebyshev_series_deriv(x_k);
                x_k -= delta;
                E = chebyshev_series(x_k, epsilon);
            }
            xgk_[2*k] = x_k;
            // copy adjacent Legendre-zero into the array:
            if (2*k+1 < n_) xgk_[2*k+1] = zeros[m_-k];
        }
    }

    Real GaussKronrod::legendre_err(int n, Real_t x, Real& err) {
        if (n == 0) {
            err = Real(0);
            return Real(1);
        }
        else if (n == 1) {
            err = Real(0);
            return x;
        }

        Real P0 = Real(1), P1 = x, P2;
        Real E0 = mpfr::machine_epsilon();
        Real E1 = abs(x) * mpfr::machine_epsilon();
        for (int k = 1; k < n; ++k) {
            P2 = ((2*k + 1) * x * P1 - k * P0) / (k + 1);
            err = ((2*k + 1) * abs(x) * E1 + k * E0) / (2*(k + 1));
            P0 = P1; P1 = P2;
            E0 = E1; E1 = err;
        }
        return P2;
    }

    Real GaussKronrod::legendre_deriv(int n, Real_t x) {
        if (n == 0) return Real(0);
        else if (n == 1) return Real(1);

        Real P0 = Real(1), P1 = x, P2;
        Real dP0 = Real(0), dP1 = Real(1), dP2;
        for (int k = 1; k < n; ++k) {
            P2 = ((2*k + 1) * x * P1 - k * P0) / (k + Real(1));
            dP2 = (2*k + 1) * P1 + dP0;
            P0 = P1; P1 = P2;
            dP0 = dP1; dP1 = dP2;
        }
        return dP2;
    }

    Real GaussKronrod::chebyshev_series(Real_t x, Real& err) {
        Real d1(0), d2(0);
        Real absc = abs(coefs[0]); // final term for truncation error
        Real y2 = 2 * x; // linear term for Clenshaw recursion

        for (int k = n_; k >= 1; --k) {
            Real temp = d1;
            d1 = y2 * d1 - d2 + coefs[k];
            d2 = temp;
            absc += abs(coefs[k]);
        }

        err = absc * mpfr::machine_epsilon();
        return x * d1 - d2 + coefs[0]/2;
    }

    Real GaussKronrod::chebyshev_series_deriv(Real_t x) {
        Real d1(0), d2(0);
        Real y2 = 2 * x; // linear term for Clenshaw recursion

        for (int k = n_; k >= 2; --k) {
            Real temp = d1;
            d1 = y2 * d1 - d2 + k * coefs[k];
            d2 = temp;
        }

        return y2 * d1 - d2 + coefs[1];
    }
    
}
