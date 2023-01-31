#pragma once
#include "mpreal.h"
#include <functional>

namespace {
    typedef mpfr::mpreal Real;
    class FtnBase {
    public:
        virtual Real operator() (const Real & x) =0;
    };
    
    class Function : public FtnBase {
    public:
        typedef std::function<Real(const Real & x, void * fdata)> f1Type;
        virtual Real operator() (const Real & x) override {
            return function_(x,fdata_);
        }
        Function(const f1Type & function, void * fdata) : function_(function), fdata_(fdata) { }
        ~Function() { }
    private:
        f1Type function_;
        void * fdata_;
    };

    class GaussKronrod {
    private:
        size_t m_;  // Gauss-Legendre degree
        size_t n_;  // size of Gauss-Kronrod arrays
        Real*  xgk_;  // Gauss-Kronrod abscissae
        Real*  wg_;   // Gauss-Legendre weights
        Real*  wgk_;  // Gauss-Kronrod weights

        Real *coefs;  // Chebyshev coefficients
        Real *zeros;  // zeros of Legendre polynomial
        Real *fv1, *fv2;  // scratch space for error estimator

        Real rescale_error(Real err, const Real result_abs, const Real result_asc);

        void legendre_zeros();
        void chebyshev_coefs();
        void gauss_kronrod_abscissae();
        void gauss_kronrod_weights();

        Real legendre_err(int deg, Real x, Real& err);
        Real legendre_deriv(int deg, Real x);
        Real chebyshev_series(Real x, Real& err);
        Real chebyshev_series_deriv(Real x);

    public:
        GaussKronrod(size_t m = 10);
        ~GaussKronrod();

        void qk(FtnBase& f, Real a, Real b, Real& result, Real& abserr, Real& resabs, Real& resasc, bool parallel=false);

        size_t size() { return n_; };

        Real xgk(int k) {
            return (0 <= k && k < n_) ? xgk_[k] : Real(0);
        }

        Real wgk(int k) {
            return (0 <= k && k < n_) ? wgk_[k] : Real(0);
        }

        Real wg(int k) {
            return (0 <= k && k < n_/2) ? wg_[k] : Real(0);
        }
    };

    GaussKronrod::GaussKronrod(size_t m) {
        m_ = m;
        n_ = m_ + 1;
        xgk_ = new Real[n_];
        wg_  = new Real[n_ / 2];
        wgk_ = new Real[n_];
        coefs = new Real[n_ + 1];
        zeros = new Real[m_ + 2];
        fv1 = new Real[n_];
        fv2 = new Real[n_];

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
        delete[] fv1;
        delete[] fv2;
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
            for (int i = 1; i < k; ++i)
                alpha[k] -= f[i] * alpha[k-i];
        }

        for (int k = 0; k <= ell; ++k) {
            coefs[m_ + 1 - 2*k] = alpha[k];
            if (m_  >= 2*k)
                coefs[m_ - 2*k] = Real(0);
        }

        delete[] alpha;
        delete[] f;
    }

    void GaussKronrod::gauss_kronrod_weights() {
        Real err;
        /* Gauss-Legendre weights:
         */
        for (int k = 0; k < n_ / 2; ++k) {
            Real x = xgk_[2*k + 1];
            wg_[k] = (Real(-2) / ((m_ + 1) * legendre_deriv(m_, x) * legendre_err(m_+1, x, err)));
        }

        /* The ratio of leading coefficients of P_n and T_{n+1} is computed
         * from the recursive formulae for the respective polynomials.
         */
        Real F_m = Real(2) / Real(2*m_ + 1);
        for (int k = 1; k <= m_; ++k) F_m *= (Real(2*k) / Real(2*k - 1));

        /* Gauss-Kronrod weights:
         */
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
            while (abs(E) > epsilon &&
                     abs(delta) > mpfr::machine_epsilon())
            {
                delta = E / chebyshev_series_deriv(x_k);
                x_k -= delta;
                E = chebyshev_series(x_k, epsilon);
            }
            xgk_[2*k] = x_k;
            // copy adjacent Legendre-zero into the array:
            if (2*k+1 < n_)
                xgk_[2*k+1] = zeros[m_-k];
        }
    }

    Real GaussKronrod::legendre_err(int n, Real x, Real& err) {
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
        for (int k = 1; k < n; ++k)
        {
            P2 = ((2*k + 1) * x * P1 - k * P0) / (k + 1);
            err = ((2*k + 1) * abs(x) * E1 + k * E0) / (2*(k + 1));
            P0 = P1; P1 = P2;
            E0 = E1; E1 = err;
        }
        return P2;
    }

    Real GaussKronrod::legendre_deriv(int n, Real x) {
        if (n == 0)
            return Real(0);
        else if (n == 1)
            return Real(1);

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

    Real GaussKronrod::chebyshev_series(Real x, Real& err) {
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

    Real GaussKronrod::chebyshev_series_deriv(Real x) {
        Real d1(0), d2(0);
        Real y2 = 2 * x; // linear term for Clenshaw recursion

        for (int k = n_; k >= 2; --k) {
            Real temp = d1;
            d1 = y2 * d1 - d2 + k * coefs[k];
          d2 = temp;
        }

        return y2 * d1 - d2 + coefs[1];
    }

    Real GaussKronrod::rescale_error (Real err, const Real result_abs, const Real result_asc) {
        err = abs(err);

        if (result_asc != Real(0) && err != Real(0)) {
            // cast 1.5 as Real number
            Real exponent = Real(3)/Real(2);
            Real scale = pow((200 * err / result_asc), exponent);

            if (scale < Real(1)) {
                err = result_asc * scale ;
            } else {
                err = result_asc ;
            }
        }

        if (result_abs > mpfr::minval() / (50 * mpfr::machine_epsilon())) {
          Real min_err = 50 * mpfr::machine_epsilon() * result_abs ;

          if (min_err > err) {
                err = min_err ;
            }
        }

        return err ;
    }

    void GaussKronrod::qk(FtnBase& f, Real a, Real b,
                                         Real& result, Real& abserr,
                                         Real& resabs, Real& resasc,
                                         bool parallel) {
        const Real center = (a + b) / 2;
        const Real half_length = (b - a) / 2;
        const Real abs_half_length = abs(half_length);
        // const Real f_center = f.function(center, f.params);
        const Real f_center = f(center);

        Real result_gauss = Real(0);
        Real result_kronrod = f_center * wgk_[n_ - 1];
        Real result_abs = abs(result_kronrod);
        Real result_asc = Real(0);
        Real mean = Real(0), err = Real(0);

        int j;

         if (n_ % 2 == 0) {
             result_gauss = f_center * wg_[n_/2 - 1];
         }

if(parallel) {
    auto prec = mpfr::mpreal::get_default_prec();
    auto rnd = mpfr::mpreal::get_default_rnd();
    #pragma omp parallel for schedule(dynamic,1)
    for(int jj = 0; jj < 2*n_ ; jj++) {
        int j = jj/2, j2 = jj % 2;
        mpfr::mpreal::set_default_prec(prec);
        mpfr::mpreal::set_default_rnd(rnd);
        Real abscissa = half_length * xgk_[j];
        if(j2==0) fv1[j] = f(center - abscissa);
        else fv2[j] = f(center + abscissa);
        mpfr_free_cache();
    }
    
    for(int j = 0; j < (n_ - 1) / 2; j++) {
        int jtw = j * 2 + 1;
        Real fval1 = fv1[jtw];
        Real fval2 = fv2[jtw];
        Real fsum = fval1 + fval2;
        result_gauss += wg_[j] * fsum;
        result_kronrod += wgk_[jtw] * fsum;
        result_abs += wgk_[jtw] * (abs(fval1) + abs(fval2));
    }

    for (int j = 0; j < n_ / 2; j++) {
        int jtwm1 = j * 2;
        Real fval1 = fv1[jtwm1];
        Real fval2 = fv2[jtwm1];
        result_kronrod += wgk_[jtwm1] * (fval1 + fval2);
        result_abs += wgk_[jtwm1] * (abs(fval1) + abs(fval2));
	}
} else {
	for (j = 0; j < (n_ - 1) / 2; j++) {
      int jtw = j * 2 + 1;        /* j=1,2,3 jtw=2,4,6 */
      Real abscissa = half_length * xgk_[jtw];
//      Real fval1 = f.function( center - abscissa , f.params);
//      Real fval2 = f.function( center + abscissa , f.params);
		Real fval1 = f(center - abscissa);
		Real fval2 = f(center + abscissa);
      Real fsum = fval1 + fval2;
      fv1[jtw] = fval1;
      fv2[jtw] = fval2;
      result_gauss += wg_[j] * fsum;
      result_kronrod += wgk_[jtw] * fsum;
      result_abs += wgk_[jtw] * (abs(fval1) + abs(fval2));
	}

        for (j = 0; j < n_ / 2; j++) {
          int jtwm1 = j * 2;
          Real abscissa = half_length * xgk_[jtwm1];
    //      Real fval1 = f.function( center - abscissa , f.params);
    //      Real fval2 = f.function( center + abscissa , f.params);
            Real fval1 = f(center - abscissa);
            Real fval2 = f(center + abscissa);
            fv1[jtwm1] = fval1;
          fv2[jtwm1] = fval2;
          result_kronrod += wgk_[jtwm1] * (fval1 + fval2);
          result_abs += wgk_[jtwm1] * (abs(fval1) + abs(fval2));
        };
}

        mean = result_kronrod / 2;

        result_asc = wgk_[n_ - 1] * abs(f_center - mean);

        for (j = 0; j < n_ - 1; j++) {
          result_asc += wgk_[j] * (abs(fv1[j] - mean) + abs(fv2[j] - mean));
        }

        /* scale by the width of the integration region */

        err = (result_kronrod - result_gauss) * half_length;

        result_kronrod *= half_length;
        result_abs *= abs_half_length;
        result_asc *= abs_half_length;

        result = result_kronrod;
        resabs = result_abs;
        resasc = result_asc;
        abserr = rescale_error (err, result_abs, result_asc);
    }
}
