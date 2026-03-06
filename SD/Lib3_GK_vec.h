// modified version from https://github.com/drjerry/quadpackpp
// Adaptive GaussKronrod
#pragma once
#include "Lib3_GaussKronrod.h"

namespace {
    
    class FtnBase {
    public:
        virtual int operator() (unsigned yn, Real *y, Real *e, Real_t x) = 0;
        const unsigned n;
        FtnBase(unsigned yn) : n(yn) { }
    };

    class Function : public FtnBase {
    public:
        typedef std::function<int(unsigned yn, Real *y, Real *e, Real_t x, void *fdata)> f1Type;
        virtual int operator() (unsigned yn, Real *y, Real *e, Real_t x) override { return function_(yn, y, e, x, fdata_); }
        Function(unsigned yn, const f1Type & function, void *fdata)
            : FtnBase(yn), function_(function), fdata_(fdata) { }
        ~Function() { }
        void * fdata() { return fdata_; }
    private:
        f1Type function_;
        void * fdata_;
    };

    class GK : public GaussKronrod {
    public:
        typedef void (*PrintHookerType) (Real*, Real*, size_t *, void *);
        GK(size_t n, size_t m) : GaussKronrod(m), limit(n) { }
        int QAG(FtnBase &f, Real_t a, Real_t b, Real_t epsabs, Real_t epsrel, Real *result, Real *abserr, PrintHookerType PrintHooker);
    private:
        size_t limit;
        Real rescale_error(Real err, Real_t result_abs, Real_t result_asc);
        int qag(FtnBase &f, Real_t a, Real_t b, Real *result, Real *abserr);
    };

    Real GK::rescale_error(Real err, Real_t result_abs, Real_t result_asc) {
        err = abs(err);
        if (result_asc != Real(0) && err != Real(0)) {
            // cast 1.5 as Real number
            Real exponent = Real(3)/Real(2);
            Real scale = pow((200 * err / result_asc), exponent);
            if (scale < Real(1)) err = result_asc * scale;
            else err = result_asc;
        }
        if (result_abs > mpfr::minval() / (50 * mpfr::machine_epsilon())) {
            Real min_err = 50 * mpfr::machine_epsilon() * result_abs;
            if (min_err > err) err = min_err;
        }
        return err;
    }

    int GK::qag(FtnBase &f, Real_t a, Real_t b, Real *result, Real *abserr) {
        const Real center = (a + b) / 2;
        const Real half_length = (b - a) / 2;
        const Real abs_half_length = abs(half_length);
        unsigned yn = f.n;

        std::vector<Real> f_center(yn), ef_center(yn);
        if (f(yn, f_center.data(), ef_center.data(), center)) return FAILURE;

        std::vector<Real> result_gauss(yn), result_kronrod(yn), result_abs(yn), result_asc(yn), mean(yn), err(yn);
        for (int j = 0; j < yn; j++) {
            result_gauss[j] = Real(0);
            result_kronrod[j] = f_center[j] * wgk_[n_ - 1];
            result_abs[j] = abs(result_kronrod[j]);
            result_asc[j] = mean[j] = err[j] = Real(0);
            if (n_ % 2 == 0) result_gauss[j] = f_center[j] * wg_[n_ / 2 - 1];
        }

        std::vector<Real> fv1(yn * n_), e_fv1(yn * n_), fv2(yn * n_), e_fv2(yn * n_);
        if (true) { // Parallel
            std::vector<int> RC1(n_), RC2(n_);
            std::fill(RC1.begin(), RC1.end(), 0);
            std::fill(RC2.begin(), RC2.end(), 0);
            auto prec = mpfr::mpreal::get_default_prec();
            auto rnd = mpfr::mpreal::get_default_rnd();
            #pragma omp parallel for
            for (int jj = 0; jj < 2*n_; jj++) {
                int j = jj / 2, j2 = jj % 2;
                if (mpfr::mpreal::get_default_prec() != prec) 
                    mpfr::mpreal::set_default_prec(prec);
                mpfr::mpreal::set_default_rnd(rnd);
                Real abscissa = half_length * xgk_[j];
                if (j2 == 0) 
                    RC1[j] = f(yn, fv1.data() + yn * j, e_fv1.data() + yn * j, center - abscissa);
                else 
                    RC2[j] = f(yn, fv2.data() + yn * j, e_fv2.data() + yn * j, center + abscissa);
                mpfr_free_cache();
            }
            for (int j = 0; j < n_; j++) 
                if (RC1[j] != 0 || RC2[j] != 0) return FAILURE;

            for (int jj = 0; jj < yn; jj++) { // jj for y index
                for (int j = 0; j < (n_ - 1) / 2; j++) {
                    int jtw = j * 2 + 1;
                    Real fval1 = fv1[yn * jtw + jj];
                    Real fval2 = fv2[yn * jtw + jj];
                    Real fsum = fval1 + fval2;
                    result_gauss[jj] += wg_[j] * fsum;
                    result_kronrod[jj] += wgk_[jtw] * fsum;
                    result_abs[jj] += wgk_[jtw] * (abs(fval1) + abs(fval2));
                }
                for (int j = 0; j < n_ / 2; j++) {
                    int jtwm1 = j * 2;
                    Real fval1 = fv1[yn * jtwm1 + jj];
                    Real fval2 = fv2[yn * jtwm1 + jj];
                    result_kronrod[jj] += wgk_[jtwm1] * (fval1 + fval2);
                    result_abs[jj] += wgk_[jtwm1] * (abs(fval1) + abs(fval2));
                }
            }
        } else { // Non-Parallel
            Real fsum, fval1, fval2;
            for (int j = 0; j < (n_ - 1) / 2; j++) {
                int jtw = j * 2 + 1;        // j=1,2,3 jtw=2,4,6
                Real abscissa = half_length * xgk_[jtw];
                f(yn, fv1.data() + yn * jtw, e_fv1.data() + yn * jtw, center - abscissa);
                f(yn, fv2.data() + yn * jtw, e_fv2.data() + yn * jtw, center + abscissa);
                for (int jj = 0; jj < yn; jj++) { // jj for y index
                    Real_t fval1 = fv1[yn * jtw + jj];
                    Real_t fval2 = fv2[yn * jtw + jj];
                    fsum = fval1 + fval2;
                    result_gauss[jj] += wg_[j] * fsum;
                    result_kronrod[jj] += wgk_[jtw] * fsum;
                    result_abs[jj] += wgk_[jtw] * (abs(fval1) + abs(fval2));
                }
            }
            for (int j = 0; j < n_ / 2; j++) {
                int jtwm1 = j * 2;
                Real abscissa = half_length * xgk_[jtwm1];
                f(yn, fv1.data() + yn * jtwm1, e_fv1.data() + yn * jtwm1, center - abscissa);
                f(yn, fv2.data() + yn * jtwm1, e_fv2.data() + yn * jtwm1, center + abscissa);
                for (int jj = 0; jj < yn; jj++) { // jj for y index
                    Real_t fval1 = fv1[yn * jtwm1 + jj];
                    Real_t fval2 = fv2[yn * jtwm1 + jj];
                    result_kronrod[jj] += wgk_[jtwm1] * (fval1 + fval2);
                    result_abs[jj] += wgk_[jtwm1] * (abs(fval1) + abs(fval2));
                }
            }
        }

        for (int jj = 0; jj < yn; jj++) { // jj for y index
            mean[jj] = result_kronrod[jj] / 2;
            result_asc[jj] = wgk_[n_ - 1] * abs(f_center[jj] - mean[jj]);
            for (int j = 0; j < n_ - 1; j++) {
                result_asc[jj] += wgk_[j] * (abs(fv1[yn * j + jj] - mean[jj]) + abs(fv2[yn * j + jj] - mean[jj]));
            }
            /* scale by the width of the integration region */
            err[jj] = (result_kronrod[jj] - result_gauss[jj]) * half_length;
            result_kronrod[jj] *= half_length;
            result_abs[jj] *= abs_half_length;
            result_asc[jj] *= abs_half_length;
            result[jj] = result_kronrod[jj];
            abserr[jj] = rescale_error(err[jj], result_abs[jj], result_asc[jj]);
        }
        
        return SUCCESS;
    }

    inline Real max(unsigned yn, Real *x) {
        Real max_value = x[0];
        for (int j = 1; j < yn; j++) if (max_value < x[j]) max_value = x[j];
        return max_value;
    }

    int GK::QAG(FtnBase &f, Real_t a, Real_t b, Real_t epsabs, Real_t epsrel, Real *result, Real *abserr, PrintHookerType PrintHooker) {
        size_t size = 0;
        size_t nrmax = 0;
        size_t i_work = 0;
        size_t maximum_level = 0;

        unsigned yn = f.n;
        size_t n = limit;

        std::vector<Real> alist(n), blist(n);
        std::vector<std::vector<Real>> rlist(n, std::vector<Real>(yn));
        std::vector<std::vector<Real>> elist(n, std::vector<Real>(yn));
        std::vector<size_t> order(n), level(n);
        
        std::vector<Real> area(yn), errsum(yn), result0(yn), abserr0(yn), tolerance(yn);
        size_t iteration = 0;

        {
            size = 0;
            nrmax = 0;
            i_work = 0;
            alist[0] = a;
            blist[0] = b;
            for (int j = 0; j < yn; j++) rlist[0][j] = elist[0][j] = Real(0);
            order[0] = 0;
            level[0] = 0;
            maximum_level = 0;
        }

        for (int j = 0; j < yn; j++) result[j] = abserr[j] = Real(0);

        if (this->qag(f, a, b, result0.data(), abserr0.data())) return FAILURE;

        {
            size = 1;
            for (int j = 0; j < yn; j++) {
                rlist[0][j] = result0[j];
                elist[0][j] = abserr0[j];
            }
        }

        /* Test on accuracy */
        for (int j = 0; j < yn; j++) tolerance[j] = max(epsabs, epsrel * abs(result0[j]));
        bool ok = true;

        for (int j = 0; j < yn; j++) {
            if (abserr0[j] > tolerance[j]) { ok = false; break; }
        }

        if (ok) {
            for (int j = 0; j < yn; j++) {
                result[j] = result0[j];
                abserr[j] = abserr0[j];
            }
            size_t kk = 1;
            auto fdata = ((Function&)f).fdata();
            if (PrintHooker) PrintHooker(result, abserr, &kk, fdata);
            return SUCCESS;
        }

        for (int j = 0; j < yn; j++) {
            area[j] = result0[j];
            errsum[j] = abserr0[j];
        }
        iteration = 1;

        do {
            Real a1, b1, a2, b2;
            Real a_i, b_i, r_i[yn], e_i[yn];
            std::vector<Real> area1(yn), area2(yn), area12(yn);
            std::vector<Real> error1(yn), error2(yn), error12(yn);
            
            for (int j = 0; j < yn; j++) {
                area1[j] = area2[j] = area12[j] = error1[j] = error2[j] = error12[j] = Real(0);
            }

            /* Bisect the subinterval with the largest error estimate */
            {
                a_i = alist[i_work];
                b_i = blist[i_work];
                for (int j = 0; j < yn; j++) {
                    r_i[j] = rlist[i_work][j];
                    e_i[j] = elist[i_work][j];
                }
            }

            a1 = a_i;
            b1 = (a_i + b_i) / Real(2);
            a2 = b1;
            b2 = b_i;

            if (this->qag(f, a1, b1, area1.data(), error1.data())) return FAILURE;
            if (this->qag(f, a2, b2, area2.data(), error2.data())) return FAILURE;

            for (int j = 0; j < yn; j++) {
                area12[j] = area1[j] + area2[j];
                error12[j] = error1[j] + error2[j];
                errsum[j] += (error12[j] - e_i[j]);
                area[j] += area12[j] - r_i[j];
                tolerance[j] = max(epsabs, epsrel * abs(area[j]));
            }

            {
                const size_t i_max = i_work;
                const size_t i_new = size;
                const size_t new_level = level[i_max] + 1;

                /* append the newly-created intervals to the list */
                if (max(yn, error2.data()) > max(yn, error1.data())) {
                    alist[i_max] = a2;        /* blist[maxerr] is already == b2 */
                    for (int j = 0; j < yn; j++) {
                        rlist[i_max][j] = area2[j];
                        elist[i_max][j] = error2[j];
                    }
                    level[i_max] = new_level;

                    alist[i_new] = a1;
                    blist[i_new] = b1;
                    for (int j = 0; j < yn; j++) {
                        rlist[i_new][j] = area1[j];
                        elist[i_new][j] = error1[j];
                    }
                    level[i_new] = new_level;
                } else {
                    blist[i_max] = b1;        /* alist[maxerr] is already == a1 */
                    for (int j = 0; j < yn; j++) {
                        rlist[i_max][j] = area1[j];
                        elist[i_max][j] = error1[j];
                    }
                    level[i_max] = new_level;

                    alist[i_new] = a2;
                    blist[i_new] = b2;
                    for (int j = 0; j < yn; j++) {
                        rlist[i_new][j] = area2[j];
                        elist[i_new][j] = error2[j];
                    }
                    level[i_new] = new_level;
                }

                size++;
                if (new_level > maximum_level) maximum_level = new_level;

                {
                    const size_t last = size - 1;
                    Real errmax, errmin;
                    int i, k, top;
                    size_t i_nrmax = nrmax;
                    size_t i_maxerr = order[i_nrmax];

                    /* Check whether the list contains more than two error estimates */
                    if (last < 2) {
                        order[0] = 0;
                        order[1] = 1;
                        i_work = i_maxerr;
                    } else {
                        errmax = max(yn, elist[i_maxerr].data());
                        while (i_nrmax > 0 && errmax > max(yn, elist[order[i_nrmax - 1]].data())) {
                            order[i_nrmax] = order[i_nrmax - 1];
                            i_nrmax--;
                        }
                        if (last < (limit / 2 + 2)) top = last;
                        else top = limit - last + 1;
                        i = i_nrmax + 1;
                        while (i < top && errmax < max(yn, elist[order[i]].data())) {
                            order[i - 1] = order[i];
                            i++;
                        }
                        order[i - 1] = i_maxerr;
                        errmin = max(yn, elist[last].data());
                        k = top - 1;
                        while (k > i - 2 && errmin >= max(yn, elist[order[k]].data())) {
                            order[k + 1] = order[k];
                            k--;
                        }
                        order[k + 1] = last;

                        /* Set i_max and e_max */
                        i_maxerr = order[i_nrmax];
                        i_work = i_maxerr;
                        nrmax = i_nrmax;
                    }
                }
            }

            {
                a_i = alist[i_work];
                b_i = blist[i_work];
                for (int j = 0; j < yn; j++) {
                    r_i[j] = rlist[i_work][j];
                    e_i[j] = elist[i_work][j];
                }
            }

            if (PrintHooker) {
                size_t kk = iteration;
                auto fdata = ((Function&)f).fdata();
                PrintHooker(area.data(), errsum.data(), &kk, fdata);
                if (kk != iteration) {
                    for (int j = 0; j < yn; j++) {
                        result[j] = area[j];
                        abserr[j] = errsum[j];
                    }
                    return SUCCESS;
                }
            }

            ok = true;
            for (int j = 0; j < yn; j++) {
                if (errsum[j] > tolerance[j]) { ok = false; break; }
            }
            iteration++;
        } while (iteration < limit && !ok);

        for (int j = 0; j < yn; j++) {
            result[j] = Real(0);
            for (size_t k = 0; k < size; k++) result[j] += rlist[k][j];
            abserr[j] = errsum[j];
        }
        return SUCCESS;
    }
}
