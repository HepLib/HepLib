#include "SD.h"
#include <cln/cln.h>

namespace HepLib {

GINAC_IMPLEMENT_PRINT_CONTEXT(CppFormat, print_csrc_cl_N)

ex CppFormat::q2ex(qREAL num) {
    char buffer[128];
    quadmath_snprintf(buffer, sizeof buffer, "%.36Qg", num);
    numeric ret(buffer);
    return ret;
}

qREAL CppFormat::ex2q(ex num) {
    qREAL ret = ex_to<numeric>(num).to_double();
    return ret;
}

void CppFormat::QPrint(qREAL num) {
    char buffer[128];
    quadmath_snprintf(buffer, sizeof buffer, "%.36Qg", num);
    cout << buffer << endl;
}

void CppFormat::print_integer(const CppFormat & c, const cln::cl_I & x) {
    const int max_cln_int = 536870911; // 2^29-1
    if (x >= cln::cl_I(-max_cln_int) && x <= cln::cl_I(max_cln_int))
        c.s << cln::cl_I_to_int(x) << ".0" << c.suffix;
    else {
        print_real(c, cln::cl_float(x));
    }
}

void CppFormat::print_real(const CppFormat & c, const cln::cl_R & x) {
    if (cln::instanceof(x, cln::cl_I_ring)) {
        print_integer(c, cln::the<cln::cl_I>(x));
    } else if (cln::instanceof(x, cln::cl_RA_ring)) {
        const cln::cl_I numer = cln::numerator(cln::the<cln::cl_RA>(x));
        const cln::cl_I denom = cln::denominator(cln::the<cln::cl_RA>(x));
        if (cln::plusp(x)) {
            c.s << "(";
            print_integer(c, numer);
        } else {
            c.s << "-(";
            print_integer(c, -numer);
        }
        c.s << "/";
        print_integer(c, denom);
        c.s << ")";
    } else {
        cln::cl_print_flags ourflags;
        ourflags.default_float_format = cln::float_format(cln::the<cln::cl_F>(x));
        cln::print_real(c.s, ourflags, x);
        c.s << c.suffix;
    }
}

void CppFormat::print_numeric(const numeric & p, const CppFormat & c, unsigned level) {
    if (p.is_real()) {
        print_real(c, cln::the<cln::cl_R>(p.to_cl_N()));
    } else {
        if(c.suffix=="L") {
            c.s << "complex<long double>(";
            print_real(c, cln::realpart(p.to_cl_N()));
            c.s << ",";
            print_real(c, cln::imagpart(p.to_cl_N()));
            c.s << ")";
        } else {
            c.s << "(";
            print_real(c, cln::realpart(p.to_cl_N()));
            c.s << "+(";
            print_real(c, cln::imagpart(p.to_cl_N()));
            c.s << ")*1.Qi)";
        }
    }
}

CppFormat::CppFormat(ostream &os, const char* s, unsigned opt) : print_csrc_cl_N(os, opt) {
    set_print_func<numeric, CppFormat>(CppFormat::print_numeric);
    suffix = s;
}

}



