/**
 * @file
 * @brief Format for C++ code 
 * @author F. Feng
 * @version 1.0.0
 * @date 2020-04-21
 */
 
#include "SD.h"
#include <cln/cln.h>

namespace HepLib::SD {

GINAC_IMPLEMENT_PRINT_CONTEXT(CppFormat, print_csrc_cl_N)

ex CppFormat::q2ex(qREAL num) {
    char buffer[128];
    quadmath_snprintf(buffer, sizeof buffer, "%.36QG", num);
    numeric ret(buffer);
    return ret;
}

qREAL CppFormat::ex2q(ex num) {
    ostringstream nss;
    auto oDigits = Digits;
    Digits = 40;
    nss << num.evalf() << endl;
    qREAL ret = strtoflt128(nss.str().c_str(), NULL);
    Digits = oDigits;
    return ret;
}

void CppFormat::QPrint(qREAL num) {
    char buffer[128];
    quadmath_snprintf(buffer, sizeof buffer, "%.36QG", num);
    cout << buffer << endl;
}

void CppFormat::print_integer(const CppFormat & c, const cln::cl_I & x) {
    const int max_cln_int = 536870911; // 2^29-1
    if (x >= cln::cl_I(-max_cln_int) && x <= cln::cl_I(max_cln_int)) {
        if(c.suffix=="MP") c.s << "mpREAL(" << c.MQuote;
        c.s << cln::cl_I_to_int(x);
        if(c.suffix=="MP") c.s << c.MQuote << ")";
        else c.s << ".0" << c.suffix;
    } else {
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
        if(c.suffix=="MP") c.s << "mpREAL(" << c.MQuote;
        cln::cl_print_flags ourflags;
        ourflags.default_float_format = cln::float_format(cln::the<cln::cl_F>(x));
        cln::print_real(c.s, ourflags, x);
        if(c.suffix=="MP") c.s << c.MQuote << ")";
        else c.s << c.suffix;
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
        } else if(c.suffix=="Q") {
            c.s << "(";
            print_real(c, cln::realpart(p.to_cl_N()));
            c.s << "+(";
            print_real(c, cln::imagpart(p.to_cl_N()));
            c.s << ")*1.Qi)";
        } else if(c.suffix=="MP") {
            c.s << "complex<mpREAL>(";
            print_real(c, cln::realpart(p.to_cl_N()));
            c.s << ",";
            print_real(c, cln::imagpart(p.to_cl_N()));
            c.s << ")";
        } else {
            throw Error("CppFormat: suffix is Wrong.");
        }
    }
}

CppFormat::CppFormat(ostream &os, const string & s, unsigned opt) : print_csrc_cl_N(os, opt), suffix(s) { }

}



