/**
 * @file
 * @brief Format for C++ code 
 */
 
#include "SD.h"
#include "cln/cln.h"

namespace HepLib::SD {

    GINAC_IMPLEMENT_PRINT_CONTEXT(CppFormat, print_csrc_cl_N)

    const CppFormat & CppFormat::operator << (const basic & v) const {
        v.print(*this);
        return *this;
    }
    const CppFormat & CppFormat::operator << (const ex & v) const {
        v.print(*this);
        return *this;
    }
    const CppFormat & CppFormat::operator << (const lst & v) const {
        v.print(*this);
        return *this;
    }
    const CppFormat & CppFormat::operator<<(std::ostream& (*v)(std::ostream&)) const {
        s << v;
        return *this;
    }

    void CppFormat::print_integer(const CppFormat & c, const cln::cl_I & x) {
        const int max_cln_int = 536870911; // 2^29-1
        if (x >= cln::cl_I(-max_cln_int) && x <= cln::cl_I(max_cln_int)) {
            if(c.suffix=="MP") c.s << "mpREAL(";
            c.s << cln::cl_I_to_int(x);
            if(c.suffix=="MP") c.s << ")";
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
                c.s << "/";
                print_integer(c, denom);
                c.s << ")";
            } else {
                c.s << "(-(";
                print_integer(c, -numer);
                c.s << "/";
                print_integer(c, denom);
                c.s << "))";
            }
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
            if(c.suffix=="D") {
                c.s << "complex<double>(";
                print_real(c, cln::realpart(p.to_cl_N()));
                c.s << ",";
                print_real(c, cln::imagpart(p.to_cl_N()));
                c.s << ")";
            } else if(c.suffix=="L") {
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



