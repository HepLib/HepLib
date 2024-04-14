/**
 * @file
 * @brief Format for C++ code 
 */
 
#include "SD.h"
#include "cln/cln.h"

namespace HepLib::SD {

    GINAC_IMPLEMENT_PRINT_CONTEXT(ExFormat, print_dflt)

    const ExFormat & ExFormat::operator << (const basic & v) const {
        v.print(*this);
        return *this;
    }
    const ExFormat & ExFormat::operator << (const ex & v) const {
        v.print(*this);
        return *this;
    }
    const ExFormat & ExFormat::operator << (const lst & v) const {
        v.print(*this);
        return *this;
    }
    const ExFormat & ExFormat::operator<<(std::ostream& (*v)(std::ostream&)) const {
        s << v;
        return *this;
    }

    void ExFormat::print_integer(const ExFormat & c, const cln::cl_I & x) {
        c.s << cln::cl_I_to_int(x);
    }

    void ExFormat::print_real(const ExFormat & c, const cln::cl_R & x) {
        if (cln::instanceof(x, cln::cl_I_ring)) {
            print_integer(c, cln::the<cln::cl_I>(x));
        } else if (cln::instanceof(x, cln::cl_RA_ring)) {
            const cln::cl_I numer = cln::numerator(cln::the<cln::cl_RA>(x));
            const cln::cl_I denom = cln::denominator(cln::the<cln::cl_RA>(x));
            if (cln::plusp(x)) {
                c.s << "(";
                print_integer(c, numer);
                c.s << "/" << c.type << "(";
                print_integer(c, denom);
                c.s << "))";
            } else {
                c.s << "(-(";
                print_integer(c, -numer);
                c.s << "/" << c.type << "(";
                print_integer(c, denom);
                c.s << ")))";
            }
        } else {
            c.s << c.type << "(" << c.MQuote;
            cln::cl_print_flags ourflags;
            ourflags.default_float_format = cln::float_format(cln::the<cln::cl_F>(x));
            cln::print_real(c.s, ourflags, x);
            c.s << c.MQuote << ")";
        }
    }

    void ExFormat::print_numeric(const numeric & p, const ExFormat & c, unsigned level) {
        if (p.is_real()) {
            print_real(c, cln::the<cln::cl_R>(p.to_cl_N()));
        } else {
            c.s << "complex<" << c.type << ">(";
            print_real(c, cln::realpart(p.to_cl_N()));
            c.s << ",";
            print_real(c, cln::imagpart(p.to_cl_N()));
            c.s << ")";
        }
    }
    
    void ex_print_power(const power & p, const ExFormat & c, unsigned level) {
        if (p.op(1).is_equal(ex(1)/2)) {
            c.s << "sqrt(";
            p.op(0).print(c);
            c.s << ')';
        } else {
            c.s << "pow(";
            p.op(0).print(c);
            c.s << ',';
            p.op(1).print(c);
            c.s << ')';
        }
    }
    
    ExFormat::ExFormat(ostream &os, const string & s, unsigned opt) : print_dflt(os, opt), suffix(s) {
        set_print_func<power, ExFormat>(ex_print_power);
    }

}



