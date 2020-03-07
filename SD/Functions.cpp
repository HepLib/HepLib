#include "SD.h"

namespace HepLib {

    namespace {
        /*-----------------------------------------------------*/
        // Functions used in GiNaC
        /*-----------------------------------------------------*/
        static void print_VEO(const ex & ex1_in, const ex & ex2_in, const print_context & c) {
            ex ex1 = ex1_in, ex2 = ex2_in;
            if(abs(ex1) < numeric("1E-30")) ex1 = 0;
            if(abs(ex2) < numeric("1E-30")) ex2 = 0;
            if(ex1==0 || ex2==0) {
                char bf1[128], bf2[128];
                quadmath_snprintf(bf1, sizeof bf1, "%.10QG", CppFormat::ex2q(ex1_in));
                quadmath_snprintf(bf2, sizeof bf2, "%.10QG", CppFormat::ex2q(ex2_in));
                c.s << "(" << bf1 << " +- " << bf2 << ")";
                return;
            }
            int digits = 30;
            if(!ex2.is_zero()) {
                auto ratio = ex_to<numeric>(abs(ex1/ex2));
                digits = floorq(logq(ratio.to_double())/logq(10.Q)) + 2;
                digits = digits > 1 ? digits : 1;
                digits = digits > 30 ? 30 : digits;
            }
            auto oDigits = Digits;
            ostringstream oss;
            try {
                Digits = digits;
                oss << "(" << ex1.evalf();
                Digits = 2;
                oss << " +- " << ex2.evalf() << ")";
                Digits = oDigits;
                c.s << oss.str();
            } catch(...) {
                Digits = oDigits;
                try {
                    char bf1[128], bf2[128];
                    quadmath_snprintf(bf1, sizeof bf1, "%.10QG", CppFormat::ex2q(ex1_in));
                    quadmath_snprintf(bf2, sizeof bf2, "%.10QG", CppFormat::ex2q(ex2_in));
                    c.s << "(" << bf1 << " +- " << bf2 << ")";
                } catch(...) {
                    c.s << Color_Error << "[-NaN-]" << RESET;
                }
            }
        }
        
        static ex NoDiff_1P(const ex & x, unsigned diff_param) {return 0;}
        static ex NoDiff_2P(const ex & x, const ex & y, unsigned diff_param) {return 0;}
        static ex VE_Conjugate(const ex & x, const ex & y) { return VE(x,y).hold(); }
        static ex Diff_ID(const ex & x, unsigned diff_param) {return 1;}
        
    }

    REGISTER_FUNCTION(CV, do_not_evalf_params()) // for use symbol
    REGISTER_FUNCTION(fabs, dummy())
    REGISTER_FUNCTION(PL, dummy())
    REGISTER_FUNCTION(FTX, derivative_func(NoDiff_2P))
    REGISTER_FUNCTION(CT, derivative_func(Diff_ID))
    REGISTER_FUNCTION(VE, conjugate_func(VE_Conjugate))
    REGISTER_FUNCTION(VEO, print_func<print_dflt>(print_VEO))
    REGISTER_FUNCTION(epsID, do_not_evalf_params().derivative_func(NoDiff_1P))
    REGISTER_FUNCTION(WF, dummy())

}
