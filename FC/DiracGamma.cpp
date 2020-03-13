#include "FC.h"

namespace HepLib::FC {
    
    //-----------------------------------------------------------
    // DiracGamma Class
    //-----------------------------------------------------------
    GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(DiracGamma, basic,
        print_func<print_dflt>(&DiracGamma::print).
        print_func<FormFormat>(&DiracGamma::form_print).
        print_func<FCFormat>(&DiracGamma::fc_print)
    )
        
    return_type_t DiracGamma::return_type_tinfo() const {
        return make_return_type_t<DiracGamma>(rl);
    }
    
    bool DiracGamma::match_same_type(const basic & other) const {
        const DiracGamma &o = static_cast<const DiracGamma &>(other);
        return rl == o.rl;
    }
    
    unsigned DiracGamma::get_rl() {
        return rl;
    }
    
    int DiracGamma::compare_same_type(const basic &other) const {
        const DiracGamma &o = static_cast<const DiracGamma &>(other);
        if (rl != o.rl) return rl < o.rl ? -1 : 1;
        return pi.compare(o.pi);
    }
    
    DiracGamma::DiracGamma(int int_1567, unsigned _rl) : pi(int_1567), rl(_rl) { }
    DiracGamma::DiracGamma(const Vector &p, unsigned _rl) : pi(p), rl(_rl) { }
    DiracGamma::DiracGamma(const Index &i, unsigned _rl) : pi(i), rl(_rl) { }
    DiracGamma::DiracGamma(const DiracGamma &g, unsigned _rl) : pi(g.pi), rl(_rl) { }
    
    size_t DiracGamma::nops() const { return 2; }
    ex DiracGamma::op(size_t i) const {
        if(i==0) return pi;
        else if(i==1) return rl;
        return 0;
    }
    
    void DiracGamma::print(const print_dflt &c, unsigned level) const {
        if(is_zero(pi-1)) {
            c.s << "𝕚";
            return;
        }
        c.s << "(" << "𝛾";
        if(is_a<numeric>(pi)) c.s << pi;
        else c.s << "." << pi;
        c.s << ")";
    }
    
    void DiracGamma::form_print(const FormFormat &c, unsigned level) const {
        c << "g_(" << rl;
        if(!is_zero(pi-1)) {
            c << "," << pi;
            if(is_a<numeric>(pi)) c << "_"; // 5_, 6_, 7_ in form
        }
        c << ")";
    }
    
    void DiracGamma::fc_print(const FCFormat &c, unsigned level) const {
        if(is_zero(pi-1)) {
            c << "1";
            return;
        }
        if(is_a<Vector>(pi)) c << "GSD";
        else c << "GAD";
        c << "[" << pi << "]";
    }
    
    //-----------------------------------------------------------
    // TR/GAS functions
    //-----------------------------------------------------------
    namespace {
        void TR_form_print(const ex &arg, const print_context &c0) {
            auto c = static_cast<const FormFormat &>(c0);
            c << "(" << arg << ")";
        }
        void TR_fc_print(const ex &arg, const print_context &c0) {
            auto c = static_cast<const FCFormat &>(c0);
            c << "DiracTrace[" << arg << "]";
        }
    }
    
    REGISTER_FUNCTION(TR, do_not_evalf_params().
        print_func<FormFormat>(&TR_form_print).
        print_func<FCFormat>(&TR_fc_print).
        set_return_type(return_types::commutative)
    );
    
    ex GAS(ex expr) {
        ex tmp = expand(expr);
        lst ex_lst;
        if(is_a<add>(tmp)) {
            for(auto item : tmp) ex_lst.append(item);
        } else ex_lst.append(tmp);
        ex res = 0;
        for(auto item : ex_lst) {
            lst mul_lst;
            if(is_a<mul>(item)) {
                for(auto ii : item) mul_lst.append(ii);
            } else mul_lst.append(item);
            ex c=1, g=1;
            for(auto ii : mul_lst) {
                if(is_a<Vector>(ii)) {
                    if(is_a<DiracGamma>(g)) throw Error("Something Wrong with GAS");
                    g = DiracGamma(ex_to<Vector>(ii));
                } else if(is_a<Index>(ii)) {
                    if(is_a<DiracGamma>(g)) throw Error("Something Wrong with GAS");
                    g = DiracGamma(ex_to<Index>(ii));
                } else {
                    c *= ii;
                }
            }
            if(!is_a<DiracGamma>(g)) throw Error("Something Wrong with GAS");
            res += c * g;
        }
        return res;
    }
    
    


}

