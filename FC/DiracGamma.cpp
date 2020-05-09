/**
 * @file
 * @brief Functions for DiracGamma
 * @author F. Feng
 * @version 1.0.0
 * @date 2020-04-21
 */
 
#include "FC.h"

namespace HepLib::FC {

    static ex expl_TR_diff(const ex & arg, const symbol & s) {
        return TR(arg.diff(s));
    }
    
    //-----------------------------------------------------------
    // DiracGamma Class
    //-----------------------------------------------------------
    GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(DiracGamma, basic,
        print_func<print_dflt>(&DiracGamma::print).
        print_func<FormFormat>(&DiracGamma::form_print).
        print_func<FCFormat>(&DiracGamma::fc_print)
    )
    
    DEFAULT_CTOR(DiracGamma)
    GINAC_BIND_UNARCHIVER(DiracGamma);
    IMPLEMENT_HAS(DiracGamma)
    IMPLEMENT_ALL(DiracGamma)
        
    return_type_t DiracGamma::return_type_tinfo() const {
        return make_return_type_t<clifford>(rl);
        //return make_return_type_t<DiracGamma>(rl);
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
    
    ex & DiracGamma::let_op(size_t i) {
        static ex ex_rl = numeric(rl);
        ensure_if_modifiable();
        if(i==0) return pi;
        else return ex_rl;
    }
    
    ex DiracGamma::eval() const {
        if(flags & status_flags::evaluated) return *this;
        else if(is_zero(pi) || is_zero(pi-1) || is_zero(pi-5) || is_zero(pi-6) || is_zero(pi-7)) return this->hold();
        else if(!is_a<Vector>(pi) && !is_a<Index>(pi)) return GAS(pi,rl);
        else return this->hold();
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
        c << "(" << pi << ")";
    }
    
    void DiracGamma::archive(archive_node & n) const {
        inherited::archive(n);
        n.add_ex("pi", pi);
        n.add_unsigned("rl", rl);
    }
    
    void DiracGamma::read_archive(const archive_node& n, lst& sym_lst) {
        inherited::read_archive(n, sym_lst);
        n.find_unsigned("rl", rl);
        n.find_ex("pi", pi, sym_lst);
    }
    
    ex DiracGamma::derivative(const symbol & s) const {
        return 0;
    }
    
    ex DiracGamma::conjugate() const {
        if(is_a<Index>(pi) || is_a<Vector>(pi)) return *this;
        else if(is_zero(pi-5)) return -1*DiracGamma(5, rl);
        else if(is_zero(pi-6)) return DiracGamma(7, rl);
        else if(is_zero(pi-7)) return DiracGamma(6, rl);
        else if(is_zero(pi-1)) return DiracGamma(1, rl);
        throw Error("invalid DiracGamma Found.");
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
            c << "DiracTrace(" << arg << ")";
        }
        
        ex tr_conj(const ex & e) {
            return TR(e.conjugate());
        }
    }
    
    REGISTER_FUNCTION(TR, do_not_evalf_params().
        conjugate_func(tr_conj).
        print_func<FormFormat>(&TR_form_print).
        print_func<FCFormat>(&TR_fc_print).
        set_return_type(return_types::commutative).
        expl_derivative_func(expl_TR_diff)
    );
    
    ex GAS(ex expr, unsigned rl) {
        if(is_zero(expr-1)) return DiracGamma(1,rl);
        else if(is_zero(expr-5)) return DiracGamma(5,rl);
        else if(is_zero(expr-6)) return DiracGamma(6,rl);
        else if(is_zero(expr-7)) return DiracGamma(7,rl);
        
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
                    g = DiracGamma(ex_to<Vector>(ii),rl);
                } else if(is_a<Index>(ii)) {
                    if(is_a<DiracGamma>(g)) throw Error("Something Wrong with GAS");
                    g = DiracGamma(ex_to<Index>(ii),rl);
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

