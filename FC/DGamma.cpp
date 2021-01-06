/**
 * @file
 * @brief Functions for DGamma
 */
 
#include "FC.h"

namespace HepLib::FC {

    static ex expl_TR_diff2(const ex & arg, const symbol & s) {
        auto ret = arg.diff(s);
        if(!is_a<add>(ret)) ret = lst{ret};
        ex res = 0;
        for(auto item : ret) {
            if(!is_a<mul>(item)) res += TR(item);
            else {
                ex c=1; ex v=1;
                for(auto it : item) {
                    if(!Index::has(it) && !DGamma::has(it)) c *= it;
                    else v *= it;
                }
                res += c * TR(v);
            }
        }
        return res;
    }
    
    static ex expl_TR_diff(const ex & arg, const symbol & s) {
        auto ret = arg.diff(s);
        auto nd = numer_denom(ret.normal());
        auto num = collect_common_factors(nd.op(0));
        if(!is_a<mul>(num)) return TR(num)/nd.op(1);
        ex c=1;
        ex v=1;
        for(auto it : num) {
            if(!Index::has(it) && !DGamma::has(it)) c *= it;
            else v *= it;
        }
        return c*TR(v)/nd.op(1);
    }
    
    //-----------------------------------------------------------
    // DGamma Class
    //-----------------------------------------------------------
    GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(DGamma, basic,
        print_func<print_dflt>(&DGamma::print).
        print_func<FormFormat>(&DGamma::form_print).
        print_func<FCFormat>(&DGamma::fc_print)
    )
    
    DEFAULT_CTOR(DGamma)
    GINAC_BIND_UNARCHIVER(DGamma);
    IMPLEMENT_HAS(DGamma)
    IMPLEMENT_ALL(DGamma)
        
    return_type_t DGamma::return_type_tinfo() const {
        return make_return_type_t<clifford>(rl);
        //return make_return_type_t<DGamma>(rl);
    }
    
    bool DGamma::match_same_type(const basic & other) const {
        const DGamma &o = static_cast<const DGamma &>(other);
        return rl == o.rl;
    }
    
    unsigned DGamma::get_rl() {
        return rl;
    }
    
    int DGamma::compare_same_type(const basic &other) const {
        const DGamma &o = static_cast<const DGamma &>(other);
        if (rl != o.rl) return rl < o.rl ? -1 : 1;
        return pi.compare(o.pi);
    }
    
    DGamma::DGamma(int int_1567, unsigned _rl) : pi(int_1567), rl(_rl) { }
    DGamma::DGamma(const Vector &p, unsigned _rl) : pi(p), rl(_rl) { }
    DGamma::DGamma(const Index &i, unsigned _rl) : pi(i), rl(_rl) { }
    DGamma::DGamma(const DGamma &g, unsigned _rl) : pi(g.pi), rl(_rl) { }
    
    size_t DGamma::nops() const { return 2; }
    ex DGamma::op(size_t i) const {
        if(i==0) return pi;
        else if(i==1) return rl;
        return 0;
    }
    
    ex & DGamma::let_op(size_t i) {
        static ex ex_rl = numeric(rl);
        ensure_if_modifiable();
        if(i==0) return pi;
        else return ex_rl;
    }
    
    ex DGamma::eval() const {
        if(flags & status_flags::evaluated) return *this;
        else if(is_zero(pi) || is_zero(pi-1) || is_zero(pi-5) || is_zero(pi-6) || is_zero(pi-7)) return this->hold();
        else if(!is_a<Vector>(pi) && !is_a<Index>(pi)) return GAS(pi,rl);
        else return this->hold();
    }
    
    void DGamma::print(const print_dflt &c, unsigned level) const {
        if(is_zero(pi-1)) {
            c.s << "𝕚";
            return;
        }
        c.s << "(" << "𝛾";
        if(is_a<numeric>(pi)) c.s << pi;
        else c.s << "." << pi;
        c.s << ")";
    }
    
    void DGamma::form_print(const FormFormat &c, unsigned level) const {
        c << "g_(" << rl;
        if(!is_zero(pi-1)) {
            c << "," << pi;
            if(is_a<numeric>(pi)) c << "_"; // 5_, 6_, 7_ in form
        }
        c << ")";
    }
    
    void DGamma::fc_print(const FCFormat &c, unsigned level) const {
        if(is_zero(pi-1)) {
            c << "1";
            return;
        }
        if(is_a<Vector>(pi)) c << "GSD";
        else c << "GAD";
        c << "(" << pi << ")";
    }
    
    void DGamma::archive(archive_node & n) const {
        inherited::archive(n);
        n.add_ex("pi", pi);
        n.add_unsigned("rl", rl);
    }
    
    void DGamma::read_archive(const archive_node& n, lst& sym_lst) {
        inherited::read_archive(n, sym_lst);
        n.find_unsigned("rl", rl);
        n.find_ex("pi", pi, sym_lst);
    }
    
    ex DGamma::derivative(const symbol & s) const {
        return 0;
    }
    
    ex DGamma::conjugate() const {
        if(is_a<Index>(pi) || is_a<Vector>(pi)) return *this;
        else if(is_zero(pi-5)) return -1*DGamma(5, rl);
        else if(is_zero(pi-6)) return DGamma(7, rl);
        else if(is_zero(pi-7)) return DGamma(6, rl);
        else if(is_zero(pi-1)) return DGamma(1, rl);
        throw Error("invalid Dirac Gamma Found.");
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
        
        void TTR_form_print(const ex &arg, const print_context &c0) {
            auto c = static_cast<const FormFormat &>(c0);
            if(!is_a<lst>(arg)) c << "TTR(" << arg << ")";
            else {
                bool first = true;
                for(auto item : arg) {
                    if(first) { first=false; c << "TTR(" << item; }
                    else c << "," << item;
                }
                c << ")";
            }
        }
        void TTR_fc_print(const ex &arg, const print_context &c0) {
            auto c = static_cast<const FCFormat &>(c0);
            if(!is_a<lst>(arg)) c << "SUNTrace(SUNT(" << arg << "))";
            else {
                bool first = true;
                for(auto item : arg) {
                    if(first) { first=false; c << "SUNTrace(SUNT(" << item; }
                    else c << "," << item;
                }
                c << "))";
            }
        }
        ex ttr_conj(const ex & e) {
            lst argv;
            if(!is_a<lst>(e)) return TTR(e);
            else argv = ex_to<lst>(e);
            lst as;
            for(auto it=argv.rbegin(); it!=argv.rend(); ++it) as.append(*it);
            return TTR(as);
        }
    }
    
    REGISTER_FUNCTION(TR, do_not_evalf_params().
        conjugate_func(tr_conj).
        print_func<FormFormat>(&TR_form_print).
        print_func<FCFormat>(&TR_fc_print).
        set_return_type(return_types::commutative).
        expl_derivative_func(expl_TR_diff)
    );
    
    REGISTER_FUNCTION(TTR, do_not_evalf_params().
        conjugate_func(ttr_conj).
        print_func<FormFormat>(&TTR_form_print).
        print_func<FCFormat>(&TTR_fc_print)
    );
    
    REGISTER_FUNCTION(HF, do_not_evalf_params());
    
    /**
     * @brief function similar to GAD/GSD in FeynClac
     * @param expr momentum/index or 1,5,6,7
     * @param rl the represent line number
     * @return expanded/translasted to Dirac Gamma objects
     */
    ex GAS(const ex &expr, unsigned rl) {
        if(is_zero(expr-1)) return DGamma(1,rl);
        else if(is_zero(expr-5)) return DGamma(5,rl);
        else if(is_zero(expr-6)) return DGamma(6,rl);
        else if(is_zero(expr-7)) return DGamma(7,rl);
        
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
                    if(is_a<DGamma>(g)) throw Error("Something Wrong with GAS");
                    g = DGamma(ex_to<Vector>(ii),rl);
                } else if(is_a<Index>(ii)) {
                    if(is_a<DGamma>(g)) throw Error("Something Wrong with GAS");
                    g = DGamma(ex_to<Index>(ii),rl);
                } else {
                    c *= ii;
                }
            }
            if(!is_a<DGamma>(g)) throw Error("Something Wrong with GAS");
            res += c * g;
        }
        return res;
    }
    
    


}
