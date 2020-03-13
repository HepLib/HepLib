#include "FC.h"

namespace HepLib::FC {

    //-----------------------------------------------------------
    // Pair Class
    //-----------------------------------------------------------
    GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(Pair, basic,
        print_func<print_dflt>(&Pair::print).
        print_func<FormFormat>(&Pair::form_print).
        print_func<FCFormat>(&Pair::fc_print)
    )

    Pair::Pair(const Vector &p1, const Vector &p2) : spL(p1), spR(p2) { }
    Pair::Pair(const Index &i1, const Index &i2) : spL(i1), spR(i2) { }
    Pair::Pair(const Vector &p, const Index &i) : spL(p), spR(i) { }

    int Pair::compare_same_type(const basic &other) const {
        const Pair &o = static_cast<const Pair &>(other);
        int c1 = spL.compare(o.spL);
        if(c1!=0) return c1;
        int c2 = spR.compare(o.spR);
        return c2;
    }
    
    void Pair::print(const print_dflt &c, unsigned level) const {
        c.s << spL << "." << spR;
    }
    
    void Pair::form_print(const FormFormat &c, unsigned level) const {
        if(is_a<Vector>(spL) && is_a<Vector>(spR)) c << spL << "." << spR;
        else if(is_a<Vector>(spL) && is_a<Index>(spR)) c << spL << "(" << spR << ")";
        else if(is_a<Index>(spL) && is_a<Index>(spR)) c << "d_(" << spL << "," << spR << ")";
    }
    
    void Pair::fc_print(const FCFormat &c, unsigned level) const {
        if(is_a<Vector>(spL) && is_a<Vector>(spR)) c << "SPD[" << spL << "," << spR << "]";
        else if(is_a<Vector>(spL) && is_a<Index>(spR)) c << "FVD[" << spL << "," << spR << "]";
        else if(is_a<Index>(spL) && is_a<Index>(spR)) c << "MTD[" << spL << "," << spR << "]";
    }
    
    size_t Pair::nops() const { return 2; }
    ex Pair::op(size_t i) const {
        if(i==0) return spL;
        else if(i==1) return spR;
        return 0;
    }
    
    //-----------------------------------------------------------
    // SP function - ScalarProduct
    //-----------------------------------------------------------
    ex SP(ex a, ex b) {
        if(is_a<Vector>(a) && is_a<Vector>(b)) return Pair(ex_to<Vector>(a), ex_to<Vector>(b));
        else if(is_a<Vector>(a) && is_a<Index>(b)) return Pair(ex_to<Vector>(a), ex_to<Index>(b));
        else if(is_a<Index>(a) && is_a<Vector>(b)) return Pair(ex_to<Vector>(b), ex_to<Index>(a));
        else if(is_a<Index>(a) && is_a<Index>(b)) return Pair(ex_to<Index>(a), ex_to<Index>(b));
        
        lst alst, blst;
        a = a.expand();
        if(is_a<add>(a)) {
            for(auto item : a) alst.append(item);
        } else alst.append(a);
        for(int i=0; i<alst.nops(); i++) {
            if(!is_a<mul>(alst.op(i))) alst.let_op(i) = lst{alst.op(i)};
            ex c=1;
            ex v=1;
            for(auto ii : alst.op(i)) {
                if(is_a<Vector>(ii) || is_a<Index>(ii)) {
                    if(!is_zero(v-1)) throw Error("Error Found in SP");
                    v = ii;
                } else c *= ii;
            }
            if(is_zero(v-1)) throw Error("Error Found in SP");
            alst.let_op(i) = lst{c,v};
        }
        
        b = b.expand();
        if(is_a<add>(b)) {
            for(auto item : b) blst.append(item);
        } else blst.append(b);
        for(int i=0; i<blst.nops(); i++) {
            if(!is_a<mul>(blst.op(i))) blst.let_op(i) = lst{blst.op(i)};
            ex c=1;
            ex v=1;
            for(auto ii : blst.op(i)) {
                if(is_a<Vector>(ii) || is_a<Index>(ii)) {
                    if(!is_zero(v-1)) throw Error("Error Found in SP");
                    v = ii;
                } else c *= ii;
            }
            if(is_zero(v-1)) throw Error("Error Found in SP");
            blst.let_op(i) = lst{c,v};
        }
        
        ex res = 0;
        for(auto ai : alst) {
            for(auto bi : blst) res += ai.op(0) * bi.op(0) * SP(ai.op(1), bi.op(1));
        }
        return res;
    }

}

