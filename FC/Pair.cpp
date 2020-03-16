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
    
    DEFAULT_CTOR(Pair)
    GINAC_BIND_UNARCHIVER(Pair);
    IMPLEMENT_HAS(Pair)

    Pair::Pair(const Vector &p1, const Vector &p2) : lr{p1,p2} { }
    Pair::Pair(const Index &i1, const Index &i2) : lr{i1,i2} { }
    Pair::Pair(const Vector &p, const Index &i) : lr{p,i} { }

    int Pair::compare_same_type(const basic &other) const {
        const Pair &o = static_cast<const Pair &>(other);
        int c1 = lr[0].compare(o.lr[0]);
        if(c1!=0) return c1;
        int c2 = lr[1].compare(o.lr[1]);
        return c2;
    }
    
    void Pair::print(const print_dflt &c, unsigned level) const {
        c.s << lr[0] << "." << lr[1];
    }
    
    void Pair::form_print(const FormFormat &c, unsigned level) const {
        if(is_a<Vector>(lr[0]) && is_a<Vector>(lr[1])) c << lr[0] << "." << lr[1];
        else if(is_a<Vector>(lr[0]) && is_a<Index>(lr[1])) c << lr[0] << "(" << lr[1] << ")";
        else if(is_a<Index>(lr[0]) && is_a<Index>(lr[1])) c << "d_(" << lr[0] << "," << lr[1] << ")";
    }
    
    void Pair::fc_print(const FCFormat &c, unsigned level) const {
        if(is_a<Vector>(lr[0]) && is_a<Vector>(lr[1])) c << "SPD[" << lr[0] << "," << lr[1] << "]";
        else if(is_a<Vector>(lr[0]) && is_a<Index>(lr[1])) c << "FVD[" << lr[0] << "," << lr[1] << "]";
        else if(is_a<Index>(lr[0]) && is_a<Index>(lr[1])) c << "MTD[" << lr[0] << "," << lr[1] << "]";
    }
    
    size_t Pair::nops() const { return 2; }
    ex Pair::op(size_t i) const {
        if(i==0) return lr[0];
        else if(i==1) return lr[1];
        return 0;
    }
    
    void Pair::archive(archive_node & n) const {
        inherited::archive(n);
        for(int i=0; i<2; i++) n.add_ex("lr"+to_string(i), lr[i]);
    }
    
    void Pair::read_archive(const archive_node& n, lst& sym_lst) {
        inherited::read_archive(n, sym_lst);
        for(int i=0; i<2; i++) {
            n.find_ex("lr"+to_string(i), lr[i], sym_lst);
        }
    }
    
    ex Pair::derivative(const symbol & s) const {
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

