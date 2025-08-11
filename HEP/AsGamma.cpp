/**
 * @file
 * @brief Functions for AsGamma
 */
 
#include "HEP.h"

namespace HepLib {

    namespace {
        inline int perm_sign(int perm[], int n) {
            int cnt = 0;
            // 计算逆序对的数量
            for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) if (perm[i] > perm[j]) cnt++;
            return (cnt % 2 == 0) ? 1 : -1;
        }
    }
    
    //-----------------------------------------------------------
    // AsGamma Class
    //-----------------------------------------------------------
    //GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(AsGamma, basic,print_func<print_dflt>(&AsGamma::print).print_func<FormFormat>(&AsGamma::form_print).print_func<FCFormat>(&AsGamma::fc_print))
    GiNaC::registered_class_info & AsGamma::get_class_info_static() { return reg_info; }
    AsGamma::visitor::~visitor() { }
    AsGamma * AsGamma::duplicate() const { AsGamma * bp = new AsGamma(*this); bp->setflag(GiNaC::status_flags::dynallocated); return bp; }
    void AsGamma::accept(GiNaC::visitor & v) const { if (visitor *p = dynamic_cast<visitor *>(&v)) p->visit(*this); else inherited::accept(v); }
    const GiNaC::registered_class_info &AsGamma::get_class_info() const { return get_class_info_static(); }
    GiNaC::registered_class_info &AsGamma::get_class_info() { return get_class_info_static(); }
    const char *AsGamma::class_name() const { return get_class_info_static().options.get_name(); }
    //GINAC_IMPLEMENT_REGISTERED_CLASS END
    
    DEFAULT_CTOR(AsGamma)
    IMPLEMENT_HAS(AsGamma)
    IMPLEMENT_ALL(AsGamma)
    
    return_type_t AsGamma::return_type_tinfo() const {
        return make_return_type_t<clifford>(rl);
        //return make_return_type_t<AsGamma>(rl);
    }
    
    bool AsGamma::match_same_type(const basic & other) const {
        const AsGamma &o = static_cast<const AsGamma &>(other);
        return rl == o.rl;
    }
    
    unsigned AsGamma::get_rl() {
        return rl;
    }
    
    int AsGamma::compare_same_type(const basic &other) const {
        if(!is_a<AsGamma>(other)) throw Error("AsGamma::compare_same_type");
        const AsGamma &o = static_cast<const AsGamma &>(other);
        if (rl != o.rl) return rl < o.rl ? -1 : 1;
        int n = pis.size();
        int no = o.pis.size();
        if(n!=no) return n>no ? 1 : -1;
        for(int i=0; i<n; i++) {
            auto c = pis[i].compare(o.pis[i]);
            if(c!=0) return c;
        }
        return 0;
    }
    
    bool AsGamma::is_equal_same_type(const basic & other) const {
        if(!is_a<AsGamma>(other)) throw Error("AsGamma::is_equal_same_type");
        const AsGamma &o = static_cast<const AsGamma &>(other);
        if (rl != o.rl) return false;
        int n = pis.size();
        int no = o.pis.size();
        if(n!=no) return false;
        for(int i=0; i<n; i++) {
            if(!pis[i].is_equal(o.pis[i])) return false;
        }
        return true;
    }
    
    AsGamma::AsGamma(const AsGamma &g, unsigned _rl) : pis(g.pis), rl(_rl) { }
    AsGamma::AsGamma(const exvector & _pis, int _rl) : pis(_pis), rl(_rl) {
        lst pis_sort = vec2lst(_pis);
        int n = pis_sort.nops();
        for(int i=0; i<n; i++) if(is_a<DGamma>(pis_sort.op(i))) {
            pis_sort.let_op(i) = pis_sort.op(i).op(0);
            pis[i] = pis[i].op(0);
        }
        pis_sort.sort();
        for(int i=0; i<n; i++) if(pis_sort[i] != pis[i]) throw Error("Error: Input list is NOT sorted.");
    }
    AsGamma::AsGamma(const lst & _pis, int _rl) : pis(lst2vec(_pis)), rl(_rl) {
        lst pis_sort = _pis;
        int n = pis_sort.nops();
        for(int i=0; i<n; i++) if(is_a<DGamma>(pis_sort.op(i))) {
            pis_sort.let_op(i) = pis_sort.op(i).op(0);
            pis[i] = pis[i].op(0);
        }
        pis_sort.sort();
        for(int i=0; i<n; i++) if(pis_sort[i] != pis[i]) throw Error("Error: Input list is NOT sorted.");
    }

    size_t AsGamma::nops() const { return pis.size(); }
    ex AsGamma::op(size_t i) const {
        return pis[i];
    }
    ex & AsGamma::let_op(size_t i) {
        ensure_if_modifiable();
        return pis[i];
    }
    
    ex AsGamma::eval() const {
        if(flags & status_flags::evaluated) return *this;
        // make sure pis are sorted
        lst pis_sort = vec2lst(pis);
        pis_sort.sort();
        pis_sort.unique();
        int n = pis.size();
        int n2 = pis_sort.nops();
        if(n!=n2) return 0;
        
        std::map<ex,int,ex_is_less> e2i;
        for(int i=0; i<n; i++) e2i[pis[i]] = i;
        int si[n];
        for(int i=0; i<n; i++) si[i] = e2i[pis_sort.op(i)];
        
        bool ok = true;
        for(int i=0; i<n; i++) if(pis_sort.op(i)!=pis[i]) {
            ok = false;
            break;
        }
        if(ok) return this->hold();
        else {
            AsGamma ag(pis_sort);
            return perm_sign(si, n) * ag;
        }
    }
    
    void AsGamma::print(const print_dflt &c, unsigned level) const {
        if(pis.size()<1) {
            c.s << "\u0130";
            return;
        }
        c.s << "\u0263[";
        bool first = true;
        for(auto pi : pis) {
            if(first) first=false;
            else c.s << ",";
            c.s << pi;
        }
        c.s << "]";
    }
    
    void AsGamma::form_print(const FormFormat &c, unsigned level) const {
        throw Error("AsGamma::form_print unexpected region.");
    }
    
    void AsGamma::fc_print(const FCFormat &c, unsigned level) const {
        throw Error("AsGamma::fc_print unexpected region.");
    }
    
    void AsGamma::archive(archive_node & n) const {
        inherited::archive(n);
        n.add_unsigned("rl", rl);
        int nn = pis.size();
        n.add_ex("size", nn);
        for(int i=0; i<nn; i++) n.add_ex("pis"+to_string(i), pis[i]);
    }
    
    void AsGamma::read_archive(const archive_node& n) {
        inherited::read_archive(n);
        n.find_unsigned("rl", rl);
        ex nex;
        n.find_ex("size", nex);
        int nn = ex_to<numeric>(nex).to_int();
        if(pis.size()!=nn) pis.resize(nn);
        for(int i=0; i<nn; i++) {
            n.find_ex("pis"+to_string(i), pis[i]);
        }
    }
    
    ex AsGamma::derivative(const symbol & s) const {
        return 0;
    }
    
    ex AsGamma::conjugate() const {
        throw Error("AsGamma::conjugate unexpected region.");
    }
    
    ex AsGamma::from(const lst & pis_lst, unsigned rl) {
        lst pis_sort = pis_lst;
        int n = pis_sort.nops();
        for(int i=0; i<n; i++) if(is_a<DGamma>(pis_sort.op(i))) pis_sort.let_op(i) = pis_sort.op(i).op(0);
        lst tmp;
        for(auto item : pis_sort) {
            if(is_a<Vector>(item) || is_a<Index>(item)) tmp.append(item);
            else if(!item.is_equal(1)) throw Error("AsGamma::from only support GAS(1/i/p)."); // 1 just to skip
        }
        pis_sort = tmp;
        pis_sort.sort();
        n = pis_sort.nops();
        
        std::map<ex,int,ex_is_less> e2i;
        for(int i=0; i<n; i++) e2i[pis_lst.op(i)] = i;
        if(e2i.size()!=n) return 0; // two indices are equal
        int si[n];
        for(int i=0; i<n; i++) si[i] = e2i[pis_sort.op(i)];
        
        AsGamma ag(lst2vec(pis_sort), rl);
        
        return perm_sign(si, n) * ag;
    }
    
    ex AsGamma::to() const {
        int n = pis.size();
        if(n==0) return GAS(1, rl);
        int pn[n];
        for(int i=0; i<n; i++) pn[i] = i;
        
        ex res = 0;
        do {
            ex gi = 1;
            for(int i=0; i<n; i++) gi = gi * GAS(pis[pn[i]], rl);
            res += perm_sign(pn, n) * gi;
        } while (std::next_permutation(pn, pn+n));
        
        return res/factorial(n);
    }
    
    ex AsGamma::mul_right(const ex & pi_in) const {
        ex pi = pi_in;
        if(is_a<DGamma>(pi)) pi = pi.op(0);
        if(!is_a<Vector>(pi) && !is_a<Index>(pi)) throw Error("AsGamma::mul_right, argument should be DGamma/Index/Vector");
        lst pis2 = vec2lst(pis);
        pis2.append(pi);
        ex res = AsGamma::from(pis2, rl);
        int n = pis.size();
        for(int r=0; r<n; r++) {
            lst pis3;
            for(int i=0; i<n; i++) {
                if(i==r) continue;
                pis3.append(pis[i]);
            }
            res += ((n-r-1)%2==0 ? 1 : -1) * SP(pis[r], pi) * AsGamma::from(pis3, rl);
        }
        return res;
    }
    
    ex AsGamma::mul_left(const ex & pi_in) const {
        ex pi = pi_in;
        if(is_a<DGamma>(pi)) pi = pi.op(0);
        if(!is_a<Vector>(pi) && !is_a<Index>(pi)) throw Error("AsGamma::mul_left, argument should be DGamma/Index/Vector");
        lst pis2 = vec2lst(pis);
        pis2.prepend(pi);
        ex res = AsGamma::from(pis2, rl);
        int n = pis.size();
        for(int r=0; r<n; r++) {
            lst pis3;
            for(int i=0; i<n; i++) {
                if(i==r) continue;
                pis3.append(pis[i]);
            }
            res += (r%2==0 ? 1 : -1) * SP(pis[r], pi) * AsGamma::from(pis3, rl);
        }
        return res;
    }
    
    ex AsGamma_mul_right(const ex & expr, const ex & pi_in) {
        ex pi = pi_in;
        if(is_a<DGamma>(pi)) pi = pi.op(0);
        if(!is_a<Vector>(pi) && !is_a<Index>(pi)) throw Error("AsGamma_mul_left: argument should be DGamma/Index/Vector");
        auto to_sub = MapFunction([pi_in](const ex & e, MapFunction &self)->ex{
            if(!AsGamma::has(e)) return e;
            else if(is_a<AsGamma>(e)) {
                return ex_to<AsGamma>(e).mul_right(pi_in);
            } else return e.map(self);
        });
        return to_sub(expr);
    }
    
    ex AsGamma_mul_left(const ex & expr, const ex & pi_in) {
        ex pi = pi_in;
        if(is_a<DGamma>(pi)) pi = pi.op(0);
        if(!is_a<Vector>(pi) && !is_a<Index>(pi)) throw Error("AsGamma_mul_left: argument should be DGamma/Index/Vector");
        auto to_sub = MapFunction([pi_in](const ex & e, MapFunction &self)->ex{
            if(!AsGamma::has(e)) return e;
            else if(is_a<AsGamma>(e)) {
                return ex_to<AsGamma>(e).mul_left(pi_in);
            } else return e.map(self);
        });
        return to_sub(expr);
    }
    
    ex AsGamma_to(const ex & expr) {
        auto to_sub = MapFunction([](const ex & e, MapFunction &self)->ex{
            if(!AsGamma::has(e)) return e;
            else if(is_a<AsGamma>(e)) {
                return ex_to<AsGamma>(e).to();
            } else return e.map(self);
        });
        return GMatExpand(to_sub(expr));
    }

}

