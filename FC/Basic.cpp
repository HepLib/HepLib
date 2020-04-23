/**
 * @file
 * @brief Basic Functions for FC
 * @author F. Feng
 * @version 1.0.0
 * @date 2020-04-21
 */
 
#include "FC.h"

namespace HepLib::FC {

    DEFAULT_CTOR(Index)
    GINAC_BIND_UNARCHIVER(Index);
    IMPLEMENT_HAS(Index)
    DEFAULT_CTOR(Vector)
    GINAC_BIND_UNARCHIVER(Vector);
    IMPLEMENT_HAS(Vector)
    DEFAULT_CTOR(SUNT)
    GINAC_BIND_UNARCHIVER(SUNT);
    IMPLEMENT_HAS(SUNT)
    DEFAULT_CTOR(SUNF)
    GINAC_BIND_UNARCHIVER(SUNF);
    IMPLEMENT_HAS(SUNF)

    //-----------------------------------------------------------
    // FormFormat Output
    //-----------------------------------------------------------
    FormFormat::FormFormat(ostream &os, unsigned opt) : print_dflt(os, opt) {}
    FormFormat::FormFormat() : print_dflt(std::cout) {}
    GINAC_IMPLEMENT_PRINT_CONTEXT(FormFormat, print_dflt)
    OUT_FORMAT_IMPLEMENT(FormFormat)
    
    void FormFormat::power_print(const power & p, const FormFormat & c, unsigned level) {
        if(p.op(1)==2 && !DiracGamma::has(p)) {
            c << "((" << p.op(0) << ")*(" << p.op(0) << "))";
        } else {
            c << "(" << p.op(0) << ")^(" << p.op(1) << ")";
        }
    }
    
    //-----------------------------------------------------------
    // FCFormat Output
    //-----------------------------------------------------------
    FCFormat::FCFormat(ostream &os, unsigned opt) : print_dflt(os, opt) {}
    FCFormat::FCFormat() : print_dflt(std::cout) {}
    GINAC_IMPLEMENT_PRINT_CONTEXT(FCFormat, print_dflt)
    OUT_FORMAT_IMPLEMENT(FCFormat)
    
    namespace {
        class ncmul_hack : public ncmul { // due to printseq is protected
        public:
            ncmul_hack(ncmul _nm) : ncmul(_nm){ }
            void print(const FCFormat & c, unsigned level) {
                printseq(c, '(', '.', ')', precedence(), level);
            }
        };
    }
    void FCFormat::ncmul_print(const ncmul & nm, const FCFormat & c, unsigned level) {
        ncmul_hack(nm).print(c, level);
    }

    //-----------------------------------------------------------
    // Index Class
    //-----------------------------------------------------------
    GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(Index, basic,
        print_func<print_context>(&Index::print)
    )
    
    Index::Index(const string &s, const Type t) : name(s,true,false), type(t) { }
    int Index::compare_same_type(const basic &other) const {
        const Index &o = static_cast<const Index &>(other);
        auto ret = name.get_name().compare(o.name.get_name());
        if(ret==0) return 0;
        else if(ret<0) return -1;
        else return 1;
    }
    
    void Index::print(const print_context &c, unsigned level) const {
        c.s << name;
    }
    
    Pair Index::operator() (const Index &i) {
        return Pair(*this, i);
    }
    
    Pair Index::operator() (const Vector &p) {
        return Pair(p, *this);
    }
    
    void Index::archive(archive_node & n) const {
        inherited::archive(n);
        n.add_string("name", name.get_name());
        n.add_unsigned("type", type);
    }
    
    void Index::read_archive(const archive_node& n, lst& sym_lst) {
        inherited::read_archive(n, sym_lst);
        string nstr;
        unsigned t;
        n.find_string("name", nstr);
        name = Symbol(nstr,true,false);
        n.find_unsigned("type", t);
        type = (Type)t;
    }
    
    ex Index::derivative(const symbol & s) const {
        return 0;
    }

    //-----------------------------------------------------------
    // Vector Class
    //-----------------------------------------------------------
    GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(Vector, basic,
        print_func<print_context>(&Vector::print)
    )
    
    Vector::Vector(const string &s) : name(s,true,false) { }
    int Vector::compare_same_type(const basic &other) const {
        const Vector &o = static_cast<const Vector &>(other);
        auto ret = name.get_name().compare(o.name.get_name());
        if(ret==0) return 0;
        else if(ret<0) return -1;
        else return 1;
    }
    
    void Vector::print(const print_context &c, unsigned level) const {
        c.s << name;
    }
    
    Pair Vector::operator() (const Vector &p) {
        return Pair(*this, p);
    }
    
    Pair Vector::operator() (const Index &i) {
        return Pair(*this, i);
    }
    
    void Vector::archive(archive_node & n) const {
        inherited::archive(n);
        n.add_string("name", name.get_name());
    }
    
    void Vector::read_archive(const archive_node& n, lst& sym_lst) {
        inherited::read_archive(n, sym_lst);
        string nstr;
        unsigned t;
        n.find_string("name", nstr);
        name = Symbol(nstr,true,false);
    }
    
    ex Vector::derivative(const symbol & s) const {
        return 0;
    }
    
    //-----------------------------------------------------------
    // SUNT/SUNF Class
    //-----------------------------------------------------------
    GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(SUNT, basic,
        print_func<print_dflt>(&SUNT::print).
        print_func<FormFormat>(&SUNT::form_print)
    )
    
    SUNT::SUNT(ex i, ex j, ex a) : ija{i,j,a} { }
    int SUNT::compare_same_type(const basic &other) const {
        const SUNT &o = static_cast<const SUNT &>(other);
        for(int i=0; i<3; i++) {
            auto c = ija[i].compare(o.ija[i]);
            if(c!=0) return c;
        }
        return 0;
    }
    
    void SUNT::form_print(const FormFormat &c, unsigned level) const {
        c << "T(" << ija[0] << "," << ija[1] << "," << ija[2] << ")";
    }
    void SUNT::print(const print_dflt &c, unsigned level) const {
        c.s << "T(" << ija[0] << "," << ija[1] << "," << ija[2] << ")";
    }
    
    size_t SUNT::nops() const { return 3; }
    ex SUNT::op(size_t i) const {
        return ija[i];
    }
    ex & SUNT::let_op(size_t i) {
        ensure_if_modifiable();
        return ija[i];
    }
    
    void SUNT::archive(archive_node & n) const {
        inherited::archive(n);
        n.add_ex("i", ija[0]);
        n.add_ex("j", ija[1]);
        n.add_ex("a", ija[2]);
    }
    
    void SUNT::read_archive(const archive_node& n, lst& sym_lst) {
        inherited::read_archive(n, sym_lst);
        ex o;
        n.find_ex("i", o, sym_lst);
        ija[0] = ex_to<Index>(o);
        n.find_ex("j", o, sym_lst);
        ija[1] = ex_to<Index>(o);
        n.find_ex("a", o, sym_lst);
        ija[2] = ex_to<Index>(o);
    }
    
    ex SUNT::derivative(const symbol & s) const {
        return 0;
    }
    
    ex SUNT::conjugate() const {
        return SUNT(ija[1], ija[0], ija[2]);
    }
    
    GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(SUNF, basic,
        print_func<print_dflt>(&SUNF::print).
        print_func<FormFormat>(&SUNF::form_print)
    )
    
    SUNF::SUNF(ex i, ex j, ex k) : ijk{i,j,k} { }
    int SUNF::compare_same_type(const basic &other) const {
        const SUNF &o = static_cast<const SUNF &>(other);
        for(int i=0; i<3; i++) {
            auto c = ijk[i].compare(o.ijk[i]);
            if(c!=0) return c;
        }
        return 0;
    }
    
    void SUNF::print(const print_dflt &c, unsigned) const {
        c.s << "f(" << ijk[0] << "," << ijk[1] << "," << ijk[2] << ")";
    }
    
    void SUNF::form_print(const FormFormat &c, unsigned) const {
        c << "f(" << ijk[0] << "," << ijk[1] << "," << ijk[2] << ")";
    }
    
    size_t SUNF::nops() const { return 3; }
    ex SUNF::op(size_t i) const {
        return ijk[i];
    }
    ex & SUNF::let_op(size_t i) {
        ensure_if_modifiable();
        return ijk[i];
    }
    
    void SUNF::archive(archive_node & n) const {
        inherited::archive(n);
        n.add_ex("i", ijk[0]);
        n.add_ex("j", ijk[1]);
        n.add_ex("k", ijk[2]);
    }
    
    void SUNF::read_archive(const archive_node& n, lst& sym_lst) {
        inherited::read_archive(n, sym_lst);
        ex o;
        n.find_ex("i", o, sym_lst);
        ijk[0] = ex_to<Index>(o);
        n.find_ex("j", o, sym_lst);
        ijk[1] = ex_to<Index>(o);
        n.find_ex("k", o, sym_lst);
        ijk[2] = ex_to<Index>(o);
    }
    
    ex SUNF::derivative(const symbol & s) const {
        return 0;
    }
    
    //-----------------------------------------------------------
    // Other Helpers
    //-----------------------------------------------------------
    
    ex SUNF4(Index i, Index j, Index k, Index l) {
        static int idx=0;
        Index e("f4idx" + to_string(idx++));
        return SUNF(i, j, e) * SUNF(e, k, l);
    }
    
    ex MatrixContract(const ex & expr_in) {
        if(expr_in.has(coVF(w))) throw Error("MatrixContract: coVF found in expr_in.");
        
        auto expr = expr_in.subs(pow(Matrix(w1,w2,w3),2)==Matrix(w1,w2,w3)*Matrix(w1,w3,w2));
        expr = mma_collect(expr, Matrix(w1, w2, w3), false, true);
        expr = MapFunction([](const ex &e, MapFunction &self)->ex {
            if(e.match(coVF(w))) {
                if(e.op(0).match(Matrix(w1, w2, w3))) return e.op(0);
                if(!is_a<mul>(e.op(0))) throw Error("MatrixContract: mma_collect error: " + ex2str(e));
                lst mats;
                std::map<ex,int,ex_is_less> to_map, from_map;
                std::set<int> todo;
                for(auto item : e.op(0)) mats.append(item);
                lst mats_idx;
                for(int i=0; i<mats.nops(); i++) {
                    auto item = mats.op(i);
                    if(item.op(0).return_type()==return_types::commutative || item.op(0).is_equal(GAS(1))) {
                        mats_idx.append(lst{item,i});
                    } else {
                        if(to_map[item.op(1)]!=0 || from_map[item.op(2)]!=0) throw Error("MatrixContract: index conflict (1).");
                        to_map[item.op(1)] = i+10; // avoid 0 in map
                        from_map[item.op(2)] = i+10; // avoid 0 in map
                    }
                    todo.insert(i);
                }
                
                // update to_map/from_map w.r.t mats_idx
                bool checked = false;
                while(true) {
                    lst mats_idx2;
                    bool ok = true; // double check
                    for(int i=0; i<mats_idx.nops(); i++) {
                        auto item = mats_idx.op(i).op(0);
                        int ii = ex_to<numeric>(mats_idx.op(i).op(1)).to_int();
                        if(!checked && 
                            to_map[item.op(1)]==0 && from_map[item.op(2)]==0 && 
                            to_map[item.op(2)]==0 && from_map[item.op(1)]==0) {
                            mats_idx2.append(mats_idx.op(i));
                            continue;
                        }
                        ok = false;
                        checked = false;
                        if(to_map[item.op(1)]==0 && from_map[item.op(2)]==0) {
                            to_map[item.op(1)] = ii+10; // avoid 0 in map
                            from_map[item.op(2)] = ii+10; // avoid 0 in map
                        } else if(to_map[item.op(2)]==0 && from_map[item.op(1)]==0) {
                            to_map[item.op(2)] = ii+10; // avoid 0 in map
                            from_map[item.op(1)] = ii+10; // avoid 0 in map
                            // need to swap the 2nd and 3rd index
                            auto li = get_op(mats, ii, 1); 
                            auto ri = get_op(mats, ii, 2);
                            let_op(mats, ii, 1, ri);
                            let_op(mats, ii, 2, li);
                        } else {
                            throw Error("MatrixContract: index conflict (2).");
                        }
                    }
                    if(mats_idx2.nops()<1) break;
                    mats_idx = mats_idx2;
                    if(ok) checked=true;
                }
                
                ex retMat = 1;
                while(todo.size()>0) {
                    int c = *(todo.begin());
                    todo.erase(c);
                    ex curMat = mats.op(c).op(0);
                    auto li=mats.op(c).op(1);
                    auto ri=mats.op(c).op(2);
                    while(true) {
                        if(is_zero(li-ri)) {
                            retMat *= TR(curMat);
                            break;
                        }
                        int ti = to_map[ri];
                        int fi = from_map[li];
                        if(ti==0 && fi==0) {
                            retMat *= Matrix(curMat, li, ri);
                            break;
                        }
                        if(ti!=0) {
                            auto mat = mats.op(ti-10).op(0);
                            if(is_zero(curMat-GAS(1))) curMat = mat;
                            else if(!is_zero(mat-GAS(1))) curMat = curMat * mat;
                            ri = mats.op(ti-10).op(2);
                            todo.erase(ti-10);
                            continue;
                        }
                        if(fi!=0) {
                            auto mat = mats.op(fi-10).op(0);
                            if(is_zero(curMat-GAS(1))) curMat = mat;
                            else if(!is_zero(mat-GAS(1))) curMat = mat * curMat;
                            li = mats.op(fi-10).op(1);
                            todo.erase(fi-10);
                            continue;
                        }
                    }
                }
                return retMat;
            } else if(!e.has(coVF(w))) return e;
            else return e.map(self);
        })(expr);
        
        return expr;
    }
    
}

