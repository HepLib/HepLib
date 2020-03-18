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
    // Error Class
    //-----------------------------------------------------------
    Error::Error(const char * _msg) : msg(_msg) { }
    const char * Error::what() const throw () {
        return msg.c_str();
    }
    
    //-----------------------------------------------------------
    // Index Class
    //-----------------------------------------------------------
    GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(Index, basic,
        print_func<print_context>(&Index::print)
    )
    
    Index::Index(const string &s, const Type t) : name(get_symbol(s)), type(t) { }
    int Index::compare_same_type(const basic &other) const {
        const Index &o = static_cast<const Index &>(other);
        auto c = name.compare(o.name);
        if(c!=0) return c;
        if(type > o.type) return 1;
        if(type < o.type) return -1;
        return 0;
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
        name = get_symbol(nstr);
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
    
    Vector::Vector(const string &s) : name(get_symbol(s)) { }
    int Vector::compare_same_type(const basic &other) const {
        const Vector &o = static_cast<const Vector &>(other);
        return name.compare(o.name);
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
        name = get_symbol(nstr);
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
    
    SUNT::SUNT(Index i, Index j, Index a) : ija{i,j,a} { }
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
    
    SUNF::SUNF(Index i, Index j, Index k) : ijk{i,j,k} { }
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
    
    ex SUNSimplify(const ex & inexpr) {
        auto ret = inexpr;
        int i = 0;
        ret = MapFunction([&](const ex &e, MapFunction &self)->ex {
            ex I2R = ex(1)/2;
            if(is_a<SUNF>(e)) {
                Index colFi1 = Index("ssCOL"+to_string(i++), Index::Type::CF);
                Index colFi2 = Index("ssCOL"+to_string(i++), Index::Type::CF);
                Index colFi3 = Index("ssCOL"+to_string(i++), Index::Type::CF);
                Index colAj1 = ex_to<Index>(e.op(0));
                Index colAj2 = ex_to<Index>(e.op(1));
                Index colAj3 = ex_to<Index>(e.op(2));
                return 1/I2R/I*SUNT(colFi1,colFi2,colAj1)*SUNT(colFi2,colFi3,colAj2)*SUNT(colFi3,colFi1,colAj3) -
                    1/I2R/I*SUNT(colFi1,colFi2,colAj3)*SUNT(colFi2,colFi3,colAj2)*SUNT(colFi3,colFi1,colAj1);
            }
            return e.map(self);
        })(ret);
        
        ret = MapFunction([&](const ex &e, MapFunction &self)->ex {
            if(is_a<SUNT>(e)) {
                Index col1 = ex_to<Index>(e.op(0));
                Index col2 = ex_to<Index>(e.op(1));
                Index col3 = ex_to<Index>(e.op(2));
                return iWF(col1, col2, col3);
            }
            return e.map(self);
        })(ret);
        ret = mma_collect(ret, iWF(w1, w2, w3));
        
        while(true) {
            auto tmp = ret.subs(iWF(w1,w2,w3)*iWF(w4,w5,w3)==iWF(w1,w2,w4,w5), subs_options::algebraic);
            if(is_zero(tmp-ret)) break;
            ret = tmp;
        }
        return ret;
    }
    
}

