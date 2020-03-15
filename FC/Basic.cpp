#include "FC.h"

namespace HepLib::FC {

    //-----------------------------------------------------------
    // FormFormat Output
    //-----------------------------------------------------------
    FormFormat::FormFormat(ostream &os, unsigned opt) : print_dflt(os, opt) {}
    FormFormat::FormFormat() : print_dflt(std::cout) {}
    GINAC_IMPLEMENT_PRINT_CONTEXT(FormFormat, print_dflt)
    OUT_FORMAT_IMPLEMENT(FormFormat)
    
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
    
    Index::Index(const string &s, const Type type) : name(get_symbol(s)), IndexType(type) { }
    int Index::compare_same_type(const basic &other) const {
        const Index &o = static_cast<const Index &>(other);
        auto c = name.compare(o.name);
        if(c!=0) return c;
        if(IndexType > o.IndexType) return 1;
        if(IndexType < o.IndexType) return -1;
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

