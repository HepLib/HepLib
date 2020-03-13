#pragma once

#include "Basic.h"
#include "Process.h"

namespace HepLib::FC {

    using namespace std;
    using namespace GiNaC;
    using namespace HepLib;
    
    extern symbol D;
    extern symbol CA;
    extern symbol CF;
    extern symbol NA;
    extern symbol NF;
    
    //-----------------------------------------------------------
    // FormFormat Output
    //-----------------------------------------------------------
    class FormFormat : public print_dflt {
        GINAC_DECLARE_PRINT_CONTEXT(FormFormat, print_dflt)
    public:
        FormFormat(ostream &os, unsigned opt=0);
        OUT_FORMAT_DECLARE(FormFormat)
    private:
    };
        
    //-----------------------------------------------------------
    // FCFormat Output
    //-----------------------------------------------------------
    class FCFormat : public print_dflt {
        GINAC_DECLARE_PRINT_CONTEXT(FCFormat, print_dflt)
    public:
        FCFormat(ostream &os, unsigned opt=0);
        static void ncmul_print(const ncmul & p, const FCFormat & c, unsigned level=0);
        OUT_FORMAT_DECLARE(FCFormat)
        // static initialization
        class _init {
            public: _init();
        };
    private:
        static _init FCFormat_init;
    };
    
    class Error : public exception {
    public:
        string msg;
        const char * what() const throw ();
        Error(const char * _msg);
    };
    
    class Index;
    class Vector;
    class Pair;
    
    //-----------------------------------------------------------
    // Index Class
    //-----------------------------------------------------------
    class Index : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(Index, basic)
    public:
        enum Type {VD, CF, CA};
        Index(const string &s, const Type type=Type::VD);
        Pair operator() (const Index & i);
        Pair operator() (const Vector & p);
        const symbol name;
        const Type IndexType;
        void print(const print_context &c, unsigned level = 0) const;
    };
    
    //-----------------------------------------------------------
    // Vector Class
    //-----------------------------------------------------------
    class Vector : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(Vector, basic)
    public:
        Vector(const string &s);
        Pair operator() (const Vector & p);
        Pair operator() (const Index & mu);
        const symbol name;
        void print(const print_context &c, unsigned level = 0) const;
    };
    
    //-----------------------------------------------------------
    // SUNT Class
    //-----------------------------------------------------------
    class SUNT : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(SUNT, basic)
    public:
        SUNT(Index a, Index i, Index j);
        const Index ija[3];
        size_t nops() const override;
        ex op(size_t i) const override;
        void form_print(const FormFormat &c, unsigned level = 0) const;
    };
    
    //-----------------------------------------------------------
    // SUNF Class
    //-----------------------------------------------------------
    class SUNF : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(SUNF, basic)
    public:
        SUNF(Index i, Index j, Index k);
        const Index ijk[3];
        size_t nops() const override;
        ex op(size_t i) const override;
        void form_print(const FormFormat &c, unsigned level = 0) const;
    };
    
    //-----------------------------------------------------------
    // Pair Class
    //-----------------------------------------------------------
    class Pair : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(Pair, basic)
    public:
        Pair(const Vector &p1, const Vector &p2);
        Pair(const Index &i1, const Index &i2);
        Pair(const Vector &p, const Index &i);
        size_t nops() const override;
        ex op(size_t i) const override;
        void print(const print_dflt &c, unsigned level = 0) const;
        void form_print(const FormFormat &c, unsigned level = 0) const;
        void fc_print(const FCFormat &c, unsigned level = 0) const;
    private:
        ex spL;
        ex spR;
    };
    
    // SP functions
    inline Pair SP(const Vector &p1, const Vector &p2) {
        if(p1.compare(p2)<0) return Pair(p2, p1);
        else return Pair(p1, p2);
    }
    inline Pair SP(const Index &i1, const Index &i2) {
        if(i1.compare(i2)<0) return Pair(i2, i1);
        else return Pair(i1, i2);
    }
    inline Pair SP(const Vector &p, const Index &i) {
        return Pair(p, i);
    }
    inline Pair SP(const Index &i, const Vector &p) {
        return Pair(p, i);
    }
    ex SP(ex a, ex b);
    
    //-----------------------------------------------------------
    // DiracGamma Class
    //-----------------------------------------------------------
    class DiracGamma : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(DiracGamma, basic)
    public:
        DiracGamma(const Vector &p, unsigned rl=0);
        DiracGamma(const Index &i, unsigned rl=0);
        DiracGamma(int int_1567, unsigned _rl=0);
        DiracGamma(const DiracGamma &g, unsigned _rl);
        void print(const print_dflt &c, unsigned level = 0) const;
        void form_print(const FormFormat &c, unsigned level = 0) const;
        void fc_print(const FCFormat &c, unsigned level = 0) const;
        return_type_t return_type_tinfo() const override;
        unsigned return_type() const override { return return_types::noncommutative; }
        bool match_same_type(const basic & other) const override;
        unsigned get_rl();
        size_t nops() const override;
        ex op(size_t i) const override;
    private:
        ex pi;
        unsigned rl;
        
    };
    
    //-----------------------------------------------------------
    // TR/GAS functions
    //-----------------------------------------------------------
    DECLARE_FUNCTION_1P(TR)
    
    inline ex GAS(const Vector &p) { return DiracGamma(p); }
    inline ex GAS(const Index &i) { return DiracGamma(i); }
    ex GAS(ex expr);
    
    // Form
    ex form(ex expr, bool verb=false);

    
    
}
