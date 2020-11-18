/* File: HepLib.i */
%module(directors="1") HepLib

%feature("director") MapFunction;
%feature("allowexcept");

%{
#define SWIG_FILE_WITH_INIT
#include "HepLibW.h"
%}

//--------------------------------------------------------------------

%include std_string.i
%include exception.i
%include std_vector.i

%exception {
    try {
        $action
    } catch(const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (const std::string& e) {
        SWIG_exception(SWIG_RuntimeError, e.c_str());
    } catch(...) {
        SWIG_exception(SWIG_RuntimeError, "HepLib unkonwn exception.");
    }
}

namespace std {
    %template(expr_vector) vector<expr>;
    %template(int_vector) vector<int>;
}

class MapFunction {
public:
    MapFunction();
    virtual ~MapFunction();
    virtual expr map(const expr &e);
    expr operator() (const expr &e);
};

class expr {
public:
    expr(int i);
    expr(const std::string &s);
    
    expr operator+(const expr &e);
    expr operator-(const expr &e);
    expr operator*(const expr &e);
    expr operator/(const expr &e);
    expr operator-();
    expr operator==(const expr &e);
    
    expr operator+(const int i);
    expr operator-(const int i);
    expr operator*(const int i);
    expr operator/(const int i);
    
    unsigned int nops();
    expr op(unsigned int i);
    void let_op(unsigned int i, expr e);
    
    expr expand();
    expr normal();
    expr factor();
    expr series(const expr &s, int o);
    expr subs(const std::vector<expr> &ev);
    expr subs(const expr &e);
    bool match(const expr &e);
    bool isSymbol();
    bool isVector();
    bool isIndex();
    bool isPair();
    bool isDiracGamma();
    bool info(std::string sflags);
    expr map(MapFunction &mf);
    
    std::string str();
    std::string __str__();
    
    %pythoncode %{
    def __radd__(self, other):
        return expr(other) + self
    def __rsub__(self, other):
        return expr(other) - self
    def __rmul__(self, other):
        return expr(other) * self
    def __rdiv__(self, other):
        return expr(other) / self
    
    %}
};

extern expr expand(const expr &e);
extern expr normal(const expr &e);
extern expr factor(const expr &e);
extern expr series(const expr &e, const expr &s, int o);
extern expr subs(const expr &e, const std::vector<expr> &ev);
extern expr subs(const expr &e1, const expr &e2);

extern expr pow(const expr &e1, const expr &e2);
extern expr pow(const expr &e, const int n);

extern expr Index(const std::string &s);
extern expr Vector(const std::string &s);
extern expr Symbol(const std::string &s);
extern expr SP(const expr &e1, const expr &e2);
extern expr GAS(const expr &e);
extern expr TR(const expr &e);
extern expr form(const expr &e);

extern void letSP(const expr &e1, const expr &e2, const expr &e12);
extern expr call(const std::string func, const std::vector<expr> &ev);
extern expr call(const std::string func, const expr &e);
extern expr w(const int wi);

extern expr x(const int i);
extern expr y(const int i);
extern expr z(const int i);

// Integral
%warnfilter(509) Integral;
class Integral {
public:
    int epN;
    int epsN;
    int verb;
    
    Integral();
    void Functions(const std::vector<expr> &ev);
    void Propagators(const std::vector<expr> &ev, const std::vector<expr> &loops={}, const std::vector<expr> &tloops={});
    void Replacements(const std::vector<expr> &ev, int lt=0);
    void Exponents(const std::vector<expr> &ev);
    void Exponents(const std::vector<int> &ev);
    void Evaluate();
    
    std::string str();
    std::string __str__();
};
