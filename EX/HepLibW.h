/**
 * @file
 * @brief Wrap Class for SWIG
 */

#include "HepLib.h"

class MapFunction;
class expr {
public:
    expr();
    expr(int i);
    expr(GiNaC::ex e);
    expr(const std::string &s);
    
    expr operator+(const expr &e);
    expr operator-(const expr &e);
    expr operator*(const expr &e);
    expr operator/(const expr &e);
    expr operator==(const expr &e);
    
    expr operator+(const int i);
    expr operator-(const int i);
    expr operator*(const int i);
    expr operator/(const int i);
    
    expr operator-();
    
    unsigned int nops();
    expr op(unsigned int i);
    void let_op(unsigned int i, expr e);
    
    expr expand();
    expr normal();
    expr factor();
    expr series(const expr &s, int o);
    expr subs(const std::vector<expr> &ev);
    expr subs(const expr &e);
    
    std::string str();
    std::string __str__();
    
    //friend functions
    friend expr expand(const expr &e);
    friend expr normal(const expr &e);
    friend expr factor(const expr &e);
    friend expr subs(const expr &e, const std::vector<expr> &ev);
    friend expr subs(const expr &e1, const expr &e2);
    bool match(const expr &e);
    bool isSymbol();
    bool isVector();
    bool isIndex();
    bool isPair();
    bool isDiracGamma();
    bool info(std::string sflags);
    expr map(MapFunction &mf);
    
    friend expr series(const expr &e, const expr &s, int o);
    friend expr pow(const expr &e1, const expr &e2);
    friend expr pow(const expr &e, const int n);
    
    friend expr SP(const expr &e1, const expr &e2);
    friend expr GAS(const expr &e);
    friend expr TR(const expr &e);
    friend expr form(const expr &e);
    
    friend void letSP(const expr &e1, const expr &e2, const expr &e12);
    friend expr call(const std::string func, const std::vector<expr> &ev);
    friend expr call(const std::string func, const expr &e);
    
    friend class MapFunction;
    friend class Integral;
    
private:
    GiNaC::ex _expr;
};

expr expand(const expr &e);
expr normal(const expr &e);
expr factor(const expr &e);
expr series(const expr &e, const expr &s, int o);
expr subs(const expr &e, const std::vector<expr> &ev);
expr subs(const expr &e1, const expr &e2);

expr pow(const expr &e1, const expr &e2);
expr pow(const expr &e1, const int n);

expr Index(const std::string &s);
expr Vector(const std::string &s);
expr Symbol(const std::string &s);

expr SP(const expr &e1, const expr &e2);
expr GAS(const expr &e);
expr TR(const expr &e);
expr form(const expr &e);

expr x(const int i);
expr y(const int i);
expr z(const int i);

void letSP(const expr &e1, const expr &e2, const expr &e12);

expr call(const std::string func, const std::vector<expr> &ev);
expr call(const std::string func, const expr &e);
expr w(const int wi);

class MapFunction {
public:
    MapFunction();
    virtual ~MapFunction();
    virtual expr map(const expr &e);
    expr operator() (const expr &e);
    friend class expr;
private:
    HepLib::MapFunction _map;
};

class Integral {
public:
    int epN = 0;
    int epsN = 0;
    int verb = 0;

    Integral();
    void Functions(const std::vector<expr> &ev);
    void Propagators(const std::vector<expr> &ev1, const std::vector<expr> &loops={}, const std::vector<expr> &tloops={});
    void Replacements(const std::vector<expr> &ev, int lt=0);
    void Exponents(const std::vector<expr> &ev);
    void Exponents(const std::vector<int> &ev);
    void Evaluate();
    std::string str();
    std::string __str__();
private:
    GiNaC::exvector _Propagators_Functions;
    GiNaC::exvector _Exponents;
    GiNaC::exvector _lReplacements;
    GiNaC::exvector _tReplacements;
    GiNaC::exvector _Loops;
    GiNaC::exvector _tLoops;
    GiNaC::ex _Result;
    bool isX = true;
};
