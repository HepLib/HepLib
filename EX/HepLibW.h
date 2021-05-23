/**
 * @file
 * @brief Wrap Class for SWIG
 */

#include "HepLib.h"

class MapFunction;
class expr {
public:
    GiNaC::ex _expr;
    
    expr();
    expr(int i);
    expr(GiNaC::ex e);
    expr(const std::string &s);
    expr(const std::string &s, std::map<std::string,expr>);
    expr(const std::vector<expr> &ev);
    
    expr operator+(const expr &e);
    expr operator-(const expr &e);
    expr operator*(const expr &e);
    expr operator/(const expr &e);
    expr operator>>(const expr &e);
    bool operator<(const expr &e) const;
    bool operator>(const expr &e) const;
    
    expr operator+(const int i);
    expr operator-(const int i);
    expr operator*(const int i);
    expr operator/(const int i);
    
    bool is_equal(const expr & e);
    unsigned gethash();
    
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
    
    bool match(const expr &e);
    bool isSymbol();
    bool isVector();
    bool isIndex();
    bool isPair();
    bool isDGamma();
    bool info(std::string sflags);
    expr map(MapFunction &mf);

};

expr expand(const expr &e);
expr normal(const expr &e);
expr factor(const expr &e);
expr series(const expr &e, const expr &s, int o);
expr subs(const expr &e, const std::vector<expr> &ev);
expr subs(const expr &e1, const expr &e2);

expr pow(const expr &e1, const expr &e2);
expr pow(const expr &e1, const int n);

expr abs(const expr &z);
expr real(const expr &z);
expr imag(const expr &z);
expr csgn(const expr &z);
expr step(const expr &z);
expr numer(const expr &z);
expr denom(const expr &z);
expr sqrt(const expr &z);
expr sin(const expr &z);
expr cos(const expr &z);
expr tan(const expr &z);
expr asin(const expr &z);
expr acos(const expr &z);
expr atan(const expr &y, const expr &x);
expr sinh(const expr &z);
expr cosh(const expr &z);
expr tanh(const expr &z);
expr asinh(const expr &z);
expr acosh(const expr &z);
expr atanh(const expr &z);
expr exp(const expr &z);
expr log(const expr &z);
expr Li2(const expr &z);
expr zeta(const expr &z);
expr tgamma(const expr &z);
expr lgamma(const expr &z);
expr psi(const expr &z);
expr psi(const expr &n, const expr &z);
expr factorial(const expr &n);
expr doublefactorial(const expr &n);
expr binomial(const expr &n, const expr &k);
expr bernoulli(const expr &n);
expr fibonacci(const expr &n);
expr mod(const expr &a, const expr &b);
expr smod(const expr &a, const expr &b);
expr irem(const expr &a, const expr &b);
expr irem(const expr &a, const expr &b, const expr &q);
expr iquo(const expr &a, const expr &b);
expr iquo(const expr &a, const expr &b, const expr &r);
expr gcd(const expr &a, const expr &b);
expr lcm(const expr &a, const expr &b);

expr Index(const std::string &s);
expr IndexCA(const std::string &s);
expr IndexCF(const std::string &s);
expr Vector(const std::string &s);
expr Symbol(const std::string &s);

expr SUNT(const expr &e, const expr &i, const expr &j);
expr SUNF(const expr &a, const expr &b, const expr &c);
expr SUNF4(const expr &a, const expr &b, const expr &c, const expr &d);
expr Eps(const expr &a, const expr &b, const expr &c, const expr &d);
expr LC(const expr &a, const expr &b, const expr &c, const expr &d);

expr SP(const expr &e1, const expr &e2);
expr GAS(const expr &e);
expr GAS(const int &i);
expr TR(const expr &e);
expr form(const expr &e, int verb=0);

expr charge_conjugate(const expr &);
expr TIR(const expr &expr_in, const std::vector<expr> &loop_ps, const std::vector<expr> &ext_ps);
expr MatrixContract(const expr & expr_in);
expr Apart(const expr &expr_in, const std::vector<expr> &vars, std::map<expr, expr> sgnmap={});
expr Apart(const expr &expr_in, const std::vector<expr> &loops, const std::vector<expr> & extmoms, std::map<expr, expr> sgnmap={});
expr ApartIR2ex(const expr & expr_in);
expr ApartIR2F(const expr & expr_in);
expr F2ex(const expr & expr_in);
expr ApartIRC(const expr & expr_in);
void ApartIBP(int IBPmethod, std::vector<expr> &io_vec, const std::vector<expr> & loops, const std::vector<expr> & exts,
    const std::vector<expr> & cut_props={});
        
class AIOption {
    expr Internal; // Internal for Apart/IBP
    expr External; // External for Apart/IBP
    expr DSP; // DSP for IBP
    std::map<expr, expr> smap; // Sign Map for Apart
    std::vector<expr> Cuts; // Cut Propagators. optional
    std::vector<expr> CSP; // SP in Cuts, to be cleared. optional
    std::vector<expr> ISP; // SP for IBP. optional
    bool CutFirst = true;
    int mcl = 1; // collect_ex level, 0-nothing, 1-exnormal, 2-exfactor
};
void ApartIBP(int IBPmethod, std::vector<expr> &io_vec, AIOption aip);

expr x(const int i);
expr y(const int i);
expr z(const int i);
expr lst(const std::vector<expr> &ev);
void set_form_using_su3(bool yn);

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

class Process {
public:
    static std::string Style;
    static void DrawPDF(std::vector<expr>, std::string);
    std::string Model;
    std::string In;
    std::string Out;
    std::string LoopPrefix = "q";
    int Loops;
    std::string Options;
    std::vector<std::string> Others;
    std::vector<expr> Amplitudes(std::map<std::string, expr>, bool debug=false);
};
void set_LineTeX(expr, std::string);
void set_InOutTeX(int, std::string);

class Function {
public:
    Function();
    virtual ~Function();
    virtual expr operator()(const expr &e);
    virtual expr operator()(const expr &e1,const expr &e2);
    virtual expr operator()(const expr &e1,const expr &e2,const expr &e3);
    virtual expr operator()(const expr &e1,const expr &e2,const expr &e3,const expr &e4);
    virtual expr operator()(const expr &e1,const expr &e2,const expr &e3,const expr &e4,const expr &e5);
};
