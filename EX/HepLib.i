/**
 * @file
 * @brief Wrap Class for SWIG
 */
 
%module(directors="1") HepLib

%feature("director") MapFunction;
%feature("director") Function;
%feature("director") ParFun;
%feature("allowexcept");

%{
#define SWIG_FILE_WITH_INIT
#include "HepLibW.h"
%}

//--------------------------------------------------------------------
%include exception.i
%include std_string.i
%include std_vector.i
%include std_map.i

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
    %template(string_expr_map) map<string, expr>;
    %template(expr_expr_map) map<expr, expr, expr_is_less>;
    %template(int_longlong_map) std::map<int, long long>;
}

/*
-----------------------------------------
    GiNaC Wrapper
-----------------------------------------
*/

class expr {
public:
    expr(int i);
    expr(const std::string &s);
    expr(const std::string &s, std::map<std::string,expr>);
    expr(const std::vector<expr> &ev);
    
    expr operator+(const expr &e);
    expr operator-(const expr &e);
    expr operator*(const expr &e);
    expr operator/(const expr &e);
    expr operator-();
    bool operator<(const expr &e) const;
    bool operator==(const expr &e) const;
    expr operator+(const int i);
    expr operator-(const int i);
    expr operator*(const int i);
    expr operator/(const int i);
    
    expr operator>>(const expr &e); // a>>b (python) v.s. a==b (GiNaC)
    expr operator>>(const int &e); // a>>b (python) v.s. a==b (GiNaC)
    
    unsigned int nops();
    expr op(unsigned int i);
    void let_op(unsigned int i, expr e);
    expr __getitem__(const int i);
    
    expr expand();
    expr normal();
    expr factor();
    expr series(const expr &s, int o);
    expr subs(const std::vector<expr> &ev);
    expr subs(const expr &e);
    expr subs(const exmap & e);
    bool match(const expr &e);
    bool info(std::string sflags);
    expr map(MapFunction &mf);
    
    unsigned gethash();
    bool isSymbol();
    bool isVector();
    bool isIndex();
    bool isPair();
    bool isDGamma();
    
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
    def __hash__(self):
        return self.gethash()
    
    %}
};

extern expr conjugate(const expr &e);
extern expr expand(const expr &e);
extern expr normal(const expr &e);
extern expr factor(const expr &e);
extern expr series(const expr &e, const expr &s, int o);
extern expr subs(const expr &e, const std::vector<expr> &ev);
extern expr subs(const expr &e1, const expr &e2);
extern expr subs(const expr &e1, const exmap &e2);
extern expr pow(const expr &e1, const expr &e2);
extern expr pow(const expr &e, const int n);
extern expr abs(const expr &z);
extern expr real(const expr &z);
extern expr imag(const expr &z);
extern expr csgn(const expr &z);
extern expr step(const expr &z);
extern expr numer(const expr &z);
extern expr denom(const expr &z);
extern expr sqrt(const expr &z);
extern expr sin(const expr &z);
extern expr cos(const expr &z);
extern expr tan(const expr &z);
extern expr asin(const expr &z);
extern expr acos(const expr &z);
extern expr atan(const expr &y, const expr &x);
extern expr sinh(const expr &z);
extern expr cosh(const expr &z);
extern expr tanh(const expr &z);
extern expr asinh(const expr &z);
extern expr acosh(const expr &z);
extern expr atanh(const expr &z);
extern expr exp(const expr &z);
extern expr log(const expr &z);
extern expr Li2(const expr &z);
extern expr zeta(const expr &z);
extern expr tgamma(const expr &z);
extern expr lgamma(const expr &z);
extern expr psi(const expr &z);
extern expr psi(const expr &n, const expr &z);
extern expr factorial(const expr &n);
extern expr doublefactorial(const expr &n);
extern expr binomial(const expr &n, const expr &k);
extern expr bernoulli(const expr &n);
extern expr fibonacci(const expr &n);
extern expr mod(const expr &a, const expr &b);
extern expr smod(const expr &a, const expr &b);
extern expr irem(const expr &a, const expr &b);
extern expr irem(const expr &a, const expr &b, const expr &q);
extern expr iquo(const expr &a, const expr &b);
extern expr iquo(const expr &a, const expr &b, const expr &r);
extern expr gcd(const expr &a, const expr &b);
extern expr lcm(const expr &a, const expr &b);

extern expr lst(const std::vector<expr> &ev);
extern bool isFunction(const expr &e, std::string sf);

class exvec {
public:
    exvec();
    exvec(const std::vector<expr> es);
    exvec(expr e);
    void push_back(expr e);
    int size();
    expr __getitem__(const int i);
    void __setitem__(const int i, expr v);
    void subs(const expr & e);
    void subs(const std::vector<expr> & e);
    void subs(const exmap  &e);
    std::string str();
    std::string __str__();
};

class exmap {
public:
    exmap();
    exmap(std::map<expr,expr,expr_is_less> es);
    expr __getitem__(expr e);
    void __setitem__(expr k, expr v);
    int size();
    std::string str();
    std::string __str__();
};

class exset {
public:
    exset();
    exset(const std::vector<expr> es);
    exset(expr e);
    void insert(expr e);
    int size();
    std::string str();
    std::string __str__();
};

class MapFunction {
public:
    MapFunction();
    virtual ~MapFunction();
    virtual expr map(const expr &e);
    expr operator() (const expr &e);
};

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

class ParFun {
public:
    ParFun();
    virtual ~ParFun();
    virtual expr __call__(const int i);
};
extern exvec Parallel(int ntotal, int nbatch,
        ParFun &f,
        const std::string & key = "",
        bool rm = true,
        const std::string & pre = "  ");
extern exvec Parallel(int ntotal,
        ParFun &f,
        const std::string & key = "",
        bool rm = true,
        const std::string & pre = "  ");
        

extern std::map<std::string,expr> garReadAll(const std::string &garfn);
extern expr garRead(const std::string &garfn, const char* key);
extern expr garRead(const std::string &garfn);
extern void garWrite(const std::string &garfn, const std::map<std::string, expr> &resMap);
extern void garWrite(const std::map<std::string, expr> &resMap, const std::string &garfn);
extern void garWrite(const std::string &garfn, const expr & res);
extern void garWrite(const expr & res, const std::string &garfn);

/*
-----------------------------------------
    HepLib Wrapper
-----------------------------------------
*/

extern void set_Verbose(int v);
extern void set_Parallel_Process(int p);
extern expr LI(const int i);
extern expr TI(const int i);
extern expr DI(const int i);
extern expr CI(const int i);
extern expr LI(const expr& i);
extern expr TI(const expr& i);
extern expr DI(const expr& i);
extern expr CI(const expr& i);
extern expr RLI(const int i);
extern expr RTI(const int i);
extern expr RDI(const int i);
extern expr RCI(const int i);
extern expr RLI(const expr& i);
extern expr RTI(const expr& i);
extern expr RDI(const expr& i);
extern expr RCI(const expr& i);
extern expr IndexL2R(const expr &e);

extern expr Index(const std::string &s);
extern expr IndexCA(const std::string &s);
extern expr IndexCF(const std::string &s);
extern expr Vector(const std::string &s);
extern expr symbol(const std::string &s);
extern expr Symbol(const std::string &s);
extern expr SP(const expr &e);
extern expr SP(const expr &e1, const expr &e2);
extern expr GAS(const expr &e);
extern expr GAS(const int &i);
extern expr TR(const expr &e);
extern expr TTR(const std::vector<expr> &ev);
extern expr SUNT(const expr &e, const expr &i, const expr &j);
extern expr SUNT(const std::vector<expr> &ev, const expr &i, const expr &j);
extern expr SUNF(const expr &a, const expr &b, const expr &c);
extern expr SUNF4(const expr &a, const expr &b, const expr &c, const expr &d);
extern expr LC(const expr &a, const expr &b, const expr &c, const expr &d);
extern expr form(const expr &e, int verb=0);

extern void letSP(const expr &e, const expr &e2);
extern void letSP(const expr &e1, const expr &e2, const expr &e12);
extern expr call(const std::string func, const std::vector<expr> &ev);
extern expr call(const std::string func, const expr &e);
extern expr w(const int wi=0);

extern expr x(const int i);
extern expr y(const int i);
extern expr z(const int i);
extern expr x(const expr & i);
extern expr y(const expr & i);
extern expr z(const expr & i);
extern expr WRA(const expr &e);

extern expr WF(const expr& e);
extern expr WF(const expr& e1, const expr& e2);
extern expr WF(const expr& e1, const expr& e2, const expr& e3);
extern expr WF(const expr& e1, const expr& e2, const expr& e3, const expr& e4);
extern expr WF(const expr& e1, const expr& e2, const expr& e3, const expr& e4, const expr& e5);

extern void set_form_using_su3(bool yn);
extern void set_form_using_dim4(bool yn);

%warnfilter(509) FeynmanParameter;
class FeynmanParameter {
public:
    exvec Propagators;
    exvec Exponents;
    exvec LoopMomenta;
    exvec tLoopMomenta;
    exmap lReplacements;
    exmap tReplacements;
    exmap nReplacements;
    expr Prefactor = expr(1);
};

%warnfilter(509) XIntegrand;
class XIntegrand {
public:
    exvec Functions;
    exvec Exponents;
    exvec Deltas;
    exmap nReplacements;
};

// SecDec
%warnfilter(509) SecDec;
class SecDec {
public:
    int epN = 0;
    int epsN = 0;
    
    // used in Contours
    bool CTMaxF = true;
    expr CTLaMax = expr(10); // CTLaMax<0 for explict REAL mode
    int CTTryPTS = 3;
    int CTSavePTS = 3;
    
    long long TryPTS = 500000;
    long long LambdaSplit = 5;
    expr IntLaMax = expr(50);
    int CTry = 1;
    int CTryLeft = 1;
    int CTryRight = 1;
    expr CTryRightRatio = expr("1.5");
    int soLimit = 10000;
    
    long long RunMAX = 20;
    long long RunPTS = 500000;
    std::map<int, long long> MinPTS;
    expr EpsAbs = expr("1E-4");
    int ReIm = 3; // 1-Re, 2-Im, 3-ReIm
    
    SecDec();
    void Initialize(FeynmanParameter fp);
    void Initialize(XIntegrand xint);
    void Normalizes();
    void Scalelesses();
    void SDPrepares();
    void EpsEpExpands();
    void RemoveDeltas();
    void XReOrders();
    void XTogethers();
    void XExpands();
    void KillPowers(int bits=1+2);
    void CIPrepares(const std::string & key = "");
    void Contours(const std::string & key = "", const std::string & pkey = "");
    void Integrates(const std::string & key="", const std::string & pkey="", int kid=0);
    void ReIntegrates(const std::string & key, const std::string & pkey, expr err);
    void Evaluate(FeynmanParameter fp, const std::string & key = "");
    void Evaluate(XIntegrand xint, const std::string & key = "");
    void Evaluate(exvec FunExp, const std::string & key = "");
    void Evaluate(const std::string & key = "");
    void MB();
    void XEnd();
    void ChengWu(const expr & ft=expr(0));
    
    exvec FunExp;
    expr ResultError;
    expr VE;
};

// Process
%warnfilter(509) Process;
class Process {
public:
    static void DrawPDF(const std::vector<expr>, std::string);
    static std::string Style;
    std::string Model;
    std::string In;
    std::string Out;
    std::string LoopPrefix;
    int Loops;
    std::string Options;
    std::vector<std::string> Others;
    std::vector<expr> Amplitudes(std::map<std::string,expr> st, bool debug=false);
};
void set_LineTeX(expr, std::string);
void set_InOutTeX(int, std::string);

extern expr charge_conjugate(const expr &);
extern expr TIR(const expr &expr_in, const std::vector<expr> &loop_ps, const std::vector<expr> &ext_ps);
extern expr MatrixContract(const expr & expr_in);
extern expr Matrix(const expr & mat, const expr &i, const expr &j);
extern expr Apart(const expr &expr_in, const std::vector<expr> &vars, std::map<expr, expr, expr_is_less> sgnmap={});
extern expr Apart(const expr &expr_in, const std::vector<expr> &loops, const std::vector<expr> & extmoms, std::map<expr, expr, expr_is_less> sgnmap={});
extern expr ApartIR2ex(const expr & expr_in);
extern expr ApartIR2F(const expr & expr_in);
extern expr F2ex(const expr & expr_in);
extern expr ApartIRC(const expr & expr_in);

// FIRE
%warnfilter(509) FIRE;
class FIRE {
public:
    static int Version;
    static int Threads;
    
    bool reCut = false;
    std::string WorkingDir;
    int ProblemNumber = 0;
    
    exvec MIntegrals;
    exvec Rules;
    exvec Internal;
    exvec External;
    exvec Replacements;
    exvec Propagators;
    exvec Integrals; // lst of index lst
    exvec PIntegrals; // lst of index lst
    exvec Cuts; // index start from 1
    exvec DSP; // { {q1,q1}, {q1,p}, ... } Diff SP
    exvec ISP; // { q1*q1, q1*p } Independent SP
    std::map<int,expr> Shift;
    void Reduce();
};
