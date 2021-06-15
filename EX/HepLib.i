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

class StopIteration { };
%exception {
    try {
        $action
    } catch(const StopIteration & e) {
        PyErr_SetString(PyExc_StopIteration, "End of iterator");
        SWIG_fail;
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
    %template(expr_expr_pair) std::pair<expr,expr>;
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
    expr collect(const expr & s);
    expr subs(const std::vector<expr> &ev);
    expr subs(const expr &e);
    expr subs(const exmap & e);
    bool match(const expr &e);
    bool has(const expr &e);
    bool info(std::string sflags);
    expr map(MapFunction &mf);
    expr evalf() const;
    
    unsigned gethash();
    bool isSymbol();
    bool isVector();
    bool isIndex();
    bool isPair();
    bool isDGamma();
    
    std::string str();
    std::string __str__();
    it4expr __iter__();
    
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
extern expr evalf(const expr &e);
extern expr normal(const expr &e);
extern expr factor(const expr &e);
extern expr series(const expr &e, const expr &s, int o);
extern expr collect(const expr & e, const expr & s);
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

class it4expr {
public:
    it4expr* __iter__();
    expr & __next__();
};

class it4vec {
public:
    it4vec* __iter__();
    expr & __next__();
};

class it4map {
public:
    it4map* __iter__();
    std::pair<expr,expr> & __next__();
};

class it4set {
public:
    it4set* __iter__();
    expr & __next__();
};

class exvec {
public:
    exvec();
    exvec(const std::vector<expr> es);
    exvec(expr e);
    void push_back(expr e);
    int size();
    int nops();
    expr __getitem__(const int i);
    expr op(const int i);
    void __setitem__(const int i, expr v);
    void subs(const expr & e);
    void subs(const std::vector<expr> & e);
    void subs(const exmap &e);
    void sort();
    std::string str();
    std::string __str__();
    it4vec __iter__();
};

class exmap {
public:
    exmap();
    exmap(std::map<expr,expr,expr_is_less> es);
    expr __getitem__(expr e);
    void __setitem__(expr k, expr v);
    int size();
    int nops();
    std::string str();
    std::string __str__();
    it4map __iter__();
};

class exset {
public:
    exset();
    exset(const std::vector<expr> es);
    exset(expr e);
    void insert(expr e);
    int size();
    int nops();
    std::string str();
    std::string __str__();
    it4set __iter__();
};

class MapFunction {
public:
    MapFunction();
    virtual ~MapFunction();
    virtual expr map(const expr &e);
    expr operator() (const expr &e);
    exvec operator() (const exvec & ev);
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

extern void set_Verbose(int v);
extern void set_Parallel_Process(int p);
extern void set_Digits(int dn);

extern expr file2expr(std::string fn);
extern std::map<std::string,expr> garReadAll(const std::string &garfn);
extern expr garRead(const std::string &garfn, const char* key);
extern expr garRead(const std::string &garfn);
extern void garWrite(const std::string &garfn, const std::map<std::string, expr> &resMap);
extern void garWrite(const std::map<std::string, expr> &resMap, const std::string &garfn);
extern void garWrite(const std::string &garfn, const expr & res);
extern void garWrite(const expr & res, const std::string &garfn);
extern std::string RunOS(const std::string & cmd);

class cout {
public:
    cout & operator << (const expr &e);
    cout & operator << (const int &e);
    cout & operator << (const std::string &e);
    cout & operator << (const char * &e);
    cout & operator<<(const exvec &ev);
    cout & operator<<(const exmap &em);
    cout & operator<<(const exset &es);
    cout & operator<<(const std::vector<expr> &ev);
    cout & operator<<(const std::map<expr,expr,expr_is_less> &ev);
};

class hout {
public:
    hout & operator << (const expr &e);
    hout & operator << (const int &e);
    hout & operator << (const std::string &e);
    hout & operator << (const char * &e);
    hout & operator<<(const exvec &ev);
    hout & operator<<(const exmap &em);
    hout & operator<<(const exset &es);
    hout & operator<<(const std::vector<expr> &ev);
    hout & operator<<(const std::map<expr,expr,expr_is_less> &ev);
};

/*
-----------------------------------------
    HepLib Wrapper
-----------------------------------------
*/


// from HepLib namespace

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
extern expr wild(const int wi=0);

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

extern expr charge_conjugate(const expr &);
extern expr TIR(const expr &expr_in, const std::vector<expr> &loop_ps, const std::vector<expr> &ext_ps);
extern expr TIR(const expr &expr_in, const exvec & loop_ps, const exvec & ext_ps);
extern expr MatrixContract(const expr & expr_in);
extern expr Matrix(const expr & mat, const expr &i, const expr &j);
extern expr Apart(const expr &expr_in, const std::vector<expr> &vars, std::map<expr, expr, expr_is_less> sgnmap={});
extern expr Apart(const expr &expr_in, const std::vector<expr> &loops, const std::vector<expr> & extmoms, std::map<expr, expr, expr_is_less> sgnmap={});
extern expr ApartIR2ex(const expr & expr_in);
extern expr ApartIR2F(const expr & expr_in);
extern expr F2ex(const expr & expr_in);
extern expr ApartIRC(const expr & expr_in);

extern exvec ApartIBP(int IBPmethod, std::vector<expr> &io_vec, const std::vector<expr> & loops, const std::vector<expr> & exts, const std::vector<expr> & cut_props={});
extern exvec ApartIBP(int IBPmethod, const exvec &io_vec, const exvec & loops, const exvec & exts, const exvec & cut_props={});

// from HepLib::QGRAF namespace

// Process
%warnfilter(509) Process;
class Process {
public:
    static void DrawPDF(const exvec &, std::string);
    static std::string Style;
    std::string Model;
    std::string In;
    std::string Out;
    std::string LoopPrefix;
    int Loops;
    std::string Options;
    std::vector<std::string> Others;
    exvec Amplitudes(std::map<std::string,expr> st, bool debug=false);
};
extern void set_LineTeX(expr, std::string);
extern void set_InOutTeX(int, std::string);
extern exvec ShrinkCut(const expr & e, exvec ev, int n);

// Feynman Rules
extern expr QuarkPropagator(expr e, expr m=expr(0), bool color=true);
extern expr GluonPropagator(expr e, bool color=true);
extern expr GhostPropagator(expr e, bool color=true);
extern expr q2gVertex(expr e, bool color=true);
extern expr g3Vertex(expr e);
extern expr g4Vertex(expr e);
extern expr gh2gVertex(expr e, bool color=true);

extern expr IndexL2R(expr e, bool all=true);
extern expr IndexCC(expr e, bool all=true);

extern expr GluonSumL(int qi, bool color=true);
extern expr QuarkSumL(int qi, expr p, expr m, bool color=true);
extern expr AntiQuarkSumL(int qi, expr p, expr m, bool color=true);
extern expr GhostSumL(int qi);
extern expr AntiGhostSumL(int qi);
extern expr J1SumL(int qi, expr p);

// RC
class RC {
public:
    static expr Z2(std::string name, expr m, int loop=2);
    static expr Zm(expr m, int loop=2);
    static expr asBare(int loop=2);
    static expr asLO();
    static expr Zas(int loop);
};

// QCD
extern expr SpinProj(std::string io, int s, expr p, expr pb, expr m, expr e, expr mu);
extern expr SpinProj(std::string io, int s, expr p, expr pb, expr m, expr e, expr mb, expr eb, expr mu);
extern expr SpinProj(std::string io, int s, expr p, expr pb, expr m, expr e, expr mu, int i, int j);
extern expr SpinProj(std::string io, int s, expr p, expr pb, expr m, expr e, expr mb, expr eb, expr mu, int i, int j);
extern expr ColorProj(int i, int j, expr a);
extern expr ColorProj(int i, int j);
extern expr S1L1Proj(expr si, expr qi, expr p);
extern expr S1L1Proj(expr si, expr qi, expr mu, expr p);
extern expr S1L1Proj(expr si, expr qi, expr mu1, expr mu2, expr p);
extern expr S1L2Proj(expr si, expr qi1, expr qi2, expr mu, expr p);
extern expr S1L2Proj(expr si, expr qi1, expr qi2, expr mu1, expr mu2, expr p);
extern expr S1L1Sum(expr si, expr siR, expr qi, expr qiR, expr p, int J);
extern expr LProj(const expr &expr_in, const exvec &pqi, std::string prefix="lpj");

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

// from HepLib::SD namespace

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


// global variables in python

%pythoncode %{
    
ep = expr("ep")
eps = expr("eps")
D = expr("D")
d = expr("d")
iEpsilon = expr("iEpsilon")

I = expr("I")
Pi = expr("Pi")
Euler = expr("Euler")

NA = expr("NA")
NF = expr("NF")
gs = expr("gs")
As = expr("as")
mu = expr("mu")
nL = expr("nL")
nH = expr("nH")

w = wild()
w0 = wild(0)
w1 = wild(1)
w2 = wild(2)
w3 = wild(3)
w4 = wild(4)
w5 = wild(5)
w6 = wild(6)
w7 = wild(7)
w8 = wild(8)
w9 = wild(9)

co = cout()
cout = co
ho = hout()
hout = ho

endl = '\n'
RESET = '\033[0m'
BLACK =  '\033[30m'
RED =  '\033[31m'
GREEN =  '\033[32m'
YELLOW =  '\033[33m'
BLUE =  '\033[34m'
MAGENTA = '\033[35m'
CYAN = '\033[36m'
WHITE = '\033[37m'
BOLDBLACK = '\033[1m\033[30m'
BOLDRED = '\033[1m\033[31m'
BOLDGREEN = '\033[1m\033[32m'
BOLDYELLOW = '\033[1m\033[33m'
BOLDBLUE = '\033[1m\033[34m'
BOLDMAGENTA = '\033[1m\033[35m'
BOLDCYAN = '\033[1m\033[36m'
BOLDWHITE = '\033[1m\033[37m'

%}
