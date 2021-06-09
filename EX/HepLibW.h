/**
 * @file
 * @brief Wrap Class for SWIG
 */

#include "HepLib.h"

class StopIteration { };
class it4expr;
class exmap;
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
    expr operator>>(const int &e);
    bool operator<(const expr &e) const;
    bool operator>(const expr &e) const;
    bool operator==(const expr &e) const;
    expr operator+(const int i);
    expr operator-(const int i);
    expr operator*(const int i);
    expr operator/(const int i);
    expr operator[](const int i);
    
    unsigned gethash();
    
    expr operator-();
    
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
    expr subs(const exmap &e);
    expr evalf();
    
    std::string str();
    std::string __str__();
    it4expr __iter__();
    
    bool match(const expr &e);
    bool has(const expr &e);
    bool isSymbol();
    bool isVector();
    bool isIndex();
    bool isPair();
    bool isDGamma();
    bool info(std::string sflags);
    expr map(MapFunction &mf);

};

struct expr_is_less {
    bool operator() (const expr &lh, const expr &rh) const { return lh._expr.compare(rh._expr) < 0; }
};

expr conjugate(const expr &e);
expr expand(const expr &e);
expr normal(const expr &e);
expr factor(const expr &e);
expr series(const expr &e, const expr &s, int o);
expr collect(const expr & e, const expr & s);
expr subs(const expr &e, const std::vector<expr> &ev);
expr subs(const expr &e1, const expr &e2);
expr subs(const expr &e1, const exmap &e2);

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

class it4expr {
public:
    it4expr() { }
    it4expr(GiNaC::const_iterator b, GiNaC::const_iterator e) : _c(b), _e(e) {  }
    it4expr* __iter__() { return this; }
    expr & __next__() {
        if (_c != _e) { e = expr(*(_c++)); return e; }
        else throw StopIteration();
    }
private:
    GiNaC::const_iterator _c;
    GiNaC::const_iterator _e;
    expr e;
};

class it4vec {
public:
    it4vec() { }
    it4vec(GiNaC::exvector::iterator b, GiNaC::exvector::iterator e) : _c(b), _e(e) {  }
    it4vec* __iter__() { return this; }
    expr & __next__() {
        if (_c != _e) { e = expr(*(_c++)); return e; }
        else throw StopIteration();
    }
private:
    std::vector<GiNaC::ex>::iterator _c;
    std::vector<GiNaC::ex>::iterator _e;
    expr e;
};

class exvec {
public:
    exvec();
    exvec(std::vector<expr> es);
    exvec(expr e);
    GiNaC::exvector _g;
    void push_back(expr e);
    expr __getitem__(const int i);
    expr op(const int i);
    void __setitem__(const int i, expr v);
    int size();
    int nops();
    void subs(const expr & e);
    void subs(const std::vector<expr> & e);
    void subs(const exmap & e);
    std::string str();
    std::string __str__();
    exvec(GiNaC::exvector es); // construct from C++
    it4vec __iter__();
    void sort();
};

class it4map {
public:
    it4map() { }
    it4map(GiNaC::exmap::iterator b, GiNaC::exmap::iterator e) : _c(b), _e(e) {  }
    it4map* __iter__() { return this; }
    std::pair<expr,expr> & __next__() {
        if (_c != _e) {
            auto e = *(_c++);
            kv = std::pair<expr,expr>(expr(e.first), expr(e.second));
            return kv;
        } else throw StopIteration();
    }
private:
    GiNaC::exmap::iterator _c;
    GiNaC::exmap::iterator _e;
    std::pair<expr,expr> kv;
};

class exmap {
public:
    exmap();
    exmap(std::map<expr,expr,expr_is_less> es);
    GiNaC::exmap _g;
    expr __getitem__(expr e);
    void __setitem__(expr k, expr v);
    int size();
    int nops();
    std::string str();
    std::string __str__();
    exmap(GiNaC::exmap es); // construct from C++
    it4map __iter__();
};

class it4set {
public:
    it4set() { }
    it4set(const GiNaC::exset::iterator b, const GiNaC::exset::iterator e) : _c(b), _e(e) {  }
    it4set* __iter__() { return this; }
    expr & __next__() {
        if (_c != _e) { e = expr(*(_c++)); return e; }
        else throw StopIteration();
    }
private:
    GiNaC::exset::iterator _c;
    GiNaC::exset::iterator _e;
    expr e;
};

class exset {
public:
    exset();
    exset(std::vector<expr> es);
    exset(expr e);
    GiNaC::exset _g;
    void insert(expr e);
    int size();
    int nops();
    std::string str();
    std::string __str__();
    exset(GiNaC::exset es); // construct from C++
    it4set __iter__();
};

class MapFunction {
public:
    MapFunction();
    virtual ~MapFunction();
    virtual expr map(const expr &e);
    expr operator() (const expr &e);
    exvec operator() (const exvec & ev);
    friend class expr;
private:
    HepLib::MapFunction _map;
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
exvec Parallel(int ntotal, int nbatch,
        ParFun &f,
        const std::string & key = "",
        bool rm = true,
        const std::string & pre = "  ");
exvec Parallel(int ntotal,
        ParFun &f,
        const std::string & key = "",
        bool rm = true,
        const std::string & pre = "  ");
        
void set_Parallel_Process(int p);
void set_Verbose(int v);
void set_Digits(int dn);
bool isFunction(const expr &e, std::string sf);

expr file2expr(std::string fn);
std::map<std::string,expr> garReadAll(const std::string &garfn);
expr garRead(const std::string &garfn, const char* key);
expr garRead(const std::string &garfn);
void garWrite(const std::string &garfn, const std::map<std::string, expr> &resMap);
void garWrite(const std::map<std::string, expr> &resMap, const std::string &garfn);
void garWrite(const std::string &garfn, const expr & res);
void garWrite(const expr & res, const std::string &garfn);

// --------------------------------------------------------

expr Index(const std::string &s);
expr IndexCA(const std::string &s);
expr IndexCF(const std::string &s);
expr Vector(const std::string &s);
expr Symbol(const std::string &s);
expr symbol(const std::string &s);
expr LI(const int i);
expr TI(const int i);
expr DI(const int i);
expr CI(const int i);
expr LI(const expr& i);
expr TI(const expr& i);
expr DI(const expr& i);
expr CI(const expr& i);
expr RLI(const int i);
expr RTI(const int i);
expr RDI(const int i);
expr RCI(const int i);
expr RLI(const expr& i);
expr RTI(const expr& i);
expr RDI(const expr& i);
expr RCI(const expr& i);

expr SUNT(const expr &e, const expr &i, const expr &j);
expr SUNT(const std::vector<expr> &ev, const expr &i, const expr &j);
expr SUNF(const expr &a, const expr &b, const expr &c);
expr SUNF4(const expr &a, const expr &b, const expr &c, const expr &d);
expr Eps(const expr &a, const expr &b, const expr &c, const expr &d);
expr LC(const expr &a, const expr &b, const expr &c, const expr &d);
expr SP(const expr &e);
expr SP(const expr &e1, const expr &e2);
expr GAS(const expr &e);
expr GAS(const int &i);
expr TR(const expr &e);
expr TTR(const std::vector<expr> &ev);
expr form(const expr &e, int verb=0);

expr charge_conjugate(const expr &);
expr TIR(const expr &expr_in, const std::vector<expr> &loop_ps, const std::vector<expr> &ext_ps);
expr TIR(const expr &expr_in, const exvec & loop_ps, const exvec & ext_ps);
expr Matrix(const expr & mat, const expr &i, const expr &j);
expr MatrixContract(const expr & expr_in);
expr Apart(const expr &expr_in, const std::vector<expr> &vars, std::map<expr, expr, expr_is_less> sgnmap={});
expr Apart(const expr &expr_in, const std::vector<expr> &loops, const std::vector<expr> & extmoms, std::map<expr, expr, expr_is_less> sgnmap={});
expr ApartIR2ex(const expr & expr_in);
expr ApartIR2F(const expr & expr_in);
expr F2ex(const expr & expr_in);
expr ApartIRC(const expr & expr_in);
exvec ApartIBP(int IBPmethod, std::vector<expr> &io_vec, const std::vector<expr> & loops, const std::vector<expr> & exts, const std::vector<expr> & cut_props={});
exvec ApartIBP(int IBPmethod, const exvec &io_vec, const exvec & loops, const exvec & exts, const exvec & cut_props={});
        
class AIOption {
    expr Internal; // Internal for Apart/IBP
    expr External; // External for Apart/IBP
    expr DSP; // DSP for IBP
    exmap smap; // Sign Map for Apart
    exvec Cuts; // Cut Propagators. optional
    exvec CSP; // SP in Cuts, to be cleared. optional
    exvec ISP; // SP for IBP. optional
    bool CutFirst = true;
    int mcl = 1; // collect_ex level, 0-nothing, 1-exnormal, 2-exfactor
};
void ApartIBP(int IBPmethod, std::vector<expr> &io_vec, AIOption aip);

expr x(const int i);
expr y(const int i);
expr z(const int i);
expr x(const expr & i);
expr y(const expr & i);
expr z(const expr & i);
expr WRA(const expr &e);
expr lst(const std::vector<expr> &ev);
void set_form_using_su3(bool yn);
void set_form_using_dim4(bool yn);
expr WF(const expr& e);
expr WF(const expr& e1, const expr& e2);
expr WF(const expr& e1, const expr& e2, const expr& e3);
expr WF(const expr& e1, const expr& e2, const expr& e3, const expr& e4);
expr WF(const expr& e1, const expr& e2, const expr& e3, const expr& e4, const expr& e5);

void letSP(const expr &e, const expr &e2);
void letSP(const expr &e1, const expr &e2, const expr &e12);

expr call(const std::string func, const std::vector<expr> &ev);
expr call(const std::string func, const expr &e);
expr wild(const int wi=0);

// --------------------------------------------------------

class Process {
public:
    static std::string Style;
    static void DrawPDF(const exvec &, std::string);
    std::string Model;
    std::string In;
    std::string Out;
    std::string LoopPrefix = "q";
    int Loops;
    std::string Options;
    std::vector<std::string> Others;
    exvec Amplitudes(std::map<std::string, expr>, bool debug=false);
};
void set_LineTeX(expr, std::string);
void set_InOutTeX(int, std::string);
exvec ShrinkCut(const expr & e, exvec ev, int n);

expr QuarkPropagator(expr e, expr m=expr(0), bool color=true);
expr GluonPropagator(expr e, bool color=true);
expr GhostPropagator(expr e, bool color=true);
expr q2gVertex(expr e, bool color=true);
expr g3Vertex(expr e);
expr g4Vertex(expr e);
expr gh2gVertex(expr e, bool color=true);

expr IndexL2R(expr e, bool all=true);
expr IndexCC(expr e, bool all=true);

expr GluonSumL(int qi, bool color=true);
expr QuarkSumL(int qi, expr p, expr m, bool color=true);
expr AntiQuarkSumL(int qi, expr p, expr m, bool color=true);
expr GhostSumL(int qi);
expr AntiGhostSumL(int qi);
expr J1SumL(int qi, expr p);

expr SpinProj(std::string io, int s, expr p, expr pb, expr m, expr e, expr mu);
expr SpinProj(std::string io, int s, expr p, expr pb, expr m, expr e, expr mb, expr eb, expr mu);
expr SpinProj(std::string io, int s, expr p, expr pb, expr m, expr e, expr mu, int i, int j);
expr SpinProj(std::string io, int s, expr p, expr pb, expr m, expr e, expr mb, expr eb, expr mu, int i, int j);
expr ColorProj(int i, int j, expr a);
expr ColorProj(int i, int j);
expr S1L1Proj(expr si, expr qi, expr p);
expr S1L1Proj(expr si, expr qi, expr mu, expr p);
expr S1L1Proj(expr si, expr qi, expr mu1, expr mu2, expr p);
expr S1L2Proj(expr si, expr qi1, expr qi2, expr mu, expr p);
expr S1L2Proj(expr si, expr qi1, expr qi2, expr mu1, expr mu2, expr p);
expr S1L1Sum(expr si, expr siR, expr qi, expr qiR, expr p, int J);
expr LProj(const expr &expr_in, const exvec &pqi, std::string prefix="lpj");

class RC {
public:
    static expr Z2(std::string name, expr m, int loop=2);
    static expr Zm(expr m, int loop=2);
    static expr asBare(int loop=2);
    static expr asLO();
    static expr Zas(int loop);
};

// --------------------------------------------------------

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

class XIntegrand {
public:
    exvec Functions;
    exvec Exponents;
    exvec Deltas;
    exmap nReplacements;
};

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
        
private:
    HepLib::SD::SecDec _SecDec;
};

// --------------------------------------------------------

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
