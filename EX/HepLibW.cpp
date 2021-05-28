/**
 * @file
 * @brief Wrap Class for SWIG
 */

#include "HepLibW.h"

expr LI(const int i) { return expr(HepLib::QGRAF::LI(i)); }
expr TI(const int i) { return expr(HepLib::QGRAF::TI(i)); }
expr DI(const int i) { return expr(HepLib::QGRAF::DI(i)); }
expr CI(const int i) { return expr(HepLib::QGRAF::CI(i)); }
expr LI(const expr& i) { return expr(HepLib::QGRAF::LI(i._expr)); }
expr TI(const expr& i) { return expr(HepLib::QGRAF::TI(i._expr)); }
expr DI(const expr& i) { return expr(HepLib::QGRAF::DI(i._expr)); }
expr CI(const expr& i) { return expr(HepLib::QGRAF::CI(i._expr)); }
expr RLI(const int i) { return expr(HepLib::QGRAF::RLI(i)); }
expr RTI(const int i) { return expr(HepLib::QGRAF::RTI(i)); }
expr RDI(const int i) { return expr(HepLib::QGRAF::RDI(i)); }
expr RCI(const int i) { return expr(HepLib::QGRAF::RCI(i)); }
expr RLI(const expr& i) { return expr(HepLib::QGRAF::RLI(i._expr)); }
expr RTI(const expr& i) { return expr(HepLib::QGRAF::RTI(i._expr)); }
expr RDI(const expr& i) { return expr(HepLib::QGRAF::RDI(i._expr)); }
expr RCI(const expr& i) { return expr(HepLib::QGRAF::RCI(i._expr)); }
expr IndexL2R(const expr &e) { return expr(HepLib::QGRAF::IndexL2R(e._expr)); }

expr Symbol(const std::string &s) {
    return expr(GiNaC::ex(HepLib::Symbol(s)));
}

expr symbol(const std::string &s) {
    return expr(GiNaC::ex(GiNaC::symbol(s)));
}

expr Index(const std::string &s) {
    return expr(GiNaC::ex(HepLib::Index(s)));
}

expr IndexCA(const std::string &s) {
    return expr(GiNaC::ex(HepLib::Index(s,HepLib::Index::Type::CA)));
}

expr IndexCF(const std::string &s) {
    return expr(GiNaC::ex(HepLib::Index(s,HepLib::Index::Type::CF)));
}

expr Vector(const std::string &s) {
    return expr(GiNaC::ex(HepLib::Vector(s)));
}

expr SUNT(const expr &e, const expr &i, const expr &j) {
    return expr(HepLib::SUNT(e._expr, i._expr, j._expr));
}

expr SUNT(const std::vector<expr> &ev, const expr &i, const expr &j) {
    GiNaC::lst es;
    for(auto item : ev) es.append(item._expr);
    return expr(HepLib::SUNT(es, i._expr, j._expr));
}

expr SUNF(const expr &a, const expr &b, const expr &c) {
    return expr(HepLib::SUNF(a._expr, b._expr, c._expr));
}

expr SUNF4(const expr &a, const expr &b, const expr &c, const expr &d) {
    return expr(HepLib::SUNF4(a._expr, b._expr, c._expr, d._expr));
}

expr LC(const expr &a, const expr &b, const expr &c, const expr &d) {
    return expr(HepLib::LC(a._expr, b._expr, c._expr, d._expr));
}

expr SP(const expr &e) {
    return expr(GiNaC::ex(HepLib::SP(e._expr)));
}

expr SP(const expr &e1, const expr &e2) {
    return expr(GiNaC::ex(HepLib::SP(e1._expr, e2._expr)));
}

expr GAS(const expr &e) {
    return expr(GiNaC::ex(HepLib::GAS(e._expr)));
}

expr GAS(const int &i) {
    return expr(GiNaC::ex(HepLib::GAS(i)));
}
        
expr TR(const expr &e) {
    return expr(GiNaC::ex(HepLib::TR(e._expr)));
}

expr TTR(const std::vector<expr> &ev) {
    GiNaC::lst alst;
    for(auto item : ev) alst.append(item._expr);
    return expr(HepLib::TTR(alst));
}

expr form(const expr &e, int verb) {
    return expr(HepLib::form(e._expr, verb));
}

void letSP(const expr &e, const expr &e2) {
    HepLib::letSP(e._expr, e._expr) = e2._expr;
}

void letSP(const expr &e1, const expr &e2, const expr &e12) {
    HepLib::letSP(e1._expr, e2._expr) = e12._expr;
}

expr w(const int wi) {
    return expr(GiNaC::wild(wi));
}

expr x(const int i) {
    return expr(GiNaC::ex(HepLib::x(i)));
}

expr y(const int i) {
    return expr(GiNaC::ex(HepLib::y(i)));
}

expr z(const int i) {
    return expr(GiNaC::ex(HepLib::z(i)));
}

expr x(const expr & i) {
    return expr(GiNaC::ex(HepLib::x(i._expr)));
}

expr y(const expr & i) {
    return expr(GiNaC::ex(HepLib::y(i._expr)));
}

expr z(const expr & i) {
    return expr(GiNaC::ex(HepLib::z(i._expr)));
}

expr WRA(const expr &e) {
    return expr(GiNaC::ex(HepLib::SD::WRA(e._expr)));
}

expr WF(const expr& e) {
    return expr(HepLib::WF(e._expr));
}
extern expr WF(const expr& e1, const expr& e2) {
    return expr(HepLib::WF(e1._expr,e2._expr));
}
extern expr WF(const expr& e1, const expr& e2, const expr& e3) {
    return expr(HepLib::WF(e1._expr,e2._expr,e3._expr));
}
extern expr WF(const expr& e1, const expr& e2, const expr& e3, const expr& e4) {
    return expr(HepLib::WF(e1._expr,e2._expr,e3._expr,e4._expr));
}
extern expr WF(const expr& e1, const expr& e2, const expr& e3, const expr& e4, const expr& e5) {
    return expr(HepLib::WF(e1._expr,e2._expr,e3._expr,e4._expr,e5._expr));
}

void set_form_using_su3(bool yn) {
    HepLib::form_using_su3 = yn;
}

void set_form_using_dim4(bool yn) {
    HepLib::form_using_dim4 = yn;
}


MapFunction::MapFunction() : _map([this](const GiNaC::ex &e, HepLib::MapFunction &self)->GiNaC::ex{ return map(expr(e))._expr; }) {}
MapFunction::~MapFunction() { }
expr MapFunction::map(const expr & e) { return e; }
expr MapFunction::operator() (const expr &e) { return _map(e._expr); }

std::string Process::Style;
std::vector<expr> Process::Amplitudes(std::map<std::string,expr> st, bool debug) {
    HepLib::QGRAF::Process proc;
    if(Style != "") HepLib::QGRAF::Process::Style = Style;
    proc.Model = Model;
    proc.In = In;
    proc.Out = Out;
    proc.LoopPrefix = LoopPrefix;
    proc.Loops = Loops;
    proc.Options = Options;
    proc.Others = Others;
    GiNaC::symtab _st;
    for(auto kv : st) _st[kv.first] = kv.second._expr;
    auto ret = proc.Amplitudes(_st, debug);
    std::vector<expr> res;
    for(auto item : ret) res.push_back(expr(item));
    return res;
}

void Process::DrawPDF(std::vector<expr> amps, std::string fn) {
    GiNaC::lst amp_lst;
    for(auto item : amps) amp_lst.append(item._expr);
    HepLib::QGRAF::DrawPDF(amp_lst, fn);
}

void set_LineTeX(expr f, std::string tex) {
    HepLib::QGRAF::LineTeX[f._expr] = tex;
}
void set_InOutTeX(int fn, std::string tex) {
    HepLib::QGRAF::InOutTeX[fn] = tex;
}

expr charge_conjugate(const expr &e) { return HepLib::charge_conjugate(e._expr); }

expr TIR(const expr &expr_in, const std::vector<expr> &loop_ps, const std::vector<expr> &ext_ps) {
    GiNaC::lst lps, eps;
    for(auto item : loop_ps) lps.append(item._expr);
    for(auto item : ext_ps) eps.append(item._expr);
    return HepLib::TIR(expr_in._expr, lps, eps);
}

expr MatrixContract(const expr & e) {
    return expr(HepLib::MatrixContract(e._expr));
}
expr Matrix(const expr & mat, const expr &i, const expr &j) {
    return expr(HepLib::Matrix(mat._expr, i._expr, j._expr));
}

expr Apart(const expr &expr_in, const std::vector<expr> &vars, std::map<expr, expr, expr_is_less> sgnmap) {
    GiNaC::lst _vars;
    for(auto item : vars) _vars.append(item._expr);
    GiNaC::exmap smap;
    for(auto kv : sgnmap) sgnmap[kv.first._expr] = kv.second._expr;
    return HepLib::Apart(expr_in._expr, _vars, smap);
}

expr Apart(const expr &expr_in, const std::vector<expr> &loops, const std::vector<expr> & extmoms, std::map<expr, expr, expr_is_less> sgnmap) {
    GiNaC::lst _loops, _extmoms;
    for(auto item : loops) _loops.append(item._expr);
    for(auto item : extmoms) _extmoms.append(item._expr);
    GiNaC::exmap smap;
    for(auto kv : sgnmap) sgnmap[kv.first._expr] = kv.second._expr;
    return HepLib::Apart(expr_in._expr, _loops, _extmoms, smap);
}

expr ApartIR2ex(const expr & e) { return HepLib::ApartIR2ex(e._expr); }
expr ApartIR2F(const expr & e) { return HepLib::ApartIR2F(e._expr); }
expr F2ex(const expr & e) { return HepLib::F2ex(e._expr); }
expr ApartIRC(const expr & e) { return HepLib::ApartIRC(e._expr); }

int FIRE::Version = 6;
int FIRE::Threads = 2;

void FIRE::Reduce() {
    HepLib::IBP::FIRE::Version = FIRE::Version;
    HepLib::IBP::FIRE::Threads = FIRE::Threads;
    HepLib::IBP::FIRE fire;
    
    fire.Internal =  HepLib::vec2lst(Internal._g);
    fire.External =  HepLib::vec2lst(External._g);
    fire.Replacements =  HepLib::vec2lst(Replacements._g);
    fire.Propagators =  HepLib::vec2lst(Propagators._g);
    fire.Integrals =  HepLib::vec2lst(Integrals._g); // lst of index lst
    fire.PIntegrals =  HepLib::vec2lst(PIntegrals._g); // lst of index lst
    fire.Cuts =  HepLib::vec2lst(Cuts._g); // index start from 1
    fire.DSP =  HepLib::vec2lst(DSP._g); // { {q1,q1}, {q1,p}, ... } Diff SP
    fire.ISP =  HepLib::vec2lst(ISP._g); // { q1*q1, q1*p } Independent SP
    for(auto kv : Shift) fire.Shift[kv.first] = kv.second._expr;

    fire.reCut = reCut;
    fire.WorkingDir = WorkingDir;
    fire.ProblemNumber = ProblemNumber;
    fire.Reduce();
    for(auto it : fire.MIntegrals) MIntegrals.push_back(expr(it));
    for(auto it : fire.Rules) Rules.push_back(expr(it));
}

// Sector Decompostion

SecDec::SecDec() { }

void SecDec::Initialize(FeynmanParameter fp) {
    HepLib::SD::FeynmanParameter _fp;
    _fp.Propagators = HepLib::vec2lst(fp.Propagators._g);
    _fp.Exponents = HepLib::vec2lst(fp.Exponents._g);
    _fp.LoopMomenta = HepLib::vec2lst(fp.LoopMomenta._g);
    _fp.tLoopMomenta = HepLib::vec2lst(fp.tLoopMomenta._g);
    _fp.tReplacements = fp.tReplacements._g;
    _fp.lReplacements = fp.lReplacements._g;
    _fp.nReplacements = fp.nReplacements._g;
    _fp.Prefactor = fp.Prefactor._expr;
    _SecDec.Initialize(_fp);
    FunExp = exvec(_SecDec.FunExp);
}

void SecDec::Initialize(XIntegrand xint) {
    HepLib::SD::XIntegrand _xint;
    _xint.Functions = HepLib::vec2lst(xint.Functions._g);
    _xint.Exponents = HepLib::vec2lst(xint.Exponents._g);
    _xint.nReplacements = xint.nReplacements._g;
    _xint.Deltas = HepLib::vec2lst(xint.Deltas._g);
    _SecDec.Initialize(_xint);
    FunExp = exvec(_SecDec.FunExp);
}

void SecDec::Normalizes() { _SecDec.Normalizes(); }
void SecDec::Scalelesses() { _SecDec.Scalelesses(); }
void SecDec::SDPrepares() { _SecDec.SDPrepares(); }
void SecDec::EpsEpExpands() {
    _SecDec.epN = epN;
    _SecDec.epsN = epsN;
    _SecDec.EpsEpExpands();
}
void SecDec::RemoveDeltas() { _SecDec.RemoveDeltas(); }
void SecDec::XReOrders() { _SecDec.XReOrders(); }
void SecDec::XTogethers() { _SecDec.XTogethers(); }
void SecDec::XExpands() { _SecDec.XExpands(); }
void SecDec::KillPowers(int bits) { _SecDec.KillPowers(bits); }
void SecDec::CIPrepares(const std::string & key) { _SecDec.CIPrepares(key); }
void SecDec::Contours(const std::string & key, const std::string & pkey) {
    _SecDec.CTMaxF = CTMaxF;
    _SecDec.CTLaMax = HepLib::ex2q(CTLaMax._expr); // CTLaMax<0 for explict REAL mode
    _SecDec.CTTryPTS = CTTryPTS;
    _SecDec.CTSavePTS = CTSavePTS;
    _SecDec.Contours(key,pkey);
}
void SecDec::Integrates(const std::string & key, const std::string & pkey, int kid) {
    _SecDec.TryPTS = 500000;
    _SecDec.LambdaSplit = 5;
    _SecDec.IntLaMax = HepLib::ex2q(IntLaMax._expr);
    _SecDec.CTry = 1;
    _SecDec.CTryLeft = 1;
    _SecDec.CTryRight = 1;
    _SecDec.CTryRightRatio = HepLib::ex2q(CTryRightRatio._expr);
    _SecDec.soLimit = soLimit;
    
    _SecDec.RunMAX = RunMAX;
    _SecDec.RunPTS = RunPTS;
    for(auto kv : MinPTS) _SecDec.MinPTS[kv.first] = kv.second;
    _SecDec.EpsAbs = HepLib::ex2q(EpsAbs._expr);
    _SecDec.ReIm = ReIm; // 1-Re, 2-Im, 3-ReIm
    _SecDec.Integrates(key,pkey,kid);
    ResultError = expr(_SecDec.ResultError);
    VE = expr(HepLib::SD::VEResult(_SecDec.ResultError));
}
void SecDec::ReIntegrates(const std::string & key, const std::string & pkey, expr err) {
    _SecDec.ReIntegrates(key,pkey,HepLib::ex2q(err._expr));
}
void SecDec::MB() { _SecDec.MB(); }
void SecDec::XEnd() { _SecDec.XEnd(); }
void SecDec::ChengWu(const expr & ft) { _SecDec.ChengWu(ft._expr); }

void SecDec::Evaluate(XIntegrand xint, const std::string & key) {
    if(HepLib::Verbose>1) std::cout << std::endl << "  Starting @ " << HepLib::now() << std::endl;
    Initialize(xint);
    Evaluate(key);
}

void SecDec::Evaluate(exvec funexp, const std::string & key) {
    if(HepLib::Verbose>1) std::cout << std::endl << "  Starting @ " << HepLib::now() << std::endl;
    FunExp = funexp;
    _SecDec.FunExp = FunExp._g;
    Evaluate(key);
}

void SecDec::Evaluate(FeynmanParameter fp, const std::string & key) {
    if(HepLib::Verbose>1) std::cout << std::endl << "  Starting @ " << HepLib::now() << std::endl;
    Initialize(fp);
    Evaluate(key);
}

void SecDec::Evaluate(const std::string & key) {
    MB();
    if(_SecDec.FunExp.size()<1) return;
    Scalelesses();
    ChengWu();
    RemoveDeltas();
    KillPowers();
    SDPrepares();
    EpsEpExpands();
    CIPrepares(key);
    auto pps = HepLib::GiNaC_Parallel_Process;
    HepLib::GiNaC_Parallel_Process = 0;
    Contours(key);
    Integrates(key);
    HepLib::GiNaC_Parallel_Process = pps;
}
