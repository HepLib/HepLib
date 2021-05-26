/**
 * @file
 * @brief Wrap Class for SWIG
 */

#include "HepLibW.h"

expr LI(const int i) { return expr(HepLib::QGRAF::LI(i)); }
expr TI(const int i) { return expr(HepLib::QGRAF::TI(i)); }
expr DI(const int i) { return expr(HepLib::QGRAF::DI(i)); }
expr CI(const int i) { return expr(HepLib::QGRAF::CI(i)); }

expr Symbol(const std::string &s) {
    return expr(GiNaC::ex(HepLib::Symbol(s)));
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

expr SUNF(const expr &a, const expr &b, const expr &c) {
    return expr(HepLib::SUNF(a._expr, b._expr, c._expr));
}

expr SUNF4(const expr &a, const expr &b, const expr &c, const expr &d) {
    return expr(HepLib::SUNF4(a._expr, b._expr, c._expr, d._expr));
}

expr LC(const expr &a, const expr &b, const expr &c, const expr &d) {
    return expr(HepLib::LC(a._expr, b._expr, c._expr, d._expr));
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

expr form(const expr &e, int verb) {
    return expr(HepLib::form(e._expr, verb));
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

void set_form_using_su3(bool yn) {
    HepLib::form_using_su3 = yn;
}

MapFunction::MapFunction() : _map([this](const GiNaC::ex &e, HepLib::MapFunction &self)->GiNaC::ex{ return map(expr(e))._expr; }) {}
MapFunction::~MapFunction() { }
expr MapFunction::map(const expr & e) { return e; }
expr MapFunction::operator() (const expr &e) { return _map(e._expr); }

Integral::Integral() { }

void Integral::Functions(const std::vector<expr> &ev) {
    _Propagators_Functions.clear();
    for(auto item : ev) _Propagators_Functions.push_back(item._expr);
    isX = true;
}

void Integral::Propagators(const std::vector<expr> &ev, const std::vector<expr> &loops, const std::vector<expr> &tloops) {
    _Propagators_Functions.clear();
    for(auto item : ev) _Propagators_Functions.push_back(item._expr);
    _Loops.clear();
    for(auto item : loops) _Loops.push_back(item._expr);
    _tLoops.clear();
    for(auto item : tloops) _tLoops.push_back(item._expr);
    isX = false;
}

void Integral::Exponents(const std::vector<expr> &ev) {
    _Exponents.clear();
    for(auto item : ev) _Exponents.push_back(item._expr);
}

void Integral::Exponents(const std::vector<int> &ev) {
    _Exponents.clear();
    for(auto item : ev) _Exponents.push_back(item);
}

void Integral::Replacements(const std::vector<expr> &ev, int lt) {
    if(lt==0) {
        _lReplacements.clear();
        for(auto item : ev) _lReplacements.push_back(item._expr);
    } else if(lt==1) {
        _tReplacements.clear();
        for(auto item : ev) _tReplacements.push_back(item._expr);
    }
}

void Integral::Evaluate() {
    if(isX) {
        HepLib::SD::XIntegrand xint;
        xint.Functions = HepLib::vec2lst(_Propagators_Functions);
        xint.Exponents = HepLib::vec2lst(_Exponents);
            
        HepLib::SD::SecDec secdec;
        HepLib::Verbose = verb;
        secdec.epN = epN;
        secdec.epsN = epsN;
        secdec.Evaluate(xint);
        _Result = secdec.ResultError;
    } else {
        HepLib::SD::FeynmanParameter fp;
        fp.Propagators = HepLib::vec2lst(_Propagators_Functions);
        fp.Exponents = HepLib::vec2lst(_Exponents);
        fp.LoopMomenta = HepLib::vec2lst(_Loops);
        fp.tLoopMomenta = HepLib::vec2lst(_tLoops);
        for(auto item : _lReplacements) fp.lReplacements[item.op(0)] = item.op(1);
        for(auto item : _lReplacements) fp.tReplacements[item.op(0)] = item.op(1);
            
        HepLib::SD::SecDec secdec;
        HepLib::Verbose = verb;
        secdec.epN = epN;
        secdec.epsN = epsN;
        secdec.Evaluate(fp);
        _Result = secdec.ResultError;
    }
}

std::string Integral::str() { return HepLib::ex2str(HepLib::SD::VEResult(_Result)); }
std::string Integral::__str__() { return HepLib::ex2str(HepLib::SD::VEResult(_Result)); }

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

Function::Function() { }
Function::~Function() { }
expr Function::operator()(const expr &e) { return expr(0); }
expr Function::operator()(const expr &e1,const expr &e2) { return expr(0); }
expr Function::operator()(const expr &e1,const expr &e2,const expr &e3) { return expr(0); }
expr Function::operator()(const expr &e1,const expr &e2,const expr &e3,const expr &e4) { return expr(0); }
expr Function::operator()(const expr &e1,const expr &e2,const expr &e3,const expr &e4,const expr &e5) { return expr(0); }

expr charge_conjugate(const expr &e) { return HepLib::charge_conjugate(e._expr); }

expr TIR(const expr &expr_in, const std::vector<expr> &loop_ps, const std::vector<expr> &ext_ps) {
    GiNaC::lst lps, eps;
    for(auto item : loop_ps) lps.append(item._expr);
    for(auto item : ext_ps) eps.append(item._expr);
    return HepLib::TIR(expr_in._expr, lps, eps);
}

expr MatrixContract(const expr & e) { return HepLib::MatrixContract(e._expr); }

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
