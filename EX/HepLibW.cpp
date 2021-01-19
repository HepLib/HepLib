/**
 * @file
 * @brief Wrap Class for SWIG
 */

#include "HepLibW.h"

expr::expr() { _expr = 0; }
expr::expr(int i) { _expr = i; }
expr::expr(GiNaC::ex e) { _expr = e; }
expr::expr(const std::string &s) { _expr = HepLib::str2ex(s); }
    
expr expr::operator+(const expr &e) { return expr(_expr + e._expr); }
expr expr::operator-(const expr &e) { return expr(_expr - e._expr); }
expr expr::operator*(const expr &e) { return expr(_expr * e._expr); }
expr expr::operator/(const expr &e) { return expr(_expr / e._expr); }
expr expr::operator==(const expr &e) { return expr(_expr == e._expr); }

expr expr::operator+(const int i) { return expr(_expr + i); }
expr expr::operator-(const int i) { return expr(_expr - i); }
expr expr::operator*(const int i) { return expr(_expr * i); }
expr expr::operator/(const int i) { return expr(_expr / GiNaC::ex(i)); }

expr expr::operator-() { return expr(-_expr); }

std::string expr::str() { return HepLib::ex2str(_expr); }
std::string expr::__str__() { return HepLib::ex2str(_expr); }

unsigned int expr::nops() { return _expr.nops(); }
expr expr::op(unsigned int i) { return expr(_expr.op(i)); }
void expr::let_op(unsigned int i, expr e) { _expr.let_op(i) = e._expr; }

expr expr::subs(const std::vector<expr> &ev) {
    GiNaC::lst repl;
    for(auto item : ev) repl.append(item._expr);
    return expr(_expr.subs(repl));
}
expr expr::subs(const expr &e) {
    return expr(_expr.subs(e._expr));
}

bool expr::match(const expr &e) {
    return _expr.match(e._expr);
}

bool expr::isSymbol() { return GiNaC::is_a<HepLib::Symbol>(_expr); }
bool expr::isVector() { return GiNaC::is_a<HepLib::Vector>(_expr); }
bool expr::isIndex() { return GiNaC::is_a<HepLib::Index>(_expr); }
bool expr::isPair() { return GiNaC::is_a<HepLib::Pair>(_expr); }
bool expr::isDGamma() { return GiNaC::is_a<HepLib::DGamma>(_expr); }

bool expr::info(std::string sflags) {
    if (sflags == "numeric") return _expr.info(GiNaC::info_flags::numeric);
    else if (sflags == "real") return _expr.info(GiNaC::info_flags::real);
    else if (sflags == "rational") return _expr.info(GiNaC::info_flags::rational);
    else if (sflags == "integer") return _expr.info(GiNaC::info_flags::integer);
    else if (sflags == "crational") return _expr.info(GiNaC::info_flags::crational);
    else if (sflags == "cinteger") return _expr.info(GiNaC::info_flags::cinteger);
    else if (sflags == "positive") return _expr.info(GiNaC::info_flags::positive);
    else if (sflags == "negative") return _expr.info(GiNaC::info_flags::negative);
    else if (sflags == "nonnegative") return _expr.info(GiNaC::info_flags::nonnegative);
    else if (sflags == "posint") return _expr.info(GiNaC::info_flags::posint);
    else if (sflags == "negint") return _expr.info(GiNaC::info_flags::negint);
    else if (sflags == "nonnegint") return _expr.info(GiNaC::info_flags::nonnegint);
    else if (sflags == "even") return _expr.info(GiNaC::info_flags::even);
    else if (sflags == "odd") return _expr.info(GiNaC::info_flags::odd);
    else if (sflags == "prime") return _expr.info(GiNaC::info_flags::prime);
    else if (sflags == "relation") return _expr.info(GiNaC::info_flags::relation);
    else if (sflags == "relation_equal") return _expr.info(GiNaC::info_flags::relation_equal);
    else if (sflags == "relation_not_equal") return _expr.info(GiNaC::info_flags::relation_not_equal);
    else if (sflags == "relation_less") return _expr.info(GiNaC::info_flags::relation_less);
    else if (sflags == "relation_less_or_equal") return _expr.info(GiNaC::info_flags::relation_less_or_equal);
    else if (sflags == "relation_greater") return _expr.info(GiNaC::info_flags::relation_greater);
    else if (sflags == "relation_greater_or_equal") return _expr.info(GiNaC::info_flags::relation_greater_or_equal);
    else if (sflags == "symbol") return _expr.info(GiNaC::info_flags::symbol);
    else if (sflags == "list") return _expr.info(GiNaC::info_flags::list);
    else if (sflags == "exprseq") return _expr.info(GiNaC::info_flags::exprseq);
    else if (sflags == "polynomial") return _expr.info(GiNaC::info_flags::polynomial);
    else if (sflags == "integer_polynomial") return _expr.info(GiNaC::info_flags::integer_polynomial);
    else if (sflags == "cinteger_polynomial") return _expr.info(GiNaC::info_flags::cinteger_polynomial);
    else if (sflags == "rational_polynomial") return _expr.info(GiNaC::info_flags::rational_polynomial);
    else if (sflags == "crational_polynomial") return _expr.info(GiNaC::info_flags::crational_polynomial);
    else if (sflags == "rational_function") return _expr.info(GiNaC::info_flags::rational_function);
    else if (sflags == "has_indices") return _expr.info(GiNaC::info_flags::has_indices);
    else if (sflags == "idx") return _expr.info(GiNaC::info_flags::idx);
    else if (sflags == "expanded") return _expr.info(GiNaC::info_flags::expanded);
    else if (sflags == "indefinite") return _expr.info(GiNaC::info_flags::indefinite);
    return false;
}

expr expr::map(MapFunction &mf) {
    return expr(_expr.map(mf._map));
}

expr expr::expand() {
    return expr(GiNaC::ex(_expr.expand()));
}

expr expr::normal() {
    return expr(GiNaC::ex(HepLib::normal_fermat(_expr)));
}

expr expr::factor() {
    return expr(GiNaC::ex(HepLib::factor_form(_expr)));
}

expr expr::series(const expr &s, int o) {
    if(!GiNaC::is_a<HepLib::Symbol>(s._expr)) throw HepLib::Error("1st argument should be a Symbol.");
    return expr(HepLib::mma_series(_expr, GiNaC::ex_to<HepLib::Symbol>(s._expr), o));
}

expr expand(const expr &e) {
    return expr(e._expr.expand());
}

expr normal(const expr &e) {
    return expr(HepLib::normal_fermat(e._expr));
}

expr factor(const expr &e) {
    return expr(HepLib::factor_form(e._expr));
}

expr subs(const expr &e, const std::vector<expr> &ev) {
    GiNaC::lst repl;
    for(auto item : ev) repl.append(item._expr);
    return expr(e._expr.subs(repl));
}
expr subs(const expr &e1, const expr &e2) {
    return expr(e1._expr.subs(e2._expr));
}

expr series(const expr &e, const expr &s, int o) {
    if(!GiNaC::is_a<HepLib::Symbol>(s._expr)) throw HepLib::Error("1st argument should be a Symbol.");
    return expr(HepLib::mma_series(e._expr, GiNaC::ex_to<HepLib::Symbol>(s._expr), o));
}

expr pow(const expr &e1, const expr &e2) {
    return expr(GiNaC::ex(GiNaC::pow(e1._expr, e2._expr)));
}

expr pow(const expr &e, const int n) {
    return expr(GiNaC::ex(GiNaC::pow(e._expr, n)));
}

expr Symbol(const std::string &s) {
    return expr(GiNaC::ex(HepLib::Symbol(s)));
}

expr Index(const std::string &s) {
    return expr(GiNaC::ex(HepLib::Index(s)));
}

expr Vector(const std::string &s) {
    return expr(GiNaC::ex(HepLib::Vector(s)));
}

expr SP(const expr &e1, const expr &e2) {
    return expr(GiNaC::ex(HepLib::SP(e1._expr, e2._expr)));
}

expr GAS(const expr &e) {
    return expr(GiNaC::ex(HepLib::GAS(e._expr)));
}
        
expr TR(const expr &e) {
    return expr(GiNaC::ex(HepLib::TR(e._expr)));
}

expr form(const expr &e) {
    return expr(HepLib::form(e._expr));
}

void letSP(const expr &e1, const expr &e2, const expr &e12) {
    HepLib::letSP(e1._expr, e2._expr) = e12._expr;
}

expr call(const std::string func, const std::vector<expr> &ev) {
    int ttl = ev.size();
    std::ostringstream ss;
    GiNaC::symtab st;
    ss << func << "(";
    for(int i=0; i<ttl; i++) {
        ss << "p" << i;
        if(i+1==ttl) ss << ")";
        else ss << ",";
        st["p"+std::to_string(i)] = ev[i]._expr;
    }
    return expr(HepLib::str2ex(ss.str(), st));
}
extern expr call(const std::string func, const expr &e) {
    GiNaC::symtab st;
    st["p1"] = e._expr;
    return expr(HepLib::str2ex(func+"(p1)", st));
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
