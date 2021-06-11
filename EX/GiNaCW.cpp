/**
 * @file
 * @brief Wrap Class for SWIG
 */

#include "HepLibW.h"

void set_Parallel_Process(int p) { HepLib::GiNaC_Parallel_Process = p; }
void set_Verbose(int v) { HepLib::Verbose = v; }
void set_Digits(int dn) { GiNaC::Digits = dn; }

expr::expr() { _expr = 0; }
expr::expr(int i) { _expr = i; }
expr::expr(GiNaC::ex e) { _expr = e; }
expr::expr(const std::string &s) { _expr = HepLib::str2ex(s); }
expr::expr(const std::vector<expr> &ev) {
    GiNaC::lst ret;
    for(auto item : ev) ret.append(item._expr);
    _expr = ret;
}
expr::expr(const std::string &s, std::map<std::string,expr> st) {
    GiNaC::symtab _st;
    for(auto kv : st) _st[kv.first] = kv.second._expr;
    _expr = HepLib::str2ex(s, _st);
}
    
expr expr::operator+(const expr &e) { return expr(_expr + e._expr); }
expr expr::operator-(const expr &e) { return expr(_expr - e._expr); }
expr expr::operator*(const expr &e) { return expr(_expr * e._expr); }
expr expr::operator/(const expr &e) { return expr(_expr / e._expr); }
expr expr::operator>>(const expr &e) { return expr(_expr == e._expr); }
expr expr::operator>>(const int &e) { return expr(_expr == e); }
bool expr::operator<(const expr &e) const { return _expr < e._expr; }
bool expr::operator>(const expr &e) const { return _expr > e._expr; }
bool expr::operator==(const expr & e) const { return _expr.is_equal(e._expr); }
unsigned expr::gethash() { return _expr.gethash(); }

expr expr::operator+(const int i) { return expr(_expr + i); }
expr expr::operator-(const int i) { return expr(_expr - i); }
expr expr::operator*(const int i) { return expr(_expr * i); }
expr expr::operator/(const int i) { return expr(_expr / GiNaC::ex(i)); }
expr expr::operator-() { return expr(-_expr); }
expr expr::operator[](const int i) { return _expr.op(i); }

std::string expr::str() { return HepLib::ex2str(_expr); }
std::string expr::__str__() { return HepLib::ex2str(_expr); }

unsigned int expr::nops() { return _expr.nops(); }
expr expr::op(unsigned int i) { return expr(_expr.op(i)); }
void expr::let_op(unsigned int i, expr e) { _expr.let_op(i) = e._expr; }
expr expr::__getitem__(const int i) { return op(i); }

expr expr::subs(const std::vector<expr> &ev) {
    GiNaC::lst repl;
    for(auto item : ev) repl.append(item._expr);
    return expr(_expr.subs(repl));
}

expr expr::subs(const expr &e) {
    return expr(_expr.subs(e._expr));
}

expr expr::subs(const exmap & e) {
    return expr(_expr.subs(e._g));
}

bool expr::match(const expr &e) {
    return _expr.match(e._expr);
}

bool expr::has(const expr &e) {
    return _expr.has(e._expr);
}

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

bool expr::isSymbol() { return GiNaC::is_a<HepLib::Symbol>(_expr); }
bool expr::isVector() { return GiNaC::is_a<HepLib::Vector>(_expr); }
bool expr::isIndex() { return GiNaC::is_a<HepLib::Index>(_expr); }
bool expr::isPair() { return GiNaC::is_a<HepLib::Pair>(_expr); }
bool expr::isDGamma() { return GiNaC::is_a<HepLib::DGamma>(_expr); }

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
    return expr(HepLib::series_ex(_expr, GiNaC::ex_to<HepLib::Symbol>(s._expr), o));
}

expr expr::collect(const expr & s) {
    return expr(HepLib::collect_o(_expr, s._expr));
}

it4expr expr::__iter__() {
    return it4expr(_expr.begin(), _expr.end());
}

expr expr::evalf() const {
    return expr(_expr.evalf());
}

expr conjugate(const expr &e) {
    return expr(conjugate(e._expr));
}

expr expand(const expr &e) {
    return expr(e._expr.expand());
}

expr evalf(const expr &e) {
    return e.evalf();
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

expr subs(const expr & e1, const exmap & e2) {
    return expr(e1._expr.subs(e2._g));
}

expr series(const expr &e, const expr &s, int o) {
    if(!GiNaC::is_a<HepLib::Symbol>(s._expr)) throw HepLib::Error("1st argument should be a Symbol.");
    return expr(HepLib::series_ex(e._expr, GiNaC::ex_to<HepLib::Symbol>(s._expr), o));
}

expr collect(const expr & e, const expr & s) {
    return expr(HepLib::collect_o(e._expr, s._expr));
}

expr pow(const expr &e1, const expr &e2) {
    return expr(GiNaC::ex(GiNaC::pow(e1._expr, e2._expr)));
}

expr abs(const expr &z) {
    return expr(GiNaC::ex(GiNaC::abs(z._expr)));
}

expr real(const expr &z) {
    return expr(GiNaC::ex(GiNaC::real(GiNaC::ex_to<GiNaC::numeric>(z._expr))));
}

expr imag(const expr &z) {
    return expr(GiNaC::ex(GiNaC::imag(GiNaC::ex_to<GiNaC::numeric>(z._expr))));
}

expr csgn(const expr &z) {
    return expr(GiNaC::ex(GiNaC::csgn(z._expr)));
}

expr step(const expr &z) {
    return expr(GiNaC::ex(GiNaC::step(z._expr)));
}

expr numer(const expr &z) {
    return expr(GiNaC::ex(GiNaC::numer(z._expr)));
}

expr denom(const expr &z) {
    return expr(GiNaC::ex(GiNaC::denom(z._expr)));
}

expr sqrt(const expr &z) {
    return expr(GiNaC::ex(GiNaC::sqrt(z._expr)));
}

expr sin(const expr &z) {
    return expr(GiNaC::ex(GiNaC::sin(z._expr)));
}

expr cos(const expr &z) {
    return expr(GiNaC::ex(GiNaC::cos(z._expr)));
}

expr tan(const expr &z) {
    return expr(GiNaC::ex(GiNaC::tan(z._expr)));
}

expr asin(const expr &z) {
    return expr(GiNaC::ex(GiNaC::asin(z._expr)));
}

expr acos(const expr &z) {
    return expr(GiNaC::ex(GiNaC::acos(z._expr)));
}

expr atan(const expr &y, const expr &x) {
    return expr(GiNaC::ex(GiNaC::atan(GiNaC::ex_to<GiNaC::numeric>(y._expr), GiNaC::ex_to<GiNaC::numeric>(x._expr))));
}

expr sinh(const expr &z) {
    return expr(GiNaC::ex(GiNaC::sinh(z._expr)));
}

expr cosh(const expr &z) {
    return expr(GiNaC::ex(GiNaC::cosh(z._expr)));
}

expr tanh(const expr &z) {
    return expr(GiNaC::ex(GiNaC::tanh(z._expr)));
}

expr asinh(const expr &z) {
    return expr(GiNaC::ex(GiNaC::asinh(z._expr)));
}

expr acosh(const expr &z) {
    return expr(GiNaC::ex(GiNaC::acosh(z._expr)));
}

expr atanh(const expr &z) {
    return expr(GiNaC::ex(GiNaC::atanh(z._expr)));
}

expr exp(const expr &z) {
    return expr(GiNaC::ex(GiNaC::exp(z._expr)));
}

expr log(const expr &z) {
    return expr(GiNaC::ex(GiNaC::log(z._expr)));
}

expr Li2(const expr &z) {
    return expr(GiNaC::ex(GiNaC::Li2(z._expr)));
}

expr zeta(const expr &z) {
    return expr(GiNaC::ex(GiNaC::zeta(z._expr)));
}

expr tgamma(const expr &z) {
    return expr(GiNaC::ex(GiNaC::tgamma(z._expr)));
}

expr lgamma(const expr &z) {
    return expr(GiNaC::ex(GiNaC::lgamma(z._expr)));
}

expr psi(const expr &z) {
    return expr(GiNaC::ex(GiNaC::psi(z._expr)));
}

expr psi(const expr &n, const expr &z) {
    return expr(GiNaC::ex(GiNaC::psi(n._expr, z._expr)));
}

expr factorial(const expr &n) {
    return expr(GiNaC::ex(GiNaC::factorial(n._expr)));
}

expr doublefactorial(const expr &n) {
    return expr(GiNaC::ex(GiNaC::doublefactorial(GiNaC::ex_to<GiNaC::numeric>(n._expr))));
}

expr binomial(const expr &n, const expr &k) {
    return expr(GiNaC::ex(GiNaC::binomial(n._expr, k._expr)));
}

expr bernoulli(const expr &n) {
    return expr(GiNaC::ex(GiNaC::bernoulli(GiNaC::ex_to<GiNaC::numeric>(n._expr))));
}

expr fibonacci(const expr &n) {
    return expr(GiNaC::ex(GiNaC::fibonacci(GiNaC::ex_to<GiNaC::numeric>(n._expr))));
}

expr mod(const expr &a, const expr &b) {
    return expr(GiNaC::ex(GiNaC::mod(GiNaC::ex_to<GiNaC::numeric>(a._expr), GiNaC::ex_to<GiNaC::numeric>(b._expr))));
}

expr smod(const expr &a, const expr &b) {
    return expr(GiNaC::ex(GiNaC::smod(GiNaC::ex_to<GiNaC::numeric>(a._expr), GiNaC::ex_to<GiNaC::numeric>(b._expr))));
}

expr irem(const expr &a, const expr &b) {
    return expr(GiNaC::ex(GiNaC::irem(GiNaC::ex_to<GiNaC::numeric>(a._expr), GiNaC::ex_to<GiNaC::numeric>(b._expr))));
}

expr irem(const expr &a, const expr &b, const expr &q) {
    auto qq = GiNaC::ex_to<GiNaC::numeric>(q._expr);
    return expr(GiNaC::ex(GiNaC::irem(GiNaC::ex_to<GiNaC::numeric>(a._expr), GiNaC::ex_to<GiNaC::numeric>(b._expr), qq)));
}

expr iquo(const expr &a, const expr &b) {
    return expr(GiNaC::ex(GiNaC::iquo(GiNaC::ex_to<GiNaC::numeric>(a._expr), GiNaC::ex_to<GiNaC::numeric>(b._expr))));
}

expr iquo(const expr &a, const expr &b, const expr &r) {
    auto rr = GiNaC::ex_to<GiNaC::numeric>(r._expr);
    return expr(GiNaC::ex(GiNaC::iquo(GiNaC::ex_to<GiNaC::numeric>(a._expr), GiNaC::ex_to<GiNaC::numeric>(b._expr), rr)));
}

expr gcd(const expr &a, const expr &b) {
    return expr(GiNaC::ex(GiNaC::gcd(a._expr, b._expr)));
}

expr lcm(const expr &a, const expr &b) {
    return expr(GiNaC::ex(GiNaC::lcm(a._expr, b._expr)));
}

expr pow(const expr &e, const int n) {
    return expr(GiNaC::ex(GiNaC::pow(e._expr, n)));
}

expr wild(const int wi) {
    return expr(GiNaC::wild(wi));
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

expr lst(const std::vector<expr> &ev) {
    return expr(ev);
}

exvec::exvec() { }

exvec::exvec(std::vector<expr> es) {
    for(auto it : es) _g.push_back(it._expr);
}

exvec::exvec(expr e) {
    for(auto it : e._expr) _g.push_back(it);
}

exvec::exvec(GiNaC::exvector es) {
    for(auto it : es) _g.push_back(it);
}

void exvec::push_back(expr e) {
    _g.push_back(e._expr);
}

expr exvec::__getitem__(const int i) {
    return expr(_g[i]);
}

expr exvec::op(const int i) {
    return expr(_g[i]);
}

void exvec::__setitem__(const int i, expr v) {
    _g[i] = v._expr;
}

int exvec::size() {
    return _g.size();
}

int exvec::nops() {
    return _g.size();
}

std::string exvec::str() {
    return HepLib::ex2str(_g);
}

std::string exvec::__str__() {
    return str();
};

void exvec::subs(const expr & e) {
    for(auto & it : _g) it = it.subs(e._expr);
}
void exvec::subs(const std::vector<expr> & e) {
    GiNaC::lst ss;
    for(auto it : e) ss.append(it._expr);
    for(auto & it : _g) it = it.subs(ss);
}
void exvec::subs(const exmap & e) {
    for(auto & it : _g) it = it.subs(e._g);
}

it4vec exvec::__iter__() {
    return it4vec(_g.begin(), _g.end());
}

void exvec::sort() {
    std::sort(_g.begin(), _g.end(), HepLib::ex_less);
}

exmap::exmap() { }
exmap::exmap(std::map<expr,expr,expr_is_less> es) {
    for(auto it : es) _g[it.first._expr] = it.second._expr;
}

expr exmap::__getitem__(expr e) {
    return _g[e._expr];
}

void exmap::__setitem__(expr k, expr v) {
    _g[k._expr] = v._expr;
}

int exmap::size() {
    return _g.size();
}

int exmap::nops() {
    return _g.size();
}

std::string exmap::str() {
    return HepLib::ex2str(_g);
}

std::string exmap::__str__() {
    return str();
}

exmap::exmap(GiNaC::exmap es) {
    for(auto kv : es) _g[kv.first] = kv.second;
}

it4map exmap::__iter__() {
    return it4map(_g.begin(), _g.end());
}

exset::exset() { }

exset::exset(std::vector<expr> es) {
    for(auto it : es) _g.insert(it._expr);
}

exset::exset(expr e) {
    for(auto it : e._expr) _g.insert(it);
}

void exset::insert(expr e) {
    _g.insert(e._expr);
}

int exset::size() {
    return _g.size();
}

int exset::nops() {
    return _g.size();
}

std::string exset::str() {
    return HepLib::ex2str(_g);
}

std::string exset::__str__() {
    return str();
}

exset::exset(GiNaC::exset es) {
    for(auto it : es) _g.insert(it);
}

it4set exset::__iter__() {
    return it4set(_g.begin(), _g.end());
}

Function::Function() { }
Function::~Function() { }
expr Function::operator()(const expr &e) { return expr(0); }
expr Function::operator()(const expr &e1,const expr &e2) { return expr(0); }
expr Function::operator()(const expr &e1,const expr &e2,const expr &e3) { return expr(0); }
expr Function::operator()(const expr &e1,const expr &e2,const expr &e3,const expr &e4) { return expr(0); }
expr Function::operator()(const expr &e1,const expr &e2,const expr &e3,const expr &e4,const expr &e5) { return expr(0); }

ParFun::ParFun() { }
ParFun::~ParFun() { }
expr ParFun::__call__(const int i) { throw HepLib::Error("Not implemented yet!"); return expr(1979); }
exvec Parallel(int ntotal, int nbatch,
        ParFun &f,
        const std::string & key,
        bool rm,
        const std::string & pre) {
    auto ret = HepLib::GiNaC_Parallel(ntotal, nbatch,
        [&f](int i)->GiNaC::ex{ return f.__call__(i)._expr; },
        key, rm, pre);
    exvec res;
    for(auto item : ret) res._g.push_back(item);
    return res;
}
exvec Parallel(int ntotal,
        ParFun &f,
        const std::string & key,
        bool rm,
        const std::string & pre) {
    auto ret = HepLib::GiNaC_Parallel(ntotal,
        [&f](int i)->GiNaC::ex{ return f.__call__(i)._expr; },
        key, rm, pre);
    exvec res;
    for(auto item : ret) res._g.push_back(item);
    return res;
}

bool isFunction(const expr &e, std::string sf) {
    return HepLib::isFunction(e._expr, sf);
}

expr file2expr(std::string fn) {
    return expr(HepLib::file2ex(fn));
}

std::map<std::string,expr> garReadAll(const std::string &garfn) {
    std::map<std::string, GiNaC::ex> _resMap;
    HepLib::garRead(garfn, _resMap);
    std::map<std::string, expr> resMap;
    for(auto kv : _resMap) resMap[kv.first] = expr(kv.second);
    return resMap;
}

expr garRead(const std::string &garfn, const char* key) {
     return expr(HepLib::garRead(garfn, key));
}

expr garRead(const std::string &garfn) {
    return expr(HepLib::garRead(garfn));
}

void garWrite(const std::string &garfn, const std::map<std::string, expr> &resMap) {
    std::map<std::string, GiNaC::ex> _resMap;
    for(auto kv : resMap) _resMap[kv.first] = kv.second._expr;
    HepLib::garWrite(garfn, _resMap);
}

void garWrite(const std::map<std::string, expr> &resMap, const std::string &garfn) {
    garWrite(garfn,resMap);
}

void garWrite(const std::string &garfn, const expr & res) {
    HepLib::garWrite(garfn, res._expr);
}

void garWrite(const expr & res, const std::string &garfn) {
    garWrite(garfn,res);
}

std::string RunOS(const std::string & cmd) {
    return HepLib::RunOS(cmd);
}

MapFunction::MapFunction() : _map([this](const GiNaC::ex &e, HepLib::MapFunction &self)->GiNaC::ex{ return map(expr(e))._expr; }) {}
MapFunction::~MapFunction() { }
expr MapFunction::map(const expr & e) { return e; }
expr MapFunction::operator() (const expr &e) { return _map(e._expr); }
exvec MapFunction::operator() (const exvec &ev) {
    exvec ret;
    for(auto item : ev._g) ret._g.push_back(_map(item));
    return ret;
}

cout & cout::operator<<(const expr &e) { std::cout << e._expr; return *this; }
cout & cout::operator<<(const int &e) { std::cout << e; return *this; }
cout & cout::operator<<(const std::string &e) { std::cout << e; return *this; }
cout & cout::operator<<(const char* &e) { std::cout << e; return *this; }
cout & cout::operator<<(const exvec &ev) { std::cout << ev._g; return *this; }
cout & cout::operator<<(const exmap &em) { std::cout << em._g; return *this; }
cout & cout::operator<<(const exset &es) { std::cout << es._g; return *this; }
cout & cout::operator<<(const std::vector<expr> &ev) { std::cout << exvec(ev)._g; return *this; }
cout & cout::operator<<(const std::map<expr,expr,expr_is_less> &em) { std::cout << exmap(em)._g; return *this; }

hout & hout::operator<<(const expr &e) { HepLib::hout << e._expr; return *this; }
hout & hout::operator<<(const int &e) { HepLib::hout << e; return *this; }
hout & hout::operator<<(const std::string &e) { HepLib::hout << e; return *this; }
hout & hout::operator<<(const char* &e) { HepLib::hout << e; return *this; }
hout & hout::operator<<(const exvec &ev) { HepLib::hout << ev._g; return *this; }
hout & hout::operator<<(const exmap &em) { HepLib::hout << em._g; return *this; }
hout & hout::operator<<(const exset &es) { HepLib::hout << es._g; return *this; }
hout & hout::operator<<(const std::vector<expr> &ev) { HepLib::hout << exvec(ev)._g; return *this; }
hout & hout::operator<<(const std::map<expr,expr,expr_is_less> &em) { HepLib::hout << exmap(em)._g; return *this; }
