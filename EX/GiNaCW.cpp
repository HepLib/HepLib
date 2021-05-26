/**
 * @file
 * @brief Wrap Class for SWIG
 */

#include "HepLibW.h"

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

bool expr::match(const expr &e) {
    return _expr.match(e._expr);
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
    return expr(HepLib::series_ex(e._expr, GiNaC::ex_to<HepLib::Symbol>(s._expr), o));
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

void exvec::push_back(expr e) {
    _g.push_back(e._expr);
}

int exvec::size() {
    return _g.size();
}

std::string exvec::str() {
    return HepLib::ex2str(_g);
}

std::string exvec::__str__() {
    return str();
};

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

std::string exmap::str() {
    return HepLib::ex2str(_g);
}

std::string exmap::__str__() {
    return str();
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

std::string exset::str() {
    return HepLib::ex2str(_g);
}

std::string exset::__str__() {
    return str();
};
