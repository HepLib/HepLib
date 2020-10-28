#include "Wrap.h"

expr::expr(int i) { _expr = i; }
expr::expr(GiNaC::ex e) { _expr = e; }
expr::expr(const std::string &s) { _expr = HepLib::str2ex(s); }
    
expr expr::operator+(const expr &e) { return _expr + e._expr; }
expr expr::operator-(const expr &e) { return _expr - e._expr; }
expr expr::operator*(const expr &e) { return _expr * e._expr; }
expr expr::operator/(const expr &e) { return _expr / e._expr; }

expr expr::operator+(const int i) { return _expr + i; }
expr expr::operator-(const int i) { return _expr - i; }
expr expr::operator*(const int i) { return _expr * i; }
expr expr::operator/(const int i) { return _expr / GiNaC::ex(i); }

expr expr::operator-() { return expr(-_expr); }
    
std::string expr::__str__() { return HepLib::ex2str(_expr); }

expr expr::expand() {
    return expr(GiNaC::ex(_expr.expand()));
}

expr expr::normal() {
    return expr(GiNaC::ex(HepLib::fermat_normal(_expr)));
}

expr expr::factor() {
    return expr(GiNaC::ex(HepLib::form_factor(_expr)));
}

expr expr::series(const expr &s, int o) {
    if(!GiNaC::is_a<HepLib::Symbol>(s._expr)) throw HepLib::Error("1st argument should be a Symbol.");
    return expr(HepLib::mma_series(_expr, GiNaC::ex_to<HepLib::Symbol>(s._expr), o));
}

expr expand(const expr &e) {
    return expr(e._expr.expand());
}

expr normal(const expr &e) {
    return expr(HepLib::fermat_normal(e._expr));
}

expr factor(const expr &e) {
    return expr(HepLib::form_factor(e._expr));
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
    return expr(GiNaC::ex(HepLib::FC::Index(s)));
}

expr Vector(const std::string &s) {
    return expr(GiNaC::ex(HepLib::FC::Vector(s)));
}

expr SP(const expr &e1, const expr &e2) {
    return expr(GiNaC::ex(HepLib::FC::SP(e1._expr, e2._expr)));
}

expr GAS(const expr &e) {
    return expr(GiNaC::ex(HepLib::FC::GAS(e._expr)));
}
        
expr TR(const expr &e) {
    return expr(GiNaC::ex(HepLib::FC::TR(e._expr)));
}

expr form(const expr &e) {
    return expr(HepLib::FC::form(e._expr));
}

void letSP(const expr &e1, const expr &e2, const expr &e12) {
    HepLib::FC::letSP(e1._expr, e2._expr) = e12._expr;
}
