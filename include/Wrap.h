/**
 * @file
 * @brief Wrap Class for SWIG
 */

#include "HepLib.h"

class expr {
public:
    expr(int i);
    expr(GiNaC::ex e);
    expr(const std::string &s);
    
    expr operator+(const expr &e);
    expr operator-(const expr &e);
    expr operator*(const expr &e);
    expr operator/(const expr &e);
    
    expr operator+(const int i);
    expr operator-(const int i);
    expr operator*(const int i);
    expr operator/(const int i);
    
    expr operator-();
    
    expr expand();
    expr normal();
    expr factor();
    expr series(const expr &s, int o);
    
    std::string __str__();
    
    //friend functions
    friend expr expand(const expr &e);
    friend expr normal(const expr &e);
    friend expr factor(const expr &e);
    
    friend expr series(const expr &e, const expr &s, int o);
    friend expr pow(const expr &e1, const expr &e2);
    friend expr pow(const expr &e, const int n);
    
    friend expr SP(const expr &e1, const expr &e2);
    friend expr GAS(const expr &e);
    friend expr TR(const expr &e);
    friend expr form(const expr &e);
    
    friend void letSP(const expr &e1, const expr &e2, const expr &e12);
    
private:
    GiNaC::ex _expr;
};

expr expand(const expr &e);
expr normal(const expr &e);
expr factor(const expr &e);
expr series(const expr &e, const expr &s, int o);

expr pow(const expr &e1, const expr &e2);
expr pow(const expr &e1, const int n);

expr Index(const std::string &s);
expr Vector(const std::string &s);
expr Symbol(const std::string &s);

expr SP(const expr &e1, const expr &e2);
expr GAS(const expr &e);
expr TR(const expr &e);
expr form(const expr &e);

void letSP(const expr &e1, const expr &e2, const expr &e12);



