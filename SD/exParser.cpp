#include "SD.h"

namespace HepLib {

ex exParser::Parse(ex expr, int level) {
    if(is_a<numeric>(expr) || is_a<symbol>(expr)) return expr;
    excount[expr]++;
    if(!ex_var_map[expr].is_zero()) return ex_var_map[expr];
    ex ret = expr;
    if(is_a<add>(expr)) {
        ret = 0;
        for(auto item : expr) ret += Parse(item, level+1);
    } else if(is_a<mul>(expr)) {
        ret = 1;
        for(auto item : expr) ret *= Parse(item, level+1);
    } else if(is_a<power>(expr)) {
        ret = power(Parse(expr.op(0), level+1), expr.op(1));
    } 
    ostringstream ss;
    ss << prefix << "[" << no << "]";
    symbol so(ss.str().c_str());
    o_ex_vec.push_back(make_pair(no, ret));
    no++;
    ex_var_map[expr] = so;
    return so;
}


}
