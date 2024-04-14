/**
 * @file
 * @brief Functions for C++ code generation for Numerical Integration
 */
 
#include "SD.h"

namespace HepLib::SD {

int cseParser::on() { return no; }
vector<pair<int, ex>> cseParser::os() { return o_ex_vec;}

ex cseParser::Parse(ex expr, bool reset) {
    if(reset) {
        no = 0;
        o_ex_vec.clear();
    }
    
    if(is_a<numeric>(expr) || is_a<symbol>(expr) || is_a<symbol>(-expr) || isFunction(expr,"x") || isFunction(-expr,"x") || isFunction(expr,"PL") || isFunction(-expr,"PL")) return expr;
    if(!is_zero(ex_var_map[expr])) {
        return ex_var_map[expr];
    }
    
    ex ret = expr;
    if(is_a<add>(expr)) {
        ret = 0;
        for(auto item : expr) ret += Parse(item,false);
    } else if(is_a<lst>(expr)) {
        lst oret;
        for(auto item : expr) oret.append(Parse(item,false));
        return oret;
    } else if(is_a<mul>(expr)) {
        ret = 1;
        for(auto item : expr) ret *= Parse(item,false);
    } else if(is_a<power>(expr)) {
        ret = power(Parse(expr.op(0),false), expr.op(1));
    } else if(expr.match(log(w))) {
        ret = log(Parse(expr.op(0),false));
    } else if(expr.match(exp(w))) {
        ret = exp(Parse(expr.op(0),false));
    } else if(expr.match(sqrt(w))) {
        ret = sqrt(Parse(expr.op(0),false));
    } // else ret is expr
    stringstream ss;
    ss << oc << "[" << no << "]";
    Symbol so(ss.str().c_str());
    o_ex_vec.push_back(make_pair(no, ret));
    no++;
    ex_var_map[expr] = so;
    return so;
}


int cse_Parser::vn() { return no; }
const vector<pair<int,ex>> & cse_Parser::vs() { return on_ex_vec;}

ex cse_Parser::Parse(ex expr, bool reset) {
    if(reset) {
        no = 0;
        exv.clear();
        on_ex_vec.clear();
    }
    
    auto itr = exn.find(expr);
    if(itr!=exn.end()) {
        used[itr->second]++;
        return exv[itr->second];
    }
    
    ex ret = expr;
    if(is_a<numeric>(expr) || is_a<symbol>(expr) || is_a<symbol>(-expr) || isFunction(expr,"x") || isFunction(-expr,"x") || isFunction(expr,"PL") || isFunction(-expr,"PL")) return expr;
    if(is_a<add>(expr)) {
        ret = 0;
        for(auto item : expr) ret += Parse(item,false);
    } else if(is_a<lst>(expr)) {
        lst oret;
        for(auto item : expr) oret.append(Parse(item,false));
        return oret;
    } else if(is_a<mul>(expr)) {
        ret = 1;
        for(auto item : expr) ret *= Parse(item,false);
    } else if(is_a<power>(expr)) {
        ret = power(Parse(expr.op(0),false), expr.op(1));
    } else if(expr.match(log(w))) {
        ret = log(Parse(expr.op(0),false));
    } else if(expr.match(exp(w))) {
        ret = exp(Parse(expr.op(0),false));
    } else if(expr.match(sqrt(w))) {
        ret = sqrt(Parse(expr.op(0),false));
    } // else ret is expr
    exv.push_back(ret);
    exn[expr] = no;
    used[no] = 1;
    ex res = iWF(no);
    no++;
    if(!reset) return res;
    
    map<int,int> s2l, l2s;
    MapFunction map([&l2s,this](const ex & e, MapFunction &self)->ex{
        if(!e.has(iWF(w))) return e;
        else if(e.match(iWF(w))) {
            int wi = ex2int(e.op(0));
            used[wi]--;
            return Symbol(v+"["+to_string(l2s[wi])+"]");
        } else return e.map(self);
    });
    
    int max = 0;
    int next = -1;
    for(int i=0; i<no; i++) {
cout << i << "/" << no << "/" << max << endl;
        int on = i;
        if(next!=-1) { on = next; }
        else { on = max; max++; }
        l2s[i] = on;
        s2l[on] = i;
        on_ex_vec.push_back(make_pair(on, map(exv[i])));
        next = -1;
        for(int j=0; j<max; j++) if(used[s2l[j]]<1) {
            next = j;
            break;
        }
    }
    ret = map(res);
    no = max;
    return ret;
}


}
