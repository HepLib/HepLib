#include "SD.h"

namespace HepLib {

/*-----------------------------------------------------*/
// Global Symbol
/*-----------------------------------------------------*/
const symbol & get_symbol(const string & s) {
    static map<string, symbol> directory;
    map<string, symbol>::iterator i = directory.find(s);
    if (i != directory.end())
        return i->second;
    else
        return directory.insert(make_pair(s, symbol(s))).first->second;
}

/*-----------------------------------------------------*/
// split
/*-----------------------------------------------------*/
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

/*-----------------------------------------------------*/
// MatHelper
/*-----------------------------------------------------*/
bool MatHelper::has_zero_row(const matrix &mat) {
    for(int r=0; r<mat.rows(); r++) {
        bool iszero = true;
        for(int c=0; c<mat.cols(); c++) {
            if(!mat(r,c).is_zero()) {
                iszero = false;
                break;
            }
        }
        if(iszero) return true;
    }
    return false;
}

bool MatHelper::is_zero_row(const matrix &mat, int r) {
    assert(r<mat.rows());
    for(int c=0; c<mat.cols(); c++) {
        if(mat(r, c) !=0 ) return false;
    }
    return true;
}

vector<int> MatHelper::zero_row_index(const matrix &mat) {
    vector<int> ret;
    for(int r=0; r<mat.rows(); r++) {
        bool iszero = true;
        for(int c=0; c<mat.cols(); c++) {
            if(!mat(r,c).is_zero()) {
                iszero = false;
                break;
            }
        }
        if(iszero) ret.push_back(r);
    }
    return ret;
}

matrix MatHelper::remove_zero_rows(const matrix &mat) {
    vector<int> zri = zero_row_index(mat);
    if(zri.size()<1) return mat;
    matrix ret(mat.rows()-zri.size(), mat.cols());

    int cr = 0;
    for(int r=0; r<mat.rows(); r++) {
        bool iszero = find(zri.begin(), zri.end(), r) != zri.end();
        if(iszero) continue;
        for(int c=0; c<mat.cols(); c++) {
            ret(cr, c) = mat(r, c);
        }
        cr++;
    }
    return ret;
}

matrix MatHelper::sub(matrix mat, int r, int nr, int c, int nc) {
    matrix ret(nr, nc);
    for(int ir=0; ir<nr; ir++) {
        for(int ic=0; ic<nc; ic++) {
            ret(ir, ic) = mat(r+ir, c+ic);
        }
    }
    return ret;
}



/*-----------------------------------------------------*/
// lstHelper
/*-----------------------------------------------------*/
lst lstHelper::sub(lst m, int s, int n) {
    lst ret;
    for(int i=0; (i<=n || n<0) && i+s<m.nops(); i++) ret.append(m.op(s+i));
    return ret;
}

ex lstHelper::sum(lst m) {
    ex ret = 0;
    for(int i=0; i<m.nops(); i++) ret += m.op(i);
    return ret;
}

lst lstHelper::subs(lst m, ex r) {
    return map(m, [&](auto &&e){return e.subs(r);});
}

lst lstHelper::subs(lst m, exmap r) {
    return map(m, [&](auto &&e){return e.subs(r);});
}

/*-----------------------------------------------------*/
// now string
/*-----------------------------------------------------*/
string now(bool use_date) {
    time_t timep;
    time (&timep);
    char tmp[64];
    if(use_date) strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&timep) );
    else strftime(tmp, sizeof(tmp), "%H:%M:%S",localtime(&timep) );
    return tmp;
}

/*-----------------------------------------------------*/
// Symbols
/*-----------------------------------------------------*/
void gather_symbols_inner(const ex & e, lst & l) {
    if (is_a<symbol>(e)) {
        l.append(e);
    } else {
        for(auto const &i : e) gather_symbols_inner(i, l);
    }
}

lst gather_symbols(const ex & e) {
    lst l;
    gather_symbols_inner(e, l);
    l.sort();
    l.unique();
    return l;
}

lst gather_symbols(const vector<ex> & ve) {
    lst l;
    for(auto const & e : ve) gather_symbols_inner(e, l);
    l.sort();
    l.unique();
    return l;
}

lst gather_symbols(const vector<pair<lst,lst>> & ve) {
    lst l;
    for(auto const & kv : ve) {
        gather_symbols_inner(kv.first, l);
        gather_symbols_inner(kv.second, l);
    }
    l.sort();
    l.unique();
    return l;
}

lst gather_symbols(const vector<pair<ex,ex>> & ve) {
    lst l;
    for(auto const & kv : ve) {
        gather_symbols_inner(kv.first, l);
        gather_symbols_inner(kv.second, l);
    }
    l.sort();
    l.unique();
    return l;
}

/*-----------------------------------------------------*/
// RunOS Command
/*-----------------------------------------------------*/
string RunOS(const char * cmd) {
    char buf[128];
    ostringstream oss;
    FILE *fp = NULL;
    if( (fp = popen(cmd, "r")) != NULL) {
        while(fgets(buf, 128, fp) != NULL) {
            oss << buf;
        }
        pclose(fp);
        fp = NULL;
    }
    return oss.str();
}


/*-----------------------------------------------------*/
// garResult Function
/*-----------------------------------------------------*/
ex garResult(const char *garfn, lst syms) {
    archive ar;
    ifstream in(garfn);
    in >> ar;
    in.close();
    auto c = ar.unarchive_ex(syms, "c");
    auto res = ar.unarchive_ex(syms, "res");
    if(c!=19790923) {
        cerr << "gar file: " << garfn << endl;
        cerr << "c=" << c << ", different from 19790923!" << endl;
        assert(false);
    }
    return res;
}


/*-----------------------------------------------------*/
// Series at s=0 similar to Mathematica
/*-----------------------------------------------------*/
ex mma_series(ex expr_in, symbol s, int sn) {
    ex expr = expr_in;
    if(!expr.has(s)) return expr;
    int exN = 1;
    ex expr_input = mma_collect(expr_in,s,true);
    
    // make sure CCF has no s
    exset cset;
    expr_input.find(CCF(wild()), cset);
    for(auto ccf : cset) {
        if(ccf.has(s)) {
            cerr << "ccf = " << ccf << endl;
            assert(false);
            break;
        }
    }
    
    while(exN<10) {
        expr = expr_input + pow(s,sn+exN+2);
        ex res = expr.series(s, sn+exN);
        res = res.subs(CCF(wild())==wild()); // remove CCF
        ex ot = 0;
        for(int i=0; i<res.nops(); i++) {
            if(is_order_function(res.op(i))) {
                ot = res.op(i);
                break;
            }
        }
        if(!is_order_function(ot)) {
            cerr << RED << "Not an Order term: " << ot << RESET << endl;
            cerr << "expr = " << expr << endl;
            cerr << "res = " << res << endl;
            assert(false);
        }
        if(ot.op(0).degree(s)>sn) {
            res = series_to_poly(res);
            res = mma_collect(res,s);
            ex ret = 0;
            for(int i=res.ldegree(s); (i<=res.degree(s) && i<=sn); i++) {
                ret += res.coeff(s,i) * pow(s, i);
            }
            return ret;
        }
        exN++;
    }
    cerr << RED << "mma_series seems not working!" << RESET << endl;
    assert(false);
    return 0;
}

/*-----------------------------------------------------*/
// mma_collect
/*-----------------------------------------------------*/
ex mma_diff(ex expr, ex xp, unsigned nth, bool expand) {
    symbol s;
    ex res = expr.subs(xp==s);
    
    if(expand) {
        res = mma_collect(res, s, true);
        
        // make sure CCF has no s
        exset cset;
        res.find(CCF(wild()), cset);
        for(auto ccf : cset) {
            if(ccf.has(s)) {
                cerr << "ccf = " << ccf << endl;
                assert(false);
                break;
            }
        }
    }
    
    res = res.diff(s, nth);
    res = res.subs(CCF(wild())==wild()); // remove CCF
    res = res.subs(s==xp);
    
    return res;
}


/*-----------------------------------------------------*/
// mma_collect
/*-----------------------------------------------------*/
ex mma_collect(ex expr_in, ex pat, bool ccf, bool cvf) {
    lst cv_repl = lst{ CCF(wild(0))==wild(0), CVF(wild(1))==wild(1) };
    lst c_repl = lst{ CCF(wild())==wild() };
    lst v_repl = lst{ CVF(wild())==wild() };
    
    auto expr = expr_in.subs(cv_repl); // remove CCF & CVF
    
    if(!expr.has(pat)) return ccf ? CCF(expr) : expr;
    ex res;
    if(is_a<add>(expr)) {
        res = 0;
        for(auto item : expr) {
            res += mma_collect(item, pat, true, true);
        }
    } else if(is_a<mul>(expr)) {
        res = 1;
        for(auto item : expr) {
            res *= mma_collect(item, pat, true, true);
        }
    } else if(is_a<power>(expr) && is_a<numeric>(expr.op(1)) && ex_to<numeric>(expr.op(1)).is_nonneg_integer()) {
        res = pow(mma_collect(expr.op(0), pat, true, true), expr.op(1));
    } else {
        return cvf ? CVF(expr) : expr;
    }
    
    res = res.expand(); // expand internally
    exset vfset;
    res.find(CVF(wild()), vfset);
    lst vflst;
    for(auto vf : vfset) vflst.append(vf);
    res = collect(res, vflst, true); // collect internally
    
    res = res.subs(cv_repl); // remove CCF & CVF
    lst items;
    if(is_a<add>(res)) {
        for(auto item : res) items.append(item);
    } else {
        items.append(res);
    }
    
    ex ret = 0, cf = 0;
    for(auto item : items) {
        if(!item.has(pat)) cf += item;
        else {
            if(is_a<mul>(item)) {
                ex tc = 1, tf = 1;
                for(auto ii : item) {
                    if(!ii.has(pat)) tc *= ii;
                    else tf *= ii;
                }
                ret += CVF(tf) * CCF(tc);
            } else {
                ret += CVF(item);
            }
        }
    }
    ret += CCF(cf);
    
    ret.find(CVF(wild()), vfset);
    vflst.remove_all();
    for(auto vf : vfset) vflst.append(vf);
    ret = collect(ret, vflst, true); // collect internally
    
    if(!ccf) ret = ret.subs(c_repl);
    if(!cvf) ret = ret.subs(v_repl);
    
    return ret;
}

/*-----------------------------------------------------*/
// Customized GiNaC Function
/*-----------------------------------------------------*/
static ex CCF_Diff(const ex & x, unsigned diff_param) {return 0;}
REGISTER_FUNCTION(VF, dummy())
REGISTER_FUNCTION(VF1, dummy())
REGISTER_FUNCTION(VF2, dummy())
REGISTER_FUNCTION(VF3, dummy())

REGISTER_FUNCTION(CCF, derivative_func(CCF_Diff))
REGISTER_FUNCTION(CVF, dummy())

}
