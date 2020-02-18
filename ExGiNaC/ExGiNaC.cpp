#include "SD.h"

namespace HepLib {

/*-----------------------------------------------------*/
// Global wildcard
/*-----------------------------------------------------*/
ex w = wild();
ex w0 = wild(0);
ex w1 = wild(1);
ex w2 = wild(2);
ex w3 = wild(3);
ex w4 = wild(4);
ex w5 = wild(5);

/*-----------------------------------------------------*/
// GiNaC_Parallel
/*-----------------------------------------------------*/
vector<ex> GiNaC_Parallel(
    int nproc, vector<ex> const &invec, std::function<ex(ex const &, int)> f,
    const char* key, int verb, bool rm, int prtlvl) {
    
    auto ppid = getpid();
    int para_max_run = nproc<0 ? omp_get_num_procs() : nproc;
    ostringstream cmd;
    cmd << "mkdir -p " << ppid;
    system(cmd.str().c_str());

    int total = invec.size();
    int batch = 1;
    if(para_max_run>0) batch = total/para_max_run/10;
    if(batch<1) batch = 1;
    int btotal = total/batch + ((total%batch)==0 ? 0 : 1);

    for(int bi=0; bi<btotal; bi++) {
        if(verb > 1) {
            cout << "\r  ";
            for(int pi=0;pi<prtlvl;pi++) cout << "   ";
            cout << "\\--Evaluating ";
            if(key != NULL) cout << WHITE << key << RESET << " ";
            cout << WHITE << batch << "x" << RESET << "[" << (bi+1) << "/" << btotal << "] ... " << flush;
        }
        
        if(para_max_run>0) {
            auto pid = fork();
            if (pid < 0) perror("fork() error");
            else if (pid != 0) {
                if(bi >= para_max_run) wait(NULL);
                continue;
            }
        }
        
        try {
            for(int ri=0; ri<batch; ri++) {
                int i = bi*batch + ri;
                if(i<total) {
                    auto item = invec[i];
                    auto res = f(item, i);
                    archive ar;
                    ar.archive_ex(res, "res");
                    ar.archive_ex(19790923, "c");
                    ostringstream garfn;
                    if(key == NULL) garfn << ppid << "/" << i << ".gar";
                    else garfn << ppid << "/" << i << "." << key << ".gar";
                    ofstream outs(garfn.str().c_str());
                    outs << ar;
                    outs.close();
                }
            }
        } catch(exception &p) {
            cout << RED << "Failed in GiNaC_Parallel!" << RESET << endl;
            cout << RED << p.what() << RESET << endl;
            if(para_max_run>0) exit(0);
            throw p;
        }
        if(para_max_run>0) exit(0);
    }
    
    auto cpid = getpid();
    if(cpid!=ppid) exit(0); // make sure
    if(para_max_run>0) while (wait(NULL) != -1) { }
    if(verb > 1 && total > 0) cout << "@" << now(false) << endl;

    auto syms = gather_symbols(invec);
    vector<ex> ovec;
    for(int i=0; i<total; i++) {
        if(verb > 1) {
            if(key == NULL) {
                cout << "\r  ";
                for(int pi=0; pi<prtlvl; pi++) cout << "   ";
                cout << "\\--Reading *.gar [" << (i+1) << "/" << total << "] ... " << flush;
            } else {
                cout << "\r  ";
                for(int pi=0;pi<prtlvl;pi++) cout << "   ";
                cout << "\\--Reading *." << WHITE << key << RESET << ".gar [" << (i+1) << "/" << total << "] ... " << flush;
            }
        }

        int oDigits = Digits;
        Digits = 50; // a fix to float overflow
        
        archive ar;
        ostringstream garfn;
        if(key == NULL) garfn << ppid << "/" << i << ".gar";
        else garfn << ppid << "/" << i << "." << key << ".gar";
        ifstream ins(garfn.str().c_str());
        ins >> ar;
        ins.close();
        remove(garfn.str().c_str());
        auto c = ar.unarchive_ex(syms, "c");
        if(c!=19790923) throw runtime_error("*.gar error!");
        auto res = ar.unarchive_ex(syms, "res");
        ovec.push_back(res);
        Digits = oDigits;
    }
    
    if(rm) {
        cmd.clear();
        cmd.str("");
        cmd << "rm -fr " << ppid;
        system(cmd.str().c_str());
        system(cmd.str().c_str());
    }
    if(verb > 1 && total > 0) cout << "@" << now(false) << endl;
    return ovec;
}

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
    if(r>=mat.rows()) {
        cerr << RED << "r>=mat.rows()" << RESET << endl;
        exit(1);
    }
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
ex lstHelper::sum(lst m) {
    ex ret = 0;
    for(int i=0; i<m.nops(); i++) ret += m.op(i);
    return ret;
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
lst gather_symbols(const ex & e) {
    lst sym_lst;
    for(const_preorder_iterator i = e.preorder_begin(); i != e.preorder_end(); ++i) {
        if(is_a<symbol>(*i)) sym_lst.append(*i);
    }
    sym_lst.sort();
    sym_lst.unique();
    return sym_lst;
}

lst gather_symbols(const vector<ex> & ve) {
    lst sym_lst;
    for(auto e : ve) {
        for(const_preorder_iterator i = e.preorder_begin(); i != e.preorder_end(); ++i) {
            if(is_a<symbol>(*i)) sym_lst.append(*i);
        }
    }
    sym_lst.sort();
    sym_lst.unique();
    return sym_lst;
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
        cerr << RED << "gar file: " << garfn << endl;
        cerr << "c=" << c << ", different from 19790923!" << RESET << endl;
        exit(1);
    }
    return res;
}

/*-----------------------------------------------------*/
// str2ex Function
/*-----------------------------------------------------*/
ex str2ex(const char *expr, symtab stab) {
    parser reader(stab);
    ex ret = reader(expr);
    return ret;
}

/*-----------------------------------------------------*/
// str2lst Function
/*-----------------------------------------------------*/
lst str2lst(const char *expr, symtab stab) {
    parser reader(stab);
    ex ret = reader(expr);
    return ex_to<lst>(ret);
}

/*-----------------------------------------------------*/
// str2lst Function
/*-----------------------------------------------------*/
lst xlst(int bi, int ei) {
    lst ret;
    for(int i=bi; i<=ei; i++) ret.append(x(i));
    return ret;
}
lst xlst(int ei) {
    return xlst(0, ei);
}

/*-----------------------------------------------------*/
// Series at s=0 similar to Mathematica
/*-----------------------------------------------------*/
ex mma_series(ex expr_in, symbol s0, int sn0) {
    ex expr = expr_in;
    if(!expr.has(s0)) return expr;
    
    exset sset;
    expr.find(pow(s0, w), sset);
    numeric sn_lcm = 1;
    for(auto pi : sset) {
        auto sn = pi.op(1);
        if(!(is_a<numeric>(sn) && ex_to<numeric>(sn).is_rational())) {
            cerr << RED << "Not rational sn = " << sn << RESET << endl;
            exit(1);
        }
        sn_lcm = lcm(sn_lcm, ex_to<numeric>(sn).denom());
    }
    symbol s;
    if(!sn_lcm.is_integer()) {
        cerr << RED << "Not integer: " << sn_lcm << RESET << endl;
        exit(1);
    }
    if(sn_lcm<0) sn_lcm = numeric(0)-sn_lcm;
    int sn = sn0 * sn_lcm.to_int();
    expr = expr.subs(pow(s0,w)==pow(s,w*sn_lcm)).subs(s0==pow(s, sn_lcm));
    
    int exN = 1;
    ex expr_input = mma_collect(expr,s,true);
        
    while(exN<10) {
        expr = expr_input + pow(s,sn+exN+2)+pow(s,sn+exN+3);
        ex res = expr.series(s, sn+exN);
        res = res.subs(CCF(w)==w); // remove CCF
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
            exit(1);
        }
        if(ot.op(0).degree(s)>sn) {
            res = series_to_poly(res);
            res = mma_collect(res,s);
            ex ret = 0;
            for(int i=res.ldegree(s); (i<=res.degree(s) && i<=sn); i++) {
                ret += res.coeff(s,i) * pow(s, i);
            }
            ret = ret.subs(pow(s,w)==pow(s0,w/sn_lcm)).subs(s==pow(s0,1/sn_lcm));
            return ret;
        }
        exN++;
    }
    cerr << RED << "mma_series seems not working!" << RESET << endl;
    exit(1);
    return 0;
}

/*-----------------------------------------------------*/
// mma_diff
/*-----------------------------------------------------*/
ex mma_diff(ex expr, ex xp, unsigned nth, bool expand) {
    symbol s;
    ex res = expr.subs(xp==s);
    if(expand) res = mma_collect(res, s, true);
    res = res.diff(s, nth);
    res = res.subs(CCF(w)==w); // remove CCF
    res = res.subs(s==xp);
    return res;
}

/*-----------------------------------------------------*/
// mma_expand
/*-----------------------------------------------------*/
struct map_CCF : public map_function {
    ex pat;
    map_CCF(const ex & pat_) : pat(pat_) {}
    
    ex operator()(const ex &e) {
        if(is_a<add>(e)) {
            ex res=0, npat=0;
            for(auto item : e) {
                if(!item.has(pat)) npat += item;
                else res += item;
            }
            if(is_zero(res)) {
                if(!is_zero(npat)) return CCF(npat);
                else return 0;
            }
            if(!is_zero(npat)) res += CCF(npat);
            if(!is_a<add>(res)) {
                cerr << RED << "res is NOT add: " << res << RESET << endl;
                exit(1);
            }
            return res.map(*this);
        } else if (is_a<mul>(e)) {
            ex res=1, npat=1;
            for(auto item : e) {
                if(!item.has(pat)) npat *= item;
                else res *= item;
            }
            if(is_zero(res-1)) {
                if(!is_zero(npat-1)) return CCF(npat);
                else return 1;
            }
            if(!is_zero(npat-1)) res *= CCF(npat);
            if(!is_a<mul>(res)) {
                cerr << RED << "res is NOT mul: " << res << RESET << endl;
                exit(1);
            }
            return res.map(*this).expand();
        } else if(is_a<power>(e) && e.op(1).info(info_flags::nonnegint)) {
            return e.map(*this).expand();
        }
        return e;
    }
};

ex mma_expand(ex expr_in, ex pat) {
    map_CCF ccf(pat);
    auto expr = ccf(expr_in);
    expr = expr.subs(CCF(w)==w);
    return expr;
}

/*-----------------------------------------------------*/
// mma_collect
/*-----------------------------------------------------*/
ex mma_collect(ex expr_in, ex pat, bool ccf, bool cvf) {
    auto res = mma_expand(expr_in, pat);
    lst items;
    if(is_a<add>(res)) {
        for(auto item : res) items.append(item);
    } else items.append(res);
    
    ex cf = 0;
    map<ex, ex, ex_is_less> vc_map;
    for(auto item : items) {
        if(!item.has(pat)) cf += item;
        else {
            if(is_a<mul>(item)) {
                ex tc = 1, tf = 1;
                for(auto ii : item) {
                    if(!ii.has(pat)) tc *= ii;
                    else tf *= ii;
                }
                vc_map[tf] += tc;
            } else {
                vc_map[item] += 1;
            }
        }
    }
    res = 0;
    if(!is_zero(cf)) res += CCF(cf);
    for(auto vc : vc_map) {
        res += CVF(vc.first) * CCF(vc.second);
    }
    
    if(!ccf) res = res.subs(CCF(w)==w);
    if(!cvf) res = res.subs(CVF(w)==w);
    return res;
}

/*-----------------------------------------------------*/
// Evalf
/*-----------------------------------------------------*/
ex Evalf(ex expr) {
    exset zs;
    //patterns needing evalf()
    expr.find(zeta(w), zs);
    expr.find(zeta(w,w), zs);
    
    lst repl;
    auto dd = Digits;
    Digits = 50;
    for(auto zi : zs) {
        repl.append(zi==zi.evalf());
    }
    Digits = dd;
    return expr.subs(repl);
}

/*-----------------------------------------------------*/
// xPositive & xSign
/*-----------------------------------------------------*/
bool xPositive(ex const expr) {
    auto tmp = expr.expand();
    if(tmp.is_zero()) return true;
    bool ret = false;
    if(is_a<add>(tmp)) {
        for(auto item : tmp) {
            auto nit = item.subs(x(w)==1).normal();
            if(!(is_a<numeric>(nit) && ex_to<numeric>(nit).is_positive())) {
                return false;
            }
        }
        ret = true;
    } else {
        auto ntmp = tmp.subs(x(w)==1).normal();
        ret = (is_a<numeric>(ntmp) && ex_to<numeric>(ntmp).is_positive());
    }
    return ret;
}

int xSign(ex const expr) {
    if(xPositive(expr)) return 1;
    else if(xPositive(ex(0)-expr)) return -1;
    else return 0;
}

/*-----------------------------------------------------*/
// let_op extension
/*-----------------------------------------------------*/
void let_op_append(ex & ex_in, const ex item) {
    auto tmp = ex_to<lst>(ex_in);
    tmp.append(item);
    ex_in = tmp;
}
void let_op_prepend(ex & ex_in, const ex item) {
    auto tmp = ex_to<lst>(ex_in);
    tmp.prepend(item);
    ex_in = tmp;
}
void let_op_remove_last(ex & ex_in) {
    auto tmp = ex_to<lst>(ex_in);
    tmp.remove_last();
    ex_in = tmp;
}
void let_op_remove_first(ex & ex_in) {
    auto tmp = ex_to<lst>(ex_in);
    tmp.remove_first();
    ex_in = tmp;
}

void let_op_append(ex & ex_in, int index, ex const item) {
    auto tmp = ex_to<lst>(ex_in.op(index));
    tmp.append(item);
    ex_in.let_op(index) = tmp;
}
void let_op_append(lst & ex_in, int index, ex const item) {
    auto tmp = ex_to<lst>(ex_in.op(index));
    tmp.append(item);
    ex_in.let_op(index) = tmp;
}
void let_op_append(ex & ex_in, int index1, int index2, ex const item) {
    auto tmp = ex_to<lst>(ex_in.op(index1).op(index2));
    tmp.append(item);
    ex_in.let_op(index1).let_op(index2) = tmp;
}
void let_op_append(lst & ex_in, int index1, int index2, ex const item) {
    auto tmp = ex_to<lst>(ex_in.op(index1).op(index2));
    tmp.append(item);
    ex_in.let_op(index1).let_op(index2) = tmp;
}

void let_op_prepend(ex & ex_in, int index, ex const item) {
    auto tmp = ex_to<lst>(ex_in.op(index));
    tmp.prepend(item);
    ex_in.let_op(index) = tmp;
}
void let_op_prepend(lst & ex_in, int index, ex const item) {
    auto tmp = ex_to<lst>(ex_in.op(index));
    tmp.prepend(item);
    ex_in.let_op(index) = tmp;
}
void let_op_prepend(ex & ex_in, int index1, int index2, ex const item) {
    auto tmp = ex_to<lst>(ex_in.op(index1).op(index2));
    tmp.prepend(item);
    ex_in.let_op(index1).let_op(index2) = tmp;
}
void let_op_prepend(lst & ex_in, int index1, int index2, ex const item) {
    auto tmp = ex_to<lst>(ex_in.op(index1).op(index2));
    tmp.prepend(item);
    ex_in.let_op(index1).let_op(index2) = tmp;
}

void let_op_remove_last(ex & ex_in, int index) {
    auto tmp = ex_to<lst>(ex_in.op(index));
    tmp.remove_last();
    ex_in.let_op(index) = tmp;
}
void let_op_remove_last(lst & ex_in, int index) {
    auto tmp = ex_to<lst>(ex_in.op(index));
    tmp.remove_last();
    ex_in.let_op(index) = tmp;
}
void let_op_remove_last(ex & ex_in, int index1, int index2) {
    auto tmp = ex_to<lst>(ex_in.op(index1).op(index2));
    tmp.remove_last();
    ex_in.let_op(index1).let_op(index2) = tmp;
}
void let_op_remove_last(lst & ex_in, int index1, int index2) {
    auto tmp = ex_to<lst>(ex_in.op(index1).op(index2));
    tmp.remove_last();
    ex_in.let_op(index1).let_op(index2) = tmp;
}

void let_op_remove_first(ex & ex_in, int index) {
    auto tmp = ex_to<lst>(ex_in.op(index));
    tmp.remove_first();
    ex_in.let_op(index) = tmp;
}
void let_op_remove_first(lst & ex_in, int index) {
    auto tmp = ex_to<lst>(ex_in.op(index));
    tmp.remove_first();
    ex_in.let_op(index) = tmp;
}
void let_op_remove_first(ex & ex_in, int index1, int index2) {
    auto tmp = ex_to<lst>(ex_in.op(index1).op(index2));
    tmp.remove_first();
    ex_in.let_op(index1).let_op(index2) = tmp;
}
void let_op_remove_first(lst & ex_in, int index1, int index2) {
    auto tmp = ex_to<lst>(ex_in.op(index1).op(index2));
    tmp.remove_first();
    ex_in.let_op(index1).let_op(index2) = tmp;
}

void let_op(ex &ex_in, int index1, int index2, const ex item) {
    ex_in.let_op(index1).let_op(index2) = item;
}
void let_op(lst &ex_in, int index1, int index2, const ex item) {
    ex_in.let_op(index1).let_op(index2) = item;
}
void let_op(ex &ex_in, int index1, int index2, int index3, const ex item) {
    ex_in.let_op(index1).let_op(index2).let_op(index3) = item;
}
void let_op(lst &ex_in, int index1, int index2, int index3, const ex item) {
    ex_in.let_op(index1).let_op(index2).let_op(index3) = item;
}

ex get_op(const ex ex_in, int index1, int index2) {
    return ex_in.op(index1).op(index2);
}
ex get_op(const lst ex_in, int index1, int index2) {
    return ex_in.op(index1).op(index2);
}
ex get_op(const ex ex_in, int index1, int index2, int index3) {
    return ex_in.op(index1).op(index2).op(index3);
}
ex get_op(const lst ex_in, int index1, int index2, int index3) {
    return ex_in.op(index1).op(index2).op(index3);
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

REGISTER_FUNCTION(FF, dummy())
REGISTER_FUNCTION(CV, do_not_evalf_params())

REGISTER_FUNCTION(x, dummy())
REGISTER_FUNCTION(y, dummy())
REGISTER_FUNCTION(z, dummy())


}

