#include "Basic.h"

namespace HepLib {

    //-----------------------------------------------------------
    // Symbol Class
    //-----------------------------------------------------------
    DEFAULT_CTOR(Symbol)
    GINAC_BIND_UNARCHIVER(Symbol);
    IMPLEMENT_HAS(Symbol)
    
    GINAC_IMPLEMENT_REGISTERED_CLASS(Symbol, symbol)
    
    Symbol::Symbol(const string &s, bool is_real) : symbol(get_symbol(s)), isReal(is_real) { }
    int Symbol::compare_same_type(const basic &other) const {
        const Symbol &o = static_cast<const Symbol &>(other);
        int ret = get_name().compare(o.get_name());
        if(ret==0) return 0;
        else if(ret<0) return -1;
        else return 1;
    }
    
    void Symbol::archive(archive_node & n) const {
        inherited::archive(n);
        n.add_bool("isReal", isReal);
    }
    
    void Symbol::read_archive(const archive_node& n, lst& sym_lst) {
        inherited::read_archive(n, sym_lst);
        bool is_real;
        n.find_bool("isReal", is_real);
        *this = Symbol(get_name(), is_real);
    }
    
    ex Symbol::eval() const { return *this; }
    ex Symbol::evalf() const { return *this; }
    ex Symbol::conjugate() const { return *this * (isReal ? 1 : -1); }
    ex Symbol::real_part() const { return (isReal ? *this : ex(0)); }
    ex Symbol::imag_part() const { return (isReal ? ex(0) : *this); }

    /*-----------------------------------------------------*/
    // Global varibales
    /*-----------------------------------------------------*/
    
    const char* Color_Error = RED;
    const char* Color_Warn = MAGENTA;
    const char* Color_HighLight = WHITE;

    /*-----------------------------------------------------*/
    // GiNaC_Parallel
    /*-----------------------------------------------------*/
    void GiNaC_archive_Symbols_from(ex expr) {
        auto syms = gather_symbols(expr);
        for(auto si : syms) GiNaC_archive_Symbols.append(si);
        GiNaC_archive_Symbols.sort();
        GiNaC_archive_Symbols.unique();
    }

    void GiNaC_archive_Symbols_from(vector<ex> invec) {
        auto syms = gather_symbols(invec);
        for(auto si : syms) GiNaC_archive_Symbols.append(si);
        GiNaC_archive_Symbols.sort();
        GiNaC_archive_Symbols.unique();
    }
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
                if(key != NULL) cout << Color_HighLight << key << RESET << " ";
                cout << Color_HighLight << batch << "x" << RESET << "[" << (bi+1) << "/" << btotal << "] ... " << flush;
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
                cout << Color_Error << "Failed in GiNaC_Parallel!" << RESET << endl;
                cout << Color_Error << p.what() << RESET << endl;
                if(para_max_run>0) exit(0);
                throw p;
            }
            if(para_max_run>0) exit(0);
        }
        
        auto cpid = getpid();
        if(cpid!=ppid) exit(0); // make sure
        if(para_max_run>0) while (wait(NULL) != -1) { }
        if(verb > 1 && total > 0) cout << "@" << now(false) << endl;

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
                    cout << "\\--Reading *." << Color_HighLight << key << RESET << ".gar [" << (i+1) << "/" << total << "] ... " << flush;
                }
            }

            int oDigits = Digits;
            Digits = 50; // a fix to float overflow
            
            archive ar;
            ostringstream garfn;
            if(key == NULL) garfn << ppid << "/" << i << ".gar";
            else garfn << ppid << "/" << i << "." << key << ".gar";
            if(!file_exists(garfn.str().c_str())) {
                cerr << Color_Error << "GiNaC_Parallel: Check the error message above." << RESET << endl;
                exit(0);
            }
            ifstream ins(garfn.str().c_str());
            ins >> ar;
            ins.close();
            remove(garfn.str().c_str());
            auto c = ar.unarchive_ex(GiNaC_archive_Symbols, "c");
            if(c!=19790923) throw runtime_error("*.gar error!");
            auto res = ar.unarchive_ex(GiNaC_archive_Symbols, "res");
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
            cerr << Color_Error << "r>=mat.rows()" << RESET << endl;
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
    void garRead(const string &garfn, map<string, ex> &resMap) {
        archive ar;
        ifstream in(garfn);
        in >> ar;
        in.close();
        for(int i=0; i<ar.num_expressions(); i++) {
            string name;
            ex res = ar.unarchive_ex(GiNaC_archive_Symbols, name, i);
            resMap[name] = res;
        }
    }
    
    void garWrite(const string &garfn, const map<string, ex> &resMap) {
        archive ar;
        for(const auto & item : resMap) {
            ar.archive_ex(item.second, item.first.c_str());
        }
        ofstream out(garfn);
        out << ar;
        out.close();
    }

    ex garRead(const string &garfn, const char* key) {
        archive ar;
        ifstream in(garfn);
        in >> ar;
        in.close();
        auto res = ar.unarchive_ex(GiNaC_archive_Symbols, key);
        return res;
    }

    ex garResult(const string &garfn) {
        archive ar;
        ifstream in(garfn);
        in >> ar;
        in.close();
        auto c = ar.unarchive_ex(GiNaC_archive_Symbols, "c");
        auto res = ar.unarchive_ex(GiNaC_archive_Symbols, "res");
        if(c!=19790923) {
            cerr << Color_Error << "gar file: " << garfn << endl;
            cerr << "c=" << c << ", different from 19790923!" << RESET << endl;
            exit(1);
        }
        return res;
    }

    /*-----------------------------------------------------*/
    // str2ex Function
    /*-----------------------------------------------------*/
    ex str2ex(const string &expr, symtab stab) {
        parser reader(stab);
        ex ret = reader(expr);
        return ret;
    }

    /*-----------------------------------------------------*/
    // str2lst Function
    /*-----------------------------------------------------*/
    lst str2lst(const string &expr, symtab stab) {
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
    ex mma_series(ex const & expr_in, symbol const &s0, int sn0) {
        ex expr = expr_in;
        if(!expr.has(s0)) return expr;
        
        exset sset;
        expr.find(pow(s0, w), sset);
        numeric sn_lcm = 1;
        for(auto pi : sset) {
            auto sn = pi.op(1);
            if(!(is_a<numeric>(sn) && ex_to<numeric>(sn).is_rational())) {
                cerr << Color_Error << "mma_series: Not rational sn = " << sn << RESET << endl;
                cerr << "s = " << s0 << endl;
                cerr << "expr_in = " << expr_in << endl;
                exit(1);
            }
            sn_lcm = lcm(sn_lcm, ex_to<numeric>(sn).denom());
        }
        symbol s;
        if(!sn_lcm.is_integer()) {
            cerr << Color_Error << "Not integer: " << sn_lcm << RESET << endl;
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
            res = res.subs(coCF(w)==w); // remove coCF
            ex ot = 0;
            for(int i=0; i<res.nops(); i++) {
                if(is_order_function(res.op(i))) {
                    ot = res.op(i);
                    break;
                }
            }
            if(!is_order_function(ot)) {
                cerr << Color_Error << "Not an Order term: " << ot << RESET << endl;
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
        cerr << Color_Error << "mma_series seems not working!" << RESET << endl;
        exit(1);
        return 0;
    }

    /*-----------------------------------------------------*/
    // mma_diff
    /*-----------------------------------------------------*/
    ex mma_diff(ex const expr, ex const xp, unsigned nth, bool expand) {
        symbol s;
        ex res = expr.subs(xp==s);
        if(expand) res = mma_collect(res, s, true);
        res = res.diff(s, nth);
        if(expand) res = res.subs(coCF(w)==w); // remove coCF
        res = res.subs(s==xp);
        return res;
    }

    /*-----------------------------------------------------*/
    // mma_expand
    /*-----------------------------------------------------*/
    ex mma_expand(ex const &expr_in, std::function<bool(const ex &)> isOK, int depth) {
        if(depth>5) return expr_in.expand();
        ex expr;
        if(is_a<add>(expr_in)) {
            expr = 0;
            for(auto item : expr_in) {
                if(isOK(item)) expr += mma_expand(item, isOK, depth+1);
                else expr += item;
            }
        } else if (is_a<mul>(expr_in)) {
            expr=1;
            for(auto item : expr_in) {
                if(isOK(item)) expr *= mma_expand(item, isOK, depth+1);
                else if(is_a<numeric>(item)) expr *= item;
                else expr *= coCF(item.subs(coCF(w)==w));
            }
        } else if(is_a<power>(expr_in) && expr_in.op(1).info(info_flags::nonnegint)) {
            auto item = expr_in.op(0);
            auto ni = expr_in.op(1);
            if(isOK(item)) expr = pow(mma_expand(item, isOK, depth+1), ni).expand();
            else if(is_a<numeric>(item)) expr = expr_in;
            else expr = coCF(expr_in.subs(coCF(w)==w));
        } else {
            if(isOK(expr_in) || is_a<numeric>(expr_in)) expr = expr_in;
            else expr = coCF(expr_in.subs(coCF(w)==w));
        }
        
        expr = expr.expand();
        ex res = expr;
        if(is_a<add>(expr)) {
            res = 0;
            ex ccf_expr=0;
            for(auto item : expr) {
                if(isOK(item)) res += item;
                else ccf_expr += item;
            }
            if(!is_zero(ccf_expr)) res += coCF(ccf_expr.subs(coCF(w)==w));
        }
        if(depth==0) res = res.subs(coCF(w)==w);
        return res;
    }
    ex mma_expand(ex const &expr_in, lst const &pats, int depth) {
        return mma_expand(expr_in, [pats](const ex & e)->bool {
            for(auto pat : pats) {
                if(e.has(pat)) return true;
            }
            return false;
        }, depth);
    }
    ex mma_expand(ex const &expr_in, const ex &pat, int depth) {
        return mma_expand(expr_in, [pat](const ex & e)->bool {
            return e.has(pat);
        }, depth);
    }

    /*-----------------------------------------------------*/
    // mma_collect
    /*-----------------------------------------------------*/
    ex mma_collect(ex const &expr_in, std::function<bool(const ex &)> isOK, bool ccf, bool cvf) {
        auto res = mma_expand(expr_in, isOK);
        lst items;
        if(is_a<add>(res)) {
            for(auto item : res) items.append(item);
        } else items.append(res);
        
        ex cf = 0;
        map<ex, ex, ex_is_less> vc_map;
        for(auto item : items) {
            if(!isOK(item)) cf += item;
            else if(is_a<mul>(item)) {
                ex tc = 1, tv = 1;
                for(auto ii : item) {
                    if(!isOK(ii)) tc *= ii;
                    else tv *= ii;
                }
                vc_map[tv] += tc;
            } else {
                vc_map[item] += 1;
            }
        }
        res = 0;
        if(!is_zero(cf)) res += coCF(cf)*coVF(1);
        for(auto vc : vc_map) {
            if(isOK(vc.second)) {
                cerr << Color_Error << "mma_collect: pats founds @ " << vc.second << RESET << endl;
                exit(1);
            }
            res += coVF(vc.first) * coCF(vc.second);
        }
        
        if(!ccf) res = res.subs(coCF(w)==w);
        if(!cvf) res = res.subs(coVF(w)==w);
        return res;
    }
    
    ex mma_collect(ex const &expr_in, lst const &pats, bool ccf, bool cvf) {
        return mma_collect(expr_in, [pats](const ex & e)->bool {
            for(auto pat : pats) {
                if(e.has(pat)) return true;
            }
            return false;
        }, ccf, cvf);
    }
    ex mma_collect(ex const &expr_in, ex const &pat, bool ccf, bool cvf) {
        return mma_collect(expr_in, [pat](const ex & e)->bool {
            return e.has(pat);
        }, ccf, cvf);
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


}

