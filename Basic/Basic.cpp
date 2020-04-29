/**
 * @file
 * @brief Basic Functions, extend GiNaC
 * @author F. Feng
 * @version 1.0.0
 * @date 2020-04-20
 */

#include "Basic.h"

namespace HepLib {

    /**
     * @brief Error constructor, 
     * @param _msg error message 
     */
    Error::Error(string _msg) : msg(_msg) { }
    const char * Error::what() const throw () {
        return msg.c_str();
    }

    DEFAULT_CTOR(Symbol)
    GINAC_BIND_UNARCHIVER(Symbol);
    IMPLEMENT_HAS(Symbol)
    
    GINAC_IMPLEMENT_REGISTERED_CLASS(Symbol, symbol)
    
    /**
     * @brief Symbol constructor
     * @param s symbol name
     * @param r true for real, false for pure imaginary 
     * @param c true to check the name exist or not, if exsit error thrown
     */
    Symbol::Symbol(const string &s, bool r, bool c) : symbol(get_symbol(s,c)), isReal(r) { Table[s]=*this; }
    int Symbol::compare_same_type(const basic &other) const {
        const Symbol &o = static_cast<const Symbol &>(other);
        int ret = get_name().compare(o.get_name());
        if(ret==0) return 0;
        else if(ret<0) return -1;
        else return 1;
    }
    
    /**
     * @brief Symbol archive
     * @param n archive node
     */
    void Symbol::archive(archive_node & n) const {
        inherited::archive(n);
        n.add_bool("isReal", isReal);
    }
    
    /**
     * @brief Symbol read_archive
     * @param n archive node
     * @param sym_lst symbol lst
     */
    void Symbol::read_archive(const archive_node& n, lst& sym_lst) {
        inherited::read_archive(n, sym_lst);
        bool is_real;
        n.find_bool("isReal", is_real);
        *this = Symbol(get_name(), is_real, false);
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

    /**
     * @brief update GiNaC_archive_Symbols from expr
     * @param expr input expression, symbol in expr will be added to GiNaC_archive_Symbols
     */
    void GiNaC_archive_Symbols_from(ex expr) {
        auto syms = gather_symbols(expr);
        for(auto si : syms) GiNaC_archive_Symbols.append(si);
        GiNaC_archive_Symbols.sort();
        GiNaC_archive_Symbols.unique();
    }

    /**
     * @brief update GiNaC_archive_Symbols from expr
     * @param invec input expression, symbol in the vector will be added to GiNaC_archive_Symbols
     */
    void GiNaC_archive_Symbols_from(vector<ex> invec) {
        auto syms = gather_symbols(invec);
        for(auto si : syms) GiNaC_archive_Symbols.append(si);
        GiNaC_archive_Symbols.sort();
        GiNaC_archive_Symbols.unique();
    }
    
    /**
     * @brief GiNaC Parallel Evaluation using fork
     * @param ntotal the number of total items, 0 for non-parallel version
     * @param nbatch the batch number for each run, if nbatch<=0, then nbatch = ntotal/para_max_run/10
     * @param f function to be applied on the index, from 0 to (ntotal-1)
     * @param key key used in archive file name and display message
     * @param rm default true, and false will keep the archive file
     * @param prtlvl the indent level in the print message, when the Verbose>1
     * @return return the ntotal-element vector, i.e., [f(0), ..., f(ntotal-1)]
     */
    vector<ex> GiNaC_Parallel(
        int ntotal, int nbatch, 
        std::function<ex(int)> f,
        const string & key, bool rm, int prtlvl) {
        int nproc = ParallelProcess;
        
        // nproc=0, non-parallel
        if(nproc==0) {
            vector<ex> ovec;
            for(int i=0; i<ntotal; i++) ovec.push_back(f(i));
            return ovec;
        }
        
        auto ppid = getpid();
        int para_max_run = nproc<0 ? omp_get_num_procs()-1 : nproc;
        if(para_max_run<1) para_max_run = 1;
        ostringstream cmd;
        cmd << "mkdir -p " << ppid;
        system(cmd.str().c_str());
        
        if(nbatch<=0) nbatch = ntotal/para_max_run/10;
        if(nbatch<1) nbatch = 1;
        int btotal = ntotal/nbatch + ((ntotal%nbatch)==0 ? 0 : 1);

        for(int bi=0; bi<btotal; bi++) {
            if(Verbose > 1) {
                cout << "\r  ";
                for(int pi=0;pi<prtlvl;pi++) cout << "   ";
                cout << "\\--Evaluating ";
                if(key != "") cout << Color_HighLight << key << RESET << " ";
                cout << Color_HighLight << nbatch << "x" << RESET << "[" << (bi+1) << "/" << btotal << "] ... " << flush;
            }
            
            auto pid = fork();
            if (pid < 0) perror("fork() error");
            else if (pid != 0) {
                if(bi >= para_max_run) wait(NULL);
                continue;
            } 
            
            // child process
            PID = getpid(); // update PID after fork
            try {
                lst res_lst;
                for(int ri=0; ri<nbatch; ri++) {
                    int i = bi*nbatch + ri;
                    if(i<ntotal) res_lst.append(f(i));
                    else break;
                }
                ostringstream garfn;
                if(key == "") garfn << ppid << "/" << bi << ".gar";
                else garfn << ppid << "/" << bi << "." << key << ".gar";
                garWrite(garfn.str(), res_lst);
            } catch(exception &p) {
                cout << Color_Error << "Failed in GiNaC_Parallel!" << RESET << endl;
                cout << Color_Error << p.what() << RESET << endl;
            }
            exit(0);
        }
        
        auto cpid = getpid();
        if(cpid!=ppid) exit(0); // make sure child exit
        while (wait(NULL) != -1) { }
        if(Verbose > 1 && ntotal > 0) cout << "@" << now(false) << endl;

        vector<ex> ovec;
        for(int bi=0; bi<btotal; bi++) {
            if(Verbose > 1) {
                if(key == "") {
                    cout << "\r  ";
                    for(int pi=0; pi<prtlvl; pi++) cout << "   ";
                    cout << "\\--Reading *.gar [" << (bi+1) << "/" << btotal << "] ... " << flush;
                } else {
                    cout << "\r  ";
                    for(int pi=0;pi<prtlvl;pi++) cout << "   ";
                    cout << "\\--Reading *." << Color_HighLight << key << RESET << ".gar [" << (bi+1) << "/" << btotal << "] ... " << flush;
                }
            }

            int oDigits = Digits;
            Digits = 50; // a fix to float overflow
            
            archive ar;
            ostringstream garfn;
            if(key == "") garfn << ppid << "/" << bi << ".gar";
            else garfn << ppid << "/" << bi << "." << key << ".gar";
            if(!file_exists(garfn.str().c_str())) {
                cerr << Color_Error << "GiNaC_Parallel: Check the error message above." << RESET << endl;
                exit(0);
            }
            auto res_lst = garRead(garfn.str());
            remove(garfn.str().c_str());
            for(auto res : res_lst) ovec.push_back(res);
            Digits = oDigits;
        }
        
        if(rm) {
            cmd.clear();
            cmd.str("");
            cmd << "rm -fr " << ppid;
            system(cmd.str().c_str());
            system(cmd.str().c_str());
        }
        if(Verbose > 1 && ntotal > 0) cout << "@" << now(false) << endl;
        return ovec;
    }

    /**
     * @brief get the symbol from symbol factory, if not exsit, a new one will be created
     * @param s the name of the symbol
     * @param check true to check exist, and if so, error will be thrown
     * @return return the found/created symbol
     */
    const symbol & get_symbol(const string & s, bool check) {
        static map<string, symbol> directory;
        map<string, symbol>::iterator i = directory.find(s);
        if (i != directory.end()) {
            if(check) throw std::runtime_error("get_symbol: check failed with name: " + s);
            return i->second;
        } else {
            return directory.insert(make_pair(s, symbol(s))).first->second;
        }
    }

    /**
     * @brief split the string into serveral part, separated by the delimiter
     * @param s the input string
     * @param delimiter the char delimiter
     * @return return separated string vector
     */
    std::vector<std::string> split(const std::string& s, char delimiter) {
        std::vector<std::string> tokens;
        std::string token;
        std::istringstream tokenStream(s);
        while (std::getline(tokenStream, token, delimiter)) {
            tokens.push_back(token);
        }
        return tokens;
    }

    /**
     * @brief chech matrix has zero row
     * @param mat the input matrix
     * @return return matrix has zero row or not
     */
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

    /**
     * @brief chech r-th row of matrix is zero or not
     * @param mat the input matrix
     * @param r refers to the (r+1)-th row
     * @return return r-th row of matrix is zero or not
     */
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

    /**
     * @brief the all zero row's index, start from 0
     * @param mat the input matrix
     * @return return index vector for all zero rows
     */
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

    /**
     * @brief remove the zero row
     * @param mat the input matrix
     * @return return the matrix will zero row removed
     */
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

    /**
     * @brief part of the matrix, submatrix, row: r->r+nr, c->c+nc
     * @param mat the input matrix
     * @param r start row index, start from 0
     * @param c start column index, start from 0
     * @param nr the number of rows
     * @param nc the number of columns
     * @return return index vector for all zero rows
     */
    matrix MatHelper::sub(matrix mat, int r, int nr, int c, int nc) {
        matrix ret(nr, nc);
        for(int ir=0; ir<nr; ir++) {
            for(int ic=0; ic<nc; ic++) {
                ret(ir, ic) = mat(r+ir, c+ic);
            }
        }
        return ret;
    }

    /**
     * @brief sum all elements in the list
     * @param m the input lst
     * @return return the sum
     */
    ex lstHelper::sum(lst m) {
        ex ret = 0;
        for(int i=0; i<m.nops(); i++) ret += m.op(i);
        return ret;
    }

    /**
     * @brief date/time string
     * @param use_date defualt true to include the date in the string
     * @return return date/time string
     */
    string now(bool use_date) {
        time_t timep;
        time (&timep);
        char tmp[64];
        if(use_date) strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&timep) );
        else strftime(tmp, sizeof(tmp), "%H:%M:%S",localtime(&timep) );
        return tmp;
    }

    /**
     * @brief get all symbols from input expression
     * @param e input expression
     * @return all symbols in the input
     */
    lst gather_symbols(const ex & e) {
        lst sym_lst;
        for(const_preorder_iterator i = e.preorder_begin(); i != e.preorder_end(); ++i) {
            if(is_a<symbol>(*i)) sym_lst.append(*i);
        }
        sym_lst.sort();
        sym_lst.unique();
        return sym_lst;
    }

    /**
     * @brief get all symbols from input expression
     * @param ve input expression vector
     * @return all symbols in the input
     */
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

    /**
     * @brief run symtem command and return the output as string
     * @param cmd symtem command
     * @return the command output
     */
    string RunOS(const string & cmd) {
        char buf[128];
        ostringstream oss;
        FILE *fp = NULL;
        if( (fp = popen(cmd.c_str(), "r")) != NULL) {
            while(fgets(buf, 128, fp) != NULL) {
                oss << buf;
            }
            pclose(fp);
            fp = NULL;
        }
        return oss.str();
    }


    /**
     * @brief garRead from file, and output in a map
     * @param garfn ginac archive filename
     * @param resMap will be update, a string key map
     */
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
    
    /**
     * @brief garRead from file, only the element w.r.t. key 
     * @param garfn ginac archive filename
     * @param key the archive key, note (const char *) type, match unarchive_ex in GiNaC
     * @return the expression w.r.t. key
     */
    ex garRead(const string &garfn, const char* key) { // use the const char *, not string
        archive ar;
        ifstream in(garfn);
        in >> ar;
        in.close();
        auto res = ar.unarchive_ex(GiNaC_archive_Symbols, key);
        return res;
    }

    /**
     * @brief garRead from file, only the element w.r.t. key "res", note inner check will be performed
     * @param garfn ginac archive filename
     * @return the expression w.r.t. key "res"
     */
    ex garRead(const string &garfn) {
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
    
    /**
     * @brief garWrite to write the string-key map to the archive
     * @param garfn ginac archive filename for output
     */
    void garWrite(const string &garfn, const map<string, ex> &resMap) {
        archive ar;
        for(const auto & item : resMap) {
            ar.archive_ex(item.second, item.first.c_str());
        }
        ofstream out(garfn);
        out << ar;
        out.close();
    }

    /**
     * @brief garWrite to write the expression to the archive, with key: "res", including a check key will be written
     * @param garfn ginac archive filename for output
     * @param res the output expression
     */
    void garWrite(const string &garfn, const ex & res) {
        archive ar;
        ar.archive_ex(res, "res");
        ar.archive_ex(19790923, "c");
        ofstream out(garfn);
        out << ar;
        out.close();
    }


    /**
     * @brief convert string to ex expression, using Parser internally
     * @param expr expression in string format
     * @param stab the symtab
     * @return the parsed expression
     */
    ex str2ex(const string &expr, symtab stab) {
        Parser par(stab);
        ex ret = par.Read(expr);
        return ret;
    }

    /**
     * @brief convert string to the lst, using Parser internally
     * @param expr a lst in string format
     * @param stab the symtab
     * @return the parsed expression
     */
    lst str2lst(const string &expr, symtab stab) {
        Parser par(stab);
        ex ret = par.Read(expr);
        return ex_to<lst>(ret);
    }
    
    /**
     * @brief convert ex to output string, the defalut printer format will be used
     * @param expr a expression
     * @return the output string
     */
    string ex2str(const ex &expr) {
        ostringstream oss;
        oss << expr << endl;
        return oss.str();
    }
    
    /**
     * @brief convert exvector to lst
     * @param exvec input exvector
     * @return lst
     */
    lst exvec2lst(const exvector & exvec) {
        lst ret;
        for(auto item : exvec) ret.append(item);
        return ret;
    }
    
    /**
     * @brief convert lst to exvector
     * @param exvec input lst
     * @return exvector
     */
    exvector lst2exvec(const lst & alst) {
        exvector ret;
        for(auto item : alst) ret.push_back(item);
        return ret;
    }   

    /**
     * @brief return a lst: x(bi), x(bi+1), ..., x(ei)
     * @param bi start index
     * @param ei end index
     * @return the x lst
     */
    lst xlst(int bi, int ei) {
        lst ret;
        for(int i=bi; i<=ei; i++) ret.append(x(i));
        return ret;
    }
    
    /**
     * @brief return a lst: x(0), x(1), ..., x(ei)
     * @param ei end index
     * @return the x lst
     */
    lst xlst(int ei) {
        return xlst(0, ei);
    }

    /**
     * @brief the series like Mathematica, include s^n
     * @param expr_in input expression
     * @param s0 the variable
     * @param sn0 expanded upto sn0 order, include s0^sn0
     * @return the corresponding series
     */
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

    /**
     * @brief the differential like Mathematica
     * @param expr input expression
     * @param xp the variable, can be an expression, will replace by a symbol and back again
     * @param nth nth-derivative
     * @param expand true to call mma_expand before diff
     * @return the corresponding nth-derivative
     */
    ex mma_diff(ex const expr, ex const xp, unsigned nth, bool expand) {
        symbol s;
        ex res = expr.subs(xp==s);
        if(expand) res = mma_collect(res, s, true);
        res = res.diff(s, nth);
        if(expand) res = res.subs(coCF(w)==w); // remove coCF
        res = res.subs(s==xp);
        return res;
    }

    /**
     * @brief the expand like Mathematica
     * @param expr_in input expression
     * @param isOK only expand the element e, when isOK(e) is true
     * @param depth will be incresed by 1 when each recursively called, when depth>5, GiNaC expand will be used
     * @return the expanded expression
     */
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
    
    /**
     * @brief the expand like Mathematica
     * @param expr_in input expression
     * @param pats only expand the element e, when e has at least one pattern in pats is true
     * @param depth will be incresed by 1 when each recursively called
     * @return the expanded expression
     */
    ex mma_expand(ex const &expr_in, lst const &pats, int depth) {
        return mma_expand(expr_in, [pats](const ex & e)->bool {
            for(auto pat : pats) {
                if(e.has(pat)) return true;
            }
            return false;
        }, depth);
    }
    
    /**
     * @brief the expand like Mathematica
     * @param expr_in input expression
     * @param pat only expand the element e, when e has the pattern: pat
     * @param depth will be incresed by 1 when each recursively called
     * @return the expanded expression
     */
    ex mma_expand(ex const &expr_in, const ex &pat, int depth) {
        return mma_expand(expr_in, [pat](const ex & e)->bool {
            return e.has(pat);
        }, depth);
    }

    /**
     * @brief the collect function like Mathematica
     * @param expr_in input expression
     * @param isOK only collect the element e, when isOK(e) is true
     * @param ccf true for wrapping coefficient in coCF
     * @param cvf true to wrapping isOK-ed element in coVF
     * @return the collected expression
     */
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
    
    /**
     * @brief the collect function like Mathematica
     * @param expr_in input expression
     * @param pats only collect the element e match at least one patter in pats, like mma_expand
     * @param ccf true for wrapping coefficient in coCF
     * @param cvf true to wrapping isOK-ed element in coVF
     * @return the collected expression
     */
    ex mma_collect(ex const &expr_in, lst const &pats, bool ccf, bool cvf) {
        return mma_collect(expr_in, [pats](const ex & e)->bool {
            for(auto pat : pats) {
                if(e.has(pat)) return true;
            }
            return false;
        }, ccf, cvf);
    }
    
    /**
     * @brief the collect function like Mathematica
     * @param expr_in input expression
     * @param pat only collect the element e match the pattern, like mma_expand
     * @param ccf true for wrapping coefficient in coCF
     * @param cvf true to wrapping isOK-ed element in coVF
     * @return the collected expression
     */
    ex mma_collect(ex const &expr_in, ex const &pat, bool ccf, bool cvf) {
        return mma_collect(expr_in, [pat](const ex & e)->bool {
            return e.has(pat);
        }, ccf, cvf);
    }

    /**
     * @brief the nuerical evaluation, Digits=50 will be used
     * @param expr input expression
     * @return the nuerical expression
     */
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

    /**
     * @brief check the expr is xPositive, i.e., each x-monomial item is postive
     * @param expr input expression
     * @return xPositive or not
     */
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

    /**
     * @brief the always sign for expr
     * @param expr input expression
     * @return 1 for xPositive exprs, -1 when -expr is xPositive, 0 for others
     */
    int xSign(ex const expr) {
        if(xPositive(expr)) return 1;
        else if(xPositive(ex(0)-expr)) return -1;
        else return 0;
    }

    /**
     * @brief append item into expression
     * @param ex_in input expression will be update
     * @param item element to be appended to ex_in
     */
    void let_op_append(ex & ex_in, const ex item) {
        auto tmp = ex_to<lst>(ex_in);
        tmp.append(item);
        ex_in = tmp;
    }
    
    /**
     * @brief preppend item into expression
     * @param ex_in input expression will be update
     * @param item element to be prepended to ex_in
     */
    void let_op_prepend(ex & ex_in, const ex item) {
        auto tmp = ex_to<lst>(ex_in);
        tmp.prepend(item);
        ex_in = tmp;
    }
    
    /**
     * @brief remove last from expression
     * @param ex_in input lst and will be update, last element will be remove
     */
    void let_op_remove_last(ex & ex_in) {
        auto tmp = ex_to<lst>(ex_in);
        tmp.remove_last();
        ex_in = tmp;
    }
    
    /**
     * @brief remove first from expression
     * @param ex_in input lst and will be update, first element will be remove
     */
    void let_op_remove_first(ex & ex_in) {
        auto tmp = ex_to<lst>(ex_in);
        tmp.remove_first();
        ex_in = tmp;
    }

    /**
     * @brief append item into index-th of expression
     * @param ex_in input expression will be update
     * @param index the index, index-th should be lst, and item will be append into it
     * @param item element to be appended to index-th of ex_in
     */
    void let_op_append(ex & ex_in, int index, ex const item) {
        auto tmp = ex_to<lst>(ex_in.op(index));
        tmp.append(item);
        ex_in.let_op(index) = tmp;
    }
    
    /**
     * @brief append item into index-th of expression
     * @param ex_in input expression will be update
     * @param index the index, index-th should be lst, and item will be append into it
     * @param item element to be appended to index-th of ex_in
     */
    void let_op_append(lst & ex_in, int index, ex const item) {
        auto tmp = ex_to<lst>(ex_in.op(index));
        tmp.append(item);
        ex_in.let_op(index) = tmp;
    }
    
    /**
     * @brief append item into index1-th.index2-th of expression
     * @param ex_in input expression will be update
     * @param index1 the first index
     * @param index2 the second index
     * @param item element to be appended to index1-th.index2-th of ex_in
     */
    void let_op_append(ex & ex_in, int index1, int index2, ex const item) {
        auto tmp = ex_to<lst>(ex_in.op(index1).op(index2));
        tmp.append(item);
        ex_in.let_op(index1).let_op(index2) = tmp;
    }
    
    /**
     * @brief append item into index1-th.index2-th of expression
     * @param ex_in input expression will be update
     * @param index1 the first index
     * @param index2 the second index
     * @param item element to be appended to index1-th.index2-th of ex_in
     */
    void let_op_append(lst & ex_in, int index1, int index2, ex const item) {
        auto tmp = ex_to<lst>(ex_in.op(index1).op(index2));
        tmp.append(item);
        ex_in.let_op(index1).let_op(index2) = tmp;
    }

    /**
     * @brief prepend item into index-th of expression
     * @param ex_in input expression will be update
     * @param index the index, index-th should be lst, and item will be append into it
     * @param item element to be prepended to index-th of ex_in
     */
    void let_op_prepend(ex & ex_in, int index, ex const item) {
        auto tmp = ex_to<lst>(ex_in.op(index));
        tmp.prepend(item);
        ex_in.let_op(index) = tmp;
    }
    
    /**
     * @brief prepend item into index-th of expression
     * @param ex_in input expression will be update
     * @param index the index, index-th should be lst, and item will be append into it
     * @param item element to be prepended to index-th of ex_in
     */
    void let_op_prepend(lst & ex_in, int index, ex const item) {
        auto tmp = ex_to<lst>(ex_in.op(index));
        tmp.prepend(item);
        ex_in.let_op(index) = tmp;
    }
    
    /**
     * @brief prepend item into index1-th.index2-th of expression
     * @param ex_in input expression will be update
     * @param index1 the first index
     * @param index2 the second index
     * @param item element to be prepend to index1-th.index2-th of ex_in
     */
    void let_op_prepend(ex & ex_in, int index1, int index2, ex const item) {
        auto tmp = ex_to<lst>(ex_in.op(index1).op(index2));
        tmp.prepend(item);
        ex_in.let_op(index1).let_op(index2) = tmp;
    }
    
    /**
     * @brief prepend item into index1-th.index2-th of expression
     * @param ex_in input expression will be update
     * @param index1 the first index
     * @param index2 the second index
     * @param item element to be prepend to index1-th.index2-th of ex_in
     */
    void let_op_prepend(lst & ex_in, int index1, int index2, ex const item) {
        auto tmp = ex_to<lst>(ex_in.op(index1).op(index2));
        tmp.prepend(item);
        ex_in.let_op(index1).let_op(index2) = tmp;
    }

    /**
     * @brief remove the last in index-th of expression
     * @param ex_in input expression will be update
     * @param index the index
     */
    void let_op_remove_last(ex & ex_in, int index) {
        auto tmp = ex_to<lst>(ex_in.op(index));
        tmp.remove_last();
        ex_in.let_op(index) = tmp;
    }
    
    /**
     * @brief remove the last in index-th of expression
     * @param ex_in input expression will be update
     * @param index the index
     */
    void let_op_remove_last(lst & ex_in, int index) {
        auto tmp = ex_to<lst>(ex_in.op(index));
        tmp.remove_last();
        ex_in.let_op(index) = tmp;
    }
    
    /**
     * @brief remove the last in index1-th.index2-th of expression
     * @param ex_in input expression will be update
     * @param index the index
     */
    void let_op_remove_last(ex & ex_in, int index1, int index2) {
        auto tmp = ex_to<lst>(ex_in.op(index1).op(index2));
        tmp.remove_last();
        ex_in.let_op(index1).let_op(index2) = tmp;
    }
    
    /**
     * @brief remove the last in index1-th.index2-th of expression
     * @param ex_in input expression will be update
     * @param index the index
     */
    void let_op_remove_last(lst & ex_in, int index1, int index2) {
        auto tmp = ex_to<lst>(ex_in.op(index1).op(index2));
        tmp.remove_last();
        ex_in.let_op(index1).let_op(index2) = tmp;
    }

    /**
     * @brief remove the first in index-th of expression
     * @param ex_in input expression will be update
     * @param index the index
     */
    void let_op_remove_first(ex & ex_in, int index) {
        auto tmp = ex_to<lst>(ex_in.op(index));
        tmp.remove_first();
        ex_in.let_op(index) = tmp;
    }
    
    /**
     * @brief remove the first in index-th of expression
     * @param ex_in input expression will be update
     * @param index the index
     */
    void let_op_remove_first(lst & ex_in, int index) {
        auto tmp = ex_to<lst>(ex_in.op(index));
        tmp.remove_first();
        ex_in.let_op(index) = tmp;
    }
    
    /**
     * @brief remove the first in index1-th.index2-th of expression
     * @param ex_in input expression will be update
     * @param index the index
     */
    void let_op_remove_first(ex & ex_in, int index1, int index2) {
        auto tmp = ex_to<lst>(ex_in.op(index1).op(index2));
        tmp.remove_first();
        ex_in.let_op(index1).let_op(index2) = tmp;
    }
    
    /**
     * @brief remove the first in index1-th.index2-th of expression
     * @param ex_in input expression will be update
     * @param index the index
     */
    void let_op_remove_first(lst & ex_in, int index1, int index2) {
        auto tmp = ex_to<lst>(ex_in.op(index1).op(index2));
        tmp.remove_first();
        ex_in.let_op(index1).let_op(index2) = tmp;
    }

    /**
     * @brief update index1-th.index2-th of expression with item
     * @param ex_in input expression will be update
     * @param index1 the index
     * @param index2 the index
     * @param the new item
     */
    void let_op(ex &ex_in, int index1, int index2, const ex item) {
        ex_in.let_op(index1).let_op(index2) = item;
    }
    
    /**
     * @brief update index1-th.index2-th of expression with item
     * @param ex_in input expression will be update
     * @param index1 the index
     * @param index2 the index
     * @param the new item
     */
    void let_op(lst &ex_in, int index1, int index2, const ex item) {
        ex_in.let_op(index1).let_op(index2) = item;
    }
    
    /**
     * @brief update index1-th.index2-th.index3-th of expression with item
     * @param ex_in input expression will be update
     * @param index1 the index
     * @param index2 the index
     * @param index3 the index
     * @param the new item
     */
    void let_op(ex &ex_in, int index1, int index2, int index3, const ex item) {
        ex_in.let_op(index1).let_op(index2).let_op(index3) = item;
    }
    
    /**
     * @brief update index1-th.index2-th.index3-th of expression with item
     * @param ex_in input expression will be update
     * @param index1 the index
     * @param index2 the index
     * @param index3 the index
     * @param the new item
     */
    void let_op(lst &ex_in, int index1, int index2, int index3, const ex item) {
        ex_in.let_op(index1).let_op(index2).let_op(index3) = item;
    }

    /**
     * @brief return index1-th.index2-th of expression
     * @param ex_in input expression will be update
     * @param index1 the index
     * @param index2 the index
     * @return the element
     */
    ex get_op(const ex ex_in, int index1, int index2) {
        return ex_in.op(index1).op(index2);
    }
    
    /**
     * @brief return index1-th.index2-th of expression
     * @param ex_in input expression will be update
     * @param index1 the index
     * @param index2 the index
     * @return the element
     */
    ex get_op(const lst ex_in, int index1, int index2) {
        return ex_in.op(index1).op(index2);
    }
    
    /**
     * @brief return index1-th.index2-th.index3-th of expression
     * @param ex_in input expression will be update
     * @param index1 the index
     * @param index2 the index
     * @param index3 the index
     * @return the element
     */
    ex get_op(const ex ex_in, int index1, int index2, int index3) {
        return ex_in.op(index1).op(index2).op(index3);
    }
    
    /**
     * @brief return index1-th.index2-th.index3-th of expression
     * @param ex_in input expression will be update
     * @param index1 the index
     * @param index2 the index
     * @param index3 the index
     * @return the element
     */
    ex get_op(const lst ex_in, int index1, int index2, int index3) {
        return ex_in.op(index1).op(index2).op(index3);
    }


}

