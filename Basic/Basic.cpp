/**
 * @file
 * @brief Basic Functions, extend GiNaC
 * @author F. Feng
 * @version 1.0.0
 * @date 2020-04-20
 */

#include "Basic.h"
#include "Process.h"
#include "cln/cln.h"

namespace HepLib {

    /**
     * @brief Error constructor, 
     * @param _msg error message 
     */
    Error::Error(string _msg) : msg(_msg) { }
    const char * Error::what() const throw () {
        return msg.c_str();
    }
    
    /*-----------------------------------------------------*/
    // Symbol
    /*-----------------------------------------------------*/

    DEFAULT_CTOR(Symbol)
    GINAC_BIND_UNARCHIVER(Symbol);
    IMPLEMENT_HAS(Symbol)
    IMPLEMENT_ALL(Symbol)
    GINAC_IMPLEMENT_REGISTERED_CLASS(Symbol, symbol)
    
    /**
     * @brief Symbol constructor
     * @param s symbol name
     */
    Symbol::Symbol(const string &s) : symbol(get_symbol(s)) { Table[s]=*this; }
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
    }
    
    /**
     * @brief Symbol read_archive
     * @param n archive node
     * @param sym_lst symbol lst
     */
    void Symbol::read_archive(const archive_node& n, lst& sym_lst) {
        inherited::read_archive(n, sym_lst);
        bool is_real;
        *this = Symbol(get_name());
    }
    
    ex Symbol::eval() const { return *this; }
    ex Symbol::evalf() const { return *this; }
    ex Symbol::conjugate() const { return *this; }
    ex Symbol::real_part() const { return *this; }
    ex Symbol::imag_part() const { return 0; }
    
    exmap Symbol::AssignMap;
    void Symbol::Assign(const Symbol & s, const ex & v) { AssignMap[s] = v; }
    void Symbol::Assign(const string & str, const ex & v) { AssignMap[Symbol(str)] = v; }
    void Symbol::clearAssign(const Symbol &s) { AssignMap.erase(s); }
    void Symbol::clearAssign(const string &str) { AssignMap.erase(Symbol(str)); }
    void Symbol::clearAssign() { AssignMap.clear(); }
    
    /*-----------------------------------------------------*/
    // iSymbol
    /*-----------------------------------------------------*/
    
    DEFAULT_CTOR(iSymbol)
    GINAC_BIND_UNARCHIVER(iSymbol);
    IMPLEMENT_HAS(iSymbol)
    IMPLEMENT_ALL(iSymbol)
    GINAC_IMPLEMENT_REGISTERED_CLASS(iSymbol, symbol)
    
    /**
     * @brief Symbol constructor
     * @param s symbol name
     */
    iSymbol::iSymbol(const string &s) : symbol(get_symbol(s)) { Table[s]=*this; }
    int iSymbol::compare_same_type(const basic &other) const {
        const iSymbol &o = static_cast<const iSymbol &>(other);
        int ret = get_name().compare(o.get_name());
        if(ret==0) return 0;
        else if(ret<0) return -1;
        else return 1;
    }
    
    /**
     * @brief Symbol archive
     * @param n archive node
     */
    void iSymbol::archive(archive_node & n) const {
        inherited::archive(n);
    }
    
    /**
     * @brief Symbol read_archive
     * @param n archive node
     * @param sym_lst symbol lst
     */
    void iSymbol::read_archive(const archive_node& n, lst& sym_lst) {
        inherited::read_archive(n, sym_lst);
        *this = iSymbol(get_name());
    }
    
    ex iSymbol::eval() const { return *this; }
    ex iSymbol::evalf() const { return *this; }
    ex iSymbol::conjugate() const { return (*this)*(-1); }
    ex iSymbol::real_part() const { return ex(0); }
    ex iSymbol::imag_part() const { return (-I)*(*this); }

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
        auto pgid = ppid;
        int para_max_run = nproc<0 ? omp_get_num_procs()-1 : nproc;
        if(para_max_run<1) para_max_run = 1;
        ostringstream cmd;
        cmd << "mkdir -p " << ppid;
        if(!dir_exists(to_string(ppid))) system(cmd.str().c_str());
        
        if(nbatch<=0) nbatch = ntotal/para_max_run/10;
        else if(nbatch > ntotal/para_max_run) nbatch = ntotal/para_max_run;
        if(nbatch<1) nbatch = 1;
        int btotal = ntotal/nbatch + ((ntotal%nbatch)==0 ? 0 : 1);

        for(int bi=0; bi<btotal; bi++) {
            if(Verbose > 1) {
                cout << "\r  ";
                for(int pi=0;pi<prtlvl;pi++) cout << "   ";
                cout << "\\--Evaluating ";
                if(key != "") cout << Color_HighLight << key << RESET << " ";
                cout << Color_HighLight << nbatch << "x" << RESET << "[" << (bi+1) << "/" << btotal << "] " << flush;
            }
            
            auto pid = fork();
            if(setpgid(pid, pgid)) pgid = 1; // so -pgid to -1
            if (pid < 0) {
                bi--; 
                perror("Error @ fork()");
            }
            if (pid != 0) {
                if(bi+1 >= para_max_run || pid < 0) waitpid(-pgid,NULL,0);
                continue;
            } 
            
            // pid = 0, child process
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
        while (waitpid(-pgid,NULL,0) != -1) { }
        if(Verbose > 1 && ntotal > 0) cout << "@" << now(false) << endl;

        vector<ex> ovec;
        for(int bi=0; bi<btotal; bi++) {
            if(Verbose > 1) {
                if(key == "") {
                    cout << "\r  ";
                    for(int pi=0; pi<prtlvl; pi++) cout << "   ";
                    cout << "\\--Reading *.gar [" << (bi+1) << "/" << btotal << "] " << flush;
                } else {
                    cout << "\r  ";
                    for(int pi=0;pi<prtlvl;pi++) cout << "   ";
                    cout << "\\--Reading *." << Color_HighLight << key << RESET << ".gar [" << (bi+1) << "/" << btotal << "] " << flush;
                }
            }

            int oDigits = Digits;
            Digits = 50; // a fix to float overflow
            
            archive ar;
            ostringstream garfn;
            if(key == "") garfn << ppid << "/" << bi << ".gar";
            else garfn << ppid << "/" << bi << "." << key << ".gar";
            if(!file_exists(garfn.str())) {
                cout << "File Not Found: " << garfn.str() << endl;
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
     * @brief get the possymbol from symbol factory, if not exsit, a new one will be created
     * @param s the name of the symbol
     * @param check true to check exist, and if so, error will be thrown
     * @return return the found/created symbol
     */
    const possymbol & get_possymbol(const string & s) {
        static map<string, possymbol> directory;
        map<string, possymbol>::iterator i = directory.find(s);
        if (i != directory.end()) return i->second;
        else return directory.insert(make_pair(s, possymbol(s))).first->second;
    }
    
    /**
     * @brief get the symbol from symbol factory, if not exsit, a new one will be created
     * @param s the name of the symbol
     * @param check true to check exist, and if so, error will be thrown
     * @return return the found/created symbol
     */
    const symbol & get_symbol(const string & s) {
        static map<string, symbol> directory;
        map<string, symbol>::iterator i = directory.find(s);
        if (i != directory.end()) return i->second;
        else return directory.insert(make_pair(s, symbol(s))).first->second;
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
     * @brief convert string to ex expression, using Parser internally
     * @param expr expression in string format
     * @return the parsed expression
     */
    ex str2ex(const string &expr) {
        Parser par;
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
     * @brief convert string to the lst, using Parser internally
     * @param expr a lst in string format
     * @return the parsed expression
     */
    lst str2lst(const string &expr) {
        Parser par;
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
     * @brief read file content to string
     * @param filename file name
     * @return the file content in string
     */
    string file2str(string filename) {
        ifstream ifs(filename);
        string ostr((istreambuf_iterator<char>(ifs)), (istreambuf_iterator<char>()));
        ifs.close();
        return ostr;
    }
    
    /**
     * @brief read file content to ex
     * @param filename file name
     * @return the file content in ex
     */
    ex file2ex(string filename) {
        return str2ex(file2str(filename));
    }
    
    /**
     * @brief export expression file 
     * @param the input expression
     * @param filename file name
     * @return the file content in ex
     */
    void ex2file(const ex & expr, string filename) {
        std::ofstream ofs;
        ofs.open(filename, ios::out);
        ofs << expr << endl;
        ofs.close();
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
     * @param expand true to call mma_collect before diff
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
     * @param has only expand the element e, when has(e) is true
     * @param depth will be incresed by 1 when each recursively called, when depth>5, GiNaC expand will be used
     * @return the expanded expression, expair, coeff is constant, rest involving pattern
     */
    pair<ex,epvector> mma_expand(ex const &expr_in, std::function<bool(const ex &)> has, int depth) {
        if(!has(expr_in)) {
            ex co;
            epvector epv;
            co = expr_in;
            return make_pair(co,epv);
        } if(is_a<add>(expr_in)) {
            ex co = 0;
            epvector epv;
            exmap pcmap;
            for(auto item : expr_in) {
                if(has(item)) {
                    auto co_epv = mma_expand(item, has, depth+1);
                    co += co_epv.first;
                    for(auto ep : co_epv.second) pcmap[ep.rest] += ep.coeff;
                } else co += item;
            }
            for(auto kv : pcmap) {
                if(is_zero(kv.first) || is_zero(kv.second)) continue;
                epv.push_back(expair(kv.first, kv.second));
            }
            return make_pair(co,epv);
        } else if(is_a<mul>(expr_in)) {
            ex co = 1;
            epvector epv;
            for(auto item : expr_in) {
                if(!has(item)) {
                    for(auto & ep : epv) ep.coeff *= item;
                    co *= item;
                } else {
                    exmap pcmap;
                    auto co_epv = mma_expand(item, has, depth+1);
                    for(auto ep2 : co_epv.second) {
                        pcmap[ep2.rest] += ep2.coeff * co;
                        for(auto ep1 : epv) pcmap[ep1.rest * ep2.rest] += ep1.coeff * ep2.coeff;
                    }
                    for(auto ep1 : epv) pcmap[ep1.rest] += ep1.coeff * co_epv.first;
                    co *= co_epv.first;
                    epv.clear();
                    for(auto kv : pcmap) {
                        if(is_zero(kv.first) || is_zero(kv.second)) continue;
                        epv.push_back(expair(kv.first, kv.second));
                    }
                }
            }
            return make_pair(co,epv);
        } else if(is_a<power>(expr_in) && expr_in.op(1).info(info_flags::nonnegint)) { 
            ex co = 1;
            epvector epv;
            int n = ex_to<numeric>(expr_in.op(1)).to_int();
            auto co_epv = mma_expand(expr_in.op(0), has, depth+1);
            for(int i=0; i<n; i++) {
                exmap pcmap;
                for(auto ep2 : co_epv.second) {
                    pcmap[ep2.rest] += ep2.coeff * co;
                    for(auto ep1 : epv) pcmap[ep1.rest * ep2.rest] += ep1.coeff * ep2.coeff;
                }
                for(auto ep1 : epv) pcmap[ep1.rest] += ep1.coeff * co_epv.first;
                co *= co_epv.first;
                epv.clear();
                for(auto kv : pcmap) {
                    if(is_zero(kv.first) || is_zero(kv.second)) continue;
                    epv.push_back(expair(kv.first, kv.second));
                }
            }
            return make_pair(co,epv);
        } else {
            ex co = 0;
            epvector epv;
            epv.push_back(expair(expr_in, 1));
            return make_pair(co,epv);
        }
        throw Error("mma_expand unexpected region reached.");
    }
    
    /**
     * @brief the expand like Mathematica
     * @param expr_in input expression
     * @param pats only expand the element e, when e has at least one pattern in pats is true
     * @return the expanded expression
     */
    ex mma_expand(ex const &expr_in, std::function<bool(const ex &)> has) {
        auto co_epv = mma_expand(expr_in, has, 0);
        ex ret = co_epv.first;
        for(auto ep : co_epv.second) ret += ep.coeff * ep.rest;
        return ret;
    }
    
    /**
     * @brief the expand like Mathematica
     * @param expr_in input expression
     * @param pats only expand the element e, when e has at least one pattern in pats is true
     * @param depth will be incresed by 1 when each recursively called
     * @return the expanded expression
     */
    ex mma_expand(ex const &expr_in, lst const &pats) {
        return mma_expand(expr_in, [pats](const ex & e)->bool {
            for(auto pat : pats) {
                if(e.has(pat)) return true;
            }
            return false;
        });
    }
    
    /**
     * @brief the expand like Mathematica
     * @param expr_in input expression
     * @param pat only expand the element e, when e has the pattern: pat
     * @param depth will be incresed by 1 when each recursively called
     * @return the expanded expression
     */
    ex mma_expand(ex const &expr_in, const ex &pat) {
        return mma_expand(expr_in, [pat](const ex & e)->bool {
            return e.has(pat);
        });
    }

    /**
     * @brief the collect function like Mathematica
     * @param expr_in input expression
     * @param has only collect the element e, when has(e) is true
     * @param ccf true for wrapping coefficient in coCF
     * @param cvf true to wrapping has-ed element in coVF
     * @return the collected expression
     */
    ex mma_collect(ex const &expr_in, std::function<bool(const ex &)> has, bool ccf, bool cvf) {
        auto items = mma_expand(expr_in, has);
        if(!is_a<add>(items)) items = lst{ items };
        
        ex cf = 0;
        map<ex, ex, ex_is_less> vc_map;
        for(auto item : items) {
            if(!has(item)) cf += item;
            else if(is_a<mul>(item)) {
                ex tc = 1, tv = 1;
                for(auto ii : item) {
                    if(!has(ii)) tc *= ii;
                    else tv *= ii;
                }
                vc_map[tv] += tc;
            } else {
                vc_map[item] += 1;
            }
        }
        
        ex res = 0;
        if(!is_zero(cf)) res += coCF(cf)*coVF(1);
        for(auto vc : vc_map) {
            if(has(vc.second)) {
                cerr << Color_Error << "mma_collect: pats founds @ " << vc.second << RESET << endl;
                throw Error("mma_collect: coefficent has pat.");
            }
            res += coVF(vc.first) * coCF(vc.second);
        }
        
        if(!ccf) res = res.subs(coCF(w)==w);
        if(!cvf) res = res.subs(coVF(w)==w);
        return res;
    }
    
    /**
     * @brief the collect function like Mathematica, reture the lst { {c1,v1}, {c2,v2}, ... }
     * @param expr_in input expression
     * @param has only collect the element e, when has(e) is true
     * @return the collected expression in lst
     */
    lst mma_collect_lst(ex const &expr_in, std::function<bool(const ex &)> has) {
        auto items = mma_expand(expr_in, has);
        if(!is_a<add>(items)) items = lst{ items };
        
        ex cf = 0;
        map<ex, ex, ex_is_less> vc_map;
        for(auto item : items) {
            if(!has(item)) cf += item;
            else if(is_a<mul>(item)) {
                ex tc = 1, tv = 1;
                for(auto ii : item) {
                    if(!has(ii)) tc *= ii;
                    else tv *= ii;
                }
                vc_map[tv] += tc;
            } else {
                vc_map[item] += 1;
            }
        }
        
        lst res_lst;
        if(!is_zero(cf)) res_lst.append(lst{cf,1});
        for(auto vc : vc_map) {
            if(has(vc.second)) {
                cerr << Color_Error << "mma_collect: pats founds @ " << vc.second << RESET << endl;
                throw Error("mma_collect: coefficent has pat.");
            }
            res_lst.append(lst{vc.second, vc.first});
        }
        
        return res_lst;
    }
    
    /**
     * @brief the collect function like Mathematica
     * @param expr_in input expression
     * @param pats only collect the element e match at least one patter in pats, like mma_expand
     * @param ccf true for wrapping coefficient in coCF
     * @param cvf true to wrapping has-ed element in coVF
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
     * @param pats only collect the element e match at least one patter in pats, like mma_expand
     * @return the collected expression
     */
    lst mma_collect_lst(ex const &expr_in, lst const &pats) {
        return mma_collect_lst(expr_in, [pats](const ex & e)->bool {
            for(auto pat : pats) {
                if(e.has(pat)) return true;
            }
            return false;
        });
    }
    
    /**
     * @brief the collect function like Mathematica
     * @param expr_in input expression
     * @param pat only collect the element e match the pattern, like mma_expand
     * @param ccf true for wrapping coefficient in coCF
     * @param cvf true to wrapping has-ed element in coVF
     * @return the collected expression
     */
    ex mma_collect(ex const &expr_in, ex const &pat, bool ccf, bool cvf) {
        return mma_collect(expr_in, [pat](const ex & e)->bool {
            return e.has(pat);
        }, ccf, cvf);
    }
    
    /**
     * @brief the collect function like Mathematica
     * @param expr_in input expression
     * @param pat only collect the element e match the pattern, like mma_expand
     * @return the collected expression
     */
    lst mma_collect_lst(ex const &expr_in, ex const &pat) {
        return mma_collect_lst(expr_in, [pat](const ex & e)->bool {
            return e.has(pat);
        });
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
        auto oDigits = Digits;
        Digits = 50;
        for(auto zi : zs) {
            repl.append(zi==zi.evalf());
        }
        Digits = oDigits;
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

    /**
     * @brief count the leaf nodes
     * @param e input expression
     * @return the total number of leaf nodes
     */
    long long int LeafCount(const ex & e) {
        if(is_a<numeric>(e)) return 1;
        else if(is_a<symbol>(e)) return 2;
        else if(e.nops()==1) return 3;
        long long c = 10;
        for(auto item : e) c += LeafCount(item);
        return c;
    }
    
    bool ex_less(const ex &a, const ex &b) {
        static ex_is_less comp;
        if(is_a<numeric>(a) && is_a<numeric>(b)) {
            if(a!=b) return a<b;
            else return false;
        } if(is_a<lst>(a) && is_a<lst>(b)) {
            if(a.nops()!=b.nops()) return a.nops()<b.nops();
            auto nn = a.nops();
            for(int i=0; i<nn; i++) {
                if(is_zero(a.op(i)-b.op(i))) continue;
                return ex_less(a.op(i), b.op(i));
            }
            return false;
        } else {
            auto la = LeafCount(a);
            auto lb = LeafCount(b);
            if(la!=lb) la<lb;
            return comp(a,b);
        }
    }
     
     /**
      * @brief sort the list in less order, or the reverse
      * @param ilst input lst, will be updated after call
      * @param less true for less order
      */
     void sort_lst(lst & ilst, bool less) {
        auto ivec = lst2exvec(ilst);
        std::sort(ivec.begin(), ivec.end(), ex_less);
        auto n = ivec.size();
        if(less) for(auto i=0; i<n; i++) ilst.let_op(i) = ivec[i];
        else for(auto i=0; i<n; i++) ilst.let_op(i) = ivec[n-1-i];
     }
     
     /**
      * @brief sort the list in less order, or the reverse
      * @param ilst input lst, will be updated after call
      * @param ki the sort key is at .op(n)
      * @param less true for less order
      */
     void sort_lst_by(lst & ilst, int ki, bool less) {
        auto ivec = lst2exvec(ilst);
        std::sort(ivec.begin(), ivec.end(), [ki](const auto &as, const auto &bs){
            auto a = as.op(ki);
            auto b = bs.op(ki);
            return ex_less(a,b);
        });
        auto n = ivec.size();
        if(less) for(auto i=0; i<n; i++) ilst.let_op(i) = ivec[i];
        else for(auto i=0; i<n; i++) ilst.let_op(i) = ivec[n-1-i];
     }
     
     /**
      * @brief sort the list in less order, or the reverse
      * @param ivec input exvector, will be updated after call
      * @param less true for less order
      */
     void sort_vec(exvector & ivec, bool less) {
        std::sort(ivec.begin(), ivec.end(), [less](const auto &a0, const auto &b0){
            auto a = a0;
            auto b = b0;
            if(!less) {
                auto tmp = a;
                a = b;
                b = tmp;
            }
            return ex_less(a,b);
        });
     }
     
     /**
      * @brief sort the list in less order, or the reverse
      * @param ivec input exvector, will be updated after call
      * @param ki the sort key is at .op(n)
      * @param less true for less order
      */
     void sort_vec_by(exvector & ivec, int ki, bool less) {
        std::sort(ivec.begin(), ivec.end(), [ki,less](const auto &as, const auto &bs){
            auto a = as.op(ki);
            auto b = bs.op(ki);
            if(!less) {
                auto tmp = a;
                a = b;
                b = tmp;
            }
            return ex_less(a,b);
        });
     }
     
     
    //-----------------------------------------------------------
    // XIntegral Class
    //-----------------------------------------------------------
    GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(XIntegral, basic,
        print_func<print_dflt>(&XIntegral::print)
    )
    
    DEFAULT_CTOR(XIntegral)
    GINAC_BIND_UNARCHIVER(XIntegral);
    IMPLEMENT_HAS(XIntegral)
    IMPLEMENT_ALL(XIntegral)

    int XIntegral::compare_same_type(const basic &other) const {
        const XIntegral &o = static_cast<const XIntegral &>(other);
        auto c = Functions.compare(o.Functions);
        if(c!=0) return c;
        c = Exponents.compare(o.Exponents);
        if(c!=0) return c;
        c = Deltas.compare(o.Deltas);
        if(c!=0) return c;
        return 0;
    }
    
    void XIntegral::print(const print_dflt &c, unsigned level) const {
        c.s << "[∬" << "(" << Functions << ")^(" << Exponents << ")";
        if(Deltas.nops()>0) c.s << "*𝜹(" << Deltas << ")";
        c.s << "]";
    }
    
    size_t XIntegral::nops() const { return 3; }
    ex XIntegral::op(size_t i) const {
        if(i==0) return Functions;
        else if(i==1) return Exponents;
        else if(i==2) return Deltas;
        else throw Error("XIntegral::op, It is required that i<3.");
    }
    ex & XIntegral::let_op(size_t i) {
        ensure_if_modifiable();
        if(i==0) return Functions;
        else if(i==1) return Exponents;
        else if(i==2) return Deltas;
        else throw Error("XIntegral::let_op, It is required that i<3.");
    }
    
    void XIntegral::archive(archive_node & n) const {
        inherited::archive(n);
        n.add_ex("Functions", Functions);
        n.add_ex("Exponents", Exponents);
        n.add_ex("Deltas", Deltas);
    }
    
    void XIntegral::read_archive(const archive_node& n, lst& sym_lst) {
        inherited::read_archive(n, sym_lst);
        n.find_ex("Functions", Functions, sym_lst);
        n.find_ex("Exponents", Exponents, sym_lst);
        n.find_ex("Deltas", Deltas, sym_lst);
    }
        
    XIntegral::XIntegral(ex fed)  { 
        Functions = fed.op(0);
        Exponents = fed.op(1);
        if(fed.nops()>2) Deltas = fed.op(2);
    }
    
    XIntegral::XIntegral(ex loops, ex ps, ex ns) {
        ex a = 0;
        ex pre = 1;
        ex rem = 0;
        lst delta;
        int xni = 0;
        lst funs, exps;
        for(int i=0; i<ps.nops(); i++) {
            if(is_zero(ns.op(i))) continue;
            auto cpi = ps.op(i).expand().subs(Symbol::AssignMap);
            bool ltQ = false; 
            for(auto li : loops) {
                if(cpi.has(li)) {
                    ltQ = true;
                    break;
                }
            }
            
            ex sgn = 0;
            if(!ltQ) {
                pre *= pow(ps.op(i).expand().subs(Symbol::AssignMap), ex(0)-ns.op(i));
                ns.let_op(i) = 0;
                ps.let_op(i) = 1;
                continue;
            } else if(is_a<numeric>(ns.op(i)) && (ns.op(i)<0)) {
                throw Error("XIntegral, negative powers not supported yet.");
            }

            // check loop^2
            for(auto li : loops) {
                if(!is_a<numeric>(cpi.coeff(li,2))) continue; // TODO: make sure positive
                numeric nm = ex_to<numeric>(cpi.coeff(li,2));
                if(nm.is_zero()) continue;
                sgn = nm>0 ? -1 : 1;
                break;
            }
            // check iEpsilon
            if(sgn.is_zero()) {
                if(!is_a<numeric>(cpi.coeff(iEpsilon))) continue;
                numeric nm = ex_to<numeric>(cpi.coeff(iEpsilon));
                if(!nm.is_zero()) sgn = nm>0 ? -1 : 1;
            }
            // others
            if(sgn.is_zero()) {
                sgn = 1;
                if(is_a<numeric>(cpi) && ex_to<numeric>(cpi)>0) sgn = -1;
                else sgn = 1;
            }
            
            cpi = (ps.op(i)*sgn).subs(iEpsilon==0);
            if(sgn==-1) pre *= exp(I * Pi * ns.op(i));
            a += ns.op(i);
            auto xi = x(xni);
            rem += xi * cpi;
            delta.append(xi);
            if(!is_zero(ns.op(i)-1)) {
                funs.append(xi);
                exps.append(ns.op(i)-1);
            }
            pre /= tgamma(ns.op(i));
            xni++;
        }
        rem = rem.expand().subs(Symbol::AssignMap);
        ex dl2 = loops.nops()*(4-2*ep)/2;
        pre *= tgamma(a-dl2) * pow(I,loops.nops()) * pow(Pi,dl2) * pow(2*Pi,loops.nops()*(2*ep-4));
        
        // Loop
        ex u=1;
        for(int i=0; i<loops.nops(); i++) {
            ex t2 = rem.coeff(loops.op(i),2);
            ex t1 = rem.coeff(loops.op(i),1);
            ex t0 = rem.coeff(loops.op(i),0);
            u *= (-t2);
            if(t2==0) {
                Functions = lst{0};
                Exponents = lst{1};
                return;
            }
            rem = (t0 - pow(t1,2)/(4*t2)).expand();
        }
        rem = (rem.subs(Symbol::AssignMap)).normal();
        u = (u.subs(Symbol::AssignMap)).normal();
        rem = (rem * u).normal();
        
        funs.prepend(rem);
        exps.prepend(dl2-a);
        funs.prepend(u);
        exps.prepend(a-dl2-(4-2*ep)/2);
        auto num_den = numer_denom(pre);
        funs.append(num_den.op(0));
        exps.append(1);
        if(!is_zero(num_den.op(1))) {
            funs.append(num_den.op(1));
            exps.append(-1);
        }
        
        Functions = funs;
        Exponents = exps;
        Deltas = lst{delta};
    }
    
    int CpuCores() {
        return omp_get_num_procs();
    }
    
    
    //-----------------------------------------------------------
    // fermat_numer_denom / fermat_normal
    //-----------------------------------------------------------
    ex fermat_numer_denom(const ex & expr) {
        static map<pid_t, Fermat> fermat_map;
        static int v_max = 0;
        static Symbol vi("vi");
        static MapFunction iMap([](const ex & e, MapFunction &self)->ex {
            if(!e.has(I)) return e;
            else if(e==I) return vi;
            else if(e.info(info_flags::numeric) && !e.info(info_flags::real)) return real_part(e)+imag_part(e)*vi;
            else return e.map(self);
        });

        auto pid = getpid();
        if((fermat_map.find(pid)==fermat_map.end())) { // init section
            fermat_map[pid].Init();
            v_max = 0;
        }
        Fermat &fermat = fermat_map[pid];
        
        auto expr_in = expr;
        exmap fun2s, s2fun;
        expr_in = MapFunction([&fun2s,&s2fun](const ex & e, MapFunction &self)->ex {
            if(is_a<GiNaC::function>(e) || e.match(sqrt(w))) {
                if(is_zero(fun2s[e])) {
                    symbol s;
                    fun2s[e] = s;
                    s2fun[s] = e;
                }
                return fun2s[e];
            } else if(!has_function(e) && !e.has(sqrt(w))) return e;
            else return e.map(self);
        })(expr_in);
        
        lst rep_vs;
        for(const_preorder_iterator i = expr_in.preorder_begin(); i != expr_in.preorder_end(); ++i) {
            auto e = (*i);
            if(is_a<symbol>(e)) rep_vs.append(e);
        }
        rep_vs.sort();
        rep_vs.unique();        
        sort_lst(rep_vs);
        
        exmap v2f, f2v;
        exmap nn_map;
        auto nn_pi = cln::nextprobprime(3);
        int fvi = 0;
        for(auto vi : rep_vs) {
            auto name = "v" + to_string(fvi);
            Symbol s(name);
            v2f[vi] = s;
            f2v[s] = vi;
            fvi++;
            nn_pi = cln::nextprobprime(nn_pi+1);
            nn_map[s] = ex(1)/numeric(nn_pi);
        }
        
        stringstream ss;
        if(fvi>111) {
            cout << rep_vs << endl;
            throw Error("Fermat: Too many variables.");
        }
        if(fvi>v_max) {
            if(v_max==0) ss << "&(J=vi);" << endl;
            for(int i=v_max; i<fvi; i++) ss << "&(J=v" << i << ");" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
            v_max = fvi;
        }
        
        ex nn_chk(1), num(1), den(1);
        if(!is_a<mul>(expr_in)) expr_in = lst{expr_in};
        for(auto item : expr_in) {
            if(!is_a<add>(item)) item = lst{item};
            if(fermat_use_array) ss << "Array m[" << item.nops() << "];" << endl;
            else ss << "res:=0;" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
            
            ex nn_chk2=0;
            for(int i=0; i<item.nops(); i++) {
                ex tt = item.op(i).subs(v2f);
                nn_chk2 += tt.subs(nn_map);
                tt = iMap(tt);
                if(fermat_use_array) ss << "m[" << (i+1) << "]:=";
                else ss << "item:=";
                ss << tt << ";" << endl;
                if(!fermat_use_array) ss << "res:=res+item;" << endl;
                fermat.Execute(ss.str());
                ss.clear();
                ss.str("");
            }
            if(fermat_use_array) {
                ss << "res:=Sumup([m]);" << endl;
                fermat.Execute(ss.str());
                ss.clear();
                ss.str("");
            }
            nn_chk *= nn_chk2;
            
            static string bstr("[-begin-]"), estr("[-end-]");
            ss << "&(U=1);" << endl; // ugly printing, the whitespace matters
            ss << "!('" <<bstr<< "','{',Numer(res),',',Denom(res),'}','" <<estr<< "')" << endl;
            auto ostr = fermat.Execute(ss.str());
            ss.clear();
            ss.str("");

            // make sure last char is 0
            if(ostr[ostr.length()-1]!='0') throw Error("fermat_together: last char is NOT 0.");
            ostr = ostr.substr(0, ostr.length()-1);
            auto cpos = ostr.find(bstr);
            if(cpos==string::npos) throw Error(bstr+" NOT Found.");
            ostr = ostr.substr(cpos+bstr.length(),ostr.length()-cpos);
            cpos = ostr.find(estr);
            if(cpos==string::npos) throw Error(estr+" NOT Found.");
            ostr = ostr.substr(0,cpos);
            string_trim(ostr);       

            symtab st;
            st["vi"] = I; 
            Parser fp(st);
            auto ret = fp.Read(ostr);
            num *= ret.op(0);
            den *= factor(ret.op(1));
            
            ss << "&(U=0);" << endl; // disable ugly printing
            if(fermat_use_array) ss << "@(res,[m]);" << endl;
            else ss << "@(res,item);" << endl;
            ss << "&_G;" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
        }
        //fermat.Exit();
                
        auto nn_ret = subs(num/den,nn_map);
        if(nn_chk-nn_ret!=0) {
            cout << nn_chk << " : " << nn_ret << endl;
            throw Error("fermat_together: N Check Failed.");
        }
        
        num = num.subs(f2v).subs(s2fun);
        den = den.subs(f2v).subs(s2fun);
        return lst{num, den};
    }
    
    ex fermat_normal(const ex & expr) {
        auto nd = fermat_numer_denom(expr);
        return nd.op(0)/nd.op(1);
    }
    
}

