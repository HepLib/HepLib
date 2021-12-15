/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "BASIC.h"
#include "cln/cln.h"

#ifdef _USE_FLOAT128
extern "C" {
    #include <quadmath.h>
}
#endif

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
    unsigned Symbol::calchash() const {
        static unsigned _hash = symbol("_").gethash();
        return _hash;
    }
    
    void Symbol::set(const ex & v) const { vmap[*this] = v; }
    void Symbol::unset() const { vmap.erase(*this); }
    
    void Symbol::set(const Symbol & s, const ex & v) { vmap[s] = v; }
    void Symbol::set(const string & str, const ex & v) { vmap[Symbol(str)] = v; }
    void Symbol::unset(const Symbol &s) { vmap.erase(s); }
    void Symbol::unset(const string &str) { vmap.erase(Symbol(str)); }
    void Symbol::unset_all() { vmap.clear(); }
    ex Symbol::set_all(const ex & expr) {
        ex v1 = expr;
        while(true) {
            ex v2 = v1.subs(vmap);
            if(v2.is_equal(v1)) break;
            v1 = v2;
        }
        return v1;
    }
    
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
    
    const char* ErrColor = RED;
    const char* WarnColor = MAGENTA;
    const char* Color_HighLight = WHITE;
    
    static int GiNaC_Parallel_Level = 0; // for the internal usage only
    /**
     * @brief GiNaC Parallel Evaluation using fork
     * @param ntotal the number of total items, 0 for non-parallel version
     * @param nbatch the batch number for each run, if nbatch<=0, then nbatch = ntotal/para_max_run/10
     * @param f function to be applied on the index, from 0 to (ntotal-1)
     * @param key key used in archive file name and display message
     * @param rm default true, and false will keep the archive file
     * @param pre the pre-string in the print message
     * @return return the ntotal-element vector, i.e., [f(0), ..., f(ntotal-1)]
     */
    vector<ex> GiNaC_Parallel(
        int ntotal, int nbatch, 
        std::function<ex(int)> f,
        const string & key,
        bool rm,
        const string &pre) {
        int nproc = GiNaC_Parallel_Process;
        
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
        if(!dir_exists(to_string(ppid))) system(cmd.str().c_str());
        
        if(nbatch<=0) nbatch = ntotal/para_max_run/5;
        else if(nbatch > ntotal/para_max_run) nbatch = ntotal/para_max_run;
        if(nbatch>GiNaC_Parallel_BatchMax) nbatch = GiNaC_Parallel_BatchMax;
        if(nbatch<1) nbatch = 1;
        int btotal = ntotal/nbatch + ((ntotal%nbatch)==0 ? 0 : 1);
        
        bool nst = (GiNaC_Parallel_Level>0);
        pid_t npid=0, pgid;
        if(getpgid(0)!=ppid) { // not group leader
            if(nst) {
                npid = fork();
                if(npid < 0) throw Error("GiNaC_Parallel: Error (1) @ fork()");
                if(npid!=0) goto wait_label;
            } // else - for case in a shell script
            if(setpgid(0,0)) {
                if(setpgid(0,0)) throw Error("GiNaC_Parallel: setpgid(0,0) Failed.");
            }
        }
        
        pgid = getpid(); // should be groud leader @ here
        for(int bi=0; bi<btotal; bi++) {
            if(Verbose > 1 && !nst) {
                cout << "\r                                                   \r" << pre;
                cout << "\\--Evaluating ";
                if(key != "") cout << Color_HighLight << key << RESET << " ";
                cout << Color_HighLight << nbatch << "x" << RESET << "[" << (bi+1) << "/" << btotal << "] @ " << now(false) << flush;
            }
            
            auto pid = fork();
            if(pid < 0) throw Error("GiNaC_Parallel: Error (2) @ fork()");
            if(pid==0 && setpgid(0, pgid)!=0) {
                if(setpgid(0, pgid)) throw Error("GiNaC_Parallel: setpgid(0, pgid) Failed.");
            }
            if(pid != 0) {
                if(bi+1 >= para_max_run || pid < 0) waitpid(-pgid,NULL,0);
                continue; // parent process goes next cycle
            }
            
            // pid = 0, child process
            GiNaC_Parallel_Level++;
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
            } catch(exception &p) { // NOT use Error
                cout << ErrColor << "Failed in GiNaC_Parallel!" << RESET << endl;
                cout << ErrColor << p.what() << RESET << endl;
            }
            exit(0);
        }
        
        if(getpid()!=pgid) exit(0); // make sure child exit
        while (waitpid(-pgid,NULL,0) != -1) { }
        if(ntotal > 0 && !nst) {
            if(Verbose > 10) cout << endl;
            else if(Verbose > 1) cout << "\r                                                   \r" << flush;
        }
        if(getpid()!=ppid) exit(0); // make the forked group leader exit
        
        wait_label:
        if(npid!=0) waitpid(npid,NULL,0); // wait nest process to exit
        
        vector<ex> ovec;
        for(int bi=0; bi<btotal; bi++) {
            if(Verbose > 1 && !nst) {
                if(key == "") {
                    cout << "\r                                                   \r" << pre;
                    cout << "\\--Reading *.gar [" << (bi+1) << "/" << btotal << "] @ " << now(false) << flush;
                } else {
                    cout << "\r                                                   \r" << pre;
                    cout << "\\--Reading *." << Color_HighLight << key << RESET << ".gar [" << (bi+1) << "/" << btotal << "] @ " << now(false) << flush;
                }
            }

            archive ar;
            ostringstream garfn;
            if(key == "") garfn << ppid << "/" << bi << ".gar";
            else garfn << ppid << "/" << bi << "." << key << ".gar";
            if(!file_exists(garfn.str())) throw Error("File Not Found: " + garfn.str());
            auto res_lst = garRead(garfn.str());
            remove(garfn.str().c_str());
            for(auto res : res_lst) ovec.push_back(res);
        }
        
        if(rm) {
            cmd.clear();
            cmd.str("");
            cmd << "rm -fr " << ppid;
            system(cmd.str().c_str());
        }
        if(ntotal > 0 && !nst) {
            if(Verbose > 10) cout << endl;
            else if(Verbose > 1) cout << "\r                                                   \r" << flush;
        }
        return ovec;
    }
    
    /**
     * @brief get the symbol from symbol factory, if not exsit, a new one will be created
     * @param s the name of the symbol
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
        auto oDigits = Digits;
        Digits = NNDigits; // a fix to float overflow
        archive ar;
        ifstream in(garfn);
        in >> ar;
        in.close();
        for(int i=0; i<ar.num_expressions(); i++) {
            string name;
            ex res = ar.unarchive_ex(GiNaC_archive_Symbols, name, i);
            resMap[name] = res;
        }
        Digits = oDigits;
    }
    
    /**
     * @brief garRead from file, only the element w.r.t. key 
     * @param garfn ginac archive filename
     * @param key the archive key, note (const char *) type, match unarchive_ex in GiNaC
     * @return the expression w.r.t. key
     */
    ex garRead(const string &garfn, const char* key) { // use the const char *, not string
        auto oDigits = Digits;
        Digits = NNDigits; // a fix to float overflow
        archive ar;
        ifstream in(garfn);
        in >> ar;
        in.close();
        auto res = ar.unarchive_ex(GiNaC_archive_Symbols, key);
        Digits = oDigits;
        return res;
    }

    /**
     * @brief garRead from file, only the element w.r.t. key "res", note inner check will be performed
     * @param garfn ginac archive filename
     * @return the expression w.r.t. key "res"
     */
    ex garRead(const string &garfn) {
        auto oDigits = Digits;
        Digits = NNDigits; // a fix to float overflow
        archive ar;
        ifstream in(garfn);
        in >> ar;
        in.close();
        auto c = ar.unarchive_ex(GiNaC_archive_Symbols, "c");
        auto res = ar.unarchive_ex(GiNaC_archive_Symbols, "res");
        if(c!=19790923) throw Error("garRead: check faild for file: " + garfn);
        Digits = oDigits;
        return res;
    }
    
    /**
     * @brief garWrite to write the string-key map to the archive
     * @param garfn ginac archive filename for output
     * @param resMap a key-valued map
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
        oss << expr;
        return oss.str();
    }
    
    /**
     * @brief convert exvec to output string, the defalut printer format will be used
     * @param expr a exvec expression
     * @return the output string
     */
    string ex2str(const exvector &expr) {
        ostringstream oss;
        oss << expr;
        return oss.str();
    }
    
    /**
     * @brief convert exmap to output string, the defalut printer format will be used
     * @param expr a exvec expression
     * @return the output string
     */
    string ex2str(const exmap &expr) {
        ostringstream oss;
        oss << expr;
        return oss.str();
    }
    
    /**
     * @brief convert exset to output string, the defalut printer format will be used
     * @param expr a exvec expression
     * @return the output string
     */
    string ex2str(const exset &expr) {
        ostringstream oss;
        oss << expr;
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
     * @brief read file content to string vector
     * @param filename file name
     * @param skip_empty skip blank lines
     * @return the file content in string vector, lines in vector
     */
    vector<string> file2strvec(string filename, bool skip_empty) {
        ifstream fs(filename);
        vector<string> ovec;
        std::string line;
        while (std::getline(fs, line)) {
            if(!skip_empty || line.size()>0) ovec.push_back(line);
        }
        fs.close();
        return ovec;
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
     * @brief read file content to ex
     * @param filename file name
     * @param st symtab
     * @return the file content in ex
     */
    ex file2ex(string filename, symtab st) {
        return str2ex(file2str(filename), st);
    }
    
    /**
     * @brief export expression file 
     * @param expr the input expression
     * @param filename file name
     * @return the file content in ex
     */
    void ex2file(const ex & expr, string filename) {
        std::ofstream ofs;
        ofs.open(filename, ios::out);
        ofs << expr << endl;
        ofs.close();
    }
    
    void ex2file(string filename, const ex & expr) {
        ex2file(expr, filename);
    }
    
#ifdef _USE_FLOAT128
    /**
     * @brief __float128 to ex
     * @param num a __float128 number
     * @return a ex object for the number
     */
    ex q2ex(__float128 num) {
        char buffer[128];
        quadmath_snprintf(buffer, sizeof buffer, "%.36QG", num);
        numeric ret(buffer);
        return ret;
    }

    /**
     * @brief ex of numeric to __float128
     * @param num a ex number
     * @return a __float128 object for the number
     */
    __float128 ex2q(ex num) {
        ostringstream nss;
        auto oDigits = Digits;
        Digits = 40;
        nss << num.evalf() << endl;
        __float128 ret = strtoflt128(nss.str().c_str(), NULL);
        Digits = oDigits;
        return ret;
    }
#else
    /**
     * @brief long double to ex
     * @param num a long double number
     * @return a ex object for the number
     */
    ex q2ex(long double num) {
        char buffer[128];
        sprintf(buffer, "%.20LG", num);
        numeric ret(buffer);
        return ret;
    }

    /**
     * @brief ex of numeric to long double
     * @param num a ex number
     * @return a __float128 object for the number
     */
    long double ex2q(ex num) {
        ostringstream nss;
        auto oDigits = Digits;
        Digits = 40;
        nss << num.evalf() << endl;
        long double ret = strtold(nss.str().c_str(), NULL);
        Digits = oDigits;
        return ret;
    }
#endif
    
    /**
     * @brief ex to integer
     * @param num an integer in ex
     * @return the integer for the number
     */
    int ex2int(ex num) {
        return ex_to<numeric>(num).to_int();
    }
    
    /**
     * @brief convert exvector to lst
     * @param exvec input exvector
     * @return lst
     */
    lst vec2lst(const exvector & exvec) {
        lst ret;
        for(auto item : exvec) ret.append(item);
        return ret;
    }
    
    /**
     * @brief convert lst to exvector
     * @param alst input lst
     * @return exvector
     */
    exvector lst2vec(const lst & alst) {
        exvector ret;
        for(auto item : alst) ret.push_back(item);
        return ret;
    }
    
    /**
     * @brief convert add to lst
     * @param expr input expression
     * @return lst
     */
    lst add2lst(const ex & expr) {
        if(!is_a<add>(expr)) return lst{expr};
        lst ret;
        for(auto item : expr) ret.append(item);
        return ret;
    }
    
    /**
     * @brief convert mul to lst
     * @param expr input expression
     * @return lst
     */
    lst mul2lst(const ex & expr) {
        if(!is_a<mul>(expr)) return lst{expr};
        lst ret;
        for(auto item : expr) ret.append(item);
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
    ex series_ex(ex const & expr_in, symbol const &s0, int sn0) {
        ex expr = expr_in;
        if(!expr.has(s0)) return expr;
        
        exset sset;
        expr.find(pow(s0, w), sset);
        numeric sn_lcm = 1;
        exmap repl1st, repl2nd;
        for(auto pi : sset) {
            auto sn = expand(pi.op(1));
            if(!is_a<numeric>(sn)) {
                ex ex1 = 0, ex2 = 0;
                for(auto item : add2lst(sn)) {
                    if(!is_a<numeric>(item)) ex1 += item;
                    else ex2 += item;
                }
                sn = ex2;
                static symbol sss;
                repl1st[pi] = sss * pow(s0, sn);
                repl2nd[sss] = pow(s0, ex1);
            }
            if(!(is_a<numeric>(sn) && ex_to<numeric>(sn).is_rational())) {
                cerr << "s = " << s0 << endl;
                cerr << "expr_in = " << expr_in << endl;
                throw Error("series_ex: Not rational sn = " + ex2str(sn));
            }
            sn_lcm = lcm(sn_lcm, ex_to<numeric>(sn).denom());
        }
        if(expr.has(sqrt(s0))) sn_lcm = lcm(sn_lcm, numeric(2));
        expr = expr.subs(repl1st);
        
        static symbol s("s");
        if(!sn_lcm.is_integer()) throw Error("series_ex: Not integer with " + ex2str(sn_lcm));
        if(sn_lcm<0) sn_lcm = numeric(0)-sn_lcm;
        int sn = sn0 * sn_lcm.to_int();
        expr = expr.subs(pow(s0,w)==pow(s,w*sn_lcm)).subs(sqrt(s0)==pow(s,sn_lcm/2)).subs(s0==pow(s,sn_lcm));
        
        ex cvs = collect_lst(expr,s);
        ex ret = 0;
        for(auto cv : cvs) {
            bool ok = false; 
            int exN = 1;
            while(exN<10) {
                expr = cv.op(1) + pow(s,sn+exN+2)+pow(s,sn+exN+3);
                expr = expr.series(s,sn+exN);
                ex ot = 0;
                for(int i=0; i<expr.nops(); i++) {
                    if(is_order_function(expr.op(i))) {
                        ot = expr.op(i);
                        break;
                    }
                }
                if(!is_order_function(ot)) {
                    cerr << "expr = " << expr << endl;
                    throw Error("series_ex: Not an Order term with " + ex2str(ot));
                }
                if(ot.op(0).degree(s)>sn) {
                    expr = series_to_poly(expr);
                    expr = collect_ex(expr,s);
                    for(int i=expr.ldegree(s); (i<=expr.degree(s) && i<=sn); i++) {
                        ret += cv.op(0) * expr.coeff(s,i) * pow(s0,ex(i)/sn_lcm);
                    }
                    ok = true;
                    break;
                }
                exN++;
            }
            if(!ok) throw Error("series_ex seems not working!");
        }
        ret = ret.subs(s==pow(s0,ex(1)/sn_lcm)); // need this for log-terms
        ret = ret.subs(repl2nd);
        ret = collect_ex(ret,s0);
        return ret;
    }

    /**
     * @brief the differential like Mathematica
     * @param expr input expression
     * @param xp the variable, can be an expression, will replace by a symbol and back again
     * @param nth nth-derivative
     * @param expand true to call collect_ex before diff
     * @return the corresponding nth-derivative
     */
    ex diff_ex(ex const expr, ex const xp, unsigned nth, bool expand) {
        symbol s;
        ex res = expr.subs(xp==s);
        if(!expand) return res.diff(s,nth).subs(s==xp);
        auto cvs = collect_lst(res, s);
        res = 0;
        for(auto cv : cvs) {
            res += cv.op(0) * cv.op(1).diff(s,nth).subs(s==xp);
        }
        return res;
    }
    
    ex expand_using_map(ex const &expr_in, std::function<bool(const ex &)> has_func) {
        auto map = MapFunction([has_func](const ex & e, MapFunction & self)->ex {
            if(is_a<add>(e)) 
                return ex_to<basic>(e.map(self)).setflag(status_flags::expanded);
            else if(is_a<mul>(e)) 
                return ex_to<basic>(expand(e.map(self))).setflag(status_flags::expanded);
            else if(is_a<power>(e) && e.op(1).info(info_flags::nonnegint)) 
                return ex_to<basic>(expand(e.map(self))).setflag(status_flags::expanded);
            else 
                return ex_to<basic>(e).setflag(status_flags::expanded);
        });
        return map(expr_in);
    }
    
    typedef vector<pair<ex,ex>> epvec_t; // first-rest, second-coeff
    pair<ex,epvec_t> inner_expand_collect(ex const &expr_in, std::function<bool(const ex &)> has_func, int depth=0) {
        if(!has_func(expr_in)) {
            ex co;
            epvec_t epv;
            co = expr_in;
            return make_pair(co,epv);
        } if(is_a<add>(expr_in)) {
            ex co = 0;
            epvec_t epv;
            exmap pcmap;
            for(auto item : expr_in) {
                if(has_func(item)) {
                    auto co_epv = inner_expand_collect(item, has_func, depth+1);
                    co += co_epv.first;
                    for(auto ep : co_epv.second) pcmap[ep.first] += ep.second;
                } else co += item;
            }
            for(auto kv : pcmap) {
                if(is_zero(kv.first) || is_zero(kv.second)) continue;
                epv.push_back(make_pair(kv.first, kv.second));
            }
            return make_pair(co,epv);
        } else if(is_a<mul>(expr_in)) {
            ex co = 1;
            epvec_t epv;
            for(auto item : expr_in) {
                if(!has_func(item)) {
                    for(auto & ep : epv) ep.second *= item;
                    co *= item;
                } else {
                    exmap pcmap;
                    auto co_epv = inner_expand_collect(item, has_func, depth+1);
                    for(auto ep2 : co_epv.second) {
                        pcmap[ep2.first] += ep2.second * co;
                        for(auto ep1 : epv) pcmap[ep1.first * ep2.first] += ep1.second * ep2.second;
                    }
                    for(auto ep1 : epv) pcmap[ep1.first] += ep1.second * co_epv.first;
                    co *= co_epv.first;
                    epv.clear();
                    for(auto kv : pcmap) {
                        if(is_zero(kv.first) || is_zero(kv.second)) continue;
                        epv.push_back(make_pair(kv.first, kv.second));
                    }
                }
            }
            return make_pair(co,epv);
        } else if(is_a<power>(expr_in) && expr_in.op(1).info(info_flags::nonnegint)) {
            ex co = 1;
            epvec_t epv;
            int n = ex_to<numeric>(expr_in.op(1)).to_int();
            auto co_epv = inner_expand_collect(expr_in.op(0), has_func, depth+1);
            for(int i=0; i<n; i++) {
                exmap pcmap;
                for(auto ep2 : co_epv.second) {
                    pcmap[ep2.first] += ep2.second * co;
                    for(auto ep1 : epv) pcmap[ep1.first * ep2.first] += ep1.second * ep2.second;
                }
                for(auto ep1 : epv) pcmap[ep1.first] += ep1.second * co_epv.first;
                co *= co_epv.first;
                epv.clear();
                for(auto kv : pcmap) {
                    if(is_zero(kv.first) || is_zero(kv.second)) continue;
                    epv.push_back(make_pair(kv.first, kv.second));
                }
            }
            return make_pair(co,epv);
        } else {
            ex co = 0;
            epvec_t epv;
            epv.push_back(make_pair(expr_in, 1));
            return make_pair(co,epv);
        }
        throw Error("inner_expand_collect unexpected region reached.");
    }
    
    /**
     * @brief the expand like Mathematica
     * @param expr_in input expression
     * @param has_func only expand the element e, when has_func(e) is true
     * @return the expanded expression
     */
    ex expand_ex(ex const &expr_in, std::function<bool(const ex &)> has_func) {
        auto co_epv = inner_expand_collect(expr_in, has_func, 0);
        ex ret = co_epv.first;
        for(auto ep : co_epv.second) ret += ep.second * ep.first;
        return ret;
    }

    /**
     * @brief the collect function like Mathematica
     * @param expr_in input expression
     * @param has_func only collect the element e, when has_func(e) is true
     * @param cf true for wrapping coefficient in coCF
     * @param vf true to wrapping has-ed element in coVF
     * @param opt 0: do nothing, 1: using exnormal, 2: using exfactor on the coefficient
     * @return the collected expression
     */
    ex collect_ex(ex const &expr_in, std::function<bool(const ex &)> has_func, bool cf, bool vf, int opt) {
        auto cvs = collect_lst(expr_in, has_func, opt);
        ex res = 0;
        for(auto cv : cvs) {
            auto cc = cv.op(0);
            auto vv = cv.op(1);
            if(cf) cc = coCF(cc);
            if(vf) vv = coVF(vv);
            res += cc * vv;
        }
        return res;
    }
    
    /**
     * @brief the collect function like Mathematica, reture the lst { {c1,v1}, {c2,v2}, ... }
     * @param expr_in input expression
     * @param has_func only collect the element e, when has_func(e) is true
     * @param opt 0: do nothing, 1: using exnormal, 2: using exfactor on the coefficient
     * @return the collected expression in lst
     */
    lst collect_lst(ex const &expr_in, std::function<bool(const ex &)> has_func, int opt) {
        auto co_epv = inner_expand_collect(expr_in, has_func);
        ex cf = co_epv.first;
        lst res_lst;
        if(!is_zero(cf)) res_lst.append(lst{cf,1});
        for(auto ep : co_epv.second) {
            ex vv = ep.first;
            ex cc = ep.second;
            if(opt==1) cc = exnormal(cc);
            else if(opt==2) cc = exfactor(cc);
            if(!is_zero(cc)) res_lst.append(lst{cc, vv});
        }
        return res_lst;
    }
        
    /**
     * @brief the nuerical evaluation, Digits=100 will be used
     * @param expr input expression
     * @return the nuerical expression
     */
    ex EvalF(ex expr) {
        exset zs;
        //patterns needing evalf()
        expr.find(zeta(w), zs);
        expr.find(zeta(w,w), zs);
        
        lst repl;
        for(auto zi : zs) repl.append(zi==NN(zi));
        return expr.subs(repl);
    }
    
    ex EvalL(ex expr) {
        static symbol sPi("dPi"), sEuler("dEuler"), siEpsilon("diEpsilon");
        lst repl = lst{Pi==sPi, Euler==sEuler,iEpsilon==siEpsilon};
        return EvalF(expr.subs(repl));
    }
    ex EvalQ(ex expr) {
        static symbol sPi("qPi"), sEuler("qEuler"), siEpsilon("qiEpsilon");
        lst repl = lst{Pi==sPi, Euler==sEuler,iEpsilon==siEpsilon};
        return EvalF(expr.subs(repl));
    }
    ex EvalMP(ex expr) {
        static symbol sPi("mpPi"), sEuler("mpEuler"), siEpsilon("mpiEpsilon");
        lst repl = lst{Pi==sPi, Euler==sEuler,iEpsilon==siEpsilon};
        return EvalF(expr.subs(repl));
    }
    ex NN(ex expr, int digits) {
        auto oDigits = Digits;
        Digits = digits;
        auto nexpr = evalf(expr);
        Digits = oDigits;
        return nexpr;
    }

    /**
     * @brief check the expr is xPositive, i.e., each x-monomial item is postive
     * @param expr input expression
     * @return xPositive or not
     */
    bool xPositive(ex const expr) {
        auto tmp = expand(expr);
        if(tmp.is_zero()) return false;
        bool ret = false;
        if(is_a<add>(tmp)) {
            for(auto item : tmp) {
                auto nit = NN(item.subs(x(w)==1).normal());
                if(!is_a<numeric>(nit) || nit<0) return false;
            }
            ret = true;
        } else {
            auto ntmp = NN(tmp.subs(x(w)==1).normal());
            ret = (is_a<numeric>(ntmp) && ntmp>0);
        }
        return ret;
    }

    /**
     * @brief the always sign for expr
     * @param expr input expression
     * @return 1 for xPositive exprs, -1 when -expr is xPositive, 0 for others
     */
    int xSign(ex const expr) {
        if(is_zero(expr)) return 0;
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
     * @param index1 the index
     * @param index2 the index
     */
    void let_op_remove_last(ex & ex_in, int index1, int index2) {
        auto tmp = ex_to<lst>(ex_in.op(index1).op(index2));
        tmp.remove_last();
        ex_in.let_op(index1).let_op(index2) = tmp;
    }
    
    /**
     * @brief remove the last in index1-th.index2-th of expression
     * @param ex_in input expression will be update
     * @param index1 the index
     * @param index2 the index
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
     * @param index1 the index
     * @param index2 the index
     */
    void let_op_remove_first(ex & ex_in, int index1, int index2) {
        auto tmp = ex_to<lst>(ex_in.op(index1).op(index2));
        tmp.remove_first();
        ex_in.let_op(index1).let_op(index2) = tmp;
    }
    
    /**
     * @brief remove the first in index1-th.index2-th of expression
     * @param ex_in input expression will be update
     * @param index1 the index
     * @param index2 the index
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
     * @param item the new item
     */
    void let_op(ex &ex_in, int index1, int index2, const ex item) {
        ex_in.let_op(index1).let_op(index2) = item;
    }
    
    /**
     * @brief update index1-th.index2-th of expression with item
     * @param ex_in input expression will be update
     * @param index1 the index
     * @param index2 the index
     * @param item the new item
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
     * @param item the new item
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
     * @param item the new item
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
            auto cpi = Symbol::set_all(ps.op(i).expand());
            bool ltQ = false; 
            for(auto li : loops) {
                if(cpi.has(li)) {
                    ltQ = true;
                    break;
                }
            }
            
            ex sgn = 0;
            if(!ltQ) {
                pre *= pow(Symbol::set_all(ps.op(i).expand()), ex(0)-ns.op(i));
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
        rem = Symbol::set_all(rem.expand());
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
        rem = Symbol::set_all(rem).normal();
        u = Symbol::set_all(u).normal();
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
    
    /**
     * @brief return the cpu cores using OpenMP
     * @return the cpu cores
     */
    int CpuCores() {
        return omp_get_num_procs();
    }
    
    /**
     * @brief return the numerator and denominator after normalization
     * @param expr the input expression
     * @return fermat evaluated expression
     */
    ex fermat_eval(const ex & expr) {
        static map<pid_t, Fermat> fermat_map;
        static int v_max = 0;

        auto pid = getpid();
        if((fermat_map.find(pid)==fermat_map.end())) { // init section
            fermat_map[pid].Init();
            v_max = 0;
        }
        Fermat &fermat = fermat_map[pid];
        
        auto expr_in = expr;
        exmap map_rat;
        expr_in = expr_in.to_rational(map_rat);
        
        lst rep_vs;
        map<ex,long long,ex_is_less> s2c;
        for(const_preorder_iterator i = expr_in.preorder_begin(); i != expr_in.preorder_end(); ++i) {
            if(is_a<symbol>(*i)) {
                rep_vs.append(*i);
                s2c[*i]++;
            }
        }
        rep_vs.sort();
        rep_vs.unique();
                
        exmap v2f, f2v;
        int fvi = 0;
        for(auto vi : rep_vs) {
            auto name = "v" + to_string(fvi);
            Symbol s(name);
            v2f[vi] = s;
            f2v[s] = vi;
            fvi++;
        }
        expr_in = expr_in.subs(v2f);
        
        stringstream ss;
        if(fvi>111) {
            cout << rep_vs << endl;
            throw Error("Fermat: Too many variables.");
        }
        if(fvi>v_max) {
            for(int i=v_max; i<fvi; i++) ss << "&(J=v" << i << ");" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
            v_max = fvi;
        }
        
        ex res;
        ss << "res:=" << expr_in << ";" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        static string bstr("[-begin-]"), estr("[-end-]");
        ss << "&(U=1);" << endl; // ugly printing, the whitespace matters
        ss << "!('" <<bstr<< "',res,'" <<estr<< "')" << endl;
        auto ostr = fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        // note the order,(normal_fermat will be called again in factor_form)
        ss << "&(U=0);" << endl; // disable ugly printing
        ss << "@(res);" << endl;
        ss << "&_G;" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");

        // make sure last char is 0
        if(ostr[ostr.length()-1]!='0') throw Error("fermat_together: last char is NOT 0.");
        ostr = ostr.substr(0, ostr.length()-1);
        auto cpos = ostr.find(bstr);
        if(cpos==string::npos) throw Error(bstr+" NOT Found.");
        ostr = ostr.substr(cpos+bstr.length(),string::npos);
        cpos = ostr.find(estr);
        if(cpos==string::npos) throw Error(estr+" NOT Found.");
        ostr = ostr.substr(0,cpos);
        string_trim(ostr);

        symtab st;
        Parser fp(st);
        res = fp.Read(ostr);
        res = res.subs(f2v);
        return res.subs(map_rat);
    }
    
    
    /**
     * @brief return the numerator and denominator after normalization
     * @param expr the input expression
     * @param dfactor true for factorize on the denominator
     * @return a list of { numer, denom }
     */
    ex numer_denom_fermat(const ex & expr, bool dfactor) {
        static map<pid_t, Fermat> fermat_map;
        static int v_max = 0;
        bool use_ncheck = false;

        auto pid = getpid();
        if((fermat_map.find(pid)==fermat_map.end())) { // init section
            fermat_map[pid].Init();
            v_max = 0;
        }
        Fermat &fermat = fermat_map[pid];
        
        auto expr_in = expr;
        exmap map_rat;
        expr_in = expr_in.to_rational(map_rat);
        
        lst rep_vs;
        map<ex,long long,ex_is_less> s2c;
        for(const_preorder_iterator i = expr_in.preorder_begin(); i != expr_in.preorder_end(); ++i) {
            if(is_a<symbol>(*i)) {
                rep_vs.append(*i);
                s2c[*i]++;
            }
        }
        rep_vs.sort();
        rep_vs.unique();
        
        exvector vs_vec_sort;
        auto rep_vs_tot = rep_vs.nops();
        for(int i=0; i<rep_vs_tot; i++) {
            auto ss = rep_vs.op(i);
            vs_vec_sort.push_back(lst{ s2c[ss], ss.subs(map_rat), ss });
        }
        sort_vec(vs_vec_sort);
        rep_vs.remove_all();
        for(auto item : vs_vec_sort) rep_vs.append(item.op(2));
                
        exmap v2f, f2v;
        exmap nn_map;
        auto nn_pi1 = cln::nextprobprime(3);
        auto nn_pi2 = cln::nextprobprime(3);
        int fvi = 0;
        for(auto vi : rep_vs) {
            auto name = "v" + to_string(fvi);
            Symbol s(name);
            v2f[vi] = s;
            f2v[s] = vi;
            fvi++;
            nn_pi1 = cln::nextprobprime(nn_pi2+1);
            nn_pi2 = cln::nextprobprime(nn_pi1+1);
            nn_map[s] = numeric(nn_pi1)/numeric(nn_pi2);
        }
        
        stringstream ss;
        if(fvi>111) {
            cout << rep_vs << endl;
            throw Error("Fermat: Too many variables.");
        }
        if(fvi>v_max) {
            for(int i=v_max; i<fvi; i++) ss << "&(J=v" << i << ");" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
            v_max = fvi;
        }
        
        ex nn_chk=0, num, den;
        lst item;
        if(!is_a<add>(expr_in)) item = lst{expr_in};
        else for(auto ii : expr_in) item.append(ii);
        sort_lst(item);
        if(fermat_using_array) ss << "Array m[" << item.nops() << "];" << endl;
        else ss << "res:=0;" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        for(int i=0; i<item.nops(); i++) {
            ex tt = item.op(i).subs(v2f);
            if(use_ncheck) nn_chk += tt.subs(nn_map);
            if(fermat_using_array) ss << "m[" << (i+1) << "]:=";
            else ss << "item:=";
            ss << tt << ";" << endl;
            if(!fermat_using_array) ss << "res:=res+item;" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
        }
        if(fermat_using_array) {
            ss << "res:=Sumup([m]);" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
        }
        
        static string bstr("[-begin-]"), estr("[-end-]");
        ss << "&(U=1);" << endl; // ugly printing, the whitespace matters
        ss << "!('" <<bstr<< "','{',Numer(res),',',Denom(res),'}','" <<estr<< "')" << endl;
        auto ostr = fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        // note the order,(normal_fermat will be called again in factor_form)
        ss << "&(U=0);" << endl; // disable ugly printing
        if(fermat_using_array) ss << "@(res,[m]);" << endl;
        else ss << "@(res,item);" << endl;
        ss << "&_G;" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");

        // make sure last char is 0
        if(ostr[ostr.length()-1]!='0') throw Error("fermat_together: last char is NOT 0.");
        ostr = ostr.substr(0, ostr.length()-1);
        auto cpos = ostr.find(bstr);
        if(cpos==string::npos) throw Error(bstr+" NOT Found.");
        ostr = ostr.substr(cpos+bstr.length(),string::npos);
        cpos = ostr.find(estr);
        if(cpos==string::npos) throw Error(estr+" NOT Found.");
        ostr = ostr.substr(0,cpos);
        string_trim(ostr);

        symtab st;
        Parser fp(st);
        auto ret = fp.Read(ostr);
        num = ret.op(0);
        if(dfactor) den = factor_form(ret.op(1));
        else den = ret.op(1);
        //fermat.Exit();
        
        if(use_ncheck) {
            auto nn_ret = subs(num/den,nn_map);
            if(nn_chk-nn_ret!=0) {
                cout << nn_chk << " : " << nn_ret << endl;
                throw Error("fermat_together: N Check Failed.");
            }
        }
        
        num = num.subs(f2v).subs(map_rat);
        den = den.subs(f2v).subs(map_rat);
        return lst{num, den};
        
    }
    
    /**
     * @brief return the normalizatied expression, using fermat_numer_denom
     * @param expr the input expression
     * @param dfactor true for factorize on the denominator
     * @return the normalized expression: numer/denom
     */
    ex normal_fermat(const ex & expr, bool dfactor) {
        auto nd = fermat_numer_denom(expr, dfactor);
        return nd.op(0)/nd.op(1);
    }
    
    /**
     * @brief a wrapper for collect_common_factors, catch errors
     * @param expr the input expression
     * @return collect_common_factors result
     */
    ex collect_factors(const ex & expr) {
        try {
            return collect_common_factors(expr);
        } catch(...) { }
        return expr;
    }
    
    /**
     * @brief factorize a expression
     * @param expr the input expression
     * @param opt 1 to use FORM, otherwise using GiNaC for factorization
     * @return factorized result
     */
    ex exfactor(const ex & expr, int opt) {
        if(opt==1) return factor_form(expr);
        else return ginac_factor(expr);
    }
    
    /**
     * @brief factorize a expression
     * @param expr the input expression
     * @param opt 1 to use FORM, otherwise using GiNaC for factorization
     * @return factorized result
     */
    ex exexpand(const ex & expr, int opt) {
        if(opt==1) return expand(expr);
        else return expand(expr);
    }
    
    /**
     * @brief normalize a expression
     * @param expr the input expression
     * @param opt 1 to use Fermat, otherwise using GiNaC for normalization
     * @return factorized result
     */
    ex exnormal(const ex & expr, int opt) {
        if(opt==1) return normal_fermat(expr,true);
        else return normal(expr);
    }
    
    /**
     * @brief num_den a expression
     * @param expr the input expression
     * @param opt  1 to use factor_form, 2 to use Fermat, otherwise using GiNaC for numer_denom
     * @return lst of { num, den }
     */
    ex exnd(const ex & expr, int opt) {
        if(opt==1) return numer_denom(factor_form(expr));
        else if(opt==2) return numer_denom_fermat(expr);
        else return numer_denom(expr);
    }
    
    ex form_eval(const ex & expr) {
        static map<pid_t, Form> form_map;
        auto pid = getpid();
        if((form_map.find(pid)==form_map.end())) { // init section
            ostringstream ss;
            ss << "AutoDeclare Symbols viff;" << endl;
            form_map[pid].Init("form");
            form_map[pid].Execute(ss.str());
        }
        Form &fprc = form_map[pid];
        
        auto expr_in = expr;
        exmap map_rat;
        expr_in = expr_in.to_polynomial(map_rat);
        
        exvector vs;
        for(const_preorder_iterator i = expr_in.preorder_begin(); i != expr_in.preorder_end(); ++i) {
            auto e = (*i);
            if(is_a<symbol>(e)) vs.push_back(e);
        }
        
        exmap s2v;
        symtab st;
        int cid = 0;
        for(auto ss : vs) {
            cid++;
            string cvv = "viff"+to_string(cid);
            s2v[ss] = Symbol(cvv);
            st[cvv] = ss.subs(map_rat);
        }
        ostringstream oss;
        expr_in = expr_in.subs(s2v);
        oss << "Local ff = " << expr_in << ";" << endl;
        oss << ".sort" << endl;
        try {
            auto ostr = fprc.Execute(oss.str(), "ff");
            Parser fp(st);
            ex ret = fp.Read(ostr);
            return ret;
        } catch(Error& err) {
            cout << oss.str() << endl;
            form_map.erase(pid);
            throw;
        }
    }
    
    ex inner_factor_form(const ex & expr) {
        static map<pid_t, Form> form_map;
        auto pid = getpid();
        if((form_map.find(pid)==form_map.end())) { // init section
            ostringstream ss;
            ss << "AutoDeclare Symbols viff;" << endl;
            form_map[pid].Init("form");
            form_map[pid].Execute(ss.str());
        }
        Form &fprc = form_map[pid];
        
        auto expr_in = expr;
        exmap map_rat;
        expr_in = expr_in.to_polynomial(map_rat);
        
        exvector vs;
        for(const_preorder_iterator i = expr_in.preorder_begin(); i != expr_in.preorder_end(); ++i) {
            auto e = (*i);
            if(is_a<symbol>(e)) vs.push_back(e);
        }
        
        exmap s2v;
        symtab st;
        int cid = 0;
        for(auto ss : vs) {
            cid++;
            string cvv = "viff"+to_string(cid);
            s2v[ss] = Symbol(cvv);
            st[cvv] = ss.subs(map_rat);
        }
        ostringstream oss;
        expr_in = expr_in.subs(s2v);
        oss << "Local ff = " << expr_in << ";" << endl;
        oss << "Factorize ff;" << endl;
        oss << ".sort" << endl;
        try {
            auto ostr = fprc.Execute(oss.str(), "ff");
            string_replace_all(ostr, "[", "(");
            string_replace_all(ostr, "]", ")");
            string_replace_all(ostr, "\\\n", "");
            string_replace_all(ostr, " ","");
            Parser fp(st);
            ex ret = fp.Read(ostr);
            ex res = 1;
            Symbol sfac("factor_");
            for(auto item : add2lst(ret)) res *= item.subs(sfac==1);
            return res;
        } catch(Error& err) {
            cout << oss.str() << endl;
            form_map.erase(pid);
            throw;
        }
    }
    
    /**
     * @brief factorize a expression using FORM
     * @param expr the input expression
     * @return factorized result
     */
    ex factor_form(const ex & expr, bool nd) {
        if(!nd) return inner_factor_form(expr);
        auto num_den = fermat_numer_denom(expr);
        if(is_zero(num_den.op(1)-1)) return inner_factor_form(num_den.op(0));
        return inner_factor_form(num_den.op(0))/inner_factor_form(num_den.op(1));
    }
    
    //-----------------------------------------------------------
    // HepFormat Output
    //-----------------------------------------------------------
    HepFormat::HepFormat(ostream &os, unsigned opt) : print_dflt(os, opt) {}
    HepFormat::HepFormat() : print_dflt(std::cout) {}
    GINAC_IMPLEMENT_PRINT_CONTEXT(HepFormat, print_dflt)
    
    const HepFormat & HepFormat::operator << (const basic & v) const {
        v.print(*this);
        return *this;
    }
    const HepFormat & HepFormat::operator << (const ex & v) const {
        v.print(*this);
        return *this;
    }
    const HepFormat & HepFormat::operator << (const lst & v) const {
        v.print(*this);
        return *this;
    }
    const HepFormat & HepFormat::operator<<(std::ostream& (*v)(std::ostream&)) const {
        s << v;
        return *this;
    }
        
    void HepFormat::add_print(const add & a, const HepFormat & c, unsigned level) {
        auto as = add2lst(a);
        sort_lst(as);
        auto cl = a.precedence();
        bool first = true;
        if(cl<=level) c.s << '(';
        for(auto item : as) {
            if(!first) c.s << "+";
            item.print(c, cl);
            first = false;
        }
        if(a.precedence()<=level) c.s << ')';
    }
    
    void HepFormat::mul_print(const mul & m, const HepFormat & c, unsigned level) {
        auto ms = mul2lst(m);
        sort_lst(ms);
        auto cl = m.precedence();
        
        // handle negative number
        int nn = ms.nops();
        auto ex0 = ms.op(0);
        if(nn>1 && ex0.info(info_flags::real) && ex0<0) {
            ex exn = ms.op(nn-1);
            if(is_a<add>(exn)) {
                exn = numeric(-1) * exn;
                if(is_a<add>(exn)) {
                    ms.let_op(0) = numeric(-1) * ms.op(0);
                    ms.let_op(nn-1) = exn;
                }
            }
        }
        
        bool first = true;
        if(cl<=level) c.s << '(';
        for(auto item : ms) {
            if(is_a<numeric>(item) && is_zero(item-1)) continue;
            if(!first) c.s << "*";
            item.print(c, cl);
            first = false;
        }
        if(cl<=level) c.s << ')';
    }
    
    const HepFormat & HepFormat::operator << (const matrix & mat) const {
        s << "[";
        int nr = mat.rows();
        int nc = mat.cols();
        for(int r=0; r<nr; r++) {
            s << "[";
            for(int c=0; c<nc; c++) {
                mat(r,c).print(*this);
                if(c+1!=nc) s << ",";
            }
            s << "]";
            if(r+1!=nr) s << ",";
        }
        s << "]";
        return *this;
    }
    
    const HepFormat & HepFormat::operator << (const exvector & e) const {
        auto i = e.begin();
        auto vend = e.end();
        if (i==vend) { s << "[]"; return *this; }
        s << "[";
        while (true) {
            i->print(*this);
            ++i;
            if(i==vend) break;
            s << ",";
        }
        s << "]";
        return *this;
    }

    const HepFormat & HepFormat::operator << (const exset & e) const {
        auto i = e.begin();
        auto send = e.end();
        if (i==send) { s << "<>"; return *this; }
        s << "<";
        while (true) {
            i->print(*this);
            ++i;
            if(i==send) break;
            s << ",";
        }
        s << ">";
        return *this;
    }

    const HepFormat & HepFormat::operator << (const exmap & e) const {
        auto i = e.begin();
        auto mend = e.end();
        if (i==mend) { s << "{}"; return *this; }
        s << "{";
        while (true) {
            i->first.print(*this);
            s << "==";
            i->second.print(*this);
            ++i;
            if(i==mend) break;
            s << ",";
        }
        s << "}";
        return *this;
    }
    
    //-----------------------------------------------------------
    // MMAFormat Output
    //-----------------------------------------------------------
    MMAFormat::MMAFormat(ostream &os, unsigned opt) : print_dflt(os, opt) {}
    MMAFormat::MMAFormat() : print_dflt(std::cout) {}
    GINAC_IMPLEMENT_PRINT_CONTEXT(MMAFormat, print_dflt)
    
    const MMAFormat & MMAFormat::operator << (const basic & v) const {
        v.print(*this);
        return *this;
    }
    const MMAFormat & MMAFormat::operator << (const ex & v) const {
        v.print(*this);
        return *this;
    }
    const MMAFormat & MMAFormat::operator << (const lst & v) const {
        v.print(*this);
        return *this;
    }
    const MMAFormat & MMAFormat::operator<<(std::ostream& (*v)(std::ostream&)) const {
        s << v;
        return *this;
    }
        
    void MMAFormat::add_print(const add & a, const MMAFormat & c, unsigned level) {
        auto as = add2lst(a);
        sort_lst(as);
        auto cl = a.precedence();
        bool first = true;
        if(cl<=level) c.s << '(';
        for(auto item : as) {
            if(!first) c.s << "+";
            item.print(c, cl);
            first = false;
        }
        if(a.precedence()<=level) c.s << ')';
    }
    
    void MMAFormat::mul_print(const mul & m, const MMAFormat & c, unsigned level) {
        auto ms = mul2lst(m);
        sort_lst(ms);
        auto cl = m.precedence();
        
        // handle negative number
        int nn = ms.nops();
        auto ex0 = ms.op(0);
        if(nn>1 && ex0.info(info_flags::real) && ex0<0) {
            ex exn = ms.op(nn-1);
            if(is_a<add>(exn)) {
                exn = numeric(-1) * exn;
                if(is_a<add>(exn)) {
                    ms.let_op(0) = numeric(-1) * ms.op(0);
                    ms.let_op(nn-1) = exn;
                }
            }
        }
        
        bool first = true;
        if(cl<=level) c.s << '(';
        for(auto item : ms) {
            if(is_a<numeric>(item) && is_zero(item-1)) continue;
            if(!first) c.s << "*";
            item.print(c, cl);
            first = false;
        }
        if(cl<=level) c.s << ')';
    }
    
    const MMAFormat & MMAFormat::operator << (const matrix & mat) const {
        s << "{";
        int nr = mat.rows();
        int nc = mat.cols();
        for(int r=0; r<nr; r++) {
            s << "{";
            for(int c=0; c<nc; c++) {
                mat(r,c).print(*this);
                if(c+1!=nc) s << ",";
            }
            s << "}";
            if(r+1!=nr) s << ",";
        }
        s << "}";
        return *this;
    }
    
    const MMAFormat & MMAFormat::operator << (const exvector & e) const {
        auto i = e.begin();
        auto vend = e.end();
        if (i==vend) { s << "{}"; return *this; }
        s << "{";
        while (true) {
            i->print(*this);
            ++i;
            if(i==vend) break;
            s << ",";
        }
        s << "}";
        return *this;
    }

    const MMAFormat & MMAFormat::operator << (const exset & e) const {
        auto i = e.begin();
        auto send = e.end();
        if (i==send) { s << "{}"; return *this; }
        s << "{";
        while (true) {
            i->print(*this);
            ++i;
            if(i==send) break;
            s << ",";
        }
        s << "}";
        return *this;
    }

    const MMAFormat & MMAFormat::operator << (const exmap & e) const {
        auto i = e.begin();
        auto mend = e.end();
        if (i==mend) { s << "{}"; return *this; }
        s << "{";
        while (true) {
            i->first.print(*this);
            s << "->";
            i->second.print(*this);
            ++i;
            if(i==mend) break;
            s << ",";
        }
        s << "}";
        return *this;
    }
    
    void garWrite(exvector &exv, string garfn) {
        archive ar;
        ar.archive_ex(exv.size(), "size");
        for(int i=0; i<exv.size(); i++) {
            ar.archive_ex(exv[i], to_string(i).c_str());
        }
        ofstream out(garfn);
        out << ar;
        out.close();
    }
        
    void garRead(exvector &exv, string garfn) {
        archive ar;
        ifstream in(garfn);
        in >> ar;
        in.close();
        auto size = ex_to<numeric>(ar.unarchive_ex(GiNaC_archive_Symbols, "size")).to_int();
        exvector res;
        for(int i=0; i<size; i++) {
            ex tmp = ar.unarchive_ex(GiNaC_archive_Symbols, to_string(i).c_str());
            exv.push_back(tmp);
        }
    } 
    
    void exVectorPut(exvector &exv, string garfn) {
        garWrite(exv.size(), garfn+".gar");
        GiNaC_Parallel(exv.size(), [&exv,garfn](int idx) {
            garWrite(exv[idx], garfn+"-"+to_string(idx)+".gar");
            return 0;
        }, "vecPut");
    }
    
    void exVectorGet(exvector &exv, string garfn) {
        auto size = ex_to<numeric>(garRead(garfn+".gar")).to_int();
        if(exv.size()>0) {
            if(size != exv.size()) throw Error("size not matched.");
            for(int i=0; i<size; i++) exv[i] = garRead(garfn+"-"+to_string(i)+".gar");
        } else {
            for(int i=0; i<size; i++) exv.push_back(garRead(garfn+"-"+to_string(i)+".gar"));
        }
    } 
    
    ex add_collect_normal(const exvector &exv, lst const &pats) {
        auto cvs_vec = GiNaC_Parallel(exv.size(), [&exv,pats](int idx)->ex {
            return collect_lst(exv[idx], pats);
        }, "ColEx");
        exmap res_map;
        for(auto cvs : cvs_vec) for(auto cv : cvs) res_map[cv.op(1)] += cv.op(0);
        exvector res_vec;
        for(auto kv : res_map) res_vec.push_back(lst{kv.second, kv.first});
        res_vec = GiNaC_Parallel(res_vec.size(), [&res_vec](int idx)->ex {
            return exnormal(res_vec[idx].op(0)) * res_vec[idx].op(1);
        }, "NorEx");
        return add(res_vec);
    }
        
}
