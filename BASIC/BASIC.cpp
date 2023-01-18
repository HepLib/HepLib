/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "BASIC.h"
#include "cln/cln.h"

inline unsigned golden_ratio_hash(uintptr_t n) {
	return n * UINT64_C(0x4f1bbcdd);
}

extern "C" {
    #include <quadmath.h>
}

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
    
    namespace {
        const symbol & get_symbol(const string & s) {
            static map<string, symbol> dict;
            string key = s;
            if (dict.find(key) == dict.end()) dict[key] = symbol(s,0);
            return dict[key];
        }
        
        const Symbol & get_Symbol(const string & s) {
            static map<string, Symbol> dict;
            string key = s;
            if (dict.find(key) == dict.end()) dict[key] = Symbol(s);
            return dict[key];
        }
        
        const iSymbol & get_iSymbol(const string & s) {
            static map<string, iSymbol> dict;
            string key = s;
            if (dict.find(key) == dict.end()) dict[key] = iSymbol(s);
            return dict[key];
        }
    }   

    //DEFAULT_CTOR(Symbol)
    IMPLEMENT_HAS(Symbol)
    IMPLEMENT_ALL(Symbol)
    //GINAC_IMPLEMENT_REGISTERED_CLASS(Symbol, symbol)
    GiNaC::registered_class_info & Symbol::get_class_info_static() { return reg_info; }
    Symbol::visitor::~visitor() { }
    Symbol * Symbol::duplicate() const { Symbol * bp = new Symbol(*this); bp->setflag(GiNaC::status_flags::dynallocated); return bp; }
    void Symbol::accept(GiNaC::visitor & v) const { if (visitor *p = dynamic_cast<visitor *>(&v)) p->visit(*this); else inherited::accept(v); }
    const GiNaC::registered_class_info &Symbol::get_class_info() const { return get_class_info_static(); }
    GiNaC::registered_class_info &Symbol::get_class_info() { return get_class_info_static(); }
    const char *Symbol::class_name() const { return get_class_info_static().options.get_name(); }
    //GINAC_IMPLEMENT_REGISTERED_CLASS END
    
    /**
     * @brief Symbol constructor
     * @param s symbol name
     */
    Symbol::Symbol(const string &s) : symbol(get_symbol(s)) { Table[s]=*this; }
    Symbol::Symbol() : symbol(get_symbol("")) { Table[""]=*this; }
    
    int Symbol::compare_same_type(const basic &other) const {
        if(!is_a<Symbol>(other)) throw Error("Symbol::compare_same_type");;
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
     */
    void Symbol::read_archive(const archive_node& n) {
        inherited::read_archive(n);
    }
    
    ex Symbol::eval() const { return *this; }
    ex Symbol::evalf() const { return *this; }
    ex Symbol::conjugate() const { return *this; }
    ex Symbol::real_part() const { return *this; }
    ex Symbol::imag_part() const { return 0; }
    unsigned Symbol::calchash() const {
        static std::hash<std::string> hs;
        unsigned seed = hs(get_name()+typeid(*this).name());
        hashvalue = golden_ratio_hash(seed);
        setflag(status_flags::hash_calculated);
        return hashvalue;
    }
    
    bool Symbol::is_equal_same_type(const basic & other) const {
        if(!is_a<Symbol>(other)) throw Error("Symbol::is_equal_same_type");
        const Symbol *o = static_cast<const Symbol *>(&other);
        return get_name()==o->get_name(); // should be the same as symbol class
    }
    
    ex Symbol::derivative(const symbol & s) const { 
        if(!is_exactly_a<Symbol>(s)) return 0;
        if(is_equal_same_type(s)) return 1;
        else return 0;
    }
    
    ex Symbol::series(const relational & r, int order, unsigned options) const {
        epvector seq;
        const ex point = r.rhs();
    
        if (is_exactly_a<Symbol>(r.lhs()) && this->is_equal_same_type(ex_to<Symbol>(r.lhs()))) {
            if (order>0 && !point.is_zero()) seq.emplace_back(expair(point, 0));
            if (order>1) seq.emplace_back(expair(1, 1));
            else seq.emplace_back(expair(Order(1), numeric(order)));
        } else seq.emplace_back(expair(*this, 0));
        return pseries(r, std::move(seq));
    }
    
    void Symbol::set_name(string n) { throw Error("Symbol can not reset the name!"); }
    
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
    
    //DEFAULT_CTOR(iSymbol)
    IMPLEMENT_ALL(iSymbol)
    //GINAC_IMPLEMENT_REGISTERED_CLASS(iSymbol, symbol)
    GiNaC::registered_class_info & iSymbol::get_class_info_static() { return reg_info; }
    iSymbol::visitor::~visitor() { }
    iSymbol * iSymbol::duplicate() const { iSymbol * bp = new iSymbol(*this); bp->setflag(GiNaC::status_flags::dynallocated); return bp; }
    void iSymbol::accept(GiNaC::visitor & v) const { if (visitor *p = dynamic_cast<visitor *>(&v)) p->visit(*this); else inherited::accept(v); }
    const GiNaC::registered_class_info &iSymbol::get_class_info() const { return get_class_info_static(); }
    GiNaC::registered_class_info &iSymbol::get_class_info() { return get_class_info_static(); }
    const char *iSymbol::class_name() const { return get_class_info_static().options.get_name(); }
    //GINAC_IMPLEMENT_REGISTERED_CLASS END
    
    /**
     * @brief Symbol constructor
     * @param s symbol name
     */
    iSymbol::iSymbol(const string &s) : symbol(get_symbol(s)) { Table[s]=*this; }
    iSymbol::iSymbol() : symbol(get_symbol("")) { Table[""]=*this; }
    int iSymbol::compare_same_type(const basic &other) const {
        if(!is_a<iSymbol>(other)) throw Error("iSymbol::compare_same_type");
        const iSymbol &o = static_cast<const iSymbol &>(other);
        int ret = get_name().compare(o.get_name());
        if(ret==0) return 0;
        else if(ret<0) return -1;
        else return 1;
    }
    
    bool iSymbol::is_equal_same_type(const basic & other) const {
        if(!is_a<iSymbol>(other)) throw Error("iSymbol::is_equal_same_type");
        const iSymbol *o = static_cast<const iSymbol *>(&other);
        return serial==o->serial; // should be the same as symbol class
    }
    
    ex iSymbol::derivative(const symbol & s) const {
        if(!is_exactly_a<iSymbol>(s)) return 0;
        if(is_equal_same_type(s)) return 1;
        else return 0;
    }
    
    ex iSymbol::series(const relational & r, int order, unsigned options) const {
        epvector seq;
        const ex point = r.rhs();
    
        if (is_exactly_a<iSymbol>(r.lhs()) && this->is_equal_same_type(ex_to<iSymbol>(r.lhs()))) {
            if (order>0 && !point.is_zero()) seq.emplace_back(expair(point, 0));
            if (order>1) seq.emplace_back(expair(1, 1));
            else seq.emplace_back(expair(Order(1), numeric(order)));
        } else seq.emplace_back(expair(*this, 0));
        return pseries(r, std::move(seq));
    }
    
    void iSymbol::set_name(string n) { throw Error("iSymbol can not reset the name!"); }
    
    unsigned iSymbol::calchash() const {
        hashvalue = get_symbol(get_name()).gethash();
        setflag(status_flags::hash_calculated);
        return hashvalue;
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
     */
    void iSymbol::read_archive(const archive_node& n) {
        inherited::read_archive(n);
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
     * @param f function to be applied on the index, from 0 to (ntotal-1)
     * @param key key used in archive file name and display message
     * @param rm default true, and false will keep the archive file
     * @param pre the pre-string in the print message
     * @return return the ntotal-element vector, i.e., [f(0), ..., f(ntotal-1)]
     */
    exvector GiNaC_Parallel(int ntotal, std::function<ex(int)> f, const string & key) {
        if(ntotal<1) return exvector();
        int ec = 0;
        int nactive = 0;
        string exstr = "";
        int fork_retried = 0;
        int nproc = GiNaC_Parallel_Process;
        if(GiNaC_Parallel_NP.find(key)!=GiNaC_Parallel_NP.end()) nproc = GiNaC_Parallel_NP[key];
        
        bool rm = true;
        if(GiNaC_Parallel_RM.find(key)!=GiNaC_Parallel_RM.end()) rm = GiNaC_Parallel_RM[key]; 
        string pre = "  ";
        if(GiNaC_Parallel_PRE.find(key)!=GiNaC_Parallel_PRE.end()) pre = GiNaC_Parallel_PRE[key]; 
        
        int verb = Verbose;
        if(GiNaC_Parallel_Verb.find(key)!=GiNaC_Parallel_Verb.end()) verb = GiNaC_Parallel_Verb[key]; 
        
        // nproc=0, non-parallel
        if(nproc==0) {
            exvector ovec;
            for(int i=0; i<ntotal; i++) ovec.push_back(f(i));
            return exvector(std::move(ovec));
        }
        
        int para_max_run = nproc<0 ? omp_get_num_procs()-1 : nproc;
        if(para_max_run<0) para_max_run = 0;
        if(para_max_run>omp_get_num_procs()-1) para_max_run = omp_get_num_procs()-1;
        
        auto ppid = getpid();
        ostringstream cmd;
        cmd << "mkdir -p " << ppid;
        if(!dir_exists(to_string(ppid))) system(cmd.str().c_str());
        
        int nbatch = GiNaC_Parallel_Batch;
        if(GiNaC_Parallel_NB.find(key)!=GiNaC_Parallel_NB.end()) nbatch = GiNaC_Parallel_NB[key];
        if(nbatch<=0) nbatch = ntotal/para_max_run/5;
        else if(nbatch > ntotal/para_max_run) nbatch = ntotal/para_max_run;
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
        while(fork_retried<5) {
            ec = 0;
            exstr = "";
            nactive = 0;
            for(int bi=0; bi<btotal; bi++) {
                restart: ;
                for(int i=0; i<btotal; i++) if(waitpid(-pgid,NULL,WNOHANG)>0) nactive--;
                if(verb > 1 && !nst) {
                    cout << "\r                                                                  \r" << pre;
                    cout << "\\--Evaluating ";
                    if(key != "") cout << Color_HighLight << key << RESET << " ";
                    cout << Color_HighLight << nbatch << "x" << RESET << "[" << (bi+1) << "/" << btotal << "] @ " << now(false) << exstr << flush;
                }
                
                if(fork_retried>0) { // skip when bi.*.gar exists
                    ostringstream garfn;
                    if(key == "") garfn << ppid << "/" << bi << ".gar";
                    else garfn << ppid << "/" << bi << "." << key << ".gar";
                    if(file_exists(garfn.str())) continue;
                }
                
                auto pid = fork();
                if(pid < 0) {
                    if(getpid()!=pgid) exit(0); // make sure exit child process
                    ec++;
                    if(ec>3*btotal) throw Error("GiNaC_Parallel: Error (2) @ fork()");
                    exstr = " [fork(" + to_string(ec) + ")]";
                    if(waitpid(-pgid,NULL,0)>0) nactive--;
                    goto restart; // parent process goes next cycle
                }
                if(pid==0 && setpgid(0, pgid)!=0) {
                    if(setpgid(0, pgid)) throw Error("GiNaC_Parallel: setpgid(0, pgid) Failed.");
                }
                if(pid>0) {
                    nactive++;
                    if(nactive >= para_max_run) if(waitpid(-pgid,NULL,0)>0) nactive--;
                    continue; // parent process goes next cycle
                }
                
                // pid = 0, child process
                In_GiNaC_Parallel = true;
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
            if(!nst && verb > 1) cout << endl;
            if(getpid()!=ppid) exit(0); // make the forked group leader exit
            
            if(ec<2) break;
            // check all *.gar and retry
            bool all_gar_exists = true;
            for(int bi=0; bi<btotal; bi++) {
                ostringstream garfn;
                if(key == "") garfn << ppid << "/" << bi << ".gar";
                else garfn << ppid << "/" << bi << "." << key << ".gar";
                if(!file_exists(garfn.str())) {
                    all_gar_exists = false;
                    break;
                }
            }
            if(all_gar_exists) break;
            para_max_run = para_max_run/2;
            if(para_max_run<2) break;
            fork_retried++;
        }
        
        wait_label:
        if(npid!=0) waitpid(npid,NULL,0); // wait nest process to exit
        
        //#define use_gar_write // use ReWR
        
        exvector ovec;
        if(true) {
            #ifdef use_gar_write
            exvector ovec_tmp;
            #endif
            for(int bi=0; bi<btotal; bi++) {
                if(verb > 1 && !nst) {
                    if(key == "") {
                        cout << "\r                                                   \r" << pre;
                        cout << "\\--Reading *.gar [" << (bi+1) << "/" << btotal << "] @ " << now(false) << flush;
                    } else {
                        cout << "\r                                                   \r" << pre;
                        cout << "\\--Reading *." << Color_HighLight << key << RESET << ".gar [" << (bi+1) << "/" << btotal << "] @ " << now(false) << flush;
                    }
                }

                ostringstream garfn;
                if(key == "") garfn << ppid << "/" << bi << ".gar";
                else garfn << ppid << "/" << bi << "." << key << ".gar";
                lst res_lst;
                try {
                    if(file_exists(garfn.str())) {
                        res_lst = ex_to<lst>(garRead(garfn.str()));
                        remove(garfn.str().c_str());
                        goto done;
                    } 
                } catch(exception &p) { }
                if(verb > 1 && !nst) cout << " - ReTry" << endl;
                try {
                    res_lst.remove_all();
                    for(int ri=0; ri<nbatch; ri++) {
                        int i = bi*nbatch + ri;
                        if(i<ntotal) res_lst.append(f(i));
                        else break;
                    }
                } catch(exception &p) { 
                    cout << ErrColor << "Failed in GiNaC_Parallel!" << RESET << endl;
                    cout << ErrColor << p.what() << RESET << endl;
                    throw Error("GiNaC_Parallel_ReTRY: "+string(p.what()));
                }
                done: ;
                #ifdef use_gar_write
                if(GiNaC_Parallel_ReWR.find(key)==GiNaC_Parallel_ReWR.end() || GiNaC_Parallel_ReWR[key]) {
                    for(auto res : res_lst) ovec_tmp.push_back(res);
                } else {
                    for(auto res : res_lst) ovec.push_back(res);
                }
                #else 
                for(auto res : res_lst) ovec.push_back(res);
                #endif
            }
            
            #ifdef use_gar_write
            if(!nst && verb>1) {
                cout << endl;
                if(key == "") {
                    cout << "\r                                                   \r" << pre;
                    cout << "\\--ReWR" << flush;
                } else {
                    cout << "\r                                                   \r" << pre;
                    cout << "\\--ReWR " << Color_HighLight << key << RESET << flush;
                }
            }
            
            if(GiNaC_Parallel_ReWR.find(key)==GiNaC_Parallel_ReWR.end() || GiNaC_Parallel_ReWR[key]) {
                ostringstream garfn;
                if(key == "") garfn << ppid << "/ReWR.gar";
                else garfn << ppid << "/ReWR." << key << ".gar";
                garWrite(garfn.str(), ovec_tmp);
                ovec_tmp.clear();
                garRead(garfn.str(), ovec);
            }
            #endif
            
            if(!nst && verb>1) cout << endl;
        }
        
        if(rm) {
            cmd.clear();
            cmd.str("");
            cmd << "rm -fr " << ppid;
            system(cmd.str().c_str());
        }
        if(ovec.size() != ntotal) {
            throw Error("GiNaC_Parallel: The output size is wrong!");
        }
        
        return exvector(std::move(ovec));
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
        exset ss;
        for(const_preorder_iterator i = e.preorder_begin(); i != e.preorder_end(); ++i) {
            if(is_a<symbol>(*i)) ss.insert(*i);
        }
        lst sym_lst;
        for(auto item : ss) sym_lst.append(item);
        return sym_lst;
    }

    /**
     * @brief get all symbols from input expression
     * @param ve input expression vector
     * @return all symbols in the input
     */
    lst gather_symbols(const exvector & ve) {
        exset ss;
        for(auto e : ve) {
            for(const_preorder_iterator i = e.preorder_begin(); i != e.preorder_end(); ++i) {
                if(is_a<symbol>(*i)) ss.insert(*i);
            }
        }
        lst sym_lst;
        for(auto item : ss) sym_lst.append(item);
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
            ex res = ar.unarchive_ex(name, i);
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
        auto res = ar.unarchive_ex(key);
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
        auto c = ar.unarchive_ex("c");
        auto res = ar.unarchive_ex("res");
        if(c!=19790923) throw Error("garRead: check faild for file: " + garfn);
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
    
    matrix lst2mat(const lst & ls) {
        int nr = ls.nops();
        int nc = ls.op(0).nops();
        matrix mat(nr,nc);
        for(int r=0; r<nr; r++) {
            auto row = ls.op(r);
            for(int c=0; c<nc; c++) mat(r,c) = row.op(c);
        }
        return mat;
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
    
    void str2file(const string & ostr, string filename) {
        std::ofstream ofs;
        ofs.open(filename, ios::out);
        ofs << ostr << endl;
        ofs.close();
    }
    
    // not finished, just for memo
    void str2file(char * buff, FILE* fh) {
        string nstr;
        int n = nstr.length();
        char nbuff[n+1];
        strcpy(nbuff, nstr.c_str());
        FILE * f = fmemopen(nbuff, n, "r");
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
        set_precision(40);
        nss << num.evalf() << endl;
        __float128 ret = strtoflt128(nss.str().c_str(), NULL);
        reset_precision();
        return ret;
    }
    
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
    lst vec2lst(const exvector & ev) {
        return GiNaC::lst(ev.begin(), ev.end());
    }
    
    /**
     * @brief convert lst to exvector
     * @param alst input lst
     * @return exvector
     */
    exvector lst2vec(const lst & alst) {
        exvector ret(alst.begin(), alst.end());
        return exvector(std::move(ret));
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
        if(!expr.has(s0)) return (sn0>=0 ? expr : 0);
        
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
        
        symbol s("s");
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
    
    typedef vector<pair<ex,ex>> epvec_t; // first-pat, second-coeff
    typedef pair<ex,epvec_t> co_epvec_t; // first - overall coeff
    co_epvec_t power_expand_2(const co_epvec_t & co_epv, int n) {
        if(n==1) return co_epv;
        else if(n==2) {
            ex co = pow(co_epv.first,2);
            exmap pcmap;
            if(true) {
                auto epv = co_epv.second;
                if(epv.size()>10000) cout << "Warning: power_expand_2: epv.size > 10000" << endl;
                for(int i=0; i<epv.size(); i++) { 
                    pcmap[epv[i].first] += 2*co_epv.first*epv[i].second;
                    pcmap[pow(epv[i].first,2)] += pow(epv[i].second,2);
                    for(int j=i+1; j<epv.size(); j++) pcmap[epv[i].first*epv[j].first] += 2*epv[i].second*epv[j].second;
                }
            }
            epvec_t epv;
            epv.reserve(pcmap.size());
            for(auto kv : pcmap) {
                if(is_zero(kv.first) || is_zero(kv.second)) continue;
                epv.push_back(make_pair(kv.first, kv.second));
            }
            return make_pair(co, epvec_t(std::move(epv)));
        }
        int n2 = n/2;
        auto co_epv2 = power_expand_2(co_epv, n2);
        co_epv2 = power_expand_2(co_epv2, 2);
        if((n % 2)==0) return co_epvec_t(std::move(co_epv2));
        
        ex co = co_epv.first * co_epv2.first;
        exmap pcmap;
        for(auto ep1 : co_epv.second) {
            pcmap[ep1.first] += ep1.second * co_epv2.first;
            for(auto ep2 : co_epv2.second) pcmap[ep2.first * ep1.first] += ep1.second * ep2.second;
        }
        for(auto ep2 : co_epv2.second) pcmap[ep2.first] += ep2.second * co_epv.first;
        
        epvec_t epv;
        epv.reserve(pcmap.size());
        for(auto kv : pcmap) {
            if(is_zero(kv.first) || is_zero(kv.second)) continue;
            epv.push_back(make_pair(kv.first, kv.second));
        }
        return make_pair(co, epvec_t(std::move(epv)));
    }
    
    pair<ex,epvec_t> inner_expand_collect(ex const &expr_in, std::function<bool(const ex &)> has_func, int depth=0) {
        if(!has_func(expr_in)) {
            ex co;
            epvec_t epv;
            co = expr_in;
            return make_pair(co,epv);
        } if(is_a<add>(expr_in)) {
            ex co = 0;
            exmap pcmap;
            for(auto item : expr_in) {
                if(has_func(item)) {
                    auto co_epv = inner_expand_collect(item, has_func, depth+1);
                    co += co_epv.first;
                    for(auto ep : co_epv.second) pcmap[ep.first] += ep.second;
                } else co += item;
            }
            epvec_t epv;
            epv.reserve(pcmap.size());
            for(auto kv : pcmap) {
                if(is_zero(kv.first) || is_zero(kv.second)) continue;
                epv.push_back(make_pair(kv.first, kv.second));
            }
            return make_pair(co, epvec_t(std::move(epv)));
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
                    epv.reserve(pcmap.size());
                    for(auto kv : pcmap) {
                        if(is_zero(kv.first) || is_zero(kv.second)) continue;
                        epv.push_back(make_pair(kv.first, kv.second));
                    }
                }
            }
            return make_pair(co,epvec_t(std::move(epv)));
        } else if(is_a<power>(expr_in) && expr_in.op(1).info(info_flags::nonnegint)) {
            int n = ex_to<numeric>(expr_in.op(1)).to_int();
            auto co_epv = inner_expand_collect(expr_in.op(0), has_func, depth+1);
            return power_expand_2(co_epv,n);
        } else {
            ex co = 0;
            epvec_t epv;
            epv.push_back(make_pair(expr_in, 1));
            return make_pair(co,epvec_t(std::move(epv)));
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
            if(cc.is_zero() || vv.is_zero()) continue;
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
        if(!is_zero(cf)) co_epv.second.push_back(make_pair(1,cf));
        for(auto ep : co_epv.second) {
            ex vv = ep.first;
            ex cc = ep.second;
            if(opt==o_normal || opt==o_fermat || opt==o_fermatfD || opt==o_fermatN || opt==o_flint || opt==o_flintf || opt==o_flintfD) cc = exnormal(cc,opt);
            else if(opt==o_factor || opt==o_form) cc = exfactor(cc,opt);
            else if(opt==o_normal_fermat) cc = exnormal(normal(cc),o_fermat);
            else if(opt==o_normal_factor) cc = ginac_factor(normal(cc),o_fermat);
            else if(opt==o_normal_form) cc = form_factor(normal(cc),o_fermat);
            else if(opt==o_fermat_factor) cc = ginac_factor(fermat_normal(cc),o_fermat);
            else if(opt==o_fermat_form) cc = form_factor(fermat_normal(cc),o_fermat);
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
        set_precision(digits);
        auto nexpr = evalf(expr);
        reset_precision();
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
    //GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(XIntegral, basic,print_func<print_dflt>(&XIntegral::print))
    GiNaC::registered_class_info & XIntegral::get_class_info_static() { return reg_info; }
    XIntegral::visitor::~visitor() { }
    XIntegral * XIntegral::duplicate() const { XIntegral * bp = new XIntegral(*this); bp->setflag(GiNaC::status_flags::dynallocated); return bp; }
    void XIntegral::accept(GiNaC::visitor & v) const { if (visitor *p = dynamic_cast<visitor *>(&v)) p->visit(*this); else inherited::accept(v); }
    const GiNaC::registered_class_info &XIntegral::get_class_info() const { return get_class_info_static(); }
    GiNaC::registered_class_info &XIntegral::get_class_info() { return get_class_info_static(); }
    const char *XIntegral::class_name() const { return get_class_info_static().options.get_name(); }
    //GINAC_IMPLEMENT_REGISTERED_CLASS END
    
    DEFAULT_CTOR(XIntegral)
    IMPLEMENT_HAS(XIntegral)
    IMPLEMENT_ALL(XIntegral)

    int XIntegral::compare_same_type(const basic &other) const {
        if(!is_a<XIntegral>(other)) throw Error("XIntegral::compare_same_type");
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
    
    void XIntegral::read_archive(const archive_node& n) {
        inherited::read_archive(n);
        n.find_ex("Functions", Functions);
        n.find_ex("Exponents", Exponents);
        n.find_ex("Deltas", Deltas);
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
        if(opt==o_none) return expr;
        else if(opt==o_form) return factor_form(expr);
        else if(opt==o_factor) return ginac_factor(expr);
        else if(opt==o_flint || opt==o_flintf) return factor_flint(expr);
        else return expr;
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
    
    namespace {
        exvector ev_fermat_pern(const exvector & ev, int np) {
            if(ev.size()<np) return ev;
            int tn = ev.size();
            exvector evs;
            int ci = 0;
            while(ci<tn) {
                ex sum = 0;
                for(int i=0; i<np; i++) {
                    sum += ev[ci];
                    ci++;
                    if(ci==tn) break;
                }
                ex tt = normal_fermat(sum);
                evs.push_back(tt);
            }
            return exvector(std::move(evs));
        }
        
        ex normal_fermat_pern(const ex & e, int np) {
            if(!is_a<add>(e)) return normal_fermat(e);
            exvector evs(e.begin(), e.end());
            while(evs.size()>np) evs = ev_fermat_pern(evs,np);
            ex res = normal_fermat(add(evs));
            return res;
        }
    }
    
    /**
     * @brief normalize a expression
     * @param expr the input expression
     * @param opt 1 to use Fermat, otherwise using GiNaC for normalization
     * @return factorized result
     */
    ex exnormal(const ex & expr, int opt) {
        if(opt<0) return normal_fermat_pern(expr,-opt);
        else if(opt==o_normal) return normal(expr);
        else if(opt==o_fermat) return normal_fermat(expr);
        else if(opt==o_flint || opt==o_flintf || opt==o_flintfD) return normal_flint(expr, opt);
        else if(opt==o_fermatfD) return normal_fermat(expr,true);
        else if(opt==o_fermatN) return numer_fermat(expr);
        return expr;
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
            ss << "AutoDeclare Symbols fv;" << endl;
            form_map[pid].Init("form");
            form_map[pid].Execute(ss.str());
        }
        Form &fprc = form_map[pid];
        
        auto expr_in = expr;
        exmap map_rat;
        expr_in = expr_in.to_polynomial(map_rat);
        
        exset vs;
        for(const_preorder_iterator i = expr_in.preorder_begin(); i != expr_in.preorder_end(); ++i) {
            auto e = (*i);
            if(is_a<symbol>(e)) vs.insert(e);
        }
        
        exmap s2v;
        symtab st;
        int cid = 0;
        for(auto ss : vs) {
            cid++;
            string cvv = "fv"+to_string(cid);
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
        if(is_a<mul>(expr)) {
            ex res = 1;
            for(auto item : expr) res *= inner_factor_form(item);
            return res;
        }
        static map<pid_t, Form> form_map;
        auto pid = getpid();
        if((form_map.find(pid)==form_map.end())) { // init section
            ostringstream ss;
            ss << "AutoDeclare Symbols fv;" << endl;
            form_map[pid].Init("form");
            form_map[pid].Execute(ss.str());
        }
        Form &fprc = form_map[pid];
        
        auto expr_in = expr;
        exmap map_rat;
        expr_in = expr_in.to_polynomial(map_rat);
        
        exset vs;
        for(const_preorder_iterator i = expr_in.preorder_begin(); i != expr_in.preorder_end(); ++i) {
            auto e = (*i);
            if(is_a<symbol>(e)) vs.insert(e);
        }
        
        exmap s2v;
        symtab st;
        int cid = 0;
        for(auto ss : vs) {
            cid++;
            string cvv = "fv"+to_string(cid);
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
    
    void garWrite(const exvector &exv, string garfn) {
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
        
        map<string, ex> dict;
        for(int i=0; i<ar.num_expressions(); i++) {
            string name;
            ex res = ar.unarchive_ex(name, i);
            dict[name] = res;
        }
        
        auto size = ex_to<numeric>(dict["size"]).to_int();
        if(exv.size()>0 && size != exv.size()) throw Error("garRead: exvector size>0 & not match!");
        if(exv.size()<1) exv.resize(size);
        for(int i=0; i<size; i++) exv[i] = dict[to_string(i)];
    } 
    
    ex add_collect_normal(const exvector &exv, ex const &pats, int opt) {
        auto cvs_vec = GiNaC_Parallel(exv.size(), [&exv,pats](int idx)->ex {
            return collect_lst(exv[idx], pats);
        }, "ColEx");
        exmap res_map;
        for(auto cvs : cvs_vec) for(auto cv : cvs) res_map[cv.op(1)] += cv.op(0);
        exvector res_vec;
        for(auto kv : res_map) res_vec.push_back(lst{kv.second, kv.first});
        res_vec = GiNaC_Parallel(res_vec.size(), [&res_vec,opt](int idx)->ex {
            return exnormal(res_vec[idx].op(0), opt) * res_vec[idx].op(1);
        }, "NorEx");
        return add(res_vec);
    }
    
    ex add_collect_normal(const exvector &exv, init_list const &pats, int opt) {
        auto cvs_vec = GiNaC_Parallel(exv.size(), [&exv,pats](int idx)->ex {
            return collect_lst(exv[idx], pats);
        }, "ColEx");
        exmap res_map;
        for(auto cvs : cvs_vec) for(auto cv : cvs) res_map[cv.op(1)] += cv.op(0);
        exvector res_vec;
        for(auto kv : res_map) res_vec.push_back(lst{kv.second, kv.first});
        res_vec = GiNaC_Parallel(res_vec.size(), [&res_vec,opt](int idx)->ex {
            return exnormal(res_vec[idx].op(0),opt) * res_vec[idx].op(1);
        }, "NorEx");
        return add(res_vec);
    }
    
    ex add_collect_normal(const exvector &exv, lst const &pats, int opt) {
        auto cvs_vec = GiNaC_Parallel(exv.size(), [&exv,pats](int idx)->ex {
            return collect_lst(exv[idx], pats);
        }, "ColEx");
        exmap res_map;
        for(auto cvs : cvs_vec) for(auto cv : cvs) res_map[cv.op(1)] += cv.op(0);
        exvector res_vec;
        for(auto kv : res_map) res_vec.push_back(lst{kv.second, kv.first});
        res_vec = GiNaC_Parallel(res_vec.size(), [&res_vec,opt](int idx)->ex {
            return exnormal(res_vec[idx].op(0),opt) * res_vec[idx].op(1);
        }, "NorEx");
        return add(res_vec);
    }
    
    ex add_collect_normal(const ex & e, ex const &pats, int opt) {
        if(!is_a<add>(e)) throw Error("add_collect_normal: input is NOT a add class.");
        exvector exv(e.begin(), e.end());
        return add_collect_normal(exv, pats, opt);
    }
    
    ex add_collect_normal(const ex & e, lst const &pats, int opt) {
        if(!is_a<add>(e)) throw Error("add_collect_normal: input is NOT a add class.");
        exvector exv(e.begin(), e.end());
        return add_collect_normal(exv, pats, opt);
    }
    
    ex add_collect_normal(const ex & e, init_list const &pats, int opt) {
        if(!is_a<add>(e)) throw Error("add_collect_normal: input is NOT a add class.");
        exvector exv(e.begin(), e.end());
        return add_collect_normal(exv, pats, opt);
    }
        
    void ReShare(const ex & e) {
        archive ar;
        ar.archive_ex(e, "e");
    }
    
    void ReShare(const lst & es) {
        archive ar;
        for(auto const & e : es) ar.archive_ex(e, "e");
    }
    
    void ReShare(const ex & e1, const ex & e2) {
        archive ar;
        ar.archive_ex(e1, "e");
        ar.archive_ex(e2, "e");
    }
    
    void ReShare(const ex & e1, const ex & e2, const ex & e3) {
        archive ar;
        ar.archive_ex(e1, "e");
        ar.archive_ex(e2, "e");
        ar.archive_ex(e3, "e");
    }
    
    void ReShare(const exvector & ev) {
        archive ar;
        for(auto & e : ev) ar.archive_ex(e, "e");
    }
    
    void ReShare(const exvector & ev1, const exvector & ev2) {
        archive ar;
        for(auto & e : ev1) ar.archive_ex(e, "e");
        for(auto & e : ev2) ar.archive_ex(e, "e");
    }
    
    ex nextprime(const ex & n) {
        auto v = cln::the<cln::cl_I>(ex_to<numeric>(n).to_cl_N());
        return numeric(cln::nextprobprime(v));
    }
    
    ex Rationalize(const ex & e, int dn) {
        
        static MapFunction R([](const ex & e, MapFunction & self)->ex{
            if(is_a<numeric>(e)) {
                auto ne = ex_to<numeric>(e);
                if(ne.is_crational()) return e;
                auto zz = ne.to_cl_N();
                auto re = cln::rationalize(cln::realpart(zz));
                auto im = cln::rationalize(cln::imagpart(zz));
                return numeric(cln::complex(re,im));
            } else return e.map(self);
        });
        
        if(dn>0) set_precision(dn);
        ex res = R(e);
        if(dn>0) reset_precision();
        return res;
    }
    
    extern std::stack<cln::float_format_t> cln_prec_stack;
    extern std::stack<long> digits_stack;
    void set_precision(long prec, bool push) {
        if(push) {
            cln_prec_stack.push(cln::default_float_format);
            digits_stack.push(Digits);
        }
        Digits = prec;
        cln::default_float_format = cln::float_format(prec);
    }
    
    void reset_precision() {
        if(cln_prec_stack.empty()) return;
        auto digits = digits_stack.top();
        auto prec = cln_prec_stack.top();
        Digits = digits;
        cln::default_float_format = prec;
        cln_prec_stack.pop();
        digits_stack.pop();
    }
    
    long get_precision() {
        return cln::default_float_format;
    }
    
    void get_opt(int & argc, char** & argv, const string & o, map<char,string> & kv) {
        for (int opt; (opt = getopt(argc, argv, o.c_str())) != -1;) kv[opt] = optarg;
        argc -= optind;
        argv += optind;
    }
    
    bool has_symbol(const ex & e) {
        for(const_preorder_iterator i = e.preorder_begin(); i != e.preorder_end(); ++i) if(is_a<symbol>(*i)) return true;
        return false;
    }
        
}

