/**
 * @file
 * @brief Basic/Useful Functions
 */
 
#include "BASIC.h"
#include "ginac/parse_context.h"

namespace HepLib {
    
    REGISTER_FUNCTION(x, do_not_evalf_params())
    REGISTER_FUNCTION(y, do_not_evalf_params())
    REGISTER_FUNCTION(z, do_not_evalf_params())
    
    #ifndef DOXYGEN_SKIP
    
    unsigned F1_SERIAL::serial = GiNaC::function::register_new(function_options("F",1).do_not_evalf_params().overloaded(2));
    unsigned F2_SERIAL::serial = GiNaC::function::register_new(function_options("F",2).do_not_evalf_params().overloaded(2));
        
    unsigned WF1_SERIAL::serial = GiNaC::function::register_new(function_options("WF",1).do_not_evalf_params().overloaded(5));
    unsigned WF2_SERIAL::serial = GiNaC::function::register_new(function_options("WF",2).do_not_evalf_params().overloaded(5));
    unsigned WF3_SERIAL::serial = GiNaC::function::register_new(function_options("WF",3).do_not_evalf_params().overloaded(5));
    unsigned WF4_SERIAL::serial = GiNaC::function::register_new(function_options("WF",4).do_not_evalf_params().overloaded(5));
    unsigned WF5_SERIAL::serial = GiNaC::function::register_new(function_options("WF",5).do_not_evalf_params().overloaded(5));
    // iWF function, up to 5 arguments
    unsigned iWF1_SERIAL::serial = GiNaC::function::register_new(function_options("iWF",1).do_not_evalf_params().overloaded(5));
    unsigned iWF2_SERIAL::serial = GiNaC::function::register_new(function_options("iWF",2).do_not_evalf_params().overloaded(5));
    unsigned iWF3_SERIAL::serial = GiNaC::function::register_new(function_options("iWF",3).do_not_evalf_params().overloaded(5));
    unsigned iWF4_SERIAL::serial = GiNaC::function::register_new(function_options("iWF",4).do_not_evalf_params().overloaded(5));
    unsigned iWF5_SERIAL::serial = GiNaC::function::register_new(function_options("iWF",5).do_not_evalf_params().overloaded(5));
    
    #endif
    
    /*-----------------------------------------------------*/
    // MapFunction Class
    /*-----------------------------------------------------*/
    MapFunction::MapFunction(std::function<ex(const ex &, MapFunction &)> func) : Function(func) { }
    ex MapFunction::operator()(const ex &e) {
        return Function(e, *this);
    }
    
    ex MapFunction::subs(const ex & expr, const ex & pat, std::function<ex(const ex &)>f) {
        MapFunction map([pat,f](const ex & e, MapFunction &self)->ex{
            if(!e.has(pat)) return e;
            else if(e.match(pat)) return f(e);
            else return e.map(self);
        });
        return map(expr);
    }
    
    /*-----------------------------------------------------*/
    // Parser Class
    /*-----------------------------------------------------*/
    // copy from GiNaC parser
    namespace {
        class functions_hack : public GiNaC::function {
        public:
            static const std::vector<function_options>& get_registered_functions() {
                return function::registered_functions();
            }
        };
        
        static ex sqrt_reader(const exvector& ev) {
            return GiNaC::sqrt(ev[0]);
        }

        static ex pow_reader(const exvector& ev) {
            return GiNaC::pow(ev[0], ev[1]);
        }

        static ex power_reader(const exvector& ev) {
            return GiNaC::power(ev[0], ev[1]);
        }

        static ex lst_reader(const exvector& ev) {
            return GiNaC::lst(ev.begin(), ev.end());
        }
        
        const prototype_table& ginac_reader() {
            using std::make_pair;
            static bool initialized = false;
            static prototype_table reader;
            if (!initialized) {
                reader.insert({{"sqrt", 1}, reader_func(sqrt_reader)});
                reader.insert({{"pow", 2}, reader_func(pow_reader)});
                reader.insert({{"power", 2}, reader_func(power_reader)});
                reader.insert({{"lst", 0}, reader_func(lst_reader)});
                enum { log, exp, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh, atan2,
                        Li2, Li3, zetaderiv, Li, S, H, lgamma, tgamma, beta, factorial, binomial, Order,
                        NFUNCTIONS
                };
                auto it = functions_hack::get_registered_functions().begin();
                unsigned serial = 0;
                for ( ; serial<NFUNCTIONS; ++it, ++serial ) {
                    reader.insert({{it->get_name(), it->get_nparams()}, reader_func(serial)});
                }
                initialized = true;
            }
            return reader;
        }
    }
    
    // note iSymbol will override Symbol if they share same name
    Parser::Parser() : STable(Symbol::Table), FTable(ginac_reader()) { 
        for(auto item : iSymbol::Table) STable[item.first] = item.second;
    }
    Parser::Parser(symtab st) : STable(Symbol::Table), FTable(ginac_reader()) { 
        for(auto item : iSymbol::Table) STable[item.first] = item.second;
        for(auto item : st) STable[item.first] = item.second;
    }
    ex Parser::Read(const string & in_str, bool s2S) {
        string instr(in_str);
        string_replace_all(instr, "==", "=");
        unsigned serial = 0;
        for (auto & it : functions_hack::get_registered_functions()) {
            FTable.insert({{it.get_name(), it.get_nparams()}, reader_func(serial)});
            ++serial;
        }
        parser ginac_parser(STable, false, FTable);
        ex ret = ginac_parser(instr);
        if(!s2S) return ret;
        // check & replace symbol to Symbol object
        auto st = ginac_parser.get_syms();
        bool redo = false;
        exmap repl;
        for(auto kv : st) {
            if(is_exactly_a<symbol>(kv.second) && STable.find(kv.first)==STable.end()) {
                string ss = kv.first;
                STable[ss] = Symbol(ss);
                redo = true;
                repl[kv.second] = Symbol(ss);
            }
        }
        if(redo) ret = ret.subs(repl,subs_options::no_pattern);
        //ReShare(ret);
        return ret;
    }
    ex Parser::ReadFile(string filename, bool s2S) {
        ifstream ifs(filename);
        string ostr((istreambuf_iterator<char>(ifs)), (istreambuf_iterator<char>()));
        ifs.close();
        return Read(ostr,s2S);
    }
    
    /*-----------------------------------------------------*/
    // string Functions
    /*-----------------------------------------------------*/
    void string_replace_all(string &str, const string &from, const string &to) {
        size_t start_pos = 0;
        while((start_pos = str.find(from, start_pos)) != string::npos) {
            str.replace(start_pos, from.length(), to);
            start_pos += to.length();
        }
    }
    
    void string_trim(string &ostr) {
        const char* WhiteSpace = " \t\v\r\n";
        if(!ostr.empty()) {
            ostr.erase(0, ostr.find_first_not_of(WhiteSpace));
            ostr.erase(ostr.find_last_not_of(WhiteSpace)+1);
        }
    }
    
    bool string_start_with(const string &fstr, const string & sstr) {
        if (fstr.length() >= sstr.length()) {
            return (0 == fstr.compare (0, sstr.length(), sstr));
        } else {
            return false;
        }
    }
    
    bool string_end_with(const string &fstr, const string & estr) {
        if (fstr.length() >= estr.length()) {
            return (0 == fstr.compare (fstr.length() - estr.length(), estr.length(), estr));
        } else {
            return false;
        }
    }
    
    bool string_contain(const string &fstr, const string & mstr) {
        return (fstr.find(mstr) != std::string::npos);
    }
    
    void Combinations(int n, int m, std::function<void(const int*)> f) {
        if(m<1 || m>n) return;
        std::string bitmask(m,1);
        bitmask.resize(n,0);
        do {
            int is[m]; int j=0;
            for (int i=0; i<n; ++i) {
                if(bitmask[i]) { is[j]=i; j++; }
            }
            f(is);
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    }
    
    void CombinationsR(int n, int m, std::function<void(const int*)> f) {
        if(m<1) return;
        Combinations(n+m-1, n-1, [n,m,f](const int* isr) {
            int is[m];
            int mi=m;
            for(int ni=0; ni<n; ni++) {
                int c=0;
                if(ni==n-1) c=(n+m-2)-isr[n-2];
                else if(ni==0) c=isr[ni];
                else c=isr[ni]-isr[ni-1]-1;
                for(int j=0; j<c; j++) is[--mi] = n-1-ni;
            }
            f(is); 
        });
    }
    
    void Permutations(int n, std::function<void(const int*)> f) {
        if(n<1) return;
        int pis[n];
        for(int i=0; i<n; i++) pis[i]=i;
        do { f(pis); } while(std::next_permutation(pis,pis+n));
    }
    
    void Permutations(int n, int m, std::function<void(const int*)> f) {
        if(m<1 || m>n) return;
        Combinations(n, m, [m,f](const int *ns1) {
            Permutations(m, [m,f,ns1](const int *ns2) {
                int ns[m];
                for(int i=0; i<m; i++) ns[i] = ns1[ns2[i]];
                f(ns);
            });
        });
    }
    
    namespace {
        struct Generator {
        public:
            int *a;
            
            Generator(int s, int v) : cSlots(s) , cValues(v) {
                a = new int[s];
                for(int i = 0; i<cSlots-1; i++) a[i]=1-1;
                a[cSlots-1]=0-1;
                nextInd = cSlots;
            }

            ~Generator() { delete [] a; }

            bool doNext() {
                for (;;) {
                    if (a[nextInd-1]==cValues-1) {
                        nextInd--;
                        if(nextInd==0) return false;
                    } else {
                        a[nextInd-1]++;
                        while (nextInd<cSlots) {
                            nextInd++;
                            a[nextInd-1]=1-1;
                        }
                        return true;
                    }
                }
            }
            
        private:
            int cSlots;
            int cValues;
            int nextInd;
        };
    }   
    
    void PermutationsR(int n, int m, std::function<void(const int*)> f) {
        Generator g(m,n);
        while (g.doNext()) f(g.a);
    }
    
    bool isSorted(const lst & exs) {
        for(int i=0; i<exs.nops()-1; i++) {
            if(exs.op(i).is_equal(exs.op(i+1))) continue;
            if(!ex_less(exs.op(i),exs.op(i+1))) return false;
        }
        return true;
    }
    
    bool isSorted(int n, const ex exs[]) {
        for(int i=0; i<n-1; i++) {
            if(exs[i].is_equal(exs[i+1])) continue;
            if(!ex_less(exs[i],exs[i+1])) return false;
        }
        return true;
    }
    
    int ACSort(lst & exs) {
        int ac = 0;
        int n = exs.nops();
        for(int i=0; i<n-1; i++) 
        for(int j=n-1; j>i; j--) 
        if(ex_less(exs.op(j),exs.op(j-1))) {
            auto tmp = exs.op(j-1);
            exs[j-1] = exs.op(j);
            exs[j] = tmp;
            ac++;
        }
        for(int i=0; i<n-1; i++) if(exs.op(i).is_equal(exs.op(i+1))) return 0;
        return (ac%2==1) ? -1 : 1;
    }
    
    int ACSort(int n, ex exs[]) {
        int ac = 0;
        for(int i=0; i<n-1; i++) 
        for(int j=n-1; j>i; j--) 
        if(ex_less(exs[j],exs[j-1])) {
            auto tmp = exs[j-1];
            exs[j-1] = exs[j];
            exs[j] = tmp;
            ac++;
        }
        for(int i=0; i<n-1; i++) if(exs[i].is_equal(exs[i+1])) return 0;
        return (ac%2==1) ? -1 : 1;
    }
    
}

