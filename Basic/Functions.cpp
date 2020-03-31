#include "Basic.h"

namespace HepLib {

    namespace {
        static ex coCF_Diff(const ex & x, unsigned diff_param) {return 0;}
        static ex coCF_Expand(const ex & x, unsigned expand_options) {return coCF(x).hold();}
    }

    REGISTER_FUNCTION(coCF, derivative_func(coCF_Diff).expand_func(coCF_Expand))
    REGISTER_FUNCTION(coVF, do_not_evalf_params())

    REGISTER_FUNCTION(x, do_not_evalf_params())
    REGISTER_FUNCTION(y, do_not_evalf_params())
    REGISTER_FUNCTION(z, do_not_evalf_params())
    
    // WF function, up to 5 arguments
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
    
    /*-----------------------------------------------------*/
    // MapFunction Class
    /*-----------------------------------------------------*/
    MapFunction::MapFunction(std::function<ex(const ex &, MapFunction &)> func) : Function(func) { }
    ex MapFunction::operator()(const ex &e) {
        return Function(e, *this);
    }
    
    /*-----------------------------------------------------*/
    // Parser Class
    /*-----------------------------------------------------*/
    // copy from GiNaC parser, note the alignas(2)
    namespace {
        alignas(2) static ex sqrt_reader(const exvector& ev) {
            return GiNaC::sqrt(ev[0]);
        }

        alignas(2) static ex pow_reader(const exvector& ev) {
            return GiNaC::pow(ev[0], ev[1]);
        }

        alignas(2) static ex power_reader(const exvector& ev) {
            return GiNaC::power(ev[0], ev[1]);
        }

        alignas(2) static ex lst_reader(const exvector& ev) {
            return GiNaC::lst(ev.begin(), ev.end());
        }
    
        class functions_hack : public GiNaC::function {
        public:
            static const std::vector<function_options>& get_registered_functions() {
                return function::registered_functions();
            }
        };
        
        static reader_func encode_serial_as_reader_func(unsigned serial) {
            uintptr_t u = (uintptr_t)serial;
            u = (u << 1) | (uintptr_t)1;
            reader_func ptr = (reader_func)((void *)u);
            return ptr;
        }
        
        const prototype_table& ginac_reader() {
            using std::make_pair;
            static bool initialized = false;
            static prototype_table reader;
            if (!initialized) {
                reader[make_pair("sqrt", 1)] = sqrt_reader;
                reader[make_pair("pow", 2)] = pow_reader;
                reader[make_pair("power", 2)] = power_reader;
                reader[make_pair("lst", 0)] = lst_reader;
                enum { log, exp, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh, atan2,
                        Li2, Li3, zetaderiv, Li, S, H, lgamma, tgamma, beta, factorial, binomial, Order,
                        NFUNCTIONS
                };
                auto it = functions_hack::get_registered_functions().begin();
                unsigned serial = 0;
                for ( ; serial<NFUNCTIONS; ++it, ++serial ) {
                    prototype proto = make_pair(it->get_name(), it->get_nparams());
                    reader[proto] = encode_serial_as_reader_func(serial);
                }
                initialized = true;
            }
            return reader;
        }
    }
    
    Parser::Parser(symtab st) : SymDict(st), FuncDict(ginac_reader()) { }
    ex Parser::Read(string instr) {
        unsigned serial = 0;
        for (auto & it : functions_hack::get_registered_functions()) {
            prototype proto = make_pair(it.get_name(), it.get_nparams());
            FuncDict[proto] = encode_serial_as_reader_func(serial);
            ++serial;
        }
        parser ginac_parser(SymDict, false, FuncDict);
        return ginac_parser(instr);
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

            ~Generator() { delete a; }

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
}

