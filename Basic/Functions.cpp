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
    Parser::Parser() : FuncDict(ginac_reader()) { }
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
    
}

