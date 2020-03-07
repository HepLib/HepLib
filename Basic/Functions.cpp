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

}

