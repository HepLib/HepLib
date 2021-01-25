/**
 * @file
 * @brief Initailization for global objects
 */
 
#include "cln/cln.h"
#include "BASIC.h"
#include "HEP.h"
#include "QGRAF.h"
#include "IBP.h"
#include "SD.h"
#include <cstdlib>

namespace HepLib {

    // Initialization Unit
    // to make sure the initialization order is correct,
    // i.e., ordered in one compilation unit

    //----------------------------------------
    // HepLib
    //----------------------------------------
    std::map<std::string, ex> Symbol::Table; // alias as symtab in parser
    std::map<std::string, ex> iSymbol::Table; // alias as symtab in parser
    
    ex w = wild();
    ex w0 = wild(0);
    ex w1 = wild(1);
    ex w2 = wild(2);
    ex w3 = wild(3);
    ex w4 = wild(4);
    ex w5 = wild(5);
    ex w6 = wild(6);
    ex w7 = wild(7);
    ex w8 = wild(8);
    ex w9 = wild(9);
    
    const Symbol ep("ep");
    const iSymbol iEpsilon("iEpsilon");
    const ex iEpsilonN = I*pow(ex(10), -50);
    int Verbose = 0;
    int GiNaC_Parallel_Process = -1;
    const Symbol D("D");
    bool fermat_using_array = false;
    int NNDigits = 100;
    
    HepFormat::_init::_init() {
        set_print_func<add, HepFormat>(HepFormat::add_print);
        set_print_func<mul, HepFormat>(HepFormat::mul_print);
    }
    HepFormat::_init HepFormat_init;
    
    HepFormat hout(cout);
    
    MMAFormat::_init::_init() {
        set_print_func<add, MMAFormat>(MMAFormat::add_print);
        set_print_func<mul, MMAFormat>(MMAFormat::mul_print);
    }
    MMAFormat::_init MMAFormat_init;
    MMAFormat mout(cout);
    
    lst GiNaC_archive_Symbols = lst{};
    string InstallPrefix = "@CMAKE_INSTALL_PREFIX@";
    string INC_FLAGS = "@INC_FLAGS@";
    string LIB_FLAGS = "@LIB_FLAGS@";
    
    int Fermat::buffer_size = 1024*128;
    int Form::buffer_size = 1024*128;

    //----------------------------------------
    // HepLib::SD
    //----------------------------------------
    const Symbol SD::eps("eps");
    const Symbol SD::vs("vs");
    const Symbol SD::vz("vz");
    const Symbol SD::epz("epz");
    const Symbol SD::NaN("NaN");
    
    bool SD::SecDec::use_dlclose = true;
    bool SD::SecDec::debug = false;
    string SD::SecDec::cpp = "g++";

    SD::SecDec::_init::_init() {
        // for later use
        GiNaC_archive_Symbols.sort();
        GiNaC_archive_Symbols.unique();
    }
    SD::SecDec::_init SD::SecDec::SD_init;
    
    SD::CppFormat::_init::_init() {
        set_print_func<numeric, CppFormat>(CppFormat::print_numeric);
    }
    SD::CppFormat::_init SD::CppFormat::CppFormat_init;
    int SD::VEO_Digits = 10;
    
    //----------------------------------------
    // HepLib HEP
    //----------------------------------------
    const Symbol NA("NA");
    const Symbol NF("NF");
    const Symbol gs("gs");
    const Symbol as("as");
    const Symbol mu("mu");;
    const Symbol nL("nL");;
    const Symbol nH("nH");;
    
    exmap SP_map;
    map<ex,string,ex_is_less> QGRAF::LineTeX; // key is the filed
    map<ex,string,ex_is_less> QGRAF::VerTeX; // key is the fileds in vertex
    map<ex,string,ex_is_less> QGRAF::InOutTeX; // key is the id, id<0
    
    FCFormat::_init::_init() {
        set_print_func<ncmul, FCFormat>(FCFormat::ncmul_print);
    }
    FCFormat::_init FCFormat::FCFormat_init;
    
    FormFormat::_init::_init() {
        set_print_func<power, FormFormat>(FormFormat::power_print);
    }
    FormFormat::_init FormFormat::FormFormat_init;
    FCFormat fcout(cout);
    
    const int form_trace_auto = 0;
    const int form_trace_all = 1;
    const int form_trace_each_all = 2;
    const int form_trace_each_each = 3;
    int form_trace_mode = form_trace_auto;
    
    const int form_expand_none = 0;
    const int form_expand_tr = 1;
    const int form_expand_ci = 2;
    const int form_expand_li = 3;
    const int form_expand_all = 4;
    int form_expand_mode = form_expand_tr;
    
    bool Apart_using_fermat = true;
    bool form_using_su3 = false;
    bool form_using_dim4 = false;
    
    //----------------------------------------
    // Process _init
    //----------------------------------------
    QGRAF::Process::_init::_init() {
        auto A = Symbol("A");
        auto q = Symbol("q");
        auto qbar = Symbol("qbar");
        auto l = Symbol("l");
        auto lbar = Symbol("lbar");
        auto gh = Symbol("gh");
        auto ghbar = Symbol("ghbar");
        auto Q = Symbol("Q");
        auto g = Symbol("g");
        auto Qbar = Symbol("Qbar");
        auto n = Symbol("n");
        auto nbar = Symbol("nbar");
        auto e = Symbol("e");
        
        LineTeX[q] = "fermion, edge label=$q$";
        LineTeX[qbar] = "anti fermion, edge label=$q$";
        LineTeX[l] = "fermion, edge label=$l$";
        LineTeX[lbar] = "anti fermion, edge label=$l$";
        LineTeX[gh] = "ghost, edge label=$\\chi$"; 
        LineTeX[ghbar] = "ghost, edge label=$\\chi$"; 
        LineTeX[Q] = "fermion, edge label=$Q$";
        LineTeX[g] = "gluon, edge label=$g$";
        LineTeX[A] = "photon, edge label=$\\gamma$";
        LineTeX[Qbar] = "anti fermion, edge label=$Q$";
        LineTeX[n] = "double distance=1.5pt";
        LineTeX[nbar] = "double distance=1.5pt";    
        LineTeX[e] = "color=white";
        
        VerTeX[lst{Qbar, e, nbar}] = "[crossed dot]";
        VerTeX[lst{nbar, e, g}] = "[crossed dot]";
    }
    QGRAF::Process::_init QGRAF::Process::Process_init;
    
    //----------------------------------------
    // HepLib::IBP
    //----------------------------------------
    const Symbol IBP::d("d");
    int IBP::FIRE::Version = 6;
    int IBP::FIRE::Threads = 8;
    
    string IBP::UKIRA::KArgs = "";
    string IBP::KIRA::KArgs = "";
    
    //----------------------------------------
    // Rationalize
    //----------------------------------------
    
    MapFunction Rationalize([](const ex & e, MapFunction & self)->ex{
        if(is_a<numeric>(e)) {
            auto ne = ex_to<numeric>(e);
            if(ne.is_crational()) return e;
            auto zz = ne.to_cl_N();
            auto re = cln::rationalize(cln::realpart(zz));
            auto im = cln::rationalize(cln::imagpart(zz));
            return numeric(cln::complex(re,im));
        } else return e.map(self);
    });
    
    
    // global init class
    _global_init::_init::_init() {
        ostringstream oss;
        string path = InstallPrefix + "/bin";
        if(dir_exists(path)) oss << path;
        
        path = InstallPrefix + "/FIRE5/bin";
        if(dir_exists(path)) oss << ":" << path;
        
        path = InstallPrefix + "/FIRE6/bin";
        if(dir_exists(path)) oss << ":" << path;
        
        auto opath = getenv("PATH");
        if(opath != NULL) oss << ":" << opath;
        setenv("PATH", oss.str().c_str(), true);
        
    }
    _global_init::_init _global_init::init_object;
    
}

std::string HepLib::QGRAF::Process::Style = R"EOF(
<prologue>

<diagram>
((<sign><symmetry_factor>)*
<in_loop>InField(<field>,<field_index>,<momentum>)*
<end>
<back><out_loop>OutField(<field>,<field_index>,<momentum>)*
<end>
<back><propagator_loop>Propagator(
<back>Field(<field>,<field_index>),
<back>Field(<dual-field>,<dual-field_index>),
<back><momentum>)*
<end>
<back><vertex_loop>Vertex(
<back><ray_loop>Field(<field>,<field_index>,<momentum>),
<back><end><back>)*
<end><back><back>)
,
<epilogue>

<exit>

)EOF";
