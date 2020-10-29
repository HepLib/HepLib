/**
 * @file
 * @brief Initailization for global objects
 */
 
#include "cln/cln.h"
#include "SD.h"
#include "FC.h"
#include "IBP.h"
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
    
    const Symbol ep("ep");
    const iSymbol iEpsilon("iEpsilon");
    const ex iEpsilonN = I*pow(ex(10), -50);
    int Verbose = 0;
    int ParallelProcess = -1;
    const Symbol D("D");
    bool fermat_use_array = true;
    int NNDigits = 100;
    
    lst GiNaC_archive_Symbols = lst{};
    string InstallPrefix = "@CMAKE_INSTALL_PREFIX@";
    string INC_FLAGS = "@INC_FLAGS@";
    string LIB_FLAGS = "@LIB_FLAGS@";

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
    // HepLib::FC
    //----------------------------------------
    const Symbol FC::NA("NA");
    const Symbol FC::NF("NF");
    const Symbol FC::gs("gs");
    const Symbol FC::as("as");
    const Symbol FC::mu("mu");;
    const Symbol FC::nL("nL");;
    const Symbol FC::nH("nH");;
        
    exmap FC::SP_map;
    map<ex,string,ex_is_less> FC::Qgraf::LineTeX; // key is the filed
    map<ex,string,ex_is_less> FC::Qgraf::VerTeX; // key is the fileds in vertex
    map<ex,string,ex_is_less> FC::Qgraf::InOutTeX; // key is the id, id<0
    
    FC::FCFormat::_init::_init() {
        set_print_func<ncmul, FCFormat>(FCFormat::ncmul_print);
    }
    FC::FCFormat::_init FC::FCFormat::FCFormat_init;
    
    FC::FormFormat::_init::_init() {
        set_print_func<power, FormFormat>(FormFormat::power_print);
    }
    FC::FormFormat::_init FC::FormFormat::FormFormat_init;
    FC::FCFormat FC::FCout(cout);
    const int FC::trace_all = 0;
    const int FC::trace_each_all = 1;
    const int FC::trace_each_each = 2;
    int FC::trace_method = FC::trace_each_each;
    
    //----------------------------------------
    // Process _init
    //----------------------------------------
    FC::Qgraf::Process::_init::_init() {
        auto q = Symbol("q");
        auto qbar = Symbol("qbar");
        auto gh = Symbol("gh");
        auto ghbar = Symbol("ghbar");
        auto Q = Symbol("Q");
        auto g = Symbol("g");
        auto Qbar = Symbol("Qbar");
        auto n = Symbol("n");
        auto nbar = Symbol("nbar");
        auto e = Symbol("e");
        
        LineTeX[q] = "fermion, edge label=q";
        LineTeX[qbar] = "anti fermion, edge label=q";
        LineTeX[gh] = "ghost, edge label=$\\chi$"; 
        LineTeX[ghbar] = "ghost, edge label=$\\chi$"; 
        LineTeX[Q] = "fermion, edge label=Q";
        LineTeX[g] = "gluon, edge label=g";
        LineTeX[Qbar] = "anti fermion, edge label=$\\bar{Q}$";
        LineTeX[n] = "double distance=1.5pt";
        LineTeX[nbar] = "double distance=1.5pt";    
        LineTeX[e] = "color=white";
        
        VerTeX[lst{Qbar, e, nbar}] = "[crossed dot]";
        VerTeX[lst{nbar, e, g}] = "[crossed dot]";
    }
    FC::Qgraf::Process::_init FC::Qgraf::Process::Process_init;
    
    //----------------------------------------
    // HepLib::IBP
    //----------------------------------------
    const Symbol IBP::d("d");
    int IBP::FIRE::Version = 6;
    int IBP::FIRE::Threads = 8;
    
    int IBP::UKIRA::Rounds = 1;
    
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

const std::string HepLib::FC::Qgraf::Style = R"EOF(
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