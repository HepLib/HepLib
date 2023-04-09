/**
 * @file
 * @brief Initailization for global objects
 */
 
#include "cln/cln.h"
#include "BASIC.h"
#include "HEP.h"
#include "QGRAF.h"
#include "QCD.h"
#include "IBP.h"
#include "SD.h"
#include "DE.h"
#include <cstdlib>

namespace HepLib {

    // Initialization Unit
    // to make sure the initialization order is correct,
    // i.e., ordered in one compilation unit
    
    // FROM class
    
    // similar to GINAC_DECLARE_UNARCHIVER/GINAC_IMPLEMENT_REGISTERED_CLASS_OPT
    GiNaC::registered_class_info Symbol::reg_info = 
        GiNaC::registered_class_info(GiNaC::registered_class_options("Symbol", "symbol", typeid(Symbol)));
    
    GiNaC::registered_class_info iSymbol::reg_info = 
        GiNaC::registered_class_info(GiNaC::registered_class_options("iSymbol", "symbol", typeid(iSymbol)));
    
    GiNaC::registered_class_info XIntegral::reg_info = 
        GiNaC::registered_class_info(GiNaC::registered_class_options("XIntegral", "basic", typeid(XIntegral)).print_func<print_dflt>(&XIntegral::print));
    
    GiNaC::registered_class_info Index::reg_info = 
        GiNaC::registered_class_info(GiNaC::registered_class_options("Index", "basic", typeid(Index)).print_func<print_context>(&Index::print));
    
    GiNaC::registered_class_info Vector::reg_info = 
        GiNaC::registered_class_info(GiNaC::registered_class_options("Vector", "basic", typeid(Vector)).print_func<print_context>(&Vector::print));
    
    GiNaC::registered_class_info SUNT::reg_info = 
        GiNaC::registered_class_info(GiNaC::registered_class_options("SUNT", "basic", typeid(SUNT)).print_func<print_dflt>(&SUNT::print).print_func<FormFormat>(&SUNT::form_print).print_func<FCFormat>(&SUNT::fc_print));
            
    GiNaC::registered_class_info SUNF::reg_info = 
        GiNaC::registered_class_info(GiNaC::registered_class_options("SUNF", "basic", typeid(SUNF)).print_func<print_dflt>(&SUNF::print).print_func<FormFormat>(&SUNF::form_print).print_func<FCFormat>(&SUNF::fc_print));
    
    GiNaC::registered_class_info SUNF4::reg_info = 
        GiNaC::registered_class_info(GiNaC::registered_class_options("SUNF4", "basic", typeid(SUNF4)).print_func<print_dflt>(&SUNF4::print).print_func<FormFormat>(&SUNF4::form_print).print_func<FCFormat>(&SUNF4::fc_print));
    
    GiNaC::registered_class_info DGamma::reg_info = 
        GiNaC::registered_class_info(GiNaC::registered_class_options("DGamma", "basic", typeid(DGamma)).print_func<print_dflt>(&DGamma::print).print_func<FormFormat>(&DGamma::form_print).print_func<FCFormat>(&DGamma::fc_print));
    
    GiNaC::registered_class_info Eps::reg_info = 
        GiNaC::registered_class_info(GiNaC::registered_class_options("Eps", "basic", typeid(Eps)).print_func<print_dflt>(&Eps::print).print_func<FormFormat>(&Eps::form_print).print_func<FCFormat>(&Eps::fc_print));
    
    GiNaC::registered_class_info Pair::reg_info = 
        GiNaC::registered_class_info(GiNaC::registered_class_options("Pair", "basic", typeid(Pair)).print_func<print_dflt>(&Pair::print).print_func<FormFormat>(&Pair::form_print).print_func<FCFormat>(&Pair::fc_print));
    
    // FROM BASIC

    string Version = "1.1 - 2022-02-14";
    exmap Symbol::vmap;
    std::map<std::string, ex> Symbol::Table; // alias as symtab in parser
    std::map<std::string, ex> iSymbol::Table; // alias as symtab in parser
    
    unsigned nopat = GiNaC::subs_options::no_pattern;
    
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
    
    // normal/factor options
    int _o_ = 0;
    const int o_none = (_o_++);
    const int o_normal = (_o_++);
    const int o_fermat = (_o_++);
    const int o_fermatfD = (_o_++);
    const int o_fermatN = (_o_++);
    const int o_form = (_o_++);
    const int o_flint = (_o_++);
    const int o_flintf = (_o_++);
    const int o_flintfD = (_o_++);
    const int o_normal_fermat = (_o_++);
    const int o_normal_form = (_o_++);
    const int o_fermat_form = (_o_++);
    
    // Symbols
    const Symbol NA("NA");
    const Symbol NF("NF");
    const Symbol TF("TF");
    const Symbol CA("CA");
    const Symbol CF("CF");
    const Symbol gs("gs");
    const Symbol as("as");
    const Symbol mu("mu");
    const Symbol nL("nL");
    const Symbol nH("nH");
    const Symbol eps("eps");
    const Symbol vs("vs");
    const Symbol vz("vz");
    const Symbol epz("epz");
    const Symbol NaN("NaN");
    const Symbol ep("ep");
    const Symbol d("d");
    const Symbol iet("iet");
    const iSymbol iEpsilon("iEpsilon");
    
    const ex iEpsilonN = I*pow(ex(10), -50);
    int Verbose = 0;
    bool Debug = false;
    bool In_GiNaC_Parallel = false;
    int GiNaC_Parallel_Process = -1;
    map<string, int> GiNaC_Parallel_NP;
    map<string, int> GiNaC_Parallel_Verb;
    int GiNaC_Parallel_Batch = 0;
    map<string, int> GiNaC_Parallel_NB;
    map<string, bool> GiNaC_Parallel_RM;
    map<string, string> GiNaC_Parallel_PRE;
    map<string, bool> GiNaC_Parallel_ReWR;
    int fermat_using_array = 0;
    map<ex,long long,ex_is_less> fermat_weight;
    bool using_cache = true;
    long long cache_limit = -1;
    int NNDigits = 100;
    
    MMAFormat mout(cout);
    
    string InstallPrefix = "@CMAKE_INSTALL_PREFIX@";
    string INC_FLAGS = "@INC_FLAGS@";
    string LIB_FLAGS = "@LIB_FLAGS@";
    
    int Fermat::buffer_size = 1024*128;
    int Form::buffer_size = 1024*128;

    bool SD::SecDec::use_dlclose = true;
    string SD::SecDec::cpp = "g++";
    
    SD::CppFormat::_init::_init() {
        set_print_func<numeric, CppFormat>(CppFormat::print_numeric);
    }
    SD::CppFormat::_init SD::CppFormat::CppFormat_init;
    int SD::VEO_Digits = 10;
    
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
    bool form_using_su3 = true;
    bool form_using_dim4 = false;
    
    QGRAF::Process::_init::_init() {
        auto A = Symbol("A");
        auto q = Symbol("q");
        auto qbar = Symbol("qbar");
        auto l = Symbol("l");
        auto lbar = Symbol("lbar");
        auto gh = Symbol("gh");
        auto ghbar = Symbol("ghbar");
        auto g = Symbol("g");
        auto Q = Symbol("Q");
        auto Qbar = Symbol("Qbar");
        auto C = Symbol("C");
        auto Cbar = Symbol("Cbar");
        auto B = Symbol("B");
        auto Bbar = Symbol("Bbar");
        auto n = Symbol("n");
        auto nbar = Symbol("nbar");
        auto e = Symbol("e");
        
        LineTeX[q] = "fermion, edge label=$q$";
        LineTeX[qbar] = "anti fermion, edge label=$q$";
        LineTeX[l] = "fermion, edge label=$l$";
        LineTeX[lbar] = "anti fermion, edge label=$l$";
        LineTeX[gh] = "ghost, edge label=$\\chi$"; 
        LineTeX[ghbar] = "ghost, edge label=$\\chi$"; 
        LineTeX[g] = "gluon, edge label=$g$";
        LineTeX[A] = "photon, edge label=$\\gamma$";
        LineTeX[Q] = "fermion, edge label=$Q$";
        LineTeX[Qbar] = "anti fermion, edge label=$Q$";
        LineTeX[C] = "fermion, edge label=$c$";
        LineTeX[Cbar] = "anti fermion, edge label=$c$";
        LineTeX[B] = "fermion, edge label=$b$";
        LineTeX[Bbar] = "anti fermion, edge label=$b$";
        LineTeX[n] = "double distance=1.5pt";
        LineTeX[nbar] = "double distance=1.5pt";    
        LineTeX[e] = "color=white";
        
        VerTeX[lst{Qbar, e, nbar}] = "[crossed dot]";
        VerTeX[lst{nbar, e, g}] = "[crossed dot]";
    }
    QGRAF::Process::_init QGRAF::Process::Process_init;
    
    // FROM IBP
    
    int FIRE::PosPref = 1;
    int FIRE::Threads = 4;
    int FIRE::fThreads = 1;
    int FIRE::lThreads = 0;
    int FIRE::sThreads = 0;
    int FIRE::Version = 6;
    string FIRE::Execute;
    exmap FIRE::NVariables;
    exmap MapPreSP;
    
    string UKIRA::KArgs = "";
    string KIRA::KArgs = "";
    
    // FROM QCD
    int QCD::FF::cur_mode = 0; // 0 - gluon, 1 - quark, 2 - anti-quark
    
    // FROM DE
    slong error_pass_dp = 100;
    
    // global init
    std::stack<cln::float_format_t> cln_prec_stack;
    std::stack<long> digits_stack;
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
        
        // similar to GINAC_DECLARE_UNARCHIVER/GINAC_BIND_UNARCHIVER
        GiNaC::unarchive_table_t table;
        table.insert(std::string("Symbol"), []()->GiNaC::basic*{ return new Symbol(); });
        table.insert(std::string("iSymbol"), []()->GiNaC::basic*{ return new iSymbol(); });
        table.insert(std::string("XIntegral"), []()->GiNaC::basic*{ return new XIntegral(); });
        table.insert(std::string("Index"), []()->GiNaC::basic*{ return new Index(); });
        table.insert(std::string("Vector"), []()->GiNaC::basic*{ return new Vector(); });
        table.insert(std::string("SUNT"), []()->GiNaC::basic*{ return new SUNT(); });
        table.insert(std::string("SUNF"), []()->GiNaC::basic*{ return new SUNF(); });
        table.insert(std::string("SUNF4"), []()->GiNaC::basic*{ return new SUNF4(); });
        table.insert(std::string("DGamma"), []()->GiNaC::basic*{ return new DGamma(); });
        table.insert(std::string("Eps"), []()->GiNaC::basic*{ return new Eps(); });
        table.insert(std::string("Pair"), []()->GiNaC::basic*{ return new Pair(); });
        
        // CLN configurations
        cln::cl_inhibit_floating_point_underflow = true; 
        Digits = 100;
        set_precision(100);
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

std::string HepLib::QCD::FF::GluonModel = R"EOF(
[ model = 'Gluon FF Model' ]
[q, qbar, -]
[Q, Qbar, -]
[gh, ghbar, -]
[g, g, +, notadpole]
[e, e, +, external]
[n, nbar, +]
[qbar, q, g; QCD='+1']
[Qbar, Q, g; QCD='+1']
[g, g, g, g; QCD='+2']
[g, g, g; QCD='+1']
[ghbar, gh, g; QCD='+1']
[nbar, n, g; QCD='+1']
[nbar, e, g; QCD='+0']
%[n, e, g; QCD='+0']
)EOF";

std::string HepLib::QCD::FF::QuarkModel = R"EOF(
[ model = 'Quark FF Model' ]
[q, qbar, -]
[Q, Qbar, -]
[e, ebar, -, external]
[gh, ghbar, -]
[g, g, +, notadpole]
[n, nbar, +]
[qbar, q, g; QCD='+1']
[Qbar, Q, g; QCD='+1']
[g, g, g, g; QCD='+2']
[g, g, g; QCD='+1']
[ghbar, gh, g; QCD='+1']
[nbar, n, g; QCD='+1']
[qbar, e, nbar; QCD='0']
[ebar, q, nbar; QCD='0']
[Qbar, e, nbar; QCD='0']
[ebar, Q, nbar; QCD='0']
)EOF";
