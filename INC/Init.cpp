/**
 * @file
 * @brief Initailization for global objects
 */

#include "BASIC.h"
#include "HEP.h"
#include "QGRAF.h"
#include "QCD.h"
#include "IBP.h"
#include "SD.h"
#include "DE.h"
#include "AMF.h"
#include <cstdlib>

#include "cln/cln.h"
#include <dlfcn.h>
#include <string.h>
#include <string>

namespace HepLib {


    namespace {
        string install_prefix() {
            static string wdir = "";
            if(wdir == "") {
                string path;
                Dl_info dl_info;
                dladdr((void*)install_prefix, &dl_info);
                path = dl_info.dli_fname;
                wdir = path.substr(0, path.find_last_of('/'));
                path = wdir;
                wdir = path.substr(0, path.find_last_of('/'));
            }
            return wdir;
        }
    }

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
        
    GiNaC::registered_class_info AsGamma::reg_info =
        GiNaC::registered_class_info(GiNaC::registered_class_options("AsGamma", "basic", typeid(AsGamma)).print_func<print_dflt>(&AsGamma::print).print_func<FormFormat>(&AsGamma::form_print).print_func<FCFormat>(&AsGamma::fc_print));
    
    GiNaC::registered_class_info Eps::reg_info = 
        GiNaC::registered_class_info(GiNaC::registered_class_options("Eps", "basic", typeid(Eps)).print_func<print_dflt>(&Eps::print).print_func<FormFormat>(&Eps::form_print).print_func<FCFormat>(&Eps::fc_print));
    
    GiNaC::registered_class_info Pair::reg_info = 
        GiNaC::registered_class_info(GiNaC::registered_class_options("Pair", "basic", typeid(Pair)).print_func<print_dflt>(&Pair::print).print_func<FormFormat>(&Pair::form_print).print_func<FCFormat>(&Pair::fc_print));
    
    // FROM BASIC

    string Version = "1.5 @2024-04-28";
    exmap Symbol::vmap;
    std::map<std::string, ex> Symbol::Table; // alias as symtab in parser
    std::map<std::string, ex> iSymbol::Table; // alias as symtab in parser
    GiNaC::exmap Index::Dimension; // default dimension is d, others specified here.
    
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
    string PRE = "  ";
    bool Debug = false;
    bool In_GiNaC_Parallel = false;
    int GiNaC_Parallel_Level = 0;
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
    
    bool GMat_using_cache = true;
    
    MMAFormat mout(cout);
    
    string InstallPrefix = install_prefix();
    string INC_FLAGS = "-I'"+InstallPrefix+"/include' " + "@INC_FLAGS@";
    string LIB_FLAGS = "-L'"+InstallPrefix+"/lib' -Wl,-rpath,'"+InstallPrefix+"'/lib " + "@LIB_FLAGS@";
    
    int Fermat::buffer_size = 1024*128;
    int Form::buffer_size = 1024*128;

    bool SD::SecDec::use_dlclose = true;
    string SD::SecDec::cpp = "g++ -w";
    
    SD::CppFormat::_init::_init() {
        set_print_func<numeric, CppFormat>(CppFormat::print_numeric);
    }
    SD::CppFormat::_init SD::CppFormat::CppFormat_init;
    SD::ExFormat::_init::_init() {
        set_print_func<numeric, ExFormat>(ExFormat::print_numeric);
    }
    SD::ExFormat::_init SD::ExFormat::ExFormat_init;
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
    bool form_using_su3 = false;
    bool form_using_dim4 = false;
    bool form_using_gamma5 = false;
    
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
        auto G = Symbol("G"); // graviton
        
        LineTeX[q] = "fermion, edge label=$q$";
        LineTeX[qbar] = "anti fermion, edge label=$q$";
        LineTeX[l] = "fermion, edge label=$l$";
        LineTeX[lbar] = "anti fermion, edge label=$l$";
        LineTeX[gh] = "ghost, edge label=$\\chi$"; 
        LineTeX[ghbar] = "ghost, edge label=$\\chi$"; 
        LineTeX[g] = "gluon, edge label=$g$";
        LineTeX[A] = "photon, edge label=$\\gamma$";
        LineTeX[G] = "graviton, edge label=$\\gamma$";
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
        
        {
            Symbol U("U"); // U-quark
            Symbol Ubar("Ubar"); // anti U-quark
            Symbol D("D"); // D-quark
            Symbol Dbar("Dbar"); // anti D-quark
            Symbol C("C"); // C-quark
            Symbol Cbar("Cbar"); // anti C-quark
            Symbol S("S"); // S-quark
            Symbol Sbar("Sbar"); // anti S-quark
            Symbol T("T"); // T-quark
            Symbol Tbar("Tbar"); // anti T-quark
            Symbol B("B"); // B-quark
            Symbol Bbar("Bbar"); // anti B-quark
        
            Symbol g("g"); // gluon
            Symbol gh("gh"); // gluon ghost
            Symbol ghbar("ghbar"); // anti gluon ghost
        
            Symbol A("A"); // photon
            Symbol Wm("Wm"); // W-
            Symbol Wp("Wp"); // W+
            Symbol Z("Z"); // Z
        
            Symbol ghA("SW"); // photon ghost
            Symbol ghAbar("ghAbar"); // anti photon ghost
            Symbol ghWm("ghWm"); // W- ghost
            Symbol ghWmbar("ghWmbar"); // anti W- ghost
            Symbol ghWp("ghWp"); // W+ ghost
            Symbol ghWpbar("ghWpbar"); // anti W+ ghost
            Symbol ghZ("ghZ"); // Z ghost
            Symbol ghZbar("ghZbar"); // anti Z ghost

            Symbol em("em"); // e-
            Symbol ep("ep"); // e+
            Symbol ne("ne"); // e-neutrino
            Symbol nebar("nebar"); // anti e-neutrino
            Symbol mum("mum"); // mu-
            Symbol mup("mup"); // mu+
            Symbol nmu("nmu"); // mu-neutrino
            Symbol nmubar("nmubar"); // anti mu-neutrino
            Symbol taum("taum"); // tau-
            Symbol taup("taup"); // tau+
            Symbol ntau("ntau"); // tau-neutrino
            Symbol ntaubar("ntaubar"); // anti tau-neutrino
        
            Symbol chi("chi"); // Z goldstone
            Symbol phim("phim"); // W- goldstone
            Symbol phip("phip"); // W+ goldstone
            Symbol H("H"); // higgs
            
            lst fs = { U, D, C, S, T, B, em, mum, taum, ne, nmu, ntau};
            lst fs2 = {Ubar, Dbar, Cbar, Sbar, Tbar, Bbar, ep, mup, taup, nebar, nmubar, ntaubar};
            string fsi[] = { "$u$", "$d$", "$c$", "$s$", "$t$", "$b$", "$e$", "$\\mu$", "$\\tau$", "$\\nu_e$", "$\\nu_\\mu$", "$\\nu_\\tau$"  };
            for(int i=0; i<fs.nops(); i++) {
                auto si = fsi[i];
                LineTeX[fs.op(i)] = "fermion ,edge label="+si;
                LineTeX[fs2.op(i)] = "anti fermion ,edge label="+si;
            }
            
            LineTeX[Wp] = "photon, edge label=$W$";
            LineTeX[Wm] = "photon, edge label=$W$";
            LineTeX[Z] = "photon, edge label=$Z$";
            
            LineTeX[H]="scalar, edge label=$H$";
            LineTeX[chi]="scalar, edge label=$\\chi$";
            LineTeX[phim]="charged scalar, edge label=$\\phi$";
            LineTeX[phip]="anti charged scalar, edge label=$\\phi$";
            
            lst ghs = {gh, ghA, ghWm, ghWp, ghZ};
            lst ghs2 = {ghbar, ghAbar, ghWmbar, ghWpbar, ghZbar};
            string ghsi[] = {"$c$", "$c_\\gamma$", "$c_-$", "$c_+$", "$c_Z$"};
            
            for(int i=0; i<ghs.nops(); i++) {
                auto si = ghsi[i];
                LineTeX[ghs.op(i)] = "ghost ,edge label="+si;
                LineTeX[ghs2.op(i)] = "anti ghost ,edge label="+si;
            }
        }
    }
    QGRAF::Process::_init QGRAF::Process::Process_init;
    
    // FROM IBP
    int FIRE::Version = 6;
    exmap MapPreSP;
    
    string UKIRA::KArgs = "";
    string KIRA::KArgs = "";
    
    // FROM QCD
    int QCD::FF::cur_mode = 0; // 0 - gluon, 1 - quark, 2 - anti-quark
    
    // FROM DE
    slong error_pass_dp = 100;
    
    // FROM AMF
    int DEX::Threads = 0;
    
    // gamma matrix
    ex DGamma::gi = GAS(1);
    ex DGamma::g5 = GAS(5);
    ex DGamma::C = GAS(int('c'));
    
    // global init
    std::stack<cln::float_format_t> cln_prec_stack;
    std::stack<long> digits_stack;
    _global_init::_init::_init() {
        ostringstream oss;
        string path = InstallPrefix + "/bin";
        if(dir_exists(path)) oss << path;
                
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
        table.insert(std::string("AsGamma"), []()->GiNaC::basic*{ return new AsGamma(); });
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

std::string HepLib::QGRAF::Models::GluonFF = R"EOF(
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

std::string HepLib::QGRAF::Models::QuarkFF = R"EOF(
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

std::string HepLib::QGRAF::Models::SM = R"EOF(
% From Boehm, Denner, Joos @ Feynman Gauge
[ model = 'Standard Model' ]

%------------------------------------
% Propagators
%------------------------------------

% quarks
[U, Ubar, -]
[D, Dbar, -]
[C, Cbar, -]
[S, Sbar, -]
[T, Tbar, -]
[B, Bbar, -]

% gluon and its ghost:
[g, g, +, notadpole]
[gh, ghbar, -]

% V: EWSM gauge bosons
[A, A, +, notadpole]
[Wm, Wp, +]
[Z, Z, +]

% leptons
[em, ep, -]
[ne, nebar, -]
[mum, mup, -]
[nmu, nmubar, -]
[taum, taup, -]
[ntau, ntaubar, -]

% G: Faddeev-Popov Ghosts
[ghA, ghAbar, -]
[ghZ, ghZbar, -]
[ghWm, ghWmbar, -]
[ghWp, ghWpbar, -]

% S: scalars
[H, H, +]
[chi, chi, +]
[phim, phip, +]

%------------------------------------
% Vertices
%------------------------------------

% QCD
[Ubar, U, g; QCD='+1', QED='0']
[Dbar, D, g; QCD='+1', QED='0']
[Sbar, S, g; QCD='+1', QED='0']
[Cbar, C, g; QCD='+1', QED='0']
[Bbar, B, g; QCD='+1', QED='0']
[Tbar, T, g; QCD='+1', QED='0']

[g, g, g, g; QCD='+2', QED='0']
[g, g, g; QCD='+1', QED='0']
[ghbar, gh, g; QCD='+1', QED='0']

% VVVV
[Wp, Wm, Z, Z; QCD='0', QED='+2']
[Wp, Wm, A, Z; QCD='0', QED='+2']
[Wp, Wm, A, A; QCD='0', QED='+2']
[Wp, Wp, Wm, Wm; QCD='0', QED='+2']

% VVV
[Wp, Wm, A; QCD='0', QED='+1']
[Wp, Wm, Z; QCD='0', QED='+1']

% SSSS
[H, H, H, H; QCD='0', QED='+2']
[chi, chi, chi, chi; QCD='0', QED='+2']
[H, H, chi, chi; QCD='0', QED='+2']
[H, H, phim, phip; QCD='0', QED='+2']
[chi, chi, phim, phip; QCD='0', QED='+2']
[phim, phip, phim, phip; QCD='0', QED='+2']

% SSS
[H, H, H; QCD='0', QED='+1']
[H, chi, chi; QCD='0', QED='+1']
[H, phim, phip; QCD='0', QED='+1']

% VVSS
[Z, Z, H, H; QCD='0', QED='+2']
[Z, Z, chi, chi; QCD='0', QED='+2']
[Wp, Wm, H, H; QCD='0', QED='+2']
[Wp, Wm, chi, chi; QCD='0', QED='+2']
[Wp, Wm, phip, phim; QCD='0', QED='+2']

[A, A, phip, phim; QCD='0', QED='+2']
[Z, A, phip, phim; QCD='0', QED='+2']
[Z, Z, phip, phim; QCD='0', QED='+2']

[Wp, A, phim, H; QCD='0', QED='+2']
[Wm, A, phip, H; QCD='0', QED='+2']
[Wp, A, phim, chi; QCD='0', QED='+2']
[Wm, A, phip, chi; QCD='0', QED='+2']
[Wp, Z, phim, H; QCD='0', QED='+2']
[Wm, Z, phip, H; QCD='0', QED='+2']
[Wp, Z, phim, chi; QCD='0', QED='+2']
[Wm, Z, phip, chi; QCD='0', QED='+2']

% VSS
[Z, chi, H; QCD='0', QED='+1']
[A, phip, phim; QCD='0', QED='+1']
[Z, phip, phim; QCD='0', QED='+1']
[Wp, phim, H; QCD='0', QED='+1']
[Wm, phip, H; QCD='0', QED='+1']
[Wp, phim, chi; QCD='0', QED='+1']
[Wm, phip, chi; QCD='0', QED='+1']

% SVV
[H, Z, Z; QCD='0', QED='+1']
[H, Wp, Wm; QCD='0', QED='+1']
[phip, Wm, A; QCD='0', QED='+1']
[phim, Wp, A; QCD='0', QED='+1']
[phip, Wm, Z; QCD='0', QED='+1']
[phim, Wp, Z; QCD='0', QED='+1']

% ffV
[Ubar, U, A; QCD='0', QED='+1']
[Dbar, D, A; QCD='0', QED='+1']
[Sbar, S, A; QCD='0', QED='+1']
[Cbar, C, A; QCD='0', QED='+1']
[Bbar, B, A; QCD='0', QED='+1']
[Tbar, T, A; QCD='0', QED='+1']

[ep, em, A; QCD='0', QED='+1']
[mup, mum, A; QCD='0', QED='+1']
[taup, taum, A; QCD='0', QED='+1']

[Ubar, U, Z; QCD='0', QED='+1']
[Dbar, D, Z; QCD='0', QED='+1']
[Sbar, S, Z; QCD='0', QED='+1']
[Cbar, C, Z; QCD='0', QED='+1']
[Bbar, B, Z; QCD='0', QED='+1']
[Tbar, T, Z; QCD='0', QED='+1']

[ep, em, Z; QCD='0', QED='+1']
[mup, mum, Z; QCD='0', QED='+1']
[taup, taum, Z; QCD='0', QED='+1']

[nebar, ne, Z; QCD='0', QED='+1']
[nmubar, nmu, Z; QCD='0', QED='+1']
[ntaubar, ntau, Z; QCD='0', QED='+1']

[Dbar, U, Wm; QCD='0', QED='+1']
[Dbar, C, Wm; QCD='0', QED='+1']
[Dbar, T, Wm; QCD='0', QED='+1']

[Sbar, U, Wm; QCD='0', QED='+1']
[Sbar, C, Wm; QCD='0', QED='+1']
[Sbar, T, Wm; QCD='0', QED='+1']

[Bbar, U, Wm; QCD='0', QED='+1']
[Bbar, C, Wm; QCD='0', QED='+1']
[Bbar, T, Wm; QCD='0', QED='+1']

[Ubar, D, Wp; QCD='0', QED='+1']
[Ubar, S, Wp; QCD='0', QED='+1']
[Ubar, B, Wp; QCD='0', QED='+1']

[Cbar, D, Wp; QCD='0', QED='+1']
[Cbar, S, Wp; QCD='0', QED='+1']
[Cbar, B, Wp; QCD='0', QED='+1']

[Tbar, D, Wp; QCD='0', QED='+1']
[Tbar, S, Wp; QCD='0', QED='+1']
[Tbar, B, Wp; QCD='0', QED='+1']

[ep, ne, Wm; QCD='0', QED='+1']
[mup, nmu, Wm; QCD='0', QED='+1']
[taup, ntau, Wm; QCD='0', QED='+1']

[nebar, em, Wp; QCD='0', QED='+1']
[nmubar, mum, Wp; QCD='0', QED='+1']
[ntaubar, taum, Wp; QCD='0', QED='+1']

% Sff
[Ubar, U, H; QCD='0', QED='+1']
[Dbar, D, H; QCD='0', QED='+1']
[Sbar, S, H; QCD='0', QED='+1']
[Cbar, C, H; QCD='0', QED='+1']
[Bbar, B, H; QCD='0', QED='+1']
[Tbar, T, H; QCD='0', QED='+1']

[ep, em, H; QCD='0', QED='+1']
[mup, mum, H; QCD='0', QED='+1']
[taup, taum, H; QCD='0', QED='+1']

[Ubar, U, chi; QCD='0', QED='+1']
[Dbar, D, chi; QCD='0', QED='+1']
[Sbar, S, chi; QCD='0', QED='+1']
[Cbar, C, chi; QCD='0', QED='+1']
[Bbar, B, chi; QCD='0', QED='+1']
[Tbar, T, chi; QCD='0', QED='+1']

[ep, em, chi; QCD='0', QED='+1']
[mup, mum, chi; QCD='0', QED='+1']
[taup, taum, chi; QCD='0', QED='+1']

[Dbar, U, phim; QCD='0', QED='+1']
[Dbar, C, phim; QCD='0', QED='+1']
[Dbar, T, phim; QCD='0', QED='+1']

[Sbar, U, phim; QCD='0', QED='+1']
[Sbar, C, phim; QCD='0', QED='+1']
[Sbar, T, phim; QCD='0', QED='+1']

[Bbar, U, phim; QCD='0', QED='+1']
[Bbar, C, phim; QCD='0', QED='+1']
[Bbar, T, phim; QCD='0', QED='+1']

[Ubar, D, phip; QCD='0', QED='+1']
[Ubar, S, phip; QCD='0', QED='+1']
[Ubar, B, phip; QCD='0', QED='+1']

[Cbar, D, phip; QCD='0', QED='+1']
[Cbar, S, phip; QCD='0', QED='+1']
[Cbar, B, phip; QCD='0', QED='+1']

[Tbar, D, phip; QCD='0', QED='+1']
[Tbar, S, phip; QCD='0', QED='+1']
[Tbar, B, phip; QCD='0', QED='+1']

[ep, ne, phim; QCD='0', QED='+1']
[mup, nmu, phim; QCD='0', QED='+1']
[taup, ntau, phim; QCD='0', QED='+1']

[nebar, em, phip; QCD='0', QED='+1']
[nmubar, mum, phip; QCD='0', QED='+1']
[ntaubar, taum, phip; QCD='0', QED='+1']

% VGG
[ghWpbar, ghWp, A; QCD='0', QED='+1']
[ghWmbar, ghWm, A; QCD='0', QED='+1']
[ghAbar, ghWm, Wp; QCD='0', QED='+1']
[ghAbar, ghWp, Wm; QCD='0', QED='+1']
[ghWmbar, ghA, Wm; QCD='0', QED='+1']
[ghWpbar, ghA, Wp; QCD='0', QED='+1']
[ghWpbar, ghWp, Z; QCD='0', QED='+1']
[ghWmbar, ghWm, Z; QCD='0', QED='+1']
[ghZbar, ghWm, Wp; QCD='0', QED='+1']
[ghZbar, ghWp, Wm; QCD='0', QED='+1']
[ghWmbar, ghZ, Wm; QCD='0', QED='+1']
[ghWpbar, ghZ, Wp; QCD='0', QED='+1'] 

% SGG
[ghZbar, ghZ, H; QCD='0', QED='+1']
[ghWpbar, ghWp, H; QCD='0', QED='+1']
[ghWmbar, ghWm, H; QCD='0', QED='+1']
[ghWpbar, ghWp, chi; QCD='0', QED='+1']
[ghWmbar, ghWm, chi; QCD='0', QED='+1']
[ghWpbar, ghA, phip; QCD='0', QED='+1']
[ghWmbar, ghA, phim; QCD='0', QED='+1']
[ghWpbar, ghZ, phip; QCD='0', QED='+1']
[ghWmbar, ghZ, phim; QCD='0', QED='+1']
[ghZbar, ghWp, phim; QCD='0', QED='+1']
[ghZbar, ghWm, phip; QCD='0', QED='+1']

)EOF";

std::string HepLib::QGRAF::Models::QCD = R"EOF(
[ model = 'qcd Model' ]
%------------------------------------
% Propagators
%------------------------------------
% quark
[q, qbar, -]
[Q, Qbar, -]
% gluon and its ghost:
[gh, ghbar, -]
[g, g, +, notadpole]
% external photon
[A, A, +, external]
%------------------------------------
% Vertices
%------------------------------------
% QCD
[qbar, q, g; QCD='+1']
[Qbar, Q, g; QCD='+1']
[g, g, g, g; QCD='+2']
[g, g, g; QCD='+1']
[ghbar, gh, g; QCD='+1']
% external
[qbar, q, A; QCD='+0']
[Qbar, Q, A; QCD='+0']
)EOF";
