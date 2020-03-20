#include "SD.h"
#include "FC.h"


namespace HepLib {

    // Initialization Unit
    // to make sure the initialization order is correct,
    // i.e., ordered in one compilation unit

    //----------------------------------------
    // HepLib
    //----------------------------------------
    ex w = wild();
    ex w0 = wild(0);
    ex w1 = wild(1);
    ex w2 = wild(2);
    ex w3 = wild(3);
    ex w4 = wild(4);
    ex w5 = wild(5);
    lst GiNaC_archive_Symbols = lst{};
    string InstallPrefix = "@PREFIX_DIR@";

    //----------------------------------------
    // HepLib::SD
    //----------------------------------------
    const symbol SD::ep("ep");
    const symbol SD::eps("eps");
    const symbol SD::vs("s");
    const symbol SD::vz("vz");
    const symbol SD::epz("epz");
    const symbol SD::iEpsilon("iEpsilon");
    const realsymbol SD::NaN("NaN");

    SD::SecDec::_init::_init() {
        GiNaC_archive_Symbols.append(ep);
        GiNaC_archive_Symbols.append(eps);
        GiNaC_archive_Symbols.append(vs);
        GiNaC_archive_Symbols.append(vz);
        GiNaC_archive_Symbols.append(epz);
        GiNaC_archive_Symbols.append(iEpsilon);
        GiNaC_archive_Symbols.append(NaN);
        GiNaC_archive_Symbols.sort();
        GiNaC_archive_Symbols.unique();
    }
    SD::SecDec::_init SD::SecDec::SD_init;
    
    SD::CppFormat::_init::_init() {
        set_print_func<numeric, CppFormat>(CppFormat::print_numeric);
    }
    SD::CppFormat::_init SD::CppFormat::CppFormat_init;
    
    
    //----------------------------------------
    // HepLib::FC
    //----------------------------------------
    realsymbol FC::D("D");
    realsymbol FC::CA("CA");
    realsymbol FC::CF("CF");
    realsymbol FC::NA("NA");
    realsymbol FC::NF("NF");
    realsymbol FC::gs("gs");
    exmap FC::sp_map;
    
    FC::FCFormat::_init::_init() {
        GiNaC_archive_Symbols.append(D);
        GiNaC_archive_Symbols.append(CA);
        GiNaC_archive_Symbols.append(CF);
        GiNaC_archive_Symbols.append(NA);
        GiNaC_archive_Symbols.append(NF);
        GiNaC_archive_Symbols.append(gs);
        set_print_func<ncmul, FCFormat>(FCFormat::ncmul_print);
    }
    FC::FCFormat::_init FC::FCFormat::FCFormat_init;
    
    FC::FormFormat::_init::_init() {
        set_print_func<power, FormFormat>(FormFormat::power_print);
    }
    FC::FormFormat::_init FC::FormFormat::FormFormat_init;
    
    FC::FCFormat FC::FCout(cout);
    
}
