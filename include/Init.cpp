#include "SD.h"
#include "FC.h"
#include "IBP.h"

namespace HepLib {

    // Initialization Unit
    // to make sure the initialization order is correct,
    // i.e., ordered in one compilation unit

    //----------------------------------------
    // HepLib
    //----------------------------------------
    std::map<std::string, ex> Symbol::Tables; // alias as symtab in parser
    
    ex w = wild();
    ex w0 = wild(0);
    ex w1 = wild(1);
    ex w2 = wild(2);
    ex w3 = wild(3);
    ex w4 = wild(4);
    ex w5 = wild(5);
    
    const Symbol ep("ep");
    const Symbol iEpsilon("iEpsilon",false);
    
    lst GiNaC_archive_Symbols = lst{};
    string InstallPrefix = "@PREFIX_DIR@";

    //----------------------------------------
    // HepLib::SD
    //----------------------------------------
    const Symbol SD::eps("eps");
    const Symbol SD::vs("vs");
    const Symbol SD::vz("vz");
    const Symbol SD::epz("epz");
    const Symbol SD::NaN("NaN");

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
    
    
    //----------------------------------------
    // HepLib::FC
    //----------------------------------------
    const Symbol FC::D("D");
    const Symbol FC::CA("CA");
    const Symbol FC::CF("CF");
    const Symbol FC::NA("NA");
    const Symbol FC::NF("NF");
    const Symbol FC::gs("gs");
    exmap FC::sp_map;
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
    
    //----------------------------------------
    // HepLib::IBP
    //----------------------------------------
    const Symbol IBP::d("d");
    
}