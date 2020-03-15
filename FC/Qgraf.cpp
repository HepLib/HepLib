#include "FC.h"

namespace HepLib::FC {

    namespace {
        string di_("di");
        string li_("li");
        string ti_("ti");
        string fi_("fi");
        string ci_("ci");
        string ai_("ai");
        inline string n2s(ex fn) {
            int n = ex_to<numeric>(fn).to_int();
            return (n<0 ? "m" : "") + to_string(abs(n));
        }
    }

    REGISTER_FUNCTION(Propagator, do_not_evalf_params())
    REGISTER_FUNCTION(InField, do_not_evalf_params())
    REGISTER_FUNCTION(OutField, do_not_evalf_params())
    REGISTER_FUNCTION(Matij, do_not_evalf_params().set_return_type(return_types::commutative))
    
    unsigned Field2_SERIAL::serial = GiNaC::function::register_new(function_options("Field",2).do_not_evalf_params().overloaded(2));
    unsigned Field3_SERIAL::serial = GiNaC::function::register_new(function_options("Field",3).do_not_evalf_params().overloaded(2));
    unsigned Vertex2_SERIAL::serial = GiNaC::function::register_new(function_options("Vertex",2).do_not_evalf_params().overloaded(5));
    unsigned Vertex3_SERIAL::serial = GiNaC::function::register_new(function_options("Vertex",3).do_not_evalf_params().overloaded(5));
    unsigned Vertex4_SERIAL::serial = GiNaC::function::register_new(function_options("Vertex",4).do_not_evalf_params().overloaded(5));
    unsigned Vertex5_SERIAL::serial = GiNaC::function::register_new(function_options("Vertex",5).do_not_evalf_params().overloaded(5));
    unsigned Vertex6_SERIAL::serial = GiNaC::function::register_new(function_options("Vertex",6).do_not_evalf_params().overloaded(5));

    Index Qgraf::DI(ex fn) { return Index(di_+n2s(fn),Index::Type::VD); }
    Index Qgraf::LI(ex fn) { return Index(li_+n2s(fn),Index::Type::VD); }
    Index Qgraf::TI(ex fn) { return Index(ti_+n2s(fn),Index::Type::CF); }
    Index Qgraf::FI(ex fn) { return Index(fi_+n2s(fn),Index::Type::CF); }
    Index Qgraf::CI(ex fn) { return Index(ci_+n2s(fn),Index::Type::CA); }
    Index Qgraf::AI(ex fn) { return Index(ai_+n2s(fn),Index::Type::CA); }

    ex Qgraf::Amps(symtab st) {
        std::ofstream ofs;
        ofs.open("qgraf.dat", ios::out);
                
        ofs << "output='" << Output << "';" << endl;
        ofs << "style='" << InstallPrefix << "/include/GiNaC.sty';" << endl;
        ofs << "model='" << Model << "';" << endl;
        ofs << "in=" << In << ";" << endl;
        ofs << "out=" << Out << ";" << endl;
        ofs << "loops=" << Loops << ";" << endl;
        ofs << "loop_momentum=q;" << endl;
        ofs << "options=" << Options << ";" << endl;
        for(auto vs : Others) ofs << vs << ";" << endl;
        ofs.close();
        system((InstallPrefix+"/bin/qgraf > /dev/null").c_str());
        ifstream ifs(Output);
        string ostr((istreambuf_iterator<char>(ifs)), (istreambuf_iterator<char>()));
        ifs.close();
        if(access("qgraf.dat",F_OK)!=-1) remove("qgraf.dat");
        if(access(Output.c_str(),F_OK)!=-1) remove(Output.c_str());
        const char* rm_chars = " \t\v\r\n,";
        if(!ostr.empty()) {
            ostr.erase(0, ostr.find_first_not_of(rm_chars));
            ostr.erase(ostr.find_last_not_of(rm_chars)+1);
        }
        ostr = "{" + ostr + "}";

        Parser amp(st);
        ex ret = amp.Read(ostr);
        
        return ret;
    }
    
    namespace {
        
    }
    
}

