
#include "HepLib.h"

using namespace HepLib;

int main(int argc, char **argv) {
    auto verb = Verbose;
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " [-d]/[-k]/[key] file ..." << endl;
        exit(1);
    }
    std::string arg, gar;
    if(argc>2) { arg = argv[1]; gar = argv[2]; }
    else gar = argv[1];
    try {
        bool dump_mode = false;
        if(argc>2 && arg == "-d") dump_mode = true;

        std::ifstream f(gar);
        archive ar;
        f >> ar;
        if (dump_mode) {
            ar.printraw(std::cout);
            std::cout << std::endl;
        } else if(argc>2 && arg == "-k") {
            std::cout << "keys: ";
            for (unsigned int i=0; i<ar.num_expressions(); ++i) {
                std::string name;
                ex e = ar.unarchive_ex(name, i);
                std::cout << name << " ";
            }
            std::cout << std::endl;
        } else if(argc>2) {
            ex e = ar.unarchive_ex(arg.c_str());
            std::cout << e << std::endl;
        } else {
            for (unsigned int i=0; i<ar.num_expressions(); ++i) {
                std::string name;
                ex e = ar.unarchive_ex(name, i);
                std::cout << name << " = " << e << std::endl;
            }
        }
    } catch (std::exception &e) {
        std::cerr << *argv << ": " << e.what() << std::endl;
    }
}
