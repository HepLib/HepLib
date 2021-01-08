
#include "HepLib.h"

using namespace HepLib;

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
using namespace std;

int main(int argc, char **argv) {
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " [-d] file..." << endl;
		exit(1);
	}
	--argc; ++argv;
 
	bool dump_mode = false;
	try {
		lst l;
		while (argc) {
			if (strcmp(*argv, "-d") == 0) {
				dump_mode = true;
				--argc; ++argv;
			}
			std::ifstream f(*argv, std::ios_base::binary);
			archive ar;
			f >> ar;
			if (dump_mode) {
				ar.printraw(std::cout);
				std::cout << std::endl;
			} else {
				for (unsigned int i=0; i<ar.num_expressions(); ++i) {
					std::string name;
					ex e = ar.unarchive_ex(l, name, i);
					std::cout << name << " = " << e << std::endl;
				}
			}
			--argc; ++argv;
		}
	} catch (std::exception &e) {
		std::cerr << *argv << ": " << e.what() << std::endl;
	}
}
