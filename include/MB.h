#pragma once

#include "ExGiNaC.h"

#include <dlfcn.h>

#include <string>
#include <fcntl.h>
#include <signal.h>
#include <sys/syscall.h>
#include <sys/wait.h>
#include <sstream>
#include <ios>
#include <regex>

extern "C" {
    #include <quadmath.h>
}

namespace HepLib {

/*-----------------------------------------------------*/
// MB Input
/*-----------------------------------------------------*/
struct FeynmanParameter {
    lst LoopMomenta;
    lst Propagators;
    lst Exponents;
    ex Prefactor = 1;
    exmap lReplacements;
    exmap nReplacements;
};

/*-----------------------------------------------------*/
// MB Class
/*-----------------------------------------------------*/
class MB {

public:
    static const symbol ep;
    vector<pair<lst, lst>> FunExp;
    vector<lst> Deltas;
    int Verbose = 0;
    bool IsZero = false;
    
    void Initialize(FeynmanParameter fp);
    
};



}

