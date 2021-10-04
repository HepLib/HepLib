/**
 * @file
 * @brief MB header file
 */

#pragma once

#include "BASIC.h"

#include <dlfcn.h>

#include <string>
#include <signal.h>
#include <sys/syscall.h>
#include <sys/wait.h>
#include <sstream>
#include <ios>
#include <regex>

#ifdef _USE_FLOAT128
extern "C" {
    #include <quadmath.h>
}
#endif

namespace HepLib {

    /**
     * @brief MB class still under development
     */
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

