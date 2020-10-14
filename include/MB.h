/**
 * @file
 * @brief MB header file
 */

#pragma once

#include "Basic.h"

#include <dlfcn.h>

#include <string>
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

