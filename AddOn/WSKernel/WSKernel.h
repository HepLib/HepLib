/**
 * @file 
 * @brief WSKernel header file
 */
 
#pragma once

#include "wstp.h"
#include <iostream>
#include <cstring>

namespace HepLib {

    using namespace std;
    
    /**
     * @brief interface to Wolfram Mathematica
     */
    class WSKernel {
    public:
        /**
         * @brief inner Error class for WSKernel
         */
        class Error : public exception {
        public:
            string msg;
            const char * what() const throw ();
            Error(WSLINK lp);
        };
        
        WSKernel(const string & open_str="");
        ~WSKernel();
        string Evaluate(const string & expr, const string & OutputForm="InputForm");
    private:
        WSENV ep;
        WSLINK lp;
    };
    
}
