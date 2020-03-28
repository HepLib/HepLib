#pragma once

#include "wstp.h"
#include <iostream>
#include <cstring>

namespace HepLib {

    using namespace std;
    
    class WSKernel {
    public:
        class Error : public exception {
        public:
            string msg;
            const char * what() const throw ();
            Error(WSLINK lp);
        };
        
        WSKernel(const string & open_str="");
        ~WSKernel();
        string Evaluate(const string & expr);
    private:
        WSENV ep;
        WSLINK lp;
    };
    
}
