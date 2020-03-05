#pragma once

#include "pstream.h"
#include <iostream>
#include <sstream>
#include <cstring>

namespace HepLib {
    
    using namespace std;
    
    class Process {
    public:
        static const redi::pstreams::pmode pm_in = redi::pstreams::pstdin;
        static const redi::pstreams::pmode pm_out = redi::pstreams::pstdout;
        static const redi::pstreams::pmode pm_err = redi::pstreams::pstderr;
        
        void Open(const char *cmds, const redi::pstreams::pmode pm=pm_in|pm_out|pm_err);
        string ReadLine();
        string ReadLines(string);
        
        redi::pstream &io();
        
    private:
        redi::pstream pio;
    };
    
    class Fermat {
    public:
        const char * Sentinial = "---EOF---";
        void Init(const char * fer_path);
        string Execute(string);
        void Exit();
        
        class Error : public exception {
        public:
            string msg;
            const char * what() const throw ();
            Error(const char * _msg);
        };
        
    private:
        Process fermat;
    };
    
}
