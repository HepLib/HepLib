/**
 * @file 
 */
 
#pragma once

#include "pstream.h"
#include <iostream>
#include <sstream>
#include <cstring>
#include <fcntl.h>
#include "Basic.h"

namespace HepLib {
    
    using namespace std;
    
    class Process {
    public:
        static const redi::pstreams::pmode pm_in = redi::pstreams::pstdin;
        static const redi::pstreams::pmode pm_out = redi::pstreams::pstdout;
        static const redi::pstreams::pmode pm_err = redi::pstreams::pstderr;
        
        void Open(string cmds, const redi::pstreams::pmode pm=pm_in|pm_out|pm_err);
        string ReadLine();
        string ReadLines(string);
        
        redi::pstream &io();
        
    private:
        redi::pstream pio;
    };
    
    class Fermat {
    public:
        string Sentinel = "---EOF---";
        void Init(string fer_path="fer64");
        string Execute(string);
        void Exit();
        ~Fermat();
                
    private:
        bool inited = false;
        int P2C[2];
        int C2P[2];
        pid_t pid = 0;
    };
    
    class Form {
    public:
        string Sentinel = "---EOF---";
        string Prompt = "***EOF***";
        void Init(string form_path="form");
        string Execute(string script, const string & out_var="[o]");
        void Exit();
        ~Form();
    
    private:
        bool inited = false;
        int io[2][2];
        int stdo[2];
        pid_t pid = 0;
    };
    
}
