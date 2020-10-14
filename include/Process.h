/**
 * @file
 * @brief Process header file
 */
 
#pragma once

#include <iostream>
#include <sstream>
#include <cstring>
#include "Basic.h"

namespace HepLib {
    
    using namespace std;

    /**
     * @brief interface to communicate with Fermat program
     */
    class Fermat {
    public:
        string Sentinel = "---EOF---";
        void Init(string fer_path="fer64");
        string Execute(string);
        void Exit();
        ~Fermat();
                
    private:
        bool inited = false;
        bool exited = false;
        int P2C[2];
        int C2P[2];
        pid_t fpid = 0;
        pid_t pid = 0;
    };
    
    /**
     * @brief interface to communicate with Form program
     */
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
        bool exited = false;
        int io[2][2];
        int stdo[2];
        pid_t fpid = 0;
        pid_t pid = 0;
    };
    
}
