#include "Process.h"
#include <fstream>

template class redi::basic_pstreambuf<char>;
template class redi::pstream_common<char>;
template class redi::basic_pstream<char>;
template class redi::basic_ipstream<char>;
template class redi::basic_opstream<char>;
template class redi::basic_rpstream<char>;

inline bool file_exists(const char* fn) {
    return (access(fn,F_OK)!=-1);
}

namespace HepLib {

    //-----------------------------------------------------------
    // Process Class
    //-----------------------------------------------------------
    void Process::Open(string cmds, const redi::pstreams::pmode pm) {
        pio.open(cmds, pm);
        return;
    }
    
    redi::pstream &Process::io() {
        return pio;
    }

    string Process::ReadLine() {
        string aLine;
        getline(pio, aLine);
        return aLine;
    }
    
    string Process::ReadLines(string endLine) {
        string aLines;
        string aLine;
        while(getline(pio, aLine)) {
            if(aLine == endLine) break;
            aLines += aLine;
            aLines += '\n' ;
        }
        return aLines;
    }
    
    //-----------------------------------------------------------
    // Fermat Class
    //-----------------------------------------------------------
    Fermat::Error::Error(const char * _msg) : msg(_msg) { }
    
    const char * Fermat::Error::what() const throw () {
        return msg.c_str();
    }
    
    void Fermat::Init(string fer_path) {
        fermat.Open(fer_path);
        fermat.io() << "&M" << endl << endl; // prompt
        fermat.io() << "&(_d=90000)" << endl << endl; // width of the display on the window
        fermat.io() << "&d" << endl << "0" << endl; // off floating point representation
        fermat.io() << "&(_t=0)" << endl; // off a certain fast probabalistic algorithm
        fermat.io() << "&(t=0)" << endl; // off timing
        fermat.io() << "&(U=1)" << endl; // ugly printing
        fermat.io() << "&(_s=0)" << endl;
        fermat.io() << "&(_o=1000)" << endl; // http://home.bway.net/lewis/fer64mono.html
        fermat.io() << "!('" << Sentinel << "')" << endl;
        fermat.ReadLines(Sentinel);
        fermat.ReadLine(); // read 0
    }
    
    void Fermat::Exit() {
        fermat.io() << "&q" << endl;
    }
    
    string Fermat::Execute(string expr) {
        fermat.io() << expr << endl;
        fermat.io() << "!('" << Sentinel << "')" << endl;
        auto ostr = fermat.ReadLines(Sentinel);
        const char* WhiteSpace = " \t\v\r\n";
        if(!ostr.empty()) {
            ostr.erase(0, ostr.find_first_not_of(WhiteSpace));
            ostr.erase(ostr.find_last_not_of(WhiteSpace)+1);
        }
        fermat.ReadLine(); // read 0
        if(ostr.find("***")!=string::npos) throw Error(ostr.c_str());
        return ostr;
    }
    
    //-----------------------------------------------------------
    // Form Class
    //-----------------------------------------------------------
    Form::Error::Error(const char * _msg) : msg(_msg) { }
    
    const char * Form::Error::what() const throw () {
        return msg.c_str();
    }
    
    void Form::Exit() {
        if(inited) {
            write(io[0][1], ".end\n\n", 6);
            char buffer[8];
            read(io[1][0], buffer, 8);
            close(io[0][0]);
            close(io[0][1]);
            close(io[1][0]);
            close(io[1][1]);
            close(stdo[0]);
            close(stdo[1]);
            kill(pid, SIGTERM);
        }
    }
    
    void Form::Init(string form_path_args) {
    
        if(inited) {
            close(io[0][0]);
            close(io[0][1]);
            close(io[1][0]);
            close(io[1][1]);
            close(stdo[0]);
            close(stdo[1]);
            kill(pid, SIGTERM);
        }
        inited = true;
        
        pipe(io[0]);
        pipe(io[1]);
        pipe(stdo);
        
        pid = fork();
        if (pid == 0) {
            close(io[0][1]);
            close(io[1][0]);
            close(stdo[0]);
            dup2(stdo[1], 1);
            
            auto pid = getpid(); // current process id
            ostringstream oss;
            oss << "init-" << pid << ".frm";
            
            std::ofstream ofs;
            ofs.open(oss.str().c_str(), ios::out);
            if (!ofs) throw runtime_error("failed to open init.frm file!");
            ofs << "Off Statistics;" << endl;
            ofs << "#ifndef `PIPES_'" << endl;
            ofs << "    #message \"No pipes found\";" << endl;
            ofs << "    .end;" << endl;
            ofs << "#endif" << endl;
            ofs << "#if (`PIPES_' <= 0)" << endl;
            ofs << "    #message \"No pipes found\";" << endl;
            ofs << "    .end;" << endl;
            ofs << "#endif" << endl;
            ofs << "#procedure put(mexp)" << endl;
            ofs << "    #toexternal \"%E\", `mexp'" << endl;
            ofs << "    #toexternal \""<<Sentinel<<"\"" << endl;
            ofs << "#endprocedure" << endl;
            ofs << "#setexternal `PIPE1_';" << endl;
            ofs << "#toexternal \"OK\"" << endl;
            ofs << "Local [o]=0;" << endl;
            ofs << ".sort" << endl;
            ofs << "Format Mathematica;" << endl;
            ofs << "#fromexternal" << endl;
            ofs << ".end" << endl;
            ofs.close();
    
            oss.clear();
            oss.str("");
            oss << "%s -pipe %d,%d -M init-" << pid << endl;
            
            char buffer[256];
            sprintf(buffer, oss.str().c_str(), form_path_args.c_str(), io[0][0], io[1][1]);
            system(buffer);
            exit(0);
        } else {
            close(io[0][0]);
            close(io[1][1]);
            close(stdo[1]);
            fcntl(stdo[0], F_SETFL, fcntl(stdo[0], F_GETFL, 0) | O_NONBLOCK);
            
            char buffer[1024];
            read(io[1][0], buffer, sizeof(buffer));
            char* p = strstr(buffer, "\n");
            if(p==NULL) return throw Error("Init Failed: Expect a Line break!");
            sprintf(p, ",%d\n\n\0", pid);
            write(io[0][1], buffer, strlen(buffer));
            read(io[1][0], buffer, sizeof(buffer));
            p = strstr(buffer, "OK");
            if(p==NULL || p!=buffer) throw Error("Init Failed: Expect OK!");
        }
        
        ostringstream oss;
        oss << "init-" << pid << ".frm";
        if(file_exists(oss.str().c_str())) remove(oss.str().c_str());
    }
    
    string Form::Execute(string script, const char * out_var) {
        
        script += "\n.sort\n";
        script += "#call put(";
        script += out_var;
        script += ")\n";
        script += "\n.sort\n";
        script += "#fromexternal";
        
        // replace blank line
        while(true) {
            auto pos = script.find("\n\n");
            if(pos==string::npos) break;
            script.replace(pos, 2, "\n");
        }
        script += "\n\n"; // blank line to prompt
        
        write(io[0][1], script.c_str(), script.length());
        
        string ostr;
        int n = 1024;
        char buffer[n+1]; // make sure the last one is '\0'
        int nio;
        while(true) {
            for(int i=0; i<n+1; i++) buffer[i] = '\0';
            nio = read(io[1][0], buffer, n);
            if(nio>0) ostr += buffer;
            else {
                string err_str;
                while(true) {
                    for(int i=0; i<n+1; i++) buffer[i] = '\0';
                    nio = read(stdo[0], buffer, n);
                    if(nio<=0) break;
                    err_str += buffer;
                }
                throw Error(err_str.c_str());
            }
            auto cpos = ostr.find(Sentinel);
            if(cpos!=string::npos) {
                ostr.replace(cpos, strlen(Sentinel), "");
                break;
            }
        }
        
        return ostr;
    }
    
    
}
