/**
 * @file
 * @brief Class/Function for Process
 */

#include "BASIC.h"
#include <fstream>

inline bool file_exists(const char* fn) {
    return (access(fn,F_OK)!=-1);
}

namespace HepLib {
    
    //-----------------------------------------------------------
    // Fermat Class
    //-----------------------------------------------------------
    #define ENTER endl<<endl<<endl
    
    Fermat::~Fermat() { Exit(); }
        
    void Fermat::Exit() {
        if(getpid()!=pid) return;
        if(inited) {
            ostringstream script;
            script << "&q;" << endl << "&x;" << ENTER;
            string istr = script.str();
            write(P2C[1], istr.c_str(), istr.length());
            int st;
            waitpid(fpid, &st, WUNTRACED);
            inited = false;
            exited = true;
        }
    }
    
    void Fermat::Init(string fer_path) {
        if(inited) return;
        inited = true;
        pid = getpid();
        
        if (pipe(P2C)==-1 || pipe(C2P)==-1) {
            throw Error("pipe failed in Fermat::Init.");
        }
        
        fpid = fork();
        if (fpid == 0) { // child process            
            setpgid(0,0);
            close(P2C[1]);
            close(C2P[0]);
            dup2(C2P[1], 1);
            close(C2P[1]);
            dup2(P2C[0], 0);
            close(P2C[0]);
            execlp(fer_path.c_str(), fer_path.c_str(), NULL);
            exit(0);
        }
        
        // parent process
        close(P2C[0]);  // P2C[1] for write
        close(C2P[1]); // C2P[0] for read
        
        ostringstream script;
        script << "&(_d=90000);" << endl; // width of the display on the window
        script << "&d" << endl << "0;" << endl; // off floating point representation
        script << "&(_t=0);" << endl; // off a certain fast probabalistic algorithm
        script << "&(t=0);" << endl; // off timing
        script << "&(_s=0);" << endl;
        script << "&(_o=1000);" << endl; // http://home.bway.net/lewis/fer64mono.html
        script << "&(M=' ');" << endl; // prompt
        script << "!('" << Sentinel << "');" << ENTER;
        string istr = script.str();
        write(P2C[1], istr.c_str(), istr.length());
        
        string ostr;
        int n = 1024;
        char buffer[n+1]; // make sure the last one is '\0'
        int nio;
        while(true) {
            for(int i=0; i<n+1; i++) buffer[i] = '\0';
            nio = read(C2P[0], buffer, n);
            if(nio>0) ostr += buffer;
            else throw Error(ostr);
            auto cpos = ostr.find(Sentinel);
            if(cpos!=string::npos) {
                const char* WhiteSpace = " \t\v\r\n";
                auto lpos = ostr.find_last_not_of(WhiteSpace);
                if(ostr[lpos]!='0') read(C2P[0], buffer, n); // last 0, due to Sentinel
                ostr.erase(cpos);
                break;
            }
        }
        
        string_replace_all(ostr, "*** entry > 30 or < 5 means turn off mono multiply.", "");
        string_replace_all(ostr, "*** Fermat Warning. Early exit from mod_multivar_Chinese.", "");
        if(ostr.find("***")!=string::npos) {
            cout << "Fermat script: " << endl << istr << endl << endl;
            throw Error(ostr.c_str());
        }
    }
    
    // out string still contains the last number
    string Fermat::Execute(string expr) {
        if(exited) throw Error("Fermat has already exited.");
        if(getpid() != pid) throw Error("Fermat: can not Execute on child process.");
        ostringstream script;
        script << expr << endl;
        script << "!('" << Sentinel << "')" << ENTER;
        string istr = script.str();
        write(P2C[1], istr.c_str(), istr.length());
        
        string ostr;
        int n = buffer_size;
        char buffer[n+1]; // make sure the last one is '\0'
        int nio;
        while(true) {
            for(int i=0; i<n+1; i++) buffer[i] = '\0';
            nio = read(C2P[0], buffer, n);
            if(nio>0) ostr += buffer;
            else throw Error(ostr);
            auto cpos = ostr.find(Sentinel);
            if(cpos!=string::npos) {
                const char* WhiteSpace = " \t\v\r\n";
                auto lpos = ostr.find_last_not_of(WhiteSpace);
                if(ostr[lpos]!='0') read(C2P[0], buffer, n); // last 0, due to Sentinel 
                ostr.erase(cpos);
                break;
            }
        }
        string_trim(ostr);

        if(ostr.find("***")!=string::npos) throw Error(ostr.c_str());
        return ostr;
    }
    
    //-----------------------------------------------------------
    // Form Class
    //-----------------------------------------------------------
    Form::~Form() { Exit(); }
    
    void Form::Exit() {
        if(getpid()!=pid) return;
        if(inited && !exited) {
            string exit_cmd = "\n.end\n" + Prompt +"\n";
            write(io[0][1], exit_cmd.c_str(), exit_cmd.length());
            char buffer[8];
            read(io[1][0], buffer, 8);
            int st;
            waitpid(fpid, &st, WUNTRACED);
        }
        inited = false;
        exited = true;
    }
    
    void Form::Init(string form_path) {
        if(inited) return;
        inited = true;
        pid = getpid();
        
        if (pipe(io[0])==-1 || pipe(io[1])==-1 || pipe(stdo)==-1) {
            inited = false;
            exited = true;
            throw Error("pipe failed in Form::Init.");
        }
        
        fpid = fork();
        if (fpid == 0) {
            setpgid(0,0);
            close(io[0][1]);
            close(io[1][0]);
            close(stdo[0]);
            dup2(stdo[1], 1);
                        
            auto cpid = getpid(); // current process id
            ostringstream oss;
            oss << "init-" << cpid << ".frm";
            
            std::ofstream ofs;
            ofs.open(oss.str().c_str(), ios::out);
            if (!ofs) {
                inited = false;
                exited = true;
                throw Error("failed to open init.frm file!");
            }
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
            ofs << "    Format Nospaces;" << endl;
            ofs << "    ON NoSpacesInNumbers;" << endl;
            ofs << "    #write \""<<Sentinel<<"\\n\"" << endl;
            ofs << "    #toexternal \"%E\", `mexp'" << endl;
            ofs << "    #toexternal \""<<Sentinel<<"\"" << endl;
            ofs << "    #toexternal \"\\n\"" << endl;
            ofs << "#endprocedure" << endl;
            ofs << "#setexternal `PIPE1_';" << endl;
            ofs << "#toexternal \"OK\"" << endl;
            ofs << "Local [o]=0;" << endl;
            ofs << ".sort" << endl;
            ofs << "Format Mathematica;" << endl;
            ofs << "#prompt \"" << Prompt << "\"" << endl;
            ofs << "#fromexternal-" << endl;
            ofs << ".end" << endl;
            ofs.close();
            
            execlp(form_path.c_str(), form_path.c_str(), 
                "-pipe", 
                (to_string(io[0][0])+","+to_string(io[1][1])).c_str(),
                "-M", ("init-"+to_string(cpid)).c_str(),
                NULL);
            exit(0);
        } 
        
        close(io[0][0]);
        close(io[1][1]);
        close(stdo[1]);
        
        char buffer[1024];
        read(io[1][0], buffer, sizeof(buffer));
        char* p = strstr(buffer, "\n");
        if(p==NULL){
            inited = false;
            exited = true;
            cout << "the return is: <|" << buffer << "|>" << endl;
            throw Error("Init Failed: Expect a Line break!");
        }
        sprintf(p, ",%d\n\n", fpid);
        write(io[0][1], buffer, strlen(buffer));
        read(io[1][0], buffer, sizeof(buffer));
        p = strstr(buffer, "OK");
        if(p==NULL || p!=buffer) {
            inited = false;
            exited = true;
            throw Error("Init Failed: Expect OK!");
        }
    
        ostringstream oss;
        oss << "init-" << fpid << ".frm";
        if(file_exists(oss.str().c_str())) remove(oss.str().c_str());
        
        if(true) { // read the terminal output
            string estr;
            int n = 1024;
            char buffer[n+1]; // make sure the last one is '\0'
            int nio;
            while(true) {
                for(int i=0; i<n+1; i++) buffer[i] = '\0';
                nio = read(stdo[0], buffer, n);
                if(nio>0) estr += buffer;
                auto cpos = estr.find(Sentinel);
                if(cpos!=string::npos) break;
                cpos = estr.find("-->");
                if(cpos!=string::npos) {
                    inited = false;
                    exited = true;
                    throw Error(estr);
                }
            }
        }
    }
    
    string Form::Execute(string script, const string & out_var) {
        if(exited) throw Error("Form has already exited.");
        if(getpid() != pid) throw Error("Form: can not Execute on child process.");
        string istr = script;
        istr += "\n.sort\n#call put(";
        istr += out_var;
        istr += ")\n.sort\n";
        istr += Prompt + "\n"; // prompt
        
        write(io[0][1], istr.c_str(), istr.length());

        string ostr;
        int n = buffer_size;
        char buffer[n+1]; // make sure the last one is '\0'
        int nio;
        
        string estr;
        while(true) {
            for(int i=0; i<n+1; i++) buffer[i] = '\0';
            nio = read(stdo[0], buffer, n);
            if(nio>0) estr += buffer;
            auto cpos = estr.find(Sentinel);
            if(cpos!=string::npos) break;
            cpos = estr.find("-->");
            if(cpos!=string::npos) {
                inited = false;
                exited = true;
                throw Error(estr);
            }
        }
        
        while(true) {
            for(int i=0; i<n+1; i++) buffer[i] = '\0';
            nio = read(io[1][0], buffer, n);
            if(nio>0) ostr += buffer;
            else {
                string estr;
                while(true) {
                    for(int i=0; i<n+1; i++) buffer[i] = '\0';
                    nio = read(stdo[0], buffer, n);
                    if(nio<=0) break;
                    estr += buffer;
                }
                inited = false;
                exited = true;
                throw Error(estr.c_str());
            }
            auto cpos = ostr.find(Sentinel);
            if(cpos!=string::npos) {
                ostr.replace(cpos, Sentinel.length(), "");
                break;
            }
        }
        
        return ostr;
    }
    
    
}
