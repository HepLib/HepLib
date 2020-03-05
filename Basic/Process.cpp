#include "Process.h"

template class redi::basic_pstreambuf<char>;
template class redi::pstream_common<char>;
template class redi::basic_pstream<char>;
template class redi::basic_ipstream<char>;
template class redi::basic_opstream<char>;
template class redi::basic_rpstream<char>;

namespace HepLib {

    // Process Class
    void Process::Open(const char *cmds, const redi::pstreams::pmode pm) {
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
    
    // Fermat Class
    Fermat::Error::Error(const char * _msg) {
        msg = _msg;
    }
    
    const char * Fermat::Error::what() const throw () {
        return msg.c_str();
    }
    
    void Fermat::Init(const char* fer_path) {
        fermat.Open(fer_path);
        fermat.io() << "&M" << endl << endl; // prompt
        fermat.io() << "&(_d=90000)" << endl << endl; // width of the display on the window
        fermat.io() << "&d" << endl << "0" << endl; // off floating point representation
        fermat.io() << "&(_t=0)" << endl; // off a certain fast probabalistic algorithm
        fermat.io() << "&(t=0)" << endl; // off timing
        fermat.io() << "&(U=1)" << endl; // ugly printing
        fermat.io() << "&(_s=0)" << endl;
        fermat.io() << "&(_o=1000)" << endl; // http://home.bway.net/lewis/fer64mono.html
        fermat.io() << "!('" << Sentinial << "')" << endl;
        fermat.ReadLines(Sentinial);
        fermat.ReadLine(); // read 0
    }
    
    void Fermat::Exit() {
        fermat.io() << "&q" << endl;
    }
    
    string Fermat::Execute(string expr) {
        fermat.io() << expr << endl;
        fermat.io() << "!('" << Sentinial << "')" << endl;
        auto ostr = fermat.ReadLines(Sentinial);
        const char* WhiteSpace = " \t\v\r\n";
        if(!ostr.empty()) {
            ostr.erase(0, ostr.find_first_not_of(WhiteSpace));
            ostr.erase(ostr.find_last_not_of(WhiteSpace)+1);
        }
        fermat.ReadLine(); // read 0
        if(ostr.find("***")!=string::npos) throw Error(ostr.c_str());
        return ostr;
    }
}
