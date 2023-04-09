#include "HepLib.h"

using namespace std;
using namespace HepLib;

int main(int argc,char *argv[]) {
    if(argc<4) {
        cout << "usage: " << argv[0] << " <host> <port> <exe path> [opt] [nolog]" << endl;
        return 0;
    }
    
    string sip = argv[1];
    string sport = argv[2];
    string exe = argv[3];
    string opt = "";
    bool nolog = false;
    if(argc>4) { opt = argv[4]; opt = " "+opt; }
    if(argc>5) nolog = true;
        
    while(true) {
          
        std::string si = Server::Next(sip, sport);
        if(!file_exists(si+".start")) continue;
        if(file_exists(si+".log")) continue;
        if(file_exists(si+".tables")) continue;
        system(("touch "+si+".log").c_str());
        system(("echo $HOSTNAME > "+si+".log").c_str());
        
        if(!file_exists(exe)) {
            cout << "exe: " << exe << endl;
            cout << "FIRE NOT FOUND on HOST: " << RunOS("hostname") << endl;
            break;
        }
        
        if(nolog) system((exe+opt+" -c "+si).c_str());
        else system((exe+opt+" -c "+si+" >> "+si+".log").c_str());
        system(("dbdir=$(cat "+si+".config | grep database | sed s/'#database '//);test -d $dbdir && rm -rf $dbdir").c_str());
        
        if(file_exists(si+".log")) remove((si+".log").c_str());
    }

    return 0;
}
