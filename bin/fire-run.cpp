#include "HepLib.h"

using namespace std;
using namespace HepLib;

int main(int argc,char *argv[]) {
    if(argc<3) {
        cout << "usage: " << argv[0] << " <host> <port>" << endl;
        return 0;
    }
    
    string sip = argv[1];
    string sport = argv[2];
    string suffix = "6";
    if(argc>3) suffix = argv[3];
        
    while(true) {
          
        std::string si = Server::Next(sip, sport);
        if(!file_exists(si+".start")) continue;
        if(file_exists(si+".log")) continue;
        if(file_exists(si+".tables")) continue;
        system(("touch "+si+".log").c_str());
        cout << "running si = " << si << endl;
        
        string ostr = RunOS("which FIRE" + suffix);
        if(ostr.find("which: no FIRE") != std::string::npos) {
            cout << "FIRE NOT FOUND on HOST: " << RunOS("hostname") << endl;
            break;
        }
                
        system(("$(which FIRE" + suffix + ") -silent -parallel -c "+si).c_str());
        system(("rm -rf $(cat "+si+".config | grep database | sed s/'#database '//)").c_str());
        
        if(file_exists(si+".log")) remove((si+".log").c_str());
    }

    return 0;
}
