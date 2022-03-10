#include "HepLib.h"

using namespace std;
using namespace HepLib;

Vector P("P"), Q("Q"), q1("q1"), q2("q2"), q("q");
Symbol m("m"), s("s"), zm("zm"), t("t");


ex Amps();
int main(int argc,char *argv[]) {
    if(argc<3) {
        cout << "usage: " << argv[0] << " <host> <port>" << endl;
        return 0;
    }
    
    string sip = argv[1];
    string sport = argv[2];
    
    while(true) {
          
        std::string si = Server::Next(sip, sport);        
        
        if(!file_exists(si+".start")) continue;
        if(file_exists(si+".log")) continue;
        if(file_exists(si+".tables")) continue;
        system(("touch "+si+".log").c_str());
        cout << "running si = " << si << endl;
                
        system(("$(which FIRE" + to_string(FIRE::Version) + ") -silent -parallel -c "+si).c_str());
        system(("rm -rf $(cat "+si+".config | grep database | sed s/'#database '//)").c_str());
        
        if(file_exists(si+".log")) remove((si+".log").c_str());
        
    }

    return 0;
}
