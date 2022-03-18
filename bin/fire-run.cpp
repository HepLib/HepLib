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
    string pq = "";
    if(argc>3) {
        string arg = argv[3];
        if(arg=="-p" || arg=="p") pq = "p";
        else if(arg=="-q" || arg=="q") pq = "q";
    }
    
    while(true) {
          
        std::string si = Server::Next(sip, sport);        
        
        if(!file_exists(si+".start")) continue;
        if(file_exists(si+".log")) continue;
        if(file_exists(si+".tables")) continue;
        system(("touch "+si+".log").c_str());
        cout << "running si = " << si << endl;
        
        string ostr = RunOS("which FIRE" + to_string(FIRE::Version) + pq);
        if(ostr.find("which: no FIRE") != std::string::npos) {
            cout << "FIRE NOT FOUND on HOST: " << RunOS("hostname") << endl;
            break;
        }
                
        system(("$(which FIRE" + to_string(FIRE::Version) + pq + ") -silent -parallel -c "+si).c_str());
        system(("rm -rf $(cat "+si+".config | grep database | sed s/'#database '//)").c_str());
        
        if(file_exists(si+".log")) remove((si+".log").c_str());
        
    }

    return 0;
}
