#include "HepLib.h"

using namespace std;
using namespace HepLib;

Vector P("P"), Q("Q"), q1("q1"), q2("q2"), q("q");
Symbol m("m"), s("s"), zm("zm"), t("t");


ex Amps();
int main(int argc,char *argv[]) {
    if(argc<4) {
        cout << "usage: " << argv[0] << " <table> <start_index> <end_index>" << endl;
        cout << "  e.g. " << argv[0] << " 0.tables 100 150 for 0-[100,101,...,150].tables" << endl;
        return 0;
    }
    
    string tname = argv[1];
    int si = stoi(argv[2]);
    int ei = stoi(argv[3]);
    FIRE::ThieleTables(tname, si, ei);
    return 0;
}
