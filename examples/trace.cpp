#include "HepLib.h"
using namespace HepLib;
int main(int argc, char** argv) {
 Index mu("mu"), nu("nu");
 Vector p1("p1"), p2("p2");
 Symbol m("m");
 //note GAS(1) in gline, corresponds to the identity matrix
 ex gline = GAS(p1)*GAS(mu)*(GAS(p2)+m*GAS(1))*GAS(mu);
 ex trace = form(TR(gline));
 cout << trace << endl;
 return 0;
}

