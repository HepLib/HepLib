#include "HepLib.h"
#include "WSKernel.h"

using namespace HepLib;

int main() {
    WSKernel ws;
    auto o = ws.Evaluate("Normal@Series[Exp[x],{x,0,100}]");
    cout << o << endl;
    return 0;
}
