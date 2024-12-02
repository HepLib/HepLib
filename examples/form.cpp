#include "HepLib.h"

using namespace HepLib;

int main() {

    Index i1("i1"), i2("i2"), i3("i3"), i4("i4"), i5("i5"), i6("i6");
    Index ti1("ti1",Index::Type::CF),  ti2("ti2",Index::Type::CF), ti3("ti3",Index::Type::CF), ti4("ti4",Index::Type::CF);
    Index a1("a1",Index::Type::CA), a2("a2",Index::Type::CA), a3("a3",Index::Type::CA), a4("a4",Index::Type::CA);
    Vector p1("p1"), p2("p2"), p3("p3"), p4("p4"), q1("q1"), q2("q2");
    Symbol m("m"), c("c");
    
    // ScalarProduct
    ex e1 = SP(p1+3*p2,c*p2+p4);
    cout << "(p1+3*p2).(c*p2+p4): " << e1 << endl;
    cout << endl;
    
    // Lorentz index contraction
    ex e2 = SP(p1+3*p2,i1) * SP(c*p2+p4, i2) * SP(i1,i2);
    cout << "by contraction: " << form(e2) << endl;
    cout << endl;
    
    // Dirac gamma object
    ex e3 = GAS(p1+3*p2) * GAS(i1) * (GAS(p2+p4) + m*GAS(1));
    cout << "gmamma chain: " << e3 << endl;
    ex e4 = TR(e3);
    cout << e4 << " = " << form(e4) << endl;
    cout << endl;
    
    // two traces in a single element
    ex e5 = GAS(p1+3*p2) * GAS(i1) * (GAS(p2+p4) + m*GAS(1)) * (GAS(p1+c*p3) + m*GAS(1));
    ex e6 = TR(e5) * e4;
    cout << "Two traces: " << form(e6) << endl;
    cout << endl;

    // form with verbose output
    cout << "form with verbose output:" << endl;
    auto e7 = TR(GAS(i1)*GAS(i2)*GAS(i3)*GAS(i4));
    cout << form(e7,1) << endl;
    cout << form(e7*e7,2) << endl;
    cout << endl;
    
    // SUNT object
    ex e8 = SUNT(ti1,ti2,a1) * SUNT(ti2,ti3,a2) * SUNT(ti3,ti4,a1) * SUNT(ti4,ti1,a2);
    cout << "SUNT trace: " << form(e8) << endl;
    cout << endl;
    
    // SUNF object
    ex e9 = SUNF(a1,a2,a3) * SUNF(a1,a2,a4) * SP(a3,a4);
    cout << "SUNF case: " << form(e9) << endl;
    cout << endl;
    
    // Levi-Civita object
    ex e10 = LC(i1,i2,i3,i4) * LC(i1,i2,i3,i4);
    cout << "LC case: " << exfactor(form(e10)) << endl;
    cout << endl;
    
    // form error and form again
    try {
        cout << form(TR(GAS(5)+m*GAS(1)),100) << endl;
    } catch(Error& err) {
        cout << "Error: " << err.what() << endl;
    }
    cout << "using form again after error:" << endl;
    cout << form(TR(GAS(i1)*GAS(i2)*GAS(i3)*GAS(i4))) << endl;
    
    return 0;
}

