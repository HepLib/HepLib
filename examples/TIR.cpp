#include "HepLib.h"
#include "cln/cln.h"

using namespace HepLib;

int main() {

    Vector p1("p1"), p2("p2"), p3("p3"), p4("p4"), q1("q1"), q2("q2"), q3("q3");
    Index i1("i1"), i2("i2"), i3("i3"), i4("i4");
    Symbol m("m"), s("s"), t("t");
    
    // case: q1.i1 with p1,p2
    {
        cout << "Check q1.i1 with p1,p2:" << endl;
        ex in = SP(q1,i1);
        ex res = TIR(in, lst{q1}, lst{p1,p2});
                
        lst chks = {
            SP(p1,i1),
        };
        
        for(auto chk : chks) {
            auto check = form((in-res) * chk,false);
            cout << "check: 0 = " << flush;
            cout << check.normal() << endl;
        }
        cout << endl;
    }
    
    // case: q1.i1 q1.i2 with p1
    {
        cout << "Check q1.i1*q1.i2 with p1:" << endl;
        ex in = SP(q1,i1) * SP(q1, i2);
        ex res = TIR(in, lst{q1}, lst{p1});
                
        lst chks = {
            SP(p1,i1)*SP(p1,i2),
            SP(i1,i2)
        };
        
        for(auto chk : chks) {
            auto check = form((in-res) * chk,false);
            cout << "check: 0 = " << flush;
            cout << exnormal(check) << endl;
        }
        cout << endl;
    }
        
    // case: q1.i1 q1.i2 with p1,p2
    {
        cout << "Check q1.i1*q1.i2 with p1,p2:" << endl;
        ex in = SP(q1,i1)*SP(q1,i2);
        ex res = TIR(in, lst{q1}, lst{p1,p2});
                
        lst chks = {
            SP(p1,i1)*SP(p1,i2),
            p1(i1)*p2(i2),
            SP(i1,i2)
        };
        
        for(auto chk : chks) {
            auto check = form((in-res) * chk,false);
            cout << "check: 0 = " << flush;
            cout << exnormal(check) << endl;
        }
        cout << endl;
    }
        
    // case: q1.i1 q1.i2 with p1,p2,p3,p4
    {
        cout << "Check q1.i1*q1.i2 with p1,p2,p3,p4:" << endl;
        ex e0 = SP(q1,i1)*SP(q1,i2);
        ex res = TIR(e0, lst{q1}, lst{p1,p2,p3,p4});
                
        lst chks = {
            p1(i1)*p1(i2),
            p1(i1)*p2(i2),
            SP(i1,i2)
        };
        
        auto pi = cln::nextprobprime(3);
        for(auto li : lst{q1,q2,p1,p2,p3,p4})
            for(auto ei : lst{q1,q2,p1,p2,p3,p4}) {
                pi = cln::nextprobprime(pi+1);
                letSP(li,ei) = numeric(pi);
            }
        
        for(auto ep : chks) {
            auto check = form((e0-res) * ep, false);
            cout << "check: 0 = " << flush;
            cout << exnormal(check) << endl;
        }
        cout << endl;
        clearSP();
    }
    
    // case: q1.i1 q2.i2 with p1,p2,p3,p4
    {
        cout << "Check q1.i1*q2.i2 with p1,p2,p3,p4:" << endl;
        ex e0 = SP(q1,i1)*SP(q2,i2);
        ex res = TIR(e0, lst{q1,q2}, lst{p1,p2,p3,p4});
                
        lst chks = {
            p1(i1)*p1(i2),
            p1(i1)*p2(i2),
            SP(i1,i2)
        };
        
        auto pi = cln::nextprobprime(3);
        for(auto li : lst{q1,q2,p1,p2,p3,p4})
            for(auto ei : lst{q1,q2,p1,p2,p3,p4}) {
                pi = cln::nextprobprime(pi+1);
                letSP(li,ei) = numeric(pi);
            }
        
        for(auto ep : chks) {
            auto check = form((e0-res) * ep, false);
            cout << "check: 0 = " << flush;
            cout << exnormal(check) << endl;
        }
        cout << endl;
        clearSP();
    }

    // case: q1.i1 q1.i2 q2.i3 q2.i4 with p1,p2
    {
        cout << "Check q1.i1*q1.i2*q2.i3*q2.i4 with p1,p2" << endl;
        ex e0 = SP(q1,i1)*SP(q1,i2)*SP(q2,i3)*SP(q2,i4);
        ex res = TIR(e0, lst{q1,q2}, lst{p1,p2});
                
        lst chks = {
            p1(i1)*p1(i2)*p1(i3)*p1(i4),
            p1(i1)*p2(i2)*p2(i3)*p2(i4),
            p1(i1)*p1(i2)*p2(i3)*p2(i4),
            SP(i1,i2)*SP(i3,i4)
        };
        
        auto pi = cln::nextprobprime(3);
        for(auto li : lst{q1,q2,p1,p2,p3,p4})
            for(auto ei : lst{q1,q2,p1,p2,p3,p4}) {
                pi = cln::nextprobprime(pi+1);
                letSP(li,ei) = numeric(pi);
            }
        
        for(auto ep : chks) {
            auto check = form((e0-res) * ep, false);
            cout << "check: 0 = " << flush;
            cout << exnormal(check) << endl;
        }
        cout << endl;
        clearSP();
    }
    
    // case: q1.i1 q1.i2 q2.i3 q2.i4 with p1,p2,p3
    {
        cout << "Note: this section may run form some long time." << endl;
        cout << "Check q1.i1*q1.i2*q2.i3*q2.i4 with p1,p2,p3" << endl;
        ex e0 = SP(q1,i1)*SP(q1,i2)*SP(q2,i3)*SP(q2,i4);
        ex res = TIR(e0, lst{q1,q2}, lst{p1,p2,p3});
                
        lst chks = {
            p1(i1)*p1(i2)*p1(i3)*p1(i4),
            p1(i1)*p2(i2)*p2(i3)*p2(i4),
            p1(i1)*p1(i2)*p2(i3)*p2(i4),
            SP(i1,i2)*SP(i3,i4)
        };
        
        auto pi = cln::nextprobprime(3);
        for(auto li : lst{q1,q2,p1,p2,p3,p4})
            for(auto ei : lst{q1,q2,p1,p2,p3,p4}) {
                pi = cln::nextprobprime(pi+1);
                letSP(li,ei) = numeric(pi);
            }
        
        for(auto ep : chks) {
            auto check = form((e0-res) * ep, false);
            cout << "check: 0 = " << flush;
            cout << exnormal(check) << endl;
        }
        cout << endl;
        clearSP();
    }
    
                
    return 0;
}
