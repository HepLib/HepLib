#include "HepLib.h"

using namespace HepLib;
using namespace QGRAF;

int main() {
    
    lst as = lst{ CI(1),CI(2),CI(3),CI(4),CI(5) };
    
    auto o1 = TTR(as);
    cout << "SUN Trace:" << endl;
    cout << "o1 = " << o1 << endl;
    fcout << "FC: " << ex(o1) << endl;
    cout << "conjugate(o1) = " << conjugate(o1) << endl;
    cout << "form(conjugate(o1)) = " << form(conjugate(o1)) << endl;
    cout << endl << endl;
    
    cout << "SUNT:" << endl;
    auto o2 = SUNT(lst{CI(1),CI(2)},TI(1),TI(2)) * SUNT(lst{CI(4),CI(5)},TI(3),TI(4)) * SUNT(CI(3),TI(2),TI(3));
    cout << "o2 = " << o2 << endl;
    fcout << "FC: " << o2 << endl;
    cout << "conjugate(o2) = " << conjugate(o2) << endl;
    cout << "form(conjugate(o2)) = " << form(conjugate(o2)) << endl;
    cout << "TI(1)->TI(4) = " << form(conjugate(o2).subs(TI(1)==TI(4))) << endl;
    cout << endl << endl;
    
    cout << "SUNF:" << endl;
    auto o3 = SUNF(CI(1),CI(2),CI(3)) * TTR(lst{CI(1),CI(2),CI(3),CI(4),CI(5)});
    cout << "o3 = " << o3 << endl;
    fcout << "FC: " << o3 << endl;
    cout << "form(o3) = " << form(o3) << endl;
    cout << endl << endl;
    
    cout << "SUNF4:" << endl;
    auto o4 = SUNF4(CI(1),CI(2),CI(3),CI(4)) * SUNF4(CI(5),CI(6),CI(7),CI(8));
    cout << "o4 = " << o4 << endl;
    fcout << "FC: " << o4 << endl;
    cout << "form(o4^2) = " << exfactor(form(o4*o4)) << endl;
    cout << endl << endl;
        
    return 0;
    
}
