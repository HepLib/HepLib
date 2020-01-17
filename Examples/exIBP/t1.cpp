#include "IBP.h"

using namespace HepLib;

int main(int argc, char** argv) {

    symbol q1("q1"), q2("q2");
    symbol p("p");
    symbol m("m");
    
    ex m2 = m*m;
    lst prop;
    prop.append(q1*p);
    prop.append(power(q1+p, 2) - m2);
    //prop.append(power(q2, 2));
    //prop.append(power(q2+p, 2));
    //prop.append(power(q1+q2+p, 2) - m2);
    
    lst repl = { p*p == m2 };
    lst loop = { q1 };
    lst ext = {p};
    
    lst ils = lst{3,2,1,0,-1,-2,-3};
    long long tot = 1;
    for(int i=0; i<prop.nops(); i++) tot *= ils.nops();
    
    vector<lst> seeds;
    for(long long i=0; i<tot; i++) {
        lst item;
        long long ci = i;
        ex s=0, r=0;
        for(int j=0; j<prop.nops(); j++) {
            int c = ci % ils.nops();
            item.append(ils.op(c));
            ci = ci / ils.nops();
            
            if(ils.op(c)>0) r += ils.op(c);
            else s -= ils.op(c);
        }
        seeds.push_back(item);
    }
    
    auto ibp = IBP();
    ibp.Verbose = 2;
    ibp.Variables = lst{ m };
    ibp.Prepare(loop, ext, prop, repl);
    ibp.Generate(seeds);
    ibp.Reduce();
    
    
    ex obj = F(lst{1,1});
    cout << endl << obj << "=";
    cout << subs(obj, ibp.FSolution) << endl;

    return 0;
}
