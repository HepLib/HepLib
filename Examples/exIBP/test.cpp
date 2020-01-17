#include "IBP.h"

using namespace HepLib;

int main(int argc, char** argv) {

    symbol k1("k1"), k2("k2");
    symbol p("p"), n("n");
    symbol m("m");
    
    ex z = ex(3)/4;
    ex m2 = 1;
    ex kp = 1;
    lst prop;
    
    auto eps = IBP::eps;
    
    prop.append(power(k1,2));
    prop.append(power(k2,2));
    prop.append(n*(k1+k2+2*p)-kp);
    
    prop.append(n*k2);
    prop.append(power(k1+k2, 2));
    prop.append(power(k1+k2+2*p, 2));
    prop.append(power(k2+2*p, 2));
    
    lst repl = { p*p == m2, n*n==0, n*p==z/2 };
    lst loop = { k1, k2 };
    lst ext = { p, n };
    
    lst ils = lst{2, 1,0};
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
        item.let_op(prop.nops()-1) = item.op(prop.nops()-1)+eps;
        if(r<8) seeds.push_back(item);
    }
    
    auto ibp = IBP();
    ibp.Verbose = 2;
    ibp.Variables = lst{ m };
    ibp.Cuts = lst{0,1,2};
    ibp.Prepare(loop, ext, prop, repl);
    ibp.Generate(seeds);
    ibp.Reduce();
    
    ex obj = F(lst{1, 1, 1, 0, 2, 1, 1+eps});
    cout << endl << obj << "=";
    obj = subs(obj, ibp.FSolution);
    cout << collect_common_factors(obj) << endl;

    return 0;
}
