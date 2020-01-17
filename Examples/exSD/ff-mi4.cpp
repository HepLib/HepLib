#include "SD.h"

using namespace HepLib;

int main(int argc, char** argv) {
    
    auto ep = SD::ep;
    auto iep = SD::iEpsilon;
    
    symbol p("p"), n("n");
    symbol k1("k1"), k2("k2"), k3("k3");
    symbol K1("K1"), K2("K2"), K3("K3");
    symbol kp("kp");
    symbol m("m"), zz("zz");
    
    ex m2=m*m;
    ex z0 = ex(1)/4;
    
    ex k1p = z(1)*(1-zz)*kp;
    ex k1m = pow(K1,2)/(2*k1p);
    ex k2p = z(2)*(1-zz)*kp;
    ex k2m = pow(K2,2)/(2*k2p);
    ex k3p = z(3)*(1-zz)*kp;
    ex k3m = pow(K3,2)/(2*k3p);
    
    ex pp = zz/2*kp;
    ex pm = m2/(2*pp);
    
    FeynmanParameter fp;
    
    fp.LoopMomenta = lst{};
    fp.tLoopMomenta = lst{K1,K2};
    fp.Propagators = lst {
        k2*p,k1*k2 + k1*p + k2*p
    };
    fp.Exponents = lst{ 1, 2 };
    
    fp.lReplacements[p*p] = m2;
    fp.lReplacements[n*n] = 0;
    fp.lReplacements[n*p] = pp;
    fp.lReplacements[k1*k1] = 0;
    fp.lReplacements[k2*k2] = 0;
    fp.lReplacements[k3*k3] = 0;
    fp.lReplacements[n*k1] = k1p;
    fp.lReplacements[n*k2] = k2p;
    fp.lReplacements[n*k3] = k3p;
    fp.lReplacements[p*k1] = k1p*pm+k1m*pp;
    fp.lReplacements[p*k2] = k2p*pm+k2m*pp;
    fp.lReplacements[p*k3] = k3p*pm+k3m*pp;
    fp.lReplacements[k1*k2] = k1p*k2m+k2p*k1m-K1*K2;
    fp.lReplacements[k2*k3] = k2p*k3m+k3p*k2m-K2*K3;
    fp.lReplacements[k1*k3] = k1p*k3m+k3p*k1m-K1*K3;
    fp.lReplacements[m] = 1;
    fp.lReplacements[kp] = 1;
    fp.lReplacements[zz] = z0;
    
    fp.tReplacements[p*p] = m2;
    fp.tReplacements[n*n] = 0;
    fp.tReplacements[n*p] = pp;
    fp.tReplacements[m] = 1;
    fp.tReplacements[kp] = 1;
    fp.tReplacements[zz] = z0;
    
    fp.nReplacements[ep] = ex(1)/11;
    fp.nReplacements[zz] = z0;
    fp.nReplacements[z(1)] = ex(1)/3;
    fp.nReplacements[z(2)] = ex(1)/3;
    fp.nReplacements[z(3)] = ex(1)/3;
    
    cout << endl << "Starting @ " << now() << endl;
    
    SD work;
    work.epN = 1;
    work.Verbose = 2;
    
    
    char *CFLAGS = getenv("SD_CFLAGS");
    work.CFLAGS = CFLAGS;
    
    work.Initialize(fp);
    
    if(work.SecDec==NULL) work.SecDec = new SecDecG();
    if(work.Integrator==NULL) work.Integrator = new HCubature();
    if(work.Minimizer==NULL) work.Minimizer = new MinUit();
    
    int xn = work.Deltas[0].nops();
    exmap z2x;
    
    lst zs;
    ex zFactor = 1;
    for(int i=1; i<=fp.tLoopMomenta.nops(); i++) {
        z2x[z(i)] = x(xn+i-1);
        zs.append(x(xn+i-1));
        zFactor /= x(xn+i-1);
    }
    
    work.Deltas.push_back(zs);
    
    for(auto &kv : work.FunExp) {
        
        kv.first.let_op(0) = kv.first.op(0) * zFactor;
        
        kv.first = lstHelper::subs(kv.first, z2x);
        
        auto tmp = collect_common_factors(kv.first.op(0));
        if(tmp.has(x(wild())) && is_a<mul>(tmp)) {
            ex rem = 1;
            for(auto item : tmp) {
                if(!item.has(x(wild()))) {
                    rem *= item;
                } else if(item.match(pow(wild(0), wild(1)))) {
                    kv.first.append(item.op(0));
                    kv.second.append(item.op(1));
                } else {
                    kv.first.append(item);
                    kv.second.append(1);
                }
            }
            kv.first.let_op(0) = rem;
        } else if(tmp.has(x(wild())) && tmp.match(pow(wild(0), wild(1)))) {
            kv.first.let_op(0) = 1;
            kv.first.append(tmp.op(0));
            kv.second.append(tmp.op(1));
        }
    }
    work.Normalizes();
    work.XReOrders();
    work.Normalizes();
    
    work.RemoveDeltas();
    work.SDPrepares();
    work.EpsEpExpands();
    work.CIPrepares();
    work.Contours();
    work.Integrates();
    
    delete work.SecDec;
    delete work.Integrator;
    delete work.Minimizer;
    cout << "Finished @ " << now() << endl << endl;
    
    cout << work.VEResult() << endl;
    cout << "check with:" << endl;
    cout << "(212.36633228802467 +- 8.9238E-11)+ep*(-1366.7797085867055 +- 5.4189E-7)" << endl;
    
    return 0;
}
