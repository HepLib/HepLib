#include "SD.h"

using namespace HepLib;

int main(int argc, char** argv) {

    //SD::debug = true;
    
    auto ep = SD::ep;
    auto iep = SD::iEpsilon;
    
    Digits = 30;
    
    symbol k1("k1"), k2("k2"), q1("q1"), q2("q2");
    ex mc = 1;
    
    FeynmanParameter fp;
    
    fp.LoopMomenta = lst {q1, q2};
    fp.Propagators = lst{-pow(q1, 2), 2*k2*q1 - pow(q1, 2), k1*q1 + k2*q1 - pow(q1, 2),
        -(k1*q1) + k2*q1 + k1*q2 - k2*q2 + 2*q1*q2 + pow(mc, 2) - pow(q1, 2)
        - pow(q2, 2), k1*q1 + k2*q1 - k1*q2 - k2*q2 + 2*q1*q2 - pow(mc, 2) -
        pow(q1, 2) - pow(q2, 2), pow(mc, 2) - pow(q2, 2), k1*q2 + k2*q2 -
        pow(mc, 2) - pow(q2, 2)};
    fp.Exponents = lst{ 1, 1, 1, 1, 1, 1, 1 };
    fp.lReplacements[k1*k1] = 0;
    fp.lReplacements[k2*k2] = 0;
    fp.lReplacements[k1*k2] = 2*mc*mc;
    fp.Prefactor = SD::PrefactorFIESTA(2);
        
    SD work;
    work.epN = 0;
    work.Verbose = 100;
    char *CFLAGS = getenv("SD_CFLAGS");
    work.CFLAGS = CFLAGS;
    
    cout << endl << "Starting @ " << now() << endl;
    if(work.SecDec==NULL) work.SecDec = new SecDecG();
    if(work.Integrator==NULL) work.Integrator = new HCubature();
    if(work.Minimizer==NULL) work.Minimizer = new MinUit();
    if(strlen(CFLAGS)<1) CFLAGS = getenv("SD_CFLAGS");
    
    //work.Initialize(fp);
    //if(work.FunExp.size()<1) return 0;
    //work.Scalelesses();
    //work.ChengWu();
    
    work.FunExp.clear();
    lst fun = lst{exp(2*ep*Euler)*tgamma(3 + 2*ep), pow(x(5), 2)*x(0) + pow(x(5),
    2)*x(1) + pow(x(5), 2)*x(2) + pow(x(2), 2)*x(3) + pow(x(5), 2)*x(3) +
    2*x(1)*x(2)*x(3) + pow(x(2), 2)*x(4) + pow(x(5), 2)*x(4) +
    2*x(1)*x(2)*x(4) + pow(x(2), 2)*x(5) + 2*x(1)*x(2)*x(5) +
    2*x(0)*x(3)*x(5) + 2*x(2)*x(3)*x(5) + 2*x(1)*x(4)*x(5) +
    2*x(2)*x(4)*x(5) + pow(x(2), 2)*x(6) + 2*x(1)*x(2)*x(6) +
    2*x(2)*x(3)*x(6) - 4*x(0)*x(4)*x(6), x(0)*x(3) + x(1)*x(3) +
    x(2)*x(3) + x(0)*x(4) + x(1)*x(4) + x(2)*x(4) + x(0)*x(5) + x(1)*x(5)
    + x(2)*x(5) + x(3)*x(5) + x(4)*x(5) + x(0)*x(6) + x(1)*x(6) +
    x(2)*x(6) + x(3)*x(6) + x(4)*x(6)};
    lst exp = lst{1,-3-2*ep,1+3*ep};
    work.FunExp.push_back(make_pair(fun, exp));
    
work.Deltas.push_back(lst{x(0),x(1),x(2),x(3),x(4),x(5),x(6)});
for(auto & kv : work.FunExp) {
    
    for(int i=0; i<kv.first.nops(); i++) {
        kv.first.let_op(i) = kv.first.op(i).subs(x(4)==y(4)*x(5)/x(6)).subs(y(4)==x(4));
    }
    kv.first.append(x(5));
    kv.first.append(x(6));
    kv.second.append(1);
    kv.second.append(-1);
    
    {//projective
        auto delta = work.Deltas[0];
        symbol s;
        lst sRepl;
        ex xsum = 0;
        for(int j=0; j<delta.nops(); j++) {
            sRepl.append(delta.op(j)==delta.op(j)*s);
            xsum += delta.op(j);
        }
        
        ex over_all_sn = 0;
        int nnn = kv.first.nops();
        for(int i=0; i<nnn; i++) {
            if(!kv.first.op(i).has(x(wild()))) continue;
            auto tmp = expand(kv.first.op(i));
            auto sn = tmp.subs(sRepl).degree(s);
            over_all_sn += sn*kv.second.op(i);
            lst items;
            if(is_a<add>(tmp)) {
                for(auto ii : tmp) items.append(ii);
            } else {
                items.append(tmp);
            }
            tmp = 0;
            for(auto ii : items) {
                auto sni = ii.subs(sRepl).degree(s);
                if(sni-sn!=0) cout << "sn-sni=" << sn-sni << endl;
                if(sni!=sn) tmp += ii * pow(xsum, sn-sni);
                else tmp += ii;
            }
            kv.first.let_op(i) = tmp;
        }

        over_all_sn = normal(over_all_sn+delta.nops());
        if(!over_all_sn.is_zero()) {
            kv.first.append(xsum);
            kv.second.append(ex(0)-over_all_sn);
        }
        cout << "osn1=" << over_all_sn << endl;
    }
    
    //其中一部分
    //x4==y4/5, x5==(4*y4+5*y5)/5
    auto nnn = kv.first.nops();
    for(int i=0; i<nnn; i++) {
        auto tmp = kv.first.op(i).subs(lst{x(4)==y(4)/5, x(5)==(4*y(4)+5*y(5))/5}).subs(y(wild())==x(wild()));
        auto nd = normal(tmp).numer_denom();
        kv.first.let_op(i) = nd.op(0);
        if(nd.op(1)!=1) {
            kv.first.append(nd.op(1));
            kv.second.append(ex(0)-kv.second.op(i));
        }
    }
    kv.first.let_op(0) = kv.first.op(0)/5;
    
    {//projective
        auto delta = work.Deltas[0];
        symbol s;
        lst sRepl;
        ex xsum = 0;
        for(int j=0; j<delta.nops(); j++) {
            sRepl.append(delta.op(j)==delta.op(j)*s);
            xsum += delta.op(j);
        }
        
        ex over_all_sn = 0;
        int nnn = kv.first.nops();
        for(int i=0; i<nnn; i++) {
            if(!kv.first.op(i).has(x(wild()))) continue;
            auto tmp = expand(kv.first.op(i));
            auto sn = tmp.subs(sRepl).degree(s);
            over_all_sn += sn*kv.second.op(i);
            lst items;
            if(is_a<add>(tmp)) {
                for(auto ii : tmp) items.append(ii);
            } else {
                items.append(tmp);
            }
            tmp = 0;
            for(auto ii : items) {
                auto sni = ii.subs(sRepl).degree(s);
                if(sni-sn!=0) cout << "sn-sni=" << sn-sni << endl;
                if(sni!=sn) tmp += ii * pow(xsum, sn-sni);
                else tmp += ii;
            }
            kv.first.let_op(i) = tmp;
        }

        over_all_sn = normal(over_all_sn+delta.nops());
        if(!over_all_sn.is_zero()) {
            kv.first.append(xsum);
            kv.second.append(ex(0)-over_all_sn);
        }
        cout << "osn2=" << over_all_sn << endl;
    }
        
//    cout << kv.first.op(1) << endl;
//    exit(0);
    
    
    
}
    work.ChengWu();
//exit(0);
    work.RemoveDeltas();
    work.KillPowers();
    work.SDPrepares();
    work.EpsEpExpands();
    work.CIPrepares();
    auto pps = work.ParallelProcess;
    work.ParallelProcess = 0;
    work.Contours();
    work.Integrates();
    work.ParallelProcess = pps;
    delete work.SecDec;
    delete work.Integrator;
    delete work.Minimizer;
    cout << "Finished @ " << now() << endl << endl;

    cout << work.VEResult() << endl;
            
    return 0;
}
