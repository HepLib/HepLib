#include "HepLib.h"

using namespace std;
using namespace GiNaC;
using namespace HepLib;
using namespace HepLib::FC;
using namespace HepLib::SD;
using namespace Quarkonium;
using namespace Qgraf;

Vector k1("k1"), k2("k2"), k3("k3"), q1("q1"), q2("q2");
Vector p("p"), n("n"), h("h");
Index pi("pi"), piR("piR");
Symbol m("m"), kp("kp"), pp("pp"), zz("z"), K1("K1"), K2("K2"), K3("K3"), ie("ie");
ex p1 = p;
ex p2 = p;
ex k12=m*m;
ex k22=0;
ex k32=0;
int epN = -2;

ex NLoops(const ex & amps, int loop, int tloop);
lst Amps();

int main() {
    Verbose = 100;
    
    auto amp2 = Amps();
    //amp2 = lst{amp2.op(0)}; // to select diagrams
    
    // scalar product
    letSP(p)=m*m;
    letSP(n)=0;
    letSP(p,n)=pp;
    letSP(k1)=k12;
    letSP(k2)=k22;
    letSP(k3)=k32;
    letSP(n,k1) = kp-2*pp-SP(n,k2);

    exvector amps;
    ex SpinColor = SpinProj(IO::Out, 1, p, p, m, m, LI(-11), -4, -6) * ColorProj(-4,-6);
    for(auto ampR : amp2) {
        auto tmpR = MatrixContract(ampR * SpinColor);
        tmpR = tmpR.subs(GAS(5)==Gamma5("gr"));
        tmpR = tmpR.conjugate();
        tmpR = IndexCC(tmpR);
        for(auto ampL : amp2) {
            auto tmpL = MatrixContract(ampL * SpinColor);
            tmpL = tmpL.subs(GAS(5)==Gamma5("gl"));
            auto tmpRL = tmpL * tmpR * 
                (-SP(LI(-11), RLI(-11))+SP(p,LI(-11))*SP(p,RLI(-11))/SP(p)) * 
                QuarkSum(-8, k1, m) * GluonSum(-10) *
                SP(TI(-2),RTI(-2)) * Matrix(GAS(n), DI(-2), RDI(-2));
            tmpRL = MatrixContract(tmpRL);
            amps.push_back(tmpRL);
        }
    }

    exvector fc_amps;
    if(!file_exists("FC.gar")) {
        fc_amps = GiNaC_Parallel(amps.size(), [amps](int idx)->ex {
            auto tmp = amps[idx];
            tmp = form(tmp, false);
            tmp = tmp.subs(lst{NF==3,NA==8});
            return tmp;
        }, "FORM");
        garWrite("FC.gar", exvec2lst(fc_amps));
    } else {
        fc_amps = lst2exvec(ex_to<lst>(garRead("FC.gar")));
    }
    
    ex ap_res = 0;
    if(!file_exists("AP.gar")) {
        lst vps = {q1,q2,k1,k2,k3,p,n};
        exmap SP2sp;
        for(auto vp1 : vps) {
            for(auto vp2 : vps) {
                SP2sp[SP(vp1,vp2)]=sp(vp1,vp2);
            }
        }
    
        auto ltps = lst {k1,k2};
        auto exps = lst{p,n};
        auto ap_amps = GiNaC_Parallel(fc_amps.size(), [&](int idx)->ex {
            auto amp = fc_amps[idx];
            amp = Apart(amp, ltps, exps); 
            amp = ApartIR2F(amp);
            amp = amp.subs(SP2sp);
            return amp;
        }, "AP");
        for(auto item : ap_amps) ap_res += item;
        ap_res = mma_collect(ap_res, F(w1,w2), true, true);
        garWrite("AP.gar", ap_res);
    } else {
        ap_res = garRead("AP.gar");
    }
            
    auto res = NLoops(ap_res, 0, 2);
    
    cout << res << endl<< "--" << endl << VEResult(res) << endl;
    
    return 0;
}

ex NLoops(const ex & ap_res, int loop, int tloop) {

    ex S=1;
    bool use_eps = false;
    
    clearSP();
    
    ex result = ap_res;
    result = result.subs(lst{gs==1, D==4-2*ep});
    if(!is_a<add>(result)) result = lst{ result };
    cout << "Total F-terms: " << result.nops() << endl;
 
if(false){
system("mkdir -p cm");
for(int j=0; j<result.nops(); j++) {
    auto item = result.op(j);
    ex vv = item.subs(lst{coCF(w)==1, coVF(w)==w});
    ex cc = item.subs(lst{coCF(w)==w, coVF(w)==1});
    ofstream out("cm/"+to_string(j+1)+".txt");
    out << "{" << vv.op(0) << "," << vv.op(1) << "," << cc << "}" << endl;
    out.close();
}
exit(0);
}

    ex NRes = 0;
    Verbose = 100;
    ParallelProcess = -1;
    auto vecNRes = GiNaC_Parallel(result.nops(), [&](int idx)->ex {
        Verbose = 0;
        ParallelProcess = 0;
        auto item = result.op(idx);
        if(is_zero(item)) return 0;
        else if(item.nops()!=2) {
            cout << item << endl;
            return NaN;
        }
        
        ex kp = 1;
        ex m = 1;
        ex m2=m*m;
        ex zz = ex(45)/100;
        
        ex k1p = z(1)*(1-zz)*kp;
        ex k1m = (k12+pow(K1,2))/(2*k1p);
        ex k2p = z(2)*(1-zz)*kp;
        ex k2m = (k22+pow(K2,2))/(2*k2p);
        ex k3p = z(3)*(1-zz)*kp;
        ex k3m = (k32+pow(K3,2))/(2*k3p);
        
        if(tloop<1) zz = zv = 1;
        ex pp = zz/2*kp;
        ex pm = m2/(2*pp);
        ex ppv = pp;
        
        item = item.subs(lst{
            Symbol("m")==mv, 
            Symbol("z")==zv, 
            Symbol("kp")==kp,
            Symbol("pp")==pp
        });
        ex vv = item.subs(lst{coCF(w)==1, coVF(w)==w});
        ex cc = item.subs(lst{coCF(w)==w, coVF(w)==1});

        FeynmanParameter fp;
        fp.LoopMomenta = lst{ };
        fp.tLoopMomenta = lst{ };
        if(loop>0) fp.LoopMomenta.append(q1.name);
        if(loop>1) fp.LoopMomenta.append(q2.name);
        if(tloop>0) fp.tLoopMomenta.append(K1);
        if(tloop>1) fp.tLoopMomenta.append(K2);
        if(tloop>2) fp.tLoopMomenta.append(K3);
        
        fp.lReplacements[sp(p,p)] = m2;
        fp.lReplacements[sp(n,n)] = 0;
        fp.lReplacements[sp(p,n)] = pp;
        fp.lReplacements[sp(k1,k1)] = k12;
        fp.lReplacements[sp(k2,k2)] = k22;
        fp.lReplacements[sp(k3,k3)] = k32;
        fp.lReplacements[sp(n,k1)] = k1p;
        fp.lReplacements[sp(n,k2)] = k2p;
        fp.lReplacements[sp(n,k3)] = k3p;
        fp.lReplacements[sp(p,k1)] = k1p*pm+k1m*pp;
        fp.lReplacements[sp(p,k2)] = k2p*pm+k2m*pp;
        fp.lReplacements[sp(p,k3)] = k3p*pm+k3m*pp;
        fp.lReplacements[sp(k1,k2)] = k1p*k2m+k2p*k1m-K1*K2;
        fp.lReplacements[sp(k1,k3)] = k1p*k3m+k3p*k1m-K1*K3;
        fp.lReplacements[sp(k2,k3)] = k3p*k2m+k2p*k3m-K3*K2;
                
        fp.tReplacements[sp(p,p)] = m2;
        fp.tReplacements[sp(n,n)] = 0;
        fp.tReplacements[sp(p,n)] = pp;
        
        if(tloop>0) fp.nReplacements[z(w)] = ex(1)/tloop;
        
        fp.Prefactor = pow(2*Pi, (2*ep-4)*loop) * cc;
        fp.Propagators = ex_to<lst>(vv.op(0));
        fp.Exponents = ex_to<lst>(vv.op(1));
        
        ex NCS = pow(zz,(1-2*ep))/(8*3*kp*Pi);
        ex M = 2*m;
        ex nts = fp.tLoopMomenta.nops();
        ex psFactor;
        if(nts>0) psFactor = NCS*M/S*(-((pow(2,(2-nts-(3-2*ep)*nts))*pow(Pi,(1-(3-2*ep)*nts)))/(kp*(-1+zz))));
        else psFactor = NCS*M*4*Pi/kp;
        fp.Prefactor = fp.Prefactor * psFactor;
        
        SecDec work;
        work.CheckEnd = true;
        work.epN = epN;
        work.Initialize(fp);
        
        if(tloop>0)
        for(auto &fe : work.FunExp) {
            int xn = fe.op(2).op(0).nops();
            exmap z2x;
            lst zs;
            ex zFactor = 1;
            for(int i=1; i<=tloop; i++) {
                auto xi = x(xn+i-1);
                if(fe.has(xi)) {
                    cout << "fe = " << fe << endl;
                    throw Error("Check failed: fe has xi.");
                }
                z2x[z(i)] = xi;
                zs.append(xi);
                zFactor /= xi;
            }
            let_op_append(fe, 2, zs);
            
            fe.let_op(0).let_op(0) = fe.op(0).op(0) * zFactor;
            fe.let_op(0) = subs(fe.op(0), z2x);
            
            auto tmp = Factor(fe.op(0).op(0));

            if(tmp.has(x(w)) && is_a<mul>(tmp)) {
                ex rem = 1;
                for(auto item : tmp) {
                    if(!item.has(x(w))) {
                        rem *= item;
                    } else if(item.match(pow(w0, w1))) {
                        let_op_append(fe, 0, item.op(0));
                        let_op_append(fe, 1, item.op(1));
                    } else {
                        if(!item.is_polynomial(zs)) {
                            cout << item << endl;
                            throw Error("Factor failed, item is not a polynormial w.r.t zs");
                        }
                        let_op_append(fe, 0, item);
                        let_op_append(fe, 1, 1);
                    }
                }
                fe.let_op(0).let_op(0) = rem;
            } else if(tmp.has(x(w)) && tmp.match(pow(w0, w1))) {
                fe.let_op(0).let_op(0) = 1;
                let_op_append(fe, 0, tmp.op(0));
                let_op_append(fe, 1, tmp.op(1));
            } 
            
            if(use_eps) {
                for(auto zi : zs) {
                    let_op_append(fe, 0, zi);
                    let_op_append(fe, 1, eps);
                }
            }

        }
                
        work.Normalizes();
        work.XReOrders();
        work.Normalizes();
        work.Scalelesses();
        work.RemoveDeltas();
        work.SDPrepares();
        work.EpsEpExpands();
        work.CIPrepares();
        work.Contours();
                
        work.EpsAbs = 5E-4;
        work.RunPTS = 1000000;
        work.RunMAX = 10;
        work.LambdaSplit = 5;
        work.TryPTS = 100000;
        work.CTryLeft = 3;
        work.CTryRight = 3;
        work.CTry = 1;
        work.CTryRightRatio = 5;
        work.Integrates();
                
        return work.ResultError;
    }, "SD");
    
    for(auto item : vecNRes) NRes += item;
    NRes = VESimplify(NRes);
    
    return NRes;
}

extern string model;
lst Amps() {
    
    Symbol ep("e"), np("n",true,false), nbar("nbar"); // note the ep, np
    Symbol Q("Q"), Qbar("Qbar"), q("q"), qbar("qbar"), g("g"), gh("gh"), ghbar("ghbar");
    Vector vn("n");
    
    Qgraf::Process proc;
    proc.Model = model;
    proc.In = "e[pe]";
    proc.Out = "n[pn],Q[p1],Qbar[p2],Q[k1],g[k2]";
    proc.Options = "onshell";
    proc.LoopPrefix = "q";
    
    symtab st;
    st["pe"] = p1+p2+k1+k2;
    st["pn"] = 0;
    st["p1"] = p1;
    st["p2"] = p2;
    st["q1"] = q1;
    st["q2"] = q2;
    st["k1"] = k1;
    st["k2"] = k2;
    st["k3"] = k3;
    
    proc.Loops = 0;
    //proc.Others.clear();
    auto amp2 = proc.Amplitudes(st);
        
    InOutTeX[-1]=" ";
    InOutTeX[-2]=" ";
    InOutTeX[-4]="$Q$";
    InOutTeX[-6]="$\\bar{Q}$";
    InOutTeX[-8]="${\\color{blue}Q}$";
    InOutTeX[-10]="$g$";
    
    // color singlet
    lst chk = lst{-4,-6,g};
    chk.sort();
    auto filter = [&](lst & ampi)->void {
        auto tmp = ampi;
        ampi.remove_all();
        for(auto amp : tmp) {
            if(HasLoop(amp, lst{nbar, n})) continue;
            auto cps = ShrinkCut(amp, lst{g, g}, 1);
            for(auto cpi : cps) {
                for(auto cpii : cpi) {
                    if(cpii==chk) goto amp_done;
                }
            }
            ampi.append(amp);
            amp_done: ;
        }
    };
    cout << "Filter:" << amp2.nops() << " :> ";
    filter(amp2);
    cout << amp2.nops() << endl;
    
    // generate diagrams to pdf
    if(false) {    
        DrawPDF(amp2, "amp2.pdf");
        exit(0);
    }

    // feynman rules here
    auto map = MapFunction([&](const ex &e, MapFunction &self)->ex {
        if(isFunction(e,"OutField") || isFunction(e,"InField")) return 1;
        else if(isFunction(e, "Propagator")) {
            if(e.op(0).op(0)==q) {
                return QuarkPropagator(e, 0);
            } else if(e.op(0).op(0)==Q) {
                return QuarkPropagator(e, m);
            } else if(e.op(0).op(0)==g) {
                return GluonPropagator(e);
            } else if(e.op(0).op(0)==gh) {
                return GhostPropagator(e);
            } else if(e.op(0).op(0)==np) { 
                auto ret = eikonalPropagator(e, n, 1); 
                if(!ret.has(q1) && !ret.has(q2)) ret = ret.subs(iEpsilon==0);
                return ret;
            } 
        } else if(isFunction(e, "Vertex")) {
            if(e.nops()==3 && e.op(0).op(0)==Qbar && e.op(1).op(0)==ep && e.op(2).op(0)==nbar) {
                // Qbar-e-nbar
                return QuarkFFV(e, n);
            } else if(e.nops()==3 && e.op(0).op(0)==nbar && e.op(1).op(0)==np && e.op(2).op(0)==g) {
                // nbar-n-g
                return eikonalVertex(e, n, 1);
            } else if(e.nops()==3 && ((e.op(0).op(0)==qbar && e.op(1).op(0)==q) || (e.op(0).op(0)==Qbar && e.op(1).op(0)==Q)) && (e.op(2).op(0)==g) ) {
                // qbar-q-g or Qbar-Q-g
                return q2gVertex(e);
            } else if(e.nops()==3 && e.op(0).op(0)==ghbar && e.op(1).op(0)==gh) {
                // ghbar-gh-g
                return gh2gVertex(e);
            } else if(e.nops()==3 && e.op(0).op(0)==g && e.op(1).op(0)==g && e.op(2).op(0)==g) {
                // g^3
                return g3Vertex(e);
            } else if(e.nops()==4 && e.op(0).op(0)==g && e.op(1).op(0)==g && e.op(2).op(0)==g && e.op(3).op(0)==g) {
                // g^4
                return g4Vertex(e);
            }
        } else return e.map(self);
        return e;
    });

    amp2 = ex_to<lst>(map(amp2));
    return amp2;
}

std::string model = R"EOF(
[ model = 'quark FF Model' ]

%------------------------------------
% Propagators
%------------------------------------

% quark
[q, qbar, -]
[Q, Qbar, -]
[e, ebar, -, external]

% gluon and its ghost:
[gh, ghbar, -]
[g, g, +, notadpole]

% Eikonal
[n, nbar, +]
[n1, n1bar, +]
[n2, n2bar, +]

%------------------------------------
% Vertices
%------------------------------------

% QCD
[qbar, q, g; QCD='+1']
[Qbar, Q, g; QCD='+1']

[g, g, g, g; QCD='+2']
[g, g, g; QCD='+1']
[ghbar, gh, g; QCD='+1']

% Eikonal n_i
[nbar, n, g; QCD='+1']
[n1bar, n1, g; QCD='+1']
[n2bar, n2, g; QCD='+1']

% Eikonal external, connected required : qbar
[qbar, e, nbar; QCD='0']
[ebar, q, nbar; QCD='0']

% Eikonal external, connected required : Qbar
[Qbar, e, nbar; QCD='0']
[ebar, Q, nbar; QCD='0']
)EOF";
