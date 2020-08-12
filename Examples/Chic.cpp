#include "HepLib.h"

using namespace std;
using namespace GiNaC;
using namespace HepLib;
using namespace HepLib::FC;
using namespace HepLib::SD;
using namespace Quarkonium;
using namespace Qgraf;

Vector k("k"), q1("q1"), q2("q2"), p("p"), q("q");
Symbol m("m"), s("s"), t("t");
Symbol Chi0("Chi0"), Chi1("Chi1"), Chi2("Chi2");
Symbol zm("zm");
ex p1 = p+t*q;
ex p2 = p-t*q;
bool diag = false;

ex Amps();
int main() {
    Verbose = 100;
    
    ex amp_lst; // {amp0, amp1, amp2}
    if(file_exists("amps.gar")) {
        amp_lst = garRead("amps.gar");
    } else {
        amp_lst = Amps();
        garWrite("amps.gar", amp_lst);
    }
    
    lst tree_res = {6*(-2 + D)*pow(-1 + D, -1)*pow(m, -4)*pow(s - 12*pow(m, 2), 2)*pow(s - 4*pow(m, 2), -2), 12*(-2 + D)*s*pow(m, -4)*((-3 + D)*s + 4*pow(m, 2))*pow(s - 4*pow(m, 2), -2), 6*(-2 + D)*pow(-1 + D, -1)*pow(m, -4)*(24*s*pow(m, 2) + 16*(-8 - 3*D + 2*pow(D, 2))*pow(m, 4) + (-2 + D)*pow(s, 2))*pow(s - 4*pow(m, 2), -2)};
        
    // scalar product
    letSP(p)=m*m;
    letSP(p,q)=0;
    letSP(q,q)=0;
    letSP(k)=0;
    letSP(p,k)=s/4-m*m;
    
    lst exps = lst{p,k};
    ex ampRs, ampLs;
    lst loops;
    auto DoFA = [&ampLs,&ampRs,&loops,exps,tree_res](int idx)->ex {
        ex SpinColor = SpinProj(IO::Out, 1, p1, p2, m, m, LI(-11), -2, -4) * ColorProj(-2,-4);
        int iL = idx % ampLs.nops();
        int iR = idx / ampLs.nops();
        auto ampL = ampLs.op(iL);
        auto ampR = ampRs.op(iR);

        if(loops.nops()==2) {
            ampL = ampL.subs(zm==0);
            ampR = ampR.subs(zm==0);
        }
                
        auto tmpL = MatrixContract(ampL * SpinColor).subs(SP_map);
        tmpL = tmpL.subs(GAS(5)==Gamma5("gl"));
        tmpL = mma_series(tmpL,t,1);
        tmpL = LProj(tmpL.coeff(t,1),lst{p,q,LI(-13)},"lp");
        
        auto tmpR = MatrixContract(ampR * SpinColor).subs(SP_map);
        tmpR = tmpR.subs(GAS(5)==Gamma5("gr"));
        tmpR = mma_series(tmpR,t,1);
        tmpR = LProj(tmpR.coeff(t,1),lst{p,q,LI(-13)},"rp");
        tmpR = IndexCC(conjugate(tmpR));
        
        auto tmpRL = tmpL*tmpR*(-SP(LI(-1),RLI(-1)))*(-SP(LI(-6),RLI(-6)));
        
        tmpRL = Chi0/tree_res.op(0)*tmpRL*S1L1Sum(LI(-11), RLI(-11), LI(-13), RLI(-13), p, 0)
            +Chi1/tree_res.op(1)*tmpRL*S1L1Sum(LI(-11), RLI(-11), LI(-13), RLI(-13), p, 1)
            +Chi2/tree_res.op(2)*tmpRL*S1L1Sum(LI(-11), RLI(-11), LI(-13), RLI(-13), p, 2);
        
        tmpRL = MatrixContract(tmpRL).subs(SP_map);
        tmpRL = form(tmpRL, false);
        tmpRL = tmpRL.subs(lst{NF==3,NA==8});
        
        if(loops.nops()==0) {
            tmpRL = mma_series(tmpRL,zm,2);
        } else if(loops.nops()==1) {
            tmpRL = mma_series(tmpRL,zm,1);
        }
                
        if(loops.nops()>0) tmpRL = Apart(tmpRL, loops, exps); 
        
        return tmpRL;
    };
    
    if(!file_exists("amp0.gar")) {
        ampLs = amp_lst.op(0);
        ampRs = amp_lst.op(0);
        loops = lst{};
        
        auto res = GiNaC_Parallel(ampLs.nops()*ampRs.nops(), DoFA, "FA0");
                
        ex NRes = 0;
        for(auto item : res) NRes += item;
        auto cv_lst = mma_collect_lst(NRes,F(w1,w2));
        NRes = 0;
        for(auto item : cv_lst) {
            auto cc = item.op(0);
            auto vv = item.op(1);
            NRes += vv * fermat_normal(cc);
        }
        
        garWrite("amp0.gar", NRes);
        if(!is_a<add>(NRes)) NRes = lst{NRes};
        system("rm -rf cm0;mkdir -p cm0");
        for(int i=0; i<NRes.nops(); i++) {
            ex2file(NRes.op(i), "cm0/"+to_string(i+1)+".txt");
        }
        
        cout << endl;
    }
    
    if(!file_exists("amp1.gar")) {
        ampLs = amp_lst.op(1);
        ampRs = amp_lst.op(0);
        loops = lst{q1};
        
        auto res = GiNaC_Parallel(ampLs.nops()*ampRs.nops(), DoFA, "FA1");
        
        Apart2FIRE(res, lst{loops,exps});
        
        ex NRes = 0;
        for(auto item : res) NRes += item;
        auto cv_lst = mma_collect_lst(NRes,F(w1,w2));
        NRes = 0;
        for(auto item : cv_lst) {
            auto cc = item.op(0);
            auto vv = item.op(1);
            NRes += vv * fermat_normal(cc);
        }

        garWrite("amp1.gar", NRes);
        if(!is_a<add>(NRes)) NRes = lst{NRes};
        system("rm -rf cm1;mkdir -p cm1");
        for(int i=0; i<NRes.nops(); i++) {
            ex2file(NRes.op(i), "cm1/"+to_string(i+1)+".txt");
        }
        
        cout << endl;
    }
    
    if(!file_exists("amp2.gar")) {
        ampLs = amp_lst.op(2);
        ampRs = amp_lst.op(0);
        loops = lst{q1,q2};
        
        auto res = GiNaC_Parallel(ampLs.nops()*ampRs.nops(), DoFA, "FA2");
        
        Apart2FIRE(res, lst{loops,exps});
        
        ex NRes = 0;
        for(auto item : res) NRes += item;
        auto cv_lst = mma_collect_lst(NRes,F(w1,w2));
        NRes = 0;
        for(auto item : cv_lst) {
            auto cc = item.op(0);
            auto vv = item.op(1);
            NRes += vv * fermat_normal(cc);
        }
        
        garWrite("amp2.gar", NRes);
        if(!is_a<add>(NRes)) NRes = lst{NRes};
        system("rm -rf cm2;mkdir -p cm2");
        for(int i=0; i<NRes.nops(); i++) {
            ex2file(NRes.op(i), "cm2/"+to_string(i+1)+".txt");
        }
        
        cout << endl;
    }
                        
    return 0;
}


extern string model;
ex Amps() {
    
    Symbol e("e"), Q("Q"), Qbar("Qbar"), q("q"), qbar("qbar"), g("g"), gh("gh"), ghbar("ghbar");
    
    Qgraf::Process proc;
    proc.Model = model;
    proc.In = "e[pA]";
    proc.Out = "Q[p1],Qbar[p2],e[k]";
    proc.Options = "notadpole,onshell";
    proc.LoopPrefix = "q";
    
    symtab st;
    st["p1"] = p1;
    st["p2"] = p2;
    st["k"] = k;
    st["q1"] = q1;
    st["q2"] = q2;
    st["pA"] = p1+p2+k;
    
    proc.Loops = 0;
    proc.Others.clear();
    proc.Others.push_back("true=vsum[QCD,0,0]");
    auto amp0 = proc.Amplitudes(st);
    
    proc.Loops = 1;
    proc.Others.clear();
    proc.Others.push_back("true=vsum[QCD,2,2]");
    auto amp1 = proc.Amplitudes(st);
    
    proc.Loops = 2;
    proc.Others.clear();
    proc.Others.push_back("true=vsum[QCD,4,4]");
    auto amp2 = proc.Amplitudes(st);

    lst chk = lst{-2,-4,g};
    chk.sort();
    auto filter = [&](lst & ampi)->void {
        auto tmp = ampi;
        ampi.remove_all();
        for(auto amp : tmp) {
            auto cps = ShrinkCut(amp, lst{g,g}, 1);
            for(auto cpi : cps) {
                for(auto cpii : cpi) {
                    if(cpii==chk) goto amp_done;
                }
            }
            ampi.append(amp);
            amp_done: ;
        }
    };
    cout << "Process Filter: [" << amp0.nops() << "," << amp1.nops() << "," << amp2.nops() << "] :> ";
    
    filter(amp0);
    filter(amp1);
    filter(amp2);
    
    cout << "[" << amp0.nops() << "," << amp1.nops() << "," << amp2.nops() << "]" << endl;
    
    // generate diagrams to pdf
    if(diag) {    
        LineTeX[e] = "photon, edge label=$\\gamma$";
        InOutTeX[-1]="$\\gamma^*$";
        InOutTeX[-2]="$Q$";
        InOutTeX[-4]="$\\bar{Q}$";
        InOutTeX[-6]="$\\gamma$";
        DrawPDF(amp0, "amp0.pdf");
        DrawPDF(amp1, "amp1.pdf");
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
            }  
        } else if(isFunction(e, "Vertex")) {
            if(e.nops()==3 && ((e.op(0).op(0)==qbar && e.op(1).op(0)==q) || (e.op(0).op(0)==Qbar && e.op(1).op(0)==Q)) && (e.op(2).op(0)==g || e.op(2).op(0)==Symbol("e")) ) {
                // qbar-q-g or Qbar-Q-g or g -> e
                if(e.op(2).op(0)==g) return q2gVertex(e);
                else {
                    auto fi1 = e.op(0).op(1);
                    auto fi2 = e.op(1).op(1);
                    auto fi3 = e.op(2).op(1);
                    return Matrix(GAS(LI(fi3)),DI(fi1),DI(fi2)) * SP(TI(fi1),TI(fi2));
                }
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

    amp0 = ex_to<lst>(map(amp0).subs(m==m*(1+zm)));
    amp1 = ex_to<lst>(map(amp1).subs(m==m*(1+zm)));
    amp2 = ex_to<lst>(map(amp2));
    
    return lst{amp0, amp1, amp2};
}

std::string model = R"EOF(
[ model = 'eqcd Model' ]

%------------------------------------
% Propagators
%------------------------------------

% quark
[q, qbar, -]
[Q, Qbar, -]

% gluon and its ghost:
[gh, ghbar, -]
[g, g, +, notadpole]

% external
[e, e, +, external]

%------------------------------------
% Vertices
%------------------------------------

% QCD
[qbar, q, g; QCD='+1']
[Qbar, Q, g; QCD='+1']

[g, g, g, g; QCD='+2']
[g, g, g; QCD='+1']
[ghbar, gh, g; QCD='+1']

% external
[qbar, q, e; QCD='+0']
[Qbar, Q, e; QCD='+0']
)EOF";
