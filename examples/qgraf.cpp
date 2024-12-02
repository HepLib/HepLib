#include "HepLib.h"

using namespace HepLib;
using namespace QGRAF;

extern string model;
int main() {

    Process proc;
    proc.Model = model;
    proc.In = "e[pA]";
    proc.Out = "Q[p1],Qbar[p2],e[k]";
    proc.Options = "notadpole,onshell";
    
    Vector k("k"), p("p"), h("h"), q1("q1"), q2("q2");
    Symbol e("e"), Q("Q"), Qbar("Qbar"), q("q"), qbar("qbar"), g("g"), gh("gh"), ghbar("ghbar");
    Symbol m("m"), s("s");
    
    symtab st;
    ex p1 = p + s*h;
    ex p2 = p - s*h;
    ex pA = p1+p2+k;
    st["p1"] = p1;
    st["p2"] = p2;
    st["k"] = k;
    st["pA"] = pA;
    st["q1"] = q1;
    st["q2"] = q2;
    
    // tree level
    proc.Loops = 0;
    proc.Others.clear();
    proc.Others.push_back("true=vsum[QCD,0,0]");
    auto amp0 = proc.Amplitudes(st);
    
    // NLO level
    proc.Loops = 1;
    proc.Others.clear();
    proc.Others.push_back("true=vsum[QCD,2,2]");
    auto amp1 = proc.Amplitudes(st);
    
    // NNLO level
    proc.Loops = 2;
    proc.Others.clear();
    proc.Others.push_back("true=vsum[QCD,4,4]");
    auto amp2 = proc.Amplitudes(st);
    
    // some adjustment for Feynman diagrams
    LineTeX[q] = "fermion, edge label=q";
    LineTeX[gh] = "ghost, edge label=$\\chi$";
    LineTeX[Q] = "fermion, edge label=Q";
    LineTeX[g] = "gluon, edge label=g";
    LineTeX[e] = "photon, edge label=$\\gamma$";
    LineTeX[Qbar] = "anti fermion, edge label=$\\bar{Q}$";
    
    // remove diagrams by color singlet
    lst chk = lst{-2,-4,g};
    chk.sort();
    
    auto filter = [g,chk](lst & ampi)->void {
        auto tmp = ampi;
        ampi.remove_all();
        for(auto amp : tmp) {
            auto cps = ShrinkCut(amp, lst{g, g}, 1);
            for(auto cpi : cps) {
                for(auto cpii : cpi) {
                    if(cpii==chk) {
                        goto amp_done;
                    }
                }
            }
            ampi.append(amp);
            amp_done: ;
        }
    };
    
    filter(amp0);
    filter(amp1);
    filter(amp2);
    
    // draw Feynman diagrams
    DrawPDF(amp1, "amp1.pdf");
    cout << "Feyman diagrams @ NLO are exported to amp1.pdf" << endl;
    //DrawPDF(amp2, "amp2.pdf");

    // apply the Feynman Rules: Propagator & Vertex
    auto map = MapFunction([&](const ex &e, MapFunction &self)->ex {
        if(isFunction(e,"OutField") || isFunction(e,"InField")) return 1;
        else if(isFunction(e, "Propagator")) {
            if(e.op(0).op(0)==q) {
                return QuarkPropagator(e);
            } else if(e.op(0).op(0)==Q) {
                return QuarkPropagator(e, m);
            } else if(e.op(0).op(0)==g) {
                return GluonPropagator(e);
            } else if(e.op(0).op(0)==gh) {
                return GluonGhostPropagator(e);
            }
        } else if(isFunction(e, "Vertex")) {
            if( (e.op(0).op(0)==qbar && e.op(1).op(0)==q) || (e.op(0).op(0)==Qbar && e.op(1).op(0)==Q) ) {
                if(e.op(2).op(0)==g) return QuarkGluonVertex(e);
                else return Matrix(GAS(LI(e.op(2).op(1))),DI(e.op(0).op(1)), DI(e.op(1).op(1))) * SP(TI(e.op(0).op(1)), TI(e.op(1).op(1)));
            } else if(e.op(0).op(0)==ghbar && e.op(1).op(0)==gh) {
                return GhostGluonVertex(e);
            } else if(e.nops()==3) {
                return Gluon3Vertex(e);
            } else if(e.nops()==4) {
                return Gluon4Vertex(e);
            }
        } else return e.map(self);
        return 0;
    });
    
    amp0 = ex_to<lst>(map(amp0));
    amp1 = ex_to<lst>(map(amp1));
    amp2 = ex_to<lst>(map(amp2));
    
    cout << "Amplitudes @ LO: " << endl;
    cout << amp0 << endl;
    cout << endl;
    
    return 0;
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
