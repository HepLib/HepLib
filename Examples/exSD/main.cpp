#include "SD.h"

using namespace HepLib;

int main(int argc, char** argv) {
    
    auto ep = SD::ep;
    auto iep = SD::iEpsilon;
    
    symbol k("k"), l("l"), r("r"), q("q"), q1("q1"), q2("q2");
    symbol p1("p1"), p2("p2"), p3("p3"), p4("p4"), p5("p5");
    symbol m("m"), s("s"), t("t"), u("u");
    
    ex m2=m*m;
    FeynmanParameter fp;
        
//    fp.LoopMomenta = lst{q1,q2};
//    fp.Propagators = lst{
//-2*p3*q2 - pow(q2,2),
//2*p3*q2 - pow(q2,2),
//2*p1*q1 - pow(m,2) - pow(q1,2),
//2*p1*q1 + 2*p3*q1 + 2*p1*q2 + 2*p3*q2 - 2*q1*q2 - s/4 + pow(m,2) - pow(q1,2) - pow(q2,2),
//4*p1*q1 + 4*p3*q1 - s + pow(m,2) - pow(q1,2),
//4*p1*q1 + 2*p3*q1 + 4*p1*q2 + 2*p3*q2 - 2*q1*q2 - s/2 - pow(m,2) - pow(q1,2) - pow(q2,2)
//    };

    
    fp.LoopMomenta = lst{q1,q2};
    fp.Propagators = lst {
-2*p3*q2 - pow(q2,2),
pow(m,2) - pow(q1,2),
2*p3*q2 - pow(q2,2),
2*p1*q1 + 2*p3*q1 + 2*p1*q2 + 2*p3*q2 - 2*q1*q2 - s/4 + pow(m,2) - pow(q1,2) - pow(q2,2),
4*p1*q1 + 4*p3*q1 - s + pow(m,2) - pow(q1,2),
4*p1*q1 + 2*p3*q1 + 4*p1*q2 + 2*p3*q2 - 2*q1*q2 - s/2 - pow(m,2) - pow(q1,2) -pow(q2,2)
    };

    fp.Exponents = lst{ 1, 1, 1, 1, 1, 2 };
    fp.lReplacements[p1*p3] = s/8-m2;
    fp.lReplacements[p3*p3] = m2;
    fp.lReplacements[p1*p1] = m2;
    fp.lReplacements[m] = 1;
    fp.lReplacements[s] = 39;






    
//    //FIESTA1 Example
//    fp.LoopMomenta = lst {k};
//    fp.Propagators = lst{ -pow(k,2),-pow(k + p1,2),-pow(k + p1 + p2,2),-pow(k + p1 + p2 + p4,2) };
//    fp.Exponents = lst{ 1, 1, 1, 1 };
//    fp.lReplacements[p1*p1] = 0;
//    fp.lReplacements[p2*p2] = 0;
//    fp.lReplacements[p4*p4] = 0;
//    fp.lReplacements[p1*p2] = -s/2;
//    fp.lReplacements[p2*p4] = -t/2;
//    fp.lReplacements[p1*p4] = (s+t)/2;
//    fp.lReplacements[s] = 3;
//    fp.lReplacements[t] = 1;
//    fp.Prefactor = pow(I*pow(Pi,2-ep)*exp(-ep*Euler), -1);
//    //1.33333/ep^2-0.732408/ep-4.38649


//
//    //SecDec Example - 7_epsprops_triangle_3L
//    fp.LoopMomenta = lst {k, r, q};
//    fp.Propagators = lst{ -pow(k,2),-pow(k + p1 + p2,2),-pow(-k + r,2),-pow(p1 + r,2),-pow(k - q,2),-pow(p1 + q,2) };
//    fp.Exponents = lst{ 1+3*ep,1,1,1,1,1 };
//    fp.lReplacements[p1*p1] = 0;
//    fp.lReplacements[p2*p2] = 0;
//    fp.lReplacements[p2*p1] = s/2;
//    fp.lReplacements[s] = 1;
//    fp.Prefactor = pow(I*pow(Pi,2-ep), -3) * pow(tgamma(1-ep), 3);
//    //125.32093561361496 + (5325.195869403933)*ep^2 + (1.83333333333175)*ep^(-2) + (889.9672394217711)*ep + (0.16666666666666616)*ep^(-3) + (18.123201259236488)*ep^(-1)
//


//    // SecDec Example - 5_pentagon_2L
//    symbol s12, s23, s34, s45, s51;
//    fp.LoopMomenta = lst {k, l};
//    fp.Propagators = lst{ pow(k + p1,2),pow(k + p1 + p2,2),pow(l + p1 + p2 + p3,2), pow(l + p1 + p2 + p3 + p4,2),pow(l + p1 + p2,2),pow(l,2),pow(k,2),pow(k - l,2) };
//    fp.Exponents = lst{ 1,1,1,1,  1,1,1,1 };
//    fp.lReplacements[p1*p1] = 0;
//    fp.lReplacements[p2*p2] = 0;
//    fp.lReplacements[p3*p3] = 0;
//    fp.lReplacements[p4*p4] = 0;
//    fp.lReplacements[p1*p2] = s12/2;
//    fp.lReplacements[p1*p3] = (s45-s12-s23)/2;
//    fp.lReplacements[p1*p4] = (s23-s51-s45)/2;
//    fp.lReplacements[p2*p3] = s23/2;
//    fp.lReplacements[p2*p4] = (-s23-s34+s51)/2;
//    fp.lReplacements[p3*p4] = s34/2;
//    fp.lReplacements[s34] = -3;
//    fp.lReplacements[s12] = numeric("-9.025");
//    fp.lReplacements[s23] = -2;
//    fp.lReplacements[s45] = -4;
//    fp.lReplacements[s51] = -5;
//    fp.Prefactor = pow(I*pow(Pi,2-ep), -2) / tgamma(2*(2+ep));


    
    fp.nReplacements[ep] = ex(1)/11;
    
    SD work;
    work.epN = 0;
    work.Verbose = 2;
    
    char *CFLAGS = getenv("SD_CFLAGS");
    work.CFLAGS = CFLAGS;
    
    work.Evaluate(fp);

    cout << work.ResultError << endl;

    
    return 0;
}
