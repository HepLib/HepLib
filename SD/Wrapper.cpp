#include "SD.h"

namespace HepLib {

/*********************************************************/
// Wrapper Members
/*********************************************************/
ex Wrapper::MinF;
vector<ex> Wrapper::Xs;
ex Wrapper::IntF;
ex Wrapper::Lambda;
int Wrapper::ReIm = 1;

/*********************************************************/
// Wrapper Function: MinFunction
/*********************************************************/
void Wrapper::InitMinFunction(ex minF, vector<ex> xs, int ri) {
    MinF = minF;
    Xs = xs;
    ReIm = ri;
    Digits = 25;
}

dREAL Wrapper::MinFunction(int nvars, dREAL* x, dREAL* pl, dREAL *las) {
    //using MinF
    lst xRepl;
    for(int i=0; i<nvars; i++) xRepl.append(Xs[i]==numeric((double)x[i]));
    ex res = MinF.subs(xRepl).evalf();
    if(ReIm==1) res = res.real_part();
    else if(ReIm==2) res = res.imag_part();
    return ex_to<numeric>(res).to_double();
}

/*********************************************************/
// Wrapper Function: IntFunction
/*********************************************************/
void Wrapper::InitIntFunction(ex intF, vector<ex> xs) {
    IntF = intF;
    Xs = xs;
    Digits = 25;
}

int Wrapper::IntFunction(const unsigned int xn, const qREAL x[], const unsigned int yn, qREAL y[], const qREAL pl[], const qREAL las[]) {
    lst xRepl;
    for(int i=0; i<xn; i++) xRepl.append(Xs[i]==CppFormat::q2ex(x[i]));
    xRepl.append(z(0)==Lambda);
    xRepl.append(SD::iEpsilon==I*numeric("1E-50"));
    ex res = IntF.subs(xRepl).evalf();
    y[0] = CppFormat::ex2q(res.real_part());
    y[1] = CppFormat::ex2q(res.imag_part());
    return 0;
}


}
