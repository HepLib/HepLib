#include "SD.h"

namespace HepLib {

/*-----------------------------------------------------*/
// Wrapper Members
/*-----------------------------------------------------*/
ex GWrapper::MinF;
vector<ex> GWrapper::Xs;
ex GWrapper::IntF;
ex GWrapper::Lambda;
int GWrapper::ReIm = 1;
int GWrapper::GiNaC_Digits = 35;

/*-----------------------------------------------------*/
// Wrapper Function: MinFunction
/*-----------------------------------------------------*/
void GWrapper::InitMinFunction(ex minF, vector<ex> xs, int ri) {
    MinF = minF;
    Xs = xs;
    ReIm = ri;
}

dREAL GWrapper::MinFunction(int nvars, dREAL* x, dREAL* pl, dREAL *las) {
    //using MinF
    lst xRepl;
    for(int i=0; i<nvars; i++) xRepl.append(Xs[i]==numeric((double)x[i]));
    ex res = MinF.subs(xRepl).evalf();
    if(ReIm==1) res = res.real_part();
    else if(ReIm==2) res = res.imag_part();
    return ex_to<numeric>(res).to_double();
}

/*-----------------------------------------------------*/
// Wrapper Function: IntFunction
/*-----------------------------------------------------*/
void GWrapper::InitIntFunction(ex intF, vector<ex> xs) {
    IntF = intF;
    Xs = xs;
}

int GWrapper::IntFunction(const unsigned int xn0, const qREAL x[], const unsigned int yn0, qREAL y[], const qREAL pl[], const qREAL las[]) {
    int digits = Digits;
    Digits = GiNaC_Digits;
    
    int npts = xn0/100;
    vector<int> vec_pts;
    for(int i=0; i<npts; i++) vec_pts.push_back(i);
    if(vec_pts.size()<1) vec_pts.push_back(0);
    
    int xn = xn0 % 100;
    int yn = yn0 % 100;
        
    vector<ex> y_res =
    GiNaC_Parallel(omp_get_num_procs(), {}, vec_pts, [&](auto pti, auto pr_id) {        
        lst xRepl;
        for(int i=0; i<xn; i++) xRepl.append(Xs[i]==CppFormat::q2ex(x[i+pti*xn]));
        xRepl.append(z(0)==Lambda);
        xRepl.append(SD::iEpsilon==I*numeric("1E-50"));
        ex res = IntF.subs(xRepl).evalf();
        return res;        
    }, "GiNaC-int", false, false);
    
    for(auto pti : vec_pts) {
        numeric res = ex_to<numeric>(y_res[pti]);
        y[0+pti*yn] = CppFormat::ex2q(res.real_part());
        y[1+pti*yn] = CppFormat::ex2q(res.imag_part());
    }
    
    Digits = digits;
    return 0;
}


}
