/**
 * @file 
 * @brief QGRAF header file
 */
 
#pragma once

#include "HEP.h"
#include "QGRAF.h"

/**
 * @brief namespace for generating Feynman diagrams or amplitudes.
 */
namespace HepLib::QCD {

    using namespace std;
    using namespace GiNaC;
    using namespace HepLib;
    using namespace QGRAF;
    
    /**
     * @brief namespace for functions helpful in Heavy Quarkonium
     */
    namespace Quarkonium {
        enum IO {In, Out};
        ex SpinProj(IO io, int s, ex p, ex pb, ex m, ex e, ex mu);
        ex SpinProj(IO io, int s, ex p, ex pb, ex m, ex e, ex mb, ex eb, ex mu);
        ex SpinProj(IO io, int s, ex p, ex pb, ex m, ex e, ex mu, int i, int j);
        ex SpinProj(IO io, int s, ex p, ex pb, ex m, ex e, ex mb, ex eb, ex mu, int i, int j);
        ex ColorProj(int i, int j, Index a);
        ex ColorProj(int i, int j);
        
        ex S1L1Proj(ex si, ex qi, ex p);
        ex S1L1Proj(ex si, ex qi, ex mu, ex p);
        ex S1L1Proj(ex si, ex qi, ex mu1, ex mu2, ex p);
        ex S1L2Proj(ex si, ex qi1, ex qi2, ex mu, ex p);
        ex S1L2Proj(ex si, ex qi1, ex qi2, ex mu1, ex mu2, ex p);
        ex S1L1Sum(ex si, ex siR, ex qi, ex qiR, ex p, int J);
        
        ex LProj(const ex &expr_in, const lst &pqi, string prefix="lpj");
        
        ex Gamma5(const string pre, int start=1);
        
        ex DoPS(lst moms, ex amp, int si=-1, ex q2=1);
        ex nPS(int n, ex q2=1);
    }
    
    /**
     * @brief namespace for Renormalization Constant
     */
    namespace RC {
        ex Z2(string name, ex m, int loop=2);
        ex Z2(Symbol n, ex m, int loop=2);
        ex Zm(ex m, int loop=2);
        ex asBare(int loop=2);
        ex asLO();
        ex Zas(int loop);
    }
        
    /**
     * @brief namespace for Fragmentation Function
     */
    namespace FF {
        ex Anti5R(const ex & expr);
        ex zIntFactor(int mode, int tls, const ex SF=1);
        ex zIntegrate(const ex & c1, const ex & c0, const ex & n, const ex k2=0, const ex& p2=pow(Symbol("m"),2));
    }
}

