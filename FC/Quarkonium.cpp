#include "FC.h"

namespace HepLib::FC {

    //-----------------------------------------------------------
    // Quarkonium Class
    //-----------------------------------------------------------
    namespace Quarkonium {
    
        // Using non-Relativistic Normalization
        // To Use Relativistic Normalization, Just Times Sqrt[2e 2eb]
        // SProj0In[pc,pcb,m,e,mu], where pc=P/2+q,pcb=P/2-q,e^2=m^2-q^2,
        // In CM, P=(M=2e,0) & q=(0,q), P.q=0, P.q=0
        // SProj0In[pc,pcb,m,mb,e,eb,mu], where p=e/(e+eb)P+q,pb=eb/(e+eb)P-q,e^2=m^2-q^2,eb^2=mb^2-qb^2,
        // In CM, P=(M=e+eb,0) & q=(0,q), P.q=0
        // arXiv:1003.0061
            
        ex SpinProj(IO io, int s, ex p, ex pb, ex m, ex e, ex mu) {
            ex ret;
            if(io==In && s==1) ret = (GAS(p)+m*GAS(1))*(GAS(p+pb)+2*e*GAS(1))*GAS(mu)*(GAS(pb)-m*GAS(1));
            if(io==Out && s==1) ret = (GAS(pb)-m*GAS(1))*GAS(mu)*(GAS(p+pb)+2*e*GAS(1))*(GAS(p)+m*GAS(1));
            if(io==In && s==0) ret = -(GAS(p)+m*GAS(1))*(GAS(p+pb)+2*e*GAS(1))*GAS(5)*(GAS(pb)-m*GAS(1));
            if(io==Out && s==0) ret = (GAS(pb)-m*GAS(1))*GAS(5)*(GAS(p+pb)+2*e*GAS(1))*(GAS(p)+m*GAS(1));
            ret = ret/(4*sqrt(ex(2))*e*(e+m));
            ret = ret/(2*e);
            return ret;
        }
        
        ex SpinProj(IO io, int s, ex p, ex pb, ex m, ex e, ex mb, ex eb, ex mu) {
            ex ret;
            if(io==In && s==1) ret = (GAS(p)+m*GAS(1))*(GAS(p+pb)+(e+eb)*GAS(1))*GAS(mu)*(GAS(pb)-mb*GAS(1));
            if(io==Out && s==1) ret = (GAS(pb)-mb*GAS(1))*GAS(mu)*(GAS(p+pb)+(e+eb)*GAS(1))*(GAS(p)+m*GAS(1));
            if(io==In && s==0) ret = -(GAS(p)+m*GAS(1))*(GAS(p+pb)+(e+eb)*GAS(1))*GAS(5)*(GAS(pb)-mb*GAS(1));
            if(io==Out && s==0) ret = (GAS(pb)-mb*GAS(1))*GAS(5)*(GAS(p+pb)+(e+eb)*GAS(1))*(GAS(p)+m*GAS(1));
            ret = ret/(2*sqrt(ex(2))*(e+eb)*sqrt((e+m)*(eb+mb)));
            ret = ret/(sqrt(4*e*eb));
            return ret;
        }
        
        namespace {
            ex ITD(ex si, ex qi , ex p) {
                return -SP(si,qi)+SP(p,si)*SP(p,qi)/SP(p);
            }
        }
        
        ex SL1Proj(ex si, ex qi, ex p) {
            return ITD(si,qi,p)/sqrt(D-1);
        }
        
        ex SL1Proj(ex si, ex qi, ex mu, ex p) {
            return -I*LC(si,qi,mu,p)/sqrt(2*SP(p));
        }
        
        ex SL1Proj(ex si, ex qi, ex mu1, ex mu2, ex p) {
            return (ITD(si,mu1,p)*ITD(qi,mu2,p)+ITD(qi,mu1,p)*ITD(si,mu2,p))/2 - ITD(si,qi,p)*ITD(mu1,mu2,p)/(D-1);
        }
        
        ex SL2Proj(ex si, ex qi1, ex qi2, ex mu, ex p) {
            return sqrt((D-1)/(D+1))*((ITD(si,qi1,p)*ITD(qi2,mu,p)+ITD(si,qi2,p)*ITD(qi1,mu,p))/2- ITD(si,mu,p)*ITD(qi1,qi2,p)/(D-1));
        }
        
        ex SL2Proj(ex si, ex qi1, ex qi2, ex mu1, ex mu2, ex p) {
            return -I/(sqrt(2*(D-1)*SP(p))) * (SP(mu1,qi1)*LC(qi2,si,mu2,p)+SP(mu1,qi2)*LC(qi1,si,mu2,p));
        }
        
        ex SLSum(ex si, ex siR, ex qi, ex qiR, ex p, int L) {
            if(L==0) return ITD(si,qi,p)*ITD(siR,qiR,p)/(D-1);
            else if(L==1) return (ITD(si,siR,p)*ITD(qi,qiR,p)-ITD(si,qiR,p)*ITD(qi,siR,p))/2;
            else if(L==2) return (ITD(si,siR,p)*ITD(qi,qiR,p)+ITD(si,qiR,p)*ITD(qi,siR,p))/2-ITD(si,qi,p)*ITD(siR,qiR,p)/(D-1);
            return 0;
        }
        
        ex ColorProj() {
            return 1/sqrt(NF);
        }
        ex ColorProj(Index i, Index j, Index a) {
            return sqrt(ex(2)) * SUNT(i,j,a);
        }
        
    }
}

