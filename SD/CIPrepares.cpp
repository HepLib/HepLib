/**
 * @file
 * @brief Functions for Contour/Integration Preparation
 * @author F. Feng
 * @version 1.0.0
 * @date 2020-04-21
 */
 
#include "SD.h"
#include <math.h>
#include <cmath>

namespace HepLib::SD {

    /**
     * @brief Prepare for the Contours and Integrates calls
     * .so will be generated, for more detailed exported function, check below
     * @param key used to save intermidete result
     */
    void SecDec::CIPrepares(const string & key) {
        if(expResult.size()<1) {
            IsZero = true;
        }
        
        if(IsZero) return;
        if(CFLAGS=="") CFLAGS = getenv("SD_CFLAGS");
        
        if(Verbose > 1) cout << Color_HighLight << "  CIPrepares @ " << now() << RESET << endl;
        auto pid = getpid();
                
        auto resf =
        GiNaC_Parallel(expResult.size(), [&](int idx)->ex {
            // return lst{ kv.op(0), kv.op(1), ft};
            auto kv = expResult[idx];
            auto expr = kv.op(1);
            auto xs = get_xy_from(expr);
            if(xs.size()<1) {
                return lst{kv.op(0), kv.op(1), 1};
            }

            exset ftxset;
            expr.find(FTX(w1,w2), ftxset);
            ex ft;
            int ftxsize = -1;
            for(auto item : ftxset) {
                auto xys = get_xy_from(item.op(0));
                if((int)xys.size() > ftxsize) {
                    ft = item.op(0);
                    ftxsize = xys.size();
                }
            }
            
            bool need_contour_deformation = false;
            if(ft.has(x(w)) && !ft.has(PL(w))) {
                auto tmp = ft.subs(nReplacements).subs(lst{
                    CV(w1,w2)==w2, ep==ex(1)/111, eps==ex(1)/1111
                }).expand();
                if(is_a<add>(tmp)) {
                    for(auto item : tmp) {
                        if(!is_a<numeric>(item.subs(x(w)==1))) {
                            cerr << Color_Error << "CIPrepares: (!is_a<numeric>(item.subs(x(w)==1)))" << RESET << endl;
                            exit(1);
                        }
                        if(item.subs(x(w)==1) < 0) {
                            need_contour_deformation = true;
                            break;
                        }
                    }
                } else {
                    if(!is_a<numeric>(tmp.subs(x(w)==1))) {
                        cerr << Color_Error << "CIPrepares: (!is_a<numeric>(tmp.subs(x(w)==1)))" << RESET << endl;
                        exit(1);
                    }
                    if(tmp.subs(x(w)==1) < 0) need_contour_deformation = true;
                }
                if(!need_contour_deformation) ft = 1; //note the difference with SDPrepare
            } else if(!ft.has(x(w))) {
                ft = 1;
            }
            
            ft = collect_common_factors(ft);
            return lst{ kv.op(0), kv.op(1), ft};
            
        }, "CI-F", false);
        

    //============================================================================================================
        lst fts;
        for(auto item : resf) {
            if(item.op(2).has(x(w))) {
                fts.append(item.op(2));
            }
        }
        fts.sort();
        fts.unique();
        
        vector<ex> ftnvec;
        map<ex,int,ex_is_less> ftnmap;
        int ft_n = 1;
        FT_N_XN.remove_all();
        
        //deleted from GiNaC 1.7.7
        //FT_N_XN.append(lst{0, 0});
        
        for(auto item : fts) {
            ftnvec.push_back(lst{item, ft_n});
            ftnmap[item] = ft_n;
            FT_N_XN.append(lst{item, ft_n, get_xy_from(item).size()});
            ft_n++;
        }
        //ftnvec item: lst { ft, ft-id }
        
        vector<ex> res_vec;
        map<ex, ex, ex_is_less> cf_int;
        
        for(auto &item : resf) {
            auto ii = ex_to<lst>(item);
            if(ii.op(2)==1) {
                ii.append(-1);
            } else {
                int ft_n = ftnmap[item.op(2)];
                if(ft_n==0) {
                    cerr << Color_Error << "CIPrepares: ft_n==0, " << item.op(2) << RESET << endl;
                    exit(1);
                }
                ii.append(ft_n);
            }
            if(ii.op(0).is_zero() || ii.op(1).subs(FTX(w1,w2)==1).is_zero()) continue;
            if(!use_IBF) res_vec.push_back(ii);
            else {
                ex key = ii;
                key.let_op(1) = 1;
                cf_int[key] += ii.op(1);
            }
        }
        
        if(use_IBF) {
            for(auto kv : cf_int) {
                lst ii = ex_to<lst>(kv.first);
                ii.let_op(1) = kv.second;
                res_vec.push_back(ii);
            }
        }
        
        //res_vec item: lst { coeff, integrand, ft, ft-id }
        
        if(res_vec.size()<1) {
            IsZero = true;
            return;
        }
    //============================================================================================================

        // Prepare FT-lambda
        GiNaC_Parallel(ftnvec.size(), [&](int idx)->ex {
            // return nothing
            auto kv = ftnvec[idx];
            ex ft = kv.op(0);
            ex ft_n = kv.op(1);
            auto fxs = get_xy_from(ft);
            lst las;
            
            auto pls = get_pl_from(ft);
            int npls = pls.size()>0 ? ex_to<numeric>(pls[pls.size()-1].subs(lst{PL(w)==w})).to_int() : -1;
            lst plRepl;
            for(int i=0; i<npls+1; i++) {
                ostringstream pl;
                pl << "pl[" << i << "]";
                plRepl.append(PL(i) == symbol(pl.str()));
            }
            
            ex DFs[fxs.size()], DDFs[fxs.size()*fxs.size()];
            for(int i=0; i<fxs.size(); i++) {
                auto df = mma_diff(ft, fxs[i], 1, false);
                DFs[i] = collect_common_factors(df);
                ostringstream ilaos;
                ilaos << "ilas[" << i << "]";
                symbol ila(ilaos.str());
                for(int j=0; j<fxs.size(); j++) {
                    auto ddf = mma_diff(DFs[i], fxs[j], 1, false);
                    DDFs[fxs.size()*i+j] = collect_common_factors(ddf);
                }
            }

            if(!dir_exists(to_string(pid))) system(("mkdir -p "+to_string(pid)).c_str());
            ostringstream cppfn, sofn;
            cppfn << pid << "/F" << ft_n << ".cpp";
            sofn << pid << "/F" << ft_n << ".o";
            std::ofstream ofs;
            ofs.open(cppfn.str(), ios::out);
            if (!ofs) throw runtime_error("failed to open *.cpp file! (1)");
            
            lst cxRepl, czRepl;
            for (int i=0; i<fxs.size(); i++) {
                ostringstream sx, sz;
                sx << "x[" << i << "]";
                cxRepl.append(fxs[i] == symbol(sx.str()));
                sz << "z[" << i << "]";
                czRepl.append(fxs[i] == symbol(sz.str()));
            }

    /*----------------------------------------------*/
    ofs << R"EOF(
    #include <stddef.h>
    #include <stdlib.h>
    #include <math.h>
    #include <complex>
    #include <iostream>
    extern "C" {
    #include <quadmath.h>
    }
    #include "mpreal.h"

    using namespace std;

    #define Pi 3.1415926535897932384626433832795028841971693993751L
    #define Euler 0.57721566490153286060651209008240243104215933593992L

    typedef __float128 qREAL;
    typedef __complex128 qCOMPLEX;
    typedef long double dREAL;
    typedef complex<long double> dCOMPLEX;
    typedef mpfr::mpreal mpREAL;
    typedef complex<mpREAL> mpCOMPLEX;

    dREAL expt(dREAL a, dREAL b);
    dCOMPLEX expt(dCOMPLEX a, dREAL b);
    dREAL recip(dREAL a);
    dCOMPLEX recip(dCOMPLEX a);

    qREAL expt(qREAL a, qREAL b);
    qCOMPLEX expt(qCOMPLEX a, qREAL b);
    qREAL recip(qREAL a);
    qCOMPLEX recip(qCOMPLEX a);

    mpREAL expt(mpREAL a, mpREAL b);
    mpCOMPLEX expt(mpCOMPLEX a, mpREAL b);
    mpREAL recip(mpREAL a);
    mpCOMPLEX recip(mpCOMPLEX a);

    qREAL pow(qREAL x, qREAL y);
    qREAL log(qREAL x);
    qCOMPLEX pow(qCOMPLEX x, qREAL y);
    qCOMPLEX log(qCOMPLEX x);

    )EOF" << endl;
    /*----------------------------------------------*/
            auto cppL = CppFormat(ofs, "L");
            auto cppQ = CppFormat(ofs, "Q");
            auto cppMP = CppFormat(ofs, "MP");
            ft = Evalf(ft);
            
            // FL_fid
            ofs << "dREAL FL_" << ft_n << "(const dREAL* x, const dREAL *pl) {" << endl;
            ofs << "dREAL yy = ";
            ft.subs(plRepl).subs(cxRepl).print(cppL);
            ofs << ";" << endl;
            ofs << "return yy;" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // FQ_fid
            ofs << "qREAL FQ_" << ft_n << "(const qREAL* x, const qREAL *pl) {" << endl;
            ofs << "qREAL yy = ";
            ft.subs(plRepl).subs(cxRepl).print(cppQ);
            ofs << ";" << endl;
            ofs << "return yy;" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // FMP_fid
            ofs << "mpREAL FMP_" << ft_n << "(const mpREAL* x, const mpREAL *pl) {" << endl;
            ofs << "mpREAL yy = ";
            ft.subs(plRepl).subs(cxRepl).print(cppMP);
            ofs << ";" << endl;
            ofs << "return yy;" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // D's FL_fid
            ofs << "dREAL FL_" << ft_n << "(const int i, const dREAL* x, const dREAL *pl) {" << endl;
            for(int i=0; i<fxs.size(); i++) {
                ofs << "if("<<i<<"==i) return ";
                DFs[i].subs(plRepl).subs(cxRepl).print(cppL);
                ofs << ";" << endl;
            }
            ofs << "return 0;" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // D's FQ_fid
            ofs << "qREAL FQ_" << ft_n << "(const int i, const qREAL* x, const qREAL *pl) {" << endl;
            for(int i=0; i<fxs.size(); i++) {
                ofs << "if("<<i<<"==i) return ";
                DFs[i].subs(plRepl).subs(cxRepl).print(cppQ);
                ofs << ";" << endl;
            }
            ofs << "return 0;" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // D's FMP_fid
            if(use_MP || fxs.size()<3) {
                ofs << "mpREAL FMP_" << ft_n << "(const int i, const mpREAL* x, const mpREAL *pl) {" << endl;
                for(int i=0; i<fxs.size(); i++) {
                    ofs << "if("<<i<<"==i) return ";
                    DFs[i].subs(plRepl).subs(cxRepl).print(cppMP);
                    ofs << ";" << endl;
                }
                ofs << "return 0;" << endl;
                ofs << "}" << endl;
                ofs << endl;
            }
            
            // DD's FL_fid
            ofs << "dREAL FL_" << ft_n << "(const int i, const int j, const dREAL* x, const dREAL *pl) {" << endl;
            for(int i=0; i<fxs.size(); i++) {
            for(int j=0; j<fxs.size(); j++) {
                ofs << "if("<<i<<"==i && "<<j<<"==j) return ";
                DDFs[i*fxs.size()+j].subs(plRepl).subs(cxRepl).print(cppL);
                ofs << ";" << endl;
            }}
            ofs << "return 0;" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // DD's FQ_fid
            ofs << "qREAL FQ_" << ft_n << "(const int i, const int j, const qREAL* x, const qREAL *pl) {" << endl;
            for(int i=0; i<fxs.size(); i++) {
            for(int j=0; j<fxs.size(); j++) {
                ofs << "if("<<i<<"==i && "<<j<<"==j) return ";
                DDFs[i*fxs.size()+j].subs(plRepl).subs(cxRepl).print(cppQ);
                ofs << ";" << endl;
            }}
            ofs << "return 0;" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // DD's FMP_fid
            if(use_MP || fxs.size()<3) {
                ofs << "mpREAL FMP_" << ft_n << "(const int i, const int j, const mpREAL* x, const mpREAL *pl) {" << endl;
                for(int i=0; i<fxs.size(); i++) {
                for(int j=0; j<fxs.size(); j++) {
                    ofs << "if("<<i<<"==i && "<<j<<"==j) return ";
                    DDFs[i*fxs.size()+j].subs(plRepl).subs(cxRepl).print(cppMP);
                    ofs << ";" << endl;
                }}
                ofs << "return 0;" << endl;
                ofs << "}" << endl;
                ofs << endl;
            }
            
            // X2ZL_fid
            ofs << "void X2ZL_" << ft_n << "(const dREAL* x, dCOMPLEX* z, dCOMPLEX* r, dREAL* dff, const dREAL* pl, const dREAL* las) {" << endl;
            ofs << "int nfxs="<<fxs.size()<<";" << endl;
            ofs << "dCOMPLEX ilas[nfxs];" << endl;
            ofs << "for(int i=0; i<nfxs; i++) ilas[i] = complex<long double>(0.L, las[i]);" << endl;
            ofs << "dff[nfxs] = FL_"<<ft_n<<"(x,pl);" << endl;
            ofs << "for(int i=0; i<nfxs; i++) dff[i] = FL_"<<ft_n<<"(i,x,pl);" << endl;
            ofs << "dREAL fscale=0.L;" << endl;
            if(CT_method==1) {
                ofs << "for(int i=0; i<nfxs; i++) fscale += dff[i]*dff[i];" << endl;
                ofs << "fscale += dff[nfxs]*dff[nfxs];" << endl;
            } else {
                ofs << "fscale=1.L;" << endl;
            }
            ofs << "for(int i=0; i<nfxs; i++) r[i] = dff[i]*ilas[i]/fscale;" << endl;
            ofs << "for(int i=0; i<nfxs; i++) z[i] = x[i]-x[i]*(1.L-x[i])*r[i];" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // X2ZQ_fid
            ofs << "void X2ZQ_" << ft_n << "(const qREAL* x, qCOMPLEX* z, qCOMPLEX* r, qREAL* dff, const qREAL* pl, const qREAL* las) {" << endl;
            ofs << "int nfxs="<<fxs.size()<<";" << endl;
            ofs << "qCOMPLEX ilas[nfxs];" << endl;
            ofs << "for(int i=0; i<nfxs; i++) ilas[i] = las[i] * 1.Qi;" << endl;
            ofs << "dff[nfxs] = FQ_"<<ft_n<<"(x,pl);" << endl;
            ofs << "for(int i=0; i<nfxs; i++) dff[i] = FQ_"<<ft_n<<"(i,x,pl);" << endl;
            ofs << "qREAL fscale=0.Q;" << endl;
            if(CT_method==1) {
                ofs << "for(int i=0; i<nfxs; i++) fscale += dff[i]*dff[i];" << endl;
                ofs << "fscale += dff[nfxs]*dff[nfxs];" << endl;
            } else {
                ofs << "fscale=1.Q;" << endl;
            }
            ofs << "for(int i=0; i<nfxs; i++) r[i] = dff[i]*ilas[i]/fscale;" << endl;
            ofs << "for(int i=0; i<nfxs; i++) z[i] = x[i]-x[i]*(1.Q-x[i])*r[i];" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // X2ZMP_fid
            if(use_MP || fxs.size()<3) {
                ofs << "void X2ZMP_" << ft_n << "(const mpREAL* x, mpCOMPLEX* z, mpCOMPLEX* r, mpREAL* dff, const mpREAL* pl, const mpREAL* las) {" << endl;
                ofs << "int nfxs="<<fxs.size()<<";" << endl;
                ofs << "mpCOMPLEX ilas[nfxs];" << endl;
                ofs << "for(int i=0; i<nfxs; i++) ilas[i] = complex<mpREAL>(mpREAL(0), las[i]);" << endl;
                ofs << "dff[nfxs] = FMP_"<<ft_n<<"(x,pl);" << endl;
                ofs << "for(int i=0; i<nfxs; i++) dff[i] = FMP_"<<ft_n<<"(i,x,pl);" << endl;
                ofs << "mpREAL fscale=0;" << endl;
                if(CT_method==1) {
                    ofs << "for(int i=0; i<nfxs; i++) fscale += dff[i]*dff[i];" << endl;
                    ofs << "fscale += dff[nfxs]*dff[nfxs];" << endl;
                } else {
                    ofs << "fscale=1;" << endl;
                }
                ofs << "for(int i=0; i<nfxs; i++) r[i] = dff[i]*ilas[i]/fscale;" << endl;
                ofs << "for(int i=0; i<nfxs; i++) z[i] = x[i]-x[i]*(1-x[i])*r[i];" << endl;
                ofs << "}" << endl;
                ofs << endl;
            }
            
            // MatL_id
            ofs << "void MatL_"<<ft_n<<"(dCOMPLEX* mat, const dREAL* x, const dREAL* dff, const dREAL* pl, const dREAL* las) {" << endl;
            ofs << "int nfxs="<<fxs.size()<<";" << endl;
            ofs << "dCOMPLEX ilas[nfxs];" << endl;
            ofs << "for(int i=0; i<nfxs; i++) ilas[i] = complex<long double>(0.L, las[i]);" << endl;
            ofs << "dREAL fscale=0.L;" << endl;
            if(CT_method==1) {
                ofs << "for(int i=0; i<nfxs; i++) fscale += dff[i]*dff[i];" << endl;
                ofs << "fscale += dff[nfxs]*dff[nfxs];" << endl;
            } else {
                ofs << "fscale=1.L;" << endl;
            }
            ofs << "dREAL dff2[nfxs][nfxs];" << endl;
            ofs << "for(int i=0; i<nfxs; i++) {" << endl;
            ofs << "for(int j=0; j<nfxs; j++) {" << endl;
            ofs << "dff2[i][j] = FL_"<<ft_n<<"(i,j,x,pl);" << endl;
            ofs << "}}" << endl;
            ofs << "for(int i=0; i<nfxs; i++) {" << endl;
            ofs << "for(int j=0; j<nfxs; j++) {" << endl;
            ofs << "int ij = i*nfxs+j;" << endl;
            ofs << "if(i!=j) mat[ij] = 0;" << endl;
            ofs << "else mat[ij] = 1.L-(1.L-2.L*x[i])*dff[i]*ilas[i]/fscale;" << endl;
            ofs << "mat[ij] = mat[ij]-x[i]*(1.L-x[i])*dff2[i][j]*ilas[i]/fscale;" << endl;
            if(CT_method==1) {
                ofs << "dREAL dfscale = 0;" << endl;
                ofs << "for(int ii=0; ii<nfxs; ii++) dfscale -= 2*dff[ii]*dff2[ii][j];" << endl;
                ofs << "dfscale -= 2*dff[nfxs]*dff[j];" << endl;
                ofs << "dfscale = dfscale / (fscale*fscale);" << endl;
                ofs << "mat[ij] = mat[ij]-x[i]*(1.L-x[i])*dff[i]*ilas[i]*dfscale;" << endl;
            }
            ofs << "}}" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // MatQ_fid
            ofs << "void MatQ_"<<ft_n<<"(qCOMPLEX *mat, const qREAL* x, const qREAL* dff, const qREAL *pl, const qREAL *las) {" << endl;
            ofs << "int nfxs="<<fxs.size()<<";" << endl;
            ofs << "qCOMPLEX ilas[nfxs];" << endl;
            ofs << "for(int i=0; i<nfxs; i++) ilas[i] = las[i] * 1.Qi;" << endl;
            ofs << "qREAL fscale=0.Q;" << endl;
            if(CT_method==1) {
                ofs << "for(int i=0; i<nfxs; i++) fscale += dff[i]*dff[i];" << endl;
                ofs << "fscale += dff[nfxs]*dff[nfxs];" << endl;
            } else {
                ofs << "fscale=1.Q;" << endl;
            }
            ofs << "qREAL dff2[nfxs][nfxs];" << endl;
            ofs << "for(int i=0; i<nfxs; i++) {" << endl;
            ofs << "for(int j=0; j<nfxs; j++) {" << endl;
            ofs << "dff2[i][j] = FQ_"<<ft_n<<"(i,j,x,pl);" << endl;
            ofs << "}}" << endl;
            ofs << "for(int i=0; i<nfxs; i++) {" << endl;
            ofs << "for(int j=0; j<nfxs; j++) {" << endl;
            ofs << "int ij = i*nfxs+j;" << endl;
            ofs << "if(i!=j) mat[ij] = 0;" << endl;
            ofs << "else mat[ij] = 1.Q-(1.Q-2.Q*x[i])*dff[i]*ilas[i]/fscale;" << endl;
            ofs << "mat[ij] = mat[ij]-x[i]*(1.Q-x[i])*dff2[i][j]*ilas[i]/fscale;" << endl;
            if(CT_method==1) {
                ofs << "qREAL dfscale = 0;" << endl;
                ofs << "for(int ii=0; ii<nfxs; ii++) dfscale -= 2*dff[ii]*dff2[ii][j];" << endl;
                ofs << "dfscale -= 2*dff[nfxs]*dff[j];" << endl;
                ofs << "dfscale = dfscale / (fscale*fscale);" << endl;
                ofs << "mat[ij] = mat[ij]-x[i]*(1.Q-x[i])*dff[i]*ilas[i]*dfscale;" << endl;
            }
            ofs << "}}" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // MatMP_fid
            if(use_MP || fxs.size()<3) {
                ofs << "void MatMP_"<<ft_n<<"(mpCOMPLEX *mat, const mpREAL* x, const mpREAL* dff, const mpREAL *pl, const mpREAL *las) {" << endl;
                ofs << "int nfxs="<<fxs.size()<<";" << endl;
                ofs << "mpCOMPLEX ilas[nfxs];" << endl;
                ofs << "for(int i=0; i<nfxs; i++) ilas[i] = complex<mpREAL>(mpREAL(0), las[i]);" << endl;
                ofs << "mpREAL fscale=0;" << endl;
                if(CT_method==1) {
                    ofs << "for(int i=0; i<nfxs; i++) fscale += dff[i]*dff[i];" << endl;
                    ofs << "fscale += dff[nfxs]*dff[nfxs];" << endl;
                } else {
                    ofs << "fscale=1;" << endl;
                }
                ofs << "mpREAL dff2[nfxs][nfxs];" << endl;
                ofs << "for(int i=0; i<nfxs; i++) {" << endl;
                ofs << "for(int j=0; j<nfxs; j++) {" << endl;
                ofs << "dff2[i][j] = FMP_"<<ft_n<<"(i,j,x,pl);" << endl;
                ofs << "}}" << endl;
                ofs << "for(int i=0; i<nfxs; i++) {" << endl;
                ofs << "for(int j=0; j<nfxs; j++) {" << endl;
                ofs << "int ij = i*nfxs+j;" << endl;
                ofs << "if(i!=j) mat[ij] = 0;" << endl;
                ofs << "else mat[ij] = mpREAL(1)-(1-2*x[i])*dff[i]*ilas[i]/fscale;" << endl;
                ofs << "mat[ij] = mat[ij]-x[i]*(1-x[i])*dff2[i][j]*ilas[i]/fscale;" << endl;
                if(CT_method==1) {
                    ofs << "mpREAL dfscale = 0;" << endl;
                    ofs << "for(int ii=0; ii<nfxs; ii++) dfscale -= 2*dff[ii]*dff2[ii][j];" << endl;
                    ofs << "dfscale -= 2*dff[nfxs]*dff[j];" << endl;
                    ofs << "dfscale = dfscale / (fscale*fscale);" << endl;
                    ofs << "mat[ij] = mat[ij]-x[i]*(1-x[i])*dff[i]*ilas[i]*dfscale;" << endl;
                }
                ofs << "}}" << endl;
                ofs << "}" << endl;
                ofs << endl;
            }
            
            // for Minimization of F(z)-image, x[xn-1] is the lambda
            ofs << "extern \"C\" " << endl;
            ofs << "dREAL imgF_"<<ft_n<<"(const int xn, const dREAL* x, const dREAL *pl, const dREAL *las_in) {" << endl;
            ofs << "dREAL las[xn-1];" <<endl;
            ofs << "for(int i=0; i<xn-1; i++) las[i] = las_in[i]*x[xn-1];" <<endl;
            ofs << "dCOMPLEX z[xn], r[xn];" << endl;
            ofs << "dREAL dff[xn+1];" << endl;
            ofs << "X2ZL_"<<ft_n<<"(x,z,r,dff,pl,las);" << endl;
            ofs << "dCOMPLEX zf = ";
            ft.subs(plRepl).subs(czRepl).print(cppL);
            ofs << ";" << endl;
            ofs << "return -zf.imag()/x[xn-1];" << endl; // find max image part, check with 0
            ofs << "}" << endl;
            ofs << endl;
            
            // for Minimization of F(z)-image, x[xn-1] is the lambda
            ofs << "extern \"C\" " << endl;
            ofs << "dREAL imgFQ_"<<ft_n<<"(const int xn, const dREAL* dx, const dREAL *dpl, const dREAL *dlas_in) {" << endl;
            ofs << "qREAL x[xn];" <<endl;
            ofs << "for(int i=0; i<xn; i++) x[i] = dx[i];" <<endl;
            ofs << "int npls = " << npls << ";" << endl;
            ofs << "qREAL pl[npls];" <<endl;
            ofs << "for(int i=0; i<npls; i++) pl[i] = dpl[i];" <<endl;
            ofs << "qREAL las[xn-1];" <<endl;
            ofs << "for(int i=0; i<xn-1; i++) las[i] = dlas_in[i]*x[xn-1];" <<endl;
            ofs << "qCOMPLEX z[xn], r[xn];" << endl;
            ofs << "qREAL dff[xn+1];" << endl;
            ofs << "X2ZQ_"<<ft_n<<"(x,z,r,dff,pl,las);" << endl;
            ofs << "qCOMPLEX zf = ";
            ft.subs(plRepl).subs(czRepl).print(cppQ);
            ofs << ";" << endl;
            ofs << "dREAL dret = -cimagq(zf/x[xn-1]);" << endl;
            ofs << "return dret;" << endl; // find max image part, check with 0
            ofs << "}" << endl;
            ofs << endl;
            
            // for Minimization of F
            ofs << "extern \"C\" " << endl;
            ofs << "dREAL minF_"<<ft_n<<"(const int xn, const dREAL* x, const dREAL *pl, const dREAL *las_in) {" << endl;
            ofs << "return FL_"<<ft_n<<"(x,pl);" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // for Minimization of DF-i
            for(int i=0; i<fxs.size(); i++) {
                ofs << "extern \"C\" " << endl;
                ofs << "dREAL dirF_"<<ft_n<<"_"<<i<<"(const int xn, const dREAL* x, const dREAL *pl, const dREAL *las) {" << endl;
                ofs << "dREAL yy = FL_"<<ft_n<<"("<<i<<", x, pl);" << endl;
                ofs << "return -fabs(yy);" << endl;
                ofs << "}" << endl;
                ofs << endl;
            }
            
            ostringstream cmd;
            cmd << cpp << " -fPIC -c " << CFLAGS << " -o " << sofn.str() << " " << cppfn.str();
            system(cmd.str().c_str());
            
            if(!file_exists(sofn.str().c_str())) {
                cmd.clear();
                cmd.str("");
                cmd << "cp " << sofn.str() << " . ";
                system(cmd.str().c_str());
            }
            
            return 0;
        
        }, "CI-C", false);
        
        bool hasF = (ftnvec.size()>0);
        if(hasF) {
            ostringstream sofn, cmd;
            if(key != "") {
                sofn << key << "F.so";
            } else {
                sofn << pid << "F.so";
            }
            cmd << cpp << " -rdynamic -fPIC -shared -lHepLib -lquadmath -lmpfr -lgmp " << CFLAGS << " -o " << sofn.str() << " " << pid << "/F*.o";
            system(cmd.str().c_str());
            
            cmd.clear();
            cmd.str("");
            cmd << "rm -rf " << pid;
            if(!debug) system(cmd.str().c_str());
        }


    //============================================================================================================
        // Compile the null.o
        if(true) {
            ostringstream cmd;
            if(!dir_exists(to_string(pid))) cmd << "mkdir -p " << pid;
            system(cmd.str().c_str());
            cmd.clear();
            cmd.str("");
            cmd << "echo ''>" << pid << "/null.cpp;";
            cmd << cpp << " -fPIC -c -o " << pid << "/null.o " << pid << "/null.cpp";
            system(cmd.str().c_str());
        }

        // Prepare Integrand
        auto res =
        GiNaC_Parallel(res_vec.size(), [&](int idx)->ex {
            // return lst{ no-x-result, xn, x-indepent prefactor, ft_n }
            // or     lst{ id(SD(D|Q)_id in .so), xn, x-indepent prefactor, ft_n }
            
            auto kvf = res_vec[idx];
            auto expr = kvf.op(1);
            auto xs = get_xy_from(expr);
            auto ft_n = kvf.op(3);
            bool hasF = (ft_n>0);
            
            if(xs.size()<1) {
                ostringstream cmd;
                cmd << "cp " << pid << "/null.o " << pid << "/" << idx << ".o";
                system(cmd.str().c_str());
                
                return lst{
                    expr.subs(FTX(w1,w2)==1).subs(iEpsilon==I*power(10,-50)),
                    xs.size(), kvf.op(0), -1
                };
            }
            
            auto ft = kvf.op(2);
            auto nft = ft.subs(lst{iEpsilon==0,x(w)==0});
            if(is_a<numeric>(nft) && nft>0) expr = expr.subs(iEpsilon==0);
            auto fxs = get_xy_from(ft);
            
            exset ftxset;
            expr.find(FTX(w1,w2), ftxset);
            lst ftxlst;
            for(auto it : ftxset) ftxlst.append(it);
            expr = mma_collect(expr, FTX(w1,w2));
            vector<pair<ex,ex>> ft_expr;
            for(auto item : ftxlst) {
                ft_expr.push_back(make_pair(item.op(1), expr.coeff(item)));
            }
            
            lst cxRepl, czRepl, czzRepl, plRepl;
            for (int i=0; i<fxs.size(); i++) {
                ostringstream xs, zs, zzs;
                xs << "x[" << i << "]";
                zs << "z[" << i << "]";
                zzs << "zz[" << i << "]";
                cxRepl.append(fxs[i] == symbol(xs.str()));
                czRepl.append(fxs[i] == symbol(zs.str()));
                czzRepl.append(fxs[i] == symbol(zzs.str()));
            }
            int count = fxs.size();
            for(auto xi : xs) {
                auto xii = xi.subs(czRepl);
                if(is_zero(xii-xi)) {
                    ostringstream xs, zs;
                    xs << "x[" << count << "]";
                    cxRepl.append(xi == symbol(xs.str()));
                    czRepl.append(xi == symbol(xs.str()));
                    count++;
                }
            }
            if(count!=xs.size()) {
                cerr << Color_Error << "CIPrepares: (count!=xs.size())" << RESET << endl;
                exit(1);
            }
            auto pls = get_pl_from(expr);
            int npls = pls.size()>0 ? ex_to<numeric>(pls[pls.size()-1].subs(lst{PL(w)==w})).to_int() : -1;
            for(int i=0; i<npls+1; i++) {
                ostringstream pl;
                pl << "pl[" << i << "]";
                plRepl.append(PL(i) == symbol(pl.str()));
            }
            
            if(!dir_exists(to_string(pid))) system(("mkdir -p "+to_string(pid)).c_str());
            ostringstream cppfn;
            cppfn << pid << "/" << idx << ".cpp";
            std::ofstream ofs;
            ofs.open(cppfn.str(), ios::out);
            if (!ofs) throw runtime_error("failed to open *.cpp file! (2)");

    /*----------------------------------------------*/
    ofs << R"EOF(
    #include <stddef.h>
    #include <stdlib.h>
    #include <math.h>
    #include <complex>
    #include <iostream>
    extern "C" {
    #include <quadmath.h>
    }
    #include "mpreal.h"

    using namespace std;

    typedef __float128 qREAL;
    typedef __complex128 qCOMPLEX;
    typedef long double dREAL;
    typedef complex<long double> dCOMPLEX;
    typedef mpfr::mpreal mpREAL;
    typedef complex<mpREAL> mpCOMPLEX;
    
    dCOMPLEX MatDetL(dCOMPLEX mat[], int n);
    qCOMPLEX MatDetQ(qCOMPLEX mat[], int n);
    mpCOMPLEX MatDetMP(mpCOMPLEX mat[], int n);

    extern int RCLog_NTry;
    dCOMPLEX RCLogL(dCOMPLEX ys[], int n);
    qCOMPLEX RCLogQ(qCOMPLEX ys[], int n);
    mpCOMPLEX RCLogMP(mpCOMPLEX ys[], int n);

    dREAL expt(dREAL a, dREAL b);
    dCOMPLEX expt(dCOMPLEX a, dREAL b);
    dREAL recip(dREAL a);
    dCOMPLEX recip(dCOMPLEX a);

    qREAL expt(qREAL a, qREAL b);
    qCOMPLEX expt(qCOMPLEX a, qREAL b);
    qREAL recip(qREAL a);
    qCOMPLEX recip(qCOMPLEX a);

    mpREAL expt(mpREAL a, mpREAL b);
    mpCOMPLEX expt(mpCOMPLEX a, mpREAL b);
    mpREAL recip(mpREAL a);
    mpCOMPLEX recip(mpCOMPLEX a);

    qREAL pow(qREAL x, qREAL y);
    qREAL log(qREAL x);
    qCOMPLEX pow(qCOMPLEX x, qREAL y);
    qCOMPLEX log(qCOMPLEX x);

    )EOF" << endl;
    /*----------------------------------------------*/

            if(hasF) {
                ofs << "qREAL FQ_"<<ft_n<<"(const qREAL*, const qREAL*);" << endl; // for FT only
                ofs << "dREAL FL_"<<ft_n<<"(const int, const dREAL*, const dREAL*);" << endl;
                ofs << "qREAL FQ_"<<ft_n<<"(const int, const qREAL*, const qREAL*);" << endl;
                if(use_MP || xs.size()<3) ofs << "qREAL FMP_"<<ft_n<<"(const int, const mpREAL*, const mpREAL*);" << endl;
                ofs << "dREAL FL_"<<ft_n<<"(const int, const int, const dREAL*, const dREAL*);" << endl;
                ofs << "qREAL FQ_"<<ft_n<<"(const int, const int, const qREAL*, const qREAL*);" << endl;
                if(use_MP || xs.size()<3) ofs << "qREAL FMP_"<<ft_n<<"(const int, const int, const mpREAL*, const mpREAL*);" << endl;
                ofs << "void X2ZL_"<<ft_n<<"(const dREAL*, dCOMPLEX*, dCOMPLEX*, dREAL*, const dREAL*, const dREAL*);" << endl;
                ofs << "void X2ZQ_"<<ft_n<<"(const qREAL*, qCOMPLEX*, qCOMPLEX*, qREAL*, const qREAL*, const qREAL*);" << endl;
                if(use_MP || xs.size()<3) ofs << "void X2ZMP_"<<ft_n<<"(const mpREAL*, mpCOMPLEX*, mpCOMPLEX*, mpREAL*, const mpREAL*, const mpREAL*);" << endl;
                ofs << "void MatL_"<<ft_n<<"(dCOMPLEX*, const dREAL*, const dREAL*, const dREAL*, const dREAL*);" << endl;
                ofs << "void MatQ_"<<ft_n<<"(qCOMPLEX*, const qREAL*, const qREAL*, const qREAL*, const qREAL*);" << endl;
                if(use_MP || xs.size()<3) ofs << "void MatMP_"<<ft_n<<"(mpCOMPLEX*, const mpREAL*, const mpREAL*, const mpREAL*, const mpREAL*);" << endl;
                ofs << endl << endl;
            }

    /*----------------------------------------------*/
    // long double
    /*----------------------------------------------*/
    ofs << R"EOF(
    #undef Pi
    #undef Euler
    #undef iEpsilon
    #define Pi 3.1415926535897932384626433832795028841971693993751L
    #define Euler 0.57721566490153286060651209008240243104215933593992L
    #define iEpsilon complex<long double>(0, 1.E-50)
    )EOF" << endl;
    /*----------------------------------------------*/
            auto cppL =  CppFormat(ofs, "L");
            // alwasy export non-complex function
            if(true) {
                ofs << "extern \"C\" " << endl;
                ofs << "int SDD_"<<idx<<"(const unsigned int xn, const qREAL qx[], const unsigned int yn, qREAL y[], const qREAL qpl[], const qREAL qlas[]) {" << endl;
                ofs << "dREAL x[xn], x0[xn];" << endl;
                ofs << "for(int i=0; i<xn; i++) x[i] = qx[i];" << endl;
                ofs << "dREAL pl["<<(npls<0 ? 1 : npls+1)<<"];" << endl;
                ofs << "for(int i=0; i<"<<(npls+1)<<"; i++) pl[i] = qpl[i];" << endl;
                
                if(SecDec::debug) {
                    auto tmp = expr.subs(FTX(w1,w2)==1).subs(cxRepl).subs(plRepl);
                    ofs << "//debug-int: " << tmp << endl;
                }
                
                auto intg = expr.subs(FTX(w1,w2)==1);
                bool hasF2 = intg.has(iEpsilon) || intg.has(I);
                cseParser cse;
                intg = cse.Parse(intg);
                if(hasF2) ofs << "dCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                else ofs << "dREAL "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                for(auto kv : cse.os()) {
                    ofs <<cse.oc<< "["<<kv.first<<"] = ";
                    Evalf(kv.second.subs(cxRepl).subs(plRepl)).print(cppL);
                    ofs << ";" << endl;
                }
                
                ofs << "dCOMPLEX yy = ";
                Evalf(intg.subs(cxRepl).subs(plRepl)).print(cppL);
                ofs << ";" << endl;
                
                ofs << "y[0] = yy.real();" << endl;
                ofs << "y[1] = yy.imag();" << endl;
                ofs << "return 0;" << endl;
                ofs << "}" << endl;
                ofs << endl;
            }
            
            if(hasF) {
                ofs << "extern \"C\" " << endl;
                ofs << "int CSDD_"<<idx<<"(const unsigned int xn, const qREAL qx[], const unsigned int yn, qREAL y[], const qREAL qpl[], const qREAL qlas[]) {" << endl;
                ofs << "dREAL x[xn], x0[xn];" << endl;
                ofs << "for(int i=0; i<xn; i++) x[i] = qx[i];" << endl;
                ofs << "dREAL pl["<<(npls<0 ? 1 : npls+1)<<"];" << endl;
                ofs << "for(int i=0; i<"<<(npls+1)<<"; i++) pl[i] = qpl[i];" << endl;
                
                ofs << "dCOMPLEX z[xn],zz[xn],r[xn];" << endl;
                ofs << "dREAL dff[xn+1];";
                ofs << "dCOMPLEX yy=0, ytmp, det;" << endl;
                ofs << "int ii, nfxs="<<fxs.size()<<";" << endl;
                ofs << "dREAL las[nfxs];" << endl;
                ofs << "for(int i=0; i<nfxs; i++) las[i] = qlas[i];" << endl;
                ofs << "dCOMPLEX mat[nfxs*nfxs];" << endl;
                for(auto &kv : ft_expr) {
                    ofs << "{" << endl;
                    lst xs0;
                    for(int ii=0; ii<fxs.size(); ii++) {
                        if(!kv.first.has(fxs[ii])) xs0.append(ii);
                    }
                    ofs << "for(int i=0; i<xn; i++) z[i] = x0[i] = x[i];" << endl;
                    for(auto x0i : xs0) ofs << "x0["<<x0i<<"]=0;" << endl;
                    ofs << "X2ZL_"<<ft_n<<"(x0,z,r,dff,pl,las);" << endl;
                    ofs << "MatL_"<<ft_n<<"(mat,x0,dff,pl,las);" << endl;
                    for(auto x0i : xs0) {
                        ofs << "ii = " << x0i << ";" << endl;
                        ofs << "z[ii] = x[ii]-x[ii]*(1.L-x[ii])*r[ii];" << endl;
                        ofs << "for(int j=0; j<nfxs;j++) mat[nfxs*ii+j] = 0;" << endl;
                        ofs << "for(int i=0; i<nfxs;i++) mat[nfxs*i+ii] = 0;" << endl;
                        ofs << "mat[ii*nfxs+ii] = 1.L-(1.L-2.L*x[ii])*r[ii];" << endl;
                    }
                    ofs  << "det = MatDetL(mat, nfxs);" << endl;
                    
                    ex intg = kv.second;
                    
                    if(use_RCLog) {
                        exset logs_set;
                        find(intg, log(w), logs_set);
                        if(logs_set.size()>0) {
                            lst logs;
                            for(auto item : logs_set) logs.append(item.op(0));
                            ofs << "int nlog = "<<logs.nops()<<";" << endl;
                            ofs << "dCOMPLEX CLog[nlog], CTry[nlog][RCLog_NTry];" << endl;
                            cseParser cse;
                            lst clogs = ex_to<lst>(cse.Parse(logs));
                            
                            ofs << "for(int ti=0; ti<RCLog_NTry; ti++) {" << endl;
                            ofs << "for(int i=0; i<xn; i++) zz[i] = complex<dREAL>(z[i].real(), (ti+1)*z[i].imag()/RCLog_NTry);" << endl;
                            ofs << "dCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                            for(auto kv : cse.os()) {
                                ofs <<cse.oc<< "["<<kv.first<<"] = ";
                                Evalf(kv.second.subs(czzRepl).subs(plRepl)).print(cppL);
                                ofs << ";" << endl;
                            }
                            
                            exmap log_subs;
                            for(int i=0; i<clogs.nops(); i++) {
                                ofs << "CTry["<<i<<"][ti] = ";
                                Evalf(clogs.op(i).subs(czzRepl).subs(plRepl)).print(cppL);
                                ofs << ";" << endl;
                                log_subs[log(logs.op(i))] = symbol("CLog["+to_string(i)+"]");
                            }
                            ofs << "}" << endl;
                            ofs << "for(int li=0; li<nlog; li++) CLog[li] = RCLogL(CTry[li],RCLog_NTry);" << endl;
                            
                            intg = intg.subs(log_subs);
                        }
                    }
                    
                    cseParser cse;
                    intg = cse.Parse(intg);
                    ofs << "dCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                    for(auto kv : cse.os()) {
                        ofs <<cse.oc<< "["<<kv.first<<"] = ";
                        Evalf(kv.second.subs(czRepl).subs(plRepl)).print(cppL);
                        ofs << ";" << endl;
                    }
                    
                    ofs << "ytmp = ";
                    Evalf(intg.subs(czRepl).subs(plRepl)).print(cppL);
                    ofs << ";" << endl;
                    ofs << "yy += det * ytmp;" << endl << endl;
                    ofs << "}" << endl;
                }
                
                ofs << "y[0] = yy.real();" << endl;
                ofs << "y[1] = yy.imag();" << endl;
                ofs << "return 0;" << endl;
                ofs << "}" << endl;
                ofs << endl;
            }
            
            
            
    /*----------------------------------------------*/
    // Quadruple
    /*----------------------------------------------*/
    ofs << R"EOF(
    #undef Pi
    #undef Euler
    #undef iEpsilon
    #define Pi 3.1415926535897932384626433832795028841971693993751Q
    #define Euler 0.57721566490153286060651209008240243104215933593992Q
    #define iEpsilon 1.E-50Qi
    )EOF" << endl;
    /*----------------------------------------------*/
            auto cppQ = CppFormat(ofs, "Q");
            
            // always export non-complex function
            if(true) {
                ofs << "extern \"C\" " << endl;
                ofs << "int SDQ_"<<idx<<"(const unsigned int xn, const qREAL x[], const int unsigned yn, qREAL y[], const qREAL pl[], const qREAL las[]) {" << endl;
                
                auto intg = expr.subs(FTX(w1,w2)==1);
                bool hasF2 = intg.has(iEpsilon) || intg.has(I);
                cseParser cse;
                intg = cse.Parse(intg);
                if(hasF2) ofs << "qCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                else ofs << "qREAL "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                for(auto kv : cse.os()) {
                    ofs <<cse.oc<< "["<<kv.first<<"] = ";
                    Evalf(kv.second.subs(cxRepl).subs(plRepl)).print(cppQ);
                    ofs << ";" << endl;
                }
                
                ofs << "qCOMPLEX yy = ";
                Evalf(intg.subs(cxRepl).subs(plRepl)).print(cppQ);
                ofs << ";" << endl;
                
                ofs << "y[0] = crealq(yy);" << endl;
                ofs << "y[1] = cimagq(yy);" << endl;
                ofs << "return 0;" << endl;
                ofs << "}" << endl;
                ofs << endl;
            }
            
            if(hasF) {
                ofs << "extern \"C\" " << endl;
                ofs << "int CSDQ_"<<idx<<"(const unsigned int xn, const qREAL x[], const int unsigned yn, qREAL y[], const qREAL pl[], const qREAL las[]) {" << endl;
                
                ofs << "qREAL x0[xn];" << endl;
                ofs << "qCOMPLEX z[xn],zz[xn],r[xn];" << endl;
                ofs << "qREAL dff[xn+1];";
                ofs << "qCOMPLEX yy=0, ytmp, det;;" << endl;
                ofs << "int ii, nfxs="<<fxs.size()<<";" << endl;
                ofs << "qCOMPLEX mat[nfxs*nfxs];" << endl;
                for(auto &kv : ft_expr) {
                    ofs << "{" << endl;
                    lst xs0;
                    for(int ii=0; ii<fxs.size(); ii++) {
                        if(!kv.first.has(fxs[ii])) xs0.append(ii);
                    }
                    ofs << "for(int i=0; i<xn; i++) z[i] = x0[i] = x[i];" << endl;
                    for(auto x0i : xs0) ofs << "x0["<<x0i<<"]=0;" << endl;
                    ofs << "X2ZQ_"<<ft_n<<"(x0,z,r,dff,pl,las);" << endl;
                    ofs << "MatQ_"<<ft_n<<"(mat,x0,dff,pl,las);" << endl;
                    for(auto x0i : xs0) {
                        ofs << "ii = " << x0i << ";" << endl;
                        ofs << "z[ii] = x[ii]-x[ii]*(1.Q-x[ii])*r[ii];" << endl;
                        ofs << "for(int j=0; j<nfxs;j++) mat[nfxs*ii+j] = 0;" << endl;
                        ofs << "for(int i=0; i<nfxs;i++) mat[nfxs*i+ii] = 0;" << endl;
                        ofs << "mat[ii*nfxs+ii] = 1.Q-(1.Q-2.Q*x[ii])*r[ii];" << endl;
                    }
                    ofs  << "det = MatDetQ(mat, nfxs);" << endl;
                    
                    ex intg = kv.second;
                    
                    if(use_RCLog) {
                        exset logs_set;
                        find(intg, log(w), logs_set);
                        if(logs_set.size()>0) {
                            lst logs;
                            for(auto item : logs_set) logs.append(item.op(0));
                            ofs << "int nlog = "<<logs.nops()<<";" << endl;
                            ofs << "qCOMPLEX CLog[nlog], CTry[nlog][RCLog_NTry];" << endl;
                            cseParser cse;
                            lst clogs = ex_to<lst>(cse.Parse(logs));
                            
                            ofs << "for(int ti=0; ti<RCLog_NTry; ti++) {" << endl;
                            ofs << "for(int i=0; i<xn; i++) zz[i] = crealq(z[i]) + 1.Qi * (ti+1)*cimagq(z[i])/RCLog_NTry;" << endl;
                            ofs << "qCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                            for(auto kv : cse.os()) {
                                ofs <<cse.oc<< "["<<kv.first<<"] = ";
                                Evalf(kv.second.subs(czzRepl).subs(plRepl)).print(cppQ);
                                ofs << ";" << endl;
                            }
                            
                            exmap log_subs;
                            for(int i=0; i<clogs.nops(); i++) {
                                ofs << "CTry["<<i<<"][ti] = ";
                                Evalf(clogs.op(i).subs(czzRepl).subs(plRepl)).print(cppQ);
                                ofs << ";" << endl;
                                log_subs[log(logs.op(i))] = symbol("CLog["+to_string(i)+"]");
                            }
                            ofs << "}" << endl;
                            ofs << "for(int li=0; li<nlog; li++) CLog[li] = RCLogQ(CTry[li],RCLog_NTry);" << endl;
                            
                            intg = intg.subs(log_subs);
                        }
                    }
                    
                    cseParser cse;
                    intg = cse.Parse(intg);
                    ofs << "qCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                    for(auto kv : cse.os()) {
                        ofs <<cse.oc<< "["<<kv.first<<"] = ";
                        Evalf(kv.second.subs(czRepl).subs(plRepl)).print(cppQ);
                        ofs << ";" << endl;
                    }
                    
                    ofs << "ytmp = ";
                    Evalf(intg.subs(czRepl).subs(plRepl)).print(cppQ);
                    ofs << ";" << endl;
                    ofs << "yy += det * ytmp;" << endl << endl;
                    ofs << "}" << endl;
                }
                ofs << "y[0] = crealq(yy);" << endl;
                ofs << "y[1] = cimagq(yy);" << endl;
                ofs << "return 0;" << endl;
                ofs << "}" << endl;
                ofs << endl;
                
                // Export the F-term, only Quadruple-type
                ofs << "extern \"C\" " << endl;
                ofs << "qREAL FT_"<<idx<<"(const qREAL x[], const qREAL pl[]) {" << endl;
                ofs << "qREAL yy = FQ_" << ft_n << "(x, pl);" << endl;
                ofs << "return yy;" << endl;
                ofs << "}" << endl;
                ofs << endl;
            }
            
            
    /*----------------------------------------------*/
    // Multiple Precision
    /*----------------------------------------------*/
    if(use_MP || xs.size()<3) {
    ofs << R"EOF(
    #undef Pi
    #undef Euler
    #undef iEpsilon
    #define Pi mpREAL("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068")
    #define Euler mpREAL("0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495")
    #define iEpsilon complex<mpREAL>(mpREAL(0), mpREAL(1.E-50))
    )EOF" << endl;
    /*----------------------------------------------*/
            auto cppMP =  CppFormat(ofs, "MP");
            
            // always export non-complex function
            if(true) {
                ofs << "extern \"C\" " << endl;
                ofs << "int SDMP_"<<idx<<"(const unsigned int xn, const qREAL qx[], const unsigned int yn, qREAL y[], const qREAL qpl[], const qREAL qlas[]) {" << endl;
                ofs << "mpREAL x[xn], x0[xn];" << endl;
                ofs << "for(int i=0; i<xn; i++) x[i] = mpREAL(qx[i]);" << endl;
                ofs << "mpREAL pl["<<(npls<0 ? 1 : npls+1)<<"];" << endl;
                ofs << "for(int i=0; i<"<<(npls+1)<<"; i++) pl[i] = mpREAL(qpl[i]);" << endl;
                
                auto intg = expr.subs(FTX(w1,w2)==1);
                bool hasF2 = intg.has(iEpsilon) || intg.has(I);
                cseParser cse;
                intg = cse.Parse(intg);
                if(hasF2) ofs << "mpCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                else ofs << "mpREAL "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                for(auto kv : cse.os()) {
                    ofs <<cse.oc<< "["<<kv.first<<"] = ";
                    Evalf(kv.second.subs(cxRepl).subs(plRepl)).print(cppMP);
                    ofs << ";" << endl;
                }
                
                ofs << "mpCOMPLEX yy = ";
                Evalf(intg.subs(cxRepl).subs(plRepl)).print(cppMP);
                ofs << ";" << endl;
                
                ofs << "y[0] = yy.real().toFloat128();" << endl;
                ofs << "y[1] = yy.imag().toFloat128();" << endl;
                ofs << "return 0;" << endl;
                ofs << "}" << endl;
                ofs << endl;
            }
            
            if(hasF) {
                ofs << "extern \"C\" " << endl;
                ofs << "int CSDMP_"<<idx<<"(const unsigned int xn, const qREAL qx[], const unsigned int yn, qREAL y[], const qREAL qpl[], const qREAL qlas[]) {" << endl;
                ofs << "mpREAL x[xn], x0[xn];" << endl;
                ofs << "for(int i=0; i<xn; i++) x[i] = mpREAL(qx[i]);" << endl;
                ofs << "mpREAL pl["<<(npls<0 ? 1 : npls+1)<<"];" << endl;
                ofs << "for(int i=0; i<"<<(npls+1)<<"; i++) pl[i] = mpREAL(qpl[i]);" << endl;
                
                ofs << "mpCOMPLEX z[xn],zz[xn],r[xn];" << endl;
                ofs << "mpREAL dff[xn+1];";
                ofs << "mpCOMPLEX yy=mpREAL(0), ytmp, det;" << endl;
                ofs << "int ii, nfxs="<<fxs.size()<<";" << endl;
                ofs << "mpREAL las[nfxs];" << endl;
                ofs << "for(int i=0; i<nfxs; i++) las[i] = mpREAL(qlas[i]);" << endl;
                ofs << "mpCOMPLEX mat[nfxs*nfxs];" << endl;
                for(auto &kv : ft_expr) {
                    ofs << "{" << endl;
                    lst xs0;
                    for(int ii=0; ii<fxs.size(); ii++) {
                        if(!kv.first.has(fxs[ii])) xs0.append(ii);
                    }
                    ofs << "for(int i=0; i<xn; i++) z[i] = x0[i] = x[i];" << endl;
                    for(auto x0i : xs0) ofs << "x0["<<x0i<<"]=0;" << endl;
                    ofs << "X2ZMP_"<<ft_n<<"(x0,z,r,dff,pl,las);" << endl;
                    ofs << "MatMP_"<<ft_n<<"(mat,x0,dff,pl,las);" << endl;
                    for(auto x0i : xs0) {
                        ofs << "ii = " << x0i << ";" << endl;
                        ofs << "z[ii] = x[ii]-x[ii]*(1-x[ii])*r[ii];" << endl;
                        ofs << "for(int j=0; j<nfxs;j++) mat[nfxs*ii+j] = mpREAL(0);" << endl;
                        ofs << "for(int i=0; i<nfxs;i++) mat[nfxs*i+ii] = mpREAL(0);" << endl;
                        ofs << "mat[ii*nfxs+ii] = mpREAL(1)-(1-2*x[ii])*r[ii];" << endl;
                    }
                    ofs  << "det = MatDetMP(mat, nfxs);" << endl;
                    
                    ex intg = kv.second;
                    
                    if(use_RCLog) {
                        exset logs_set;
                        find(intg, log(w), logs_set);
                        if(logs_set.size()>0) {
                            lst logs;
                            for(auto item : logs_set) logs.append(item.op(0));
                            ofs << "int nlog = "<<logs.nops()<<";" << endl;
                            ofs << "mpCOMPLEX CLog[nlog], CTry[nlog][RCLog_NTry];" << endl;
                            cseParser cse;
                            lst clogs = ex_to<lst>(cse.Parse(logs));
                            
                            ofs << "for(int ti=0; ti<RCLog_NTry; ti++) {" << endl;
                            ofs << "for(int i=0; i<xn; i++) zz[i] = complex<mpREAL>(z[i].real(), (ti+1)*z[i].imag()/RCLog_NTry);" << endl;
                            ofs << "mpCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                            for(auto kv : cse.os()) {
                                ofs <<cse.oc<< "["<<kv.first<<"] = ";
                                Evalf(kv.second.subs(czzRepl).subs(plRepl)).print(cppMP);
                                ofs << ";" << endl;
                            }
                            
                            exmap log_subs;
                            for(int i=0; i<clogs.nops(); i++) {
                                ofs << "CTry["<<i<<"][ti] = ";
                                Evalf(clogs.op(i).subs(czzRepl).subs(plRepl)).print(cppMP);
                                ofs << ";" << endl;
                                log_subs[log(logs.op(i))] = symbol("CLog["+to_string(i)+"]");
                            }
                            ofs << "}" << endl;
                            ofs << "for(int li=0; li<nlog; li++) CLog[li] = RCLogMP(CTry[li],RCLog_NTry);" << endl;
                            
                            intg = intg.subs(log_subs);
                        }
                    }
                    
                    cseParser cse;
                    intg = cse.Parse(intg);
                    ofs << "mpCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                    for(auto kv : cse.os()) {
                        ofs <<cse.oc<< "["<<kv.first<<"] = ";
                        Evalf(kv.second.subs(czRepl).subs(plRepl)).print(cppMP);
                        ofs << ";" << endl;
                    }
                    
                    ofs << "ytmp = ";
                    Evalf(intg.subs(czRepl).subs(plRepl)).print(cppMP);
                    ofs << ";" << endl;
                    ofs << "yy += det * ytmp;" << endl << endl;
                    ofs << "}" << endl;
                }
                
                ofs << "y[0] = yy.real().toFloat128();" << endl;
                ofs << "y[1] = yy.imag().toFloat128();" << endl;
                ofs << "return 0;" << endl;
                ofs << "}" << endl;
                ofs << endl;
            }
        // -----------------------
        } // end of if(use_MP)
        // -----------------------
            
            ofs.close();
            ostringstream ofn, cmd;
            ofn << pid << "/" << idx << ".o";
            cmd << cpp << " -fPIC " << CFLAGS << " -c -o " << ofn.str() << " " << cppfn.str();
            system(cmd.str().c_str());
            if(!debug) remove(cppfn.str().c_str());
            return lst{ idx, xs.size(), kvf.op(0), ft_n };
        }, "CI-I", false);
        

    //============================================================================================================

        ostringstream fsofn, sofn, garfn, cmd;
        if(key != "") {
            fsofn << key << "F.so";
            sofn << key << ".so";
            garfn << key << ".ci.gar";
            lst gar_res;
            for(auto &item : res) gar_res.append(item);
            archive ar;
            ar.archive_ex(gar_res, "res");
            ar.archive_ex(19790923, "c");
            ar.archive_ex(FT_N_XN, "ftnxn");
            ofstream out(garfn.str());
            out << ar;
            out.close();
            cmd.clear();
            cmd.str("");
            cmd << "rm -f " << key << "X*.so";
            system(cmd.str().c_str());
        } else {
            fsofn << pid << "F.so";
            sofn << pid << ".so";
            for(auto &item : res) ciResult.push_back(ex_to<lst>(item));
            cmd.clear();
            cmd.str("");
            cmd << "rm -f " << pid << "X*.so";
            system(cmd.str().c_str());
        }
        
        int res_size = res.size();
        if(GccLimit<100) GccLimit = 100;
        if(res_size>GccLimit) {
            cmd.clear();
            cmd.str("");
            cmd << cpp << " -rdynamic -fPIC -shared -lHepLib -lquadmath -lmpfr -lgmp " << CFLAGS;
            if(hasF) cmd << " " << fsofn.str();
            cmd << " -o " << sofn.str() << " $(seq -f '" << pid << "/%g.o' 0 " << (GccLimit-1) << ")";
            system(cmd.str().c_str());
            
            for(int n=1; true; n++) {
                int start = n*GccLimit;
                int end = (n+1)*GccLimit-1;
                if(end>res_size-1) end = res_size-1;
                sofn.clear();
                sofn.str("");
                if(key != "") sofn << key << "X" << n << ".so";
                else sofn << pid << "X" << n << ".so";
                cmd.clear();
                cmd.str("");
                cmd << cpp << " -rdynamic -fPIC -shared -lHepLib -lquadmath -lmpfr -lgmp " << CFLAGS;
                if(hasF) cmd << " " << fsofn.str();
                cmd << " -o " << sofn.str() << " $(seq -f '" << pid << "/%g.o' " << start << " " << end << ")";
                system(cmd.str().c_str());
                if(end>=res_size-1) break;
            }
        } else {
            cmd.clear();
            cmd.str("");
            cmd << cpp << " -rdynamic -fPIC -shared -lHepLib -lquadmath -lmpfr -lgmp " << CFLAGS;
            if(hasF) cmd << " " << fsofn.str();
            cmd << " -o " << sofn.str() << " " << pid << "/*.o";
            system(cmd.str().c_str());
        }
        cmd.clear();
        cmd.str("");
        cmd << "rm -rf " << pid;
        if(!debug) system(cmd.str().c_str());

    }
    

}
