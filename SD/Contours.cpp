/**
 * @file
 * @brief Functions for Contour Deformation
 */
 
#include "SD.h"
#include <math.h>
#include <cmath>

namespace HepLib::SD {
    
    /**
     * @brief Contours, note that here we need to provide the specific Parameter
     * @param key the to retrive data from CIPrepares call
     * @param pkey secondary key, final key will "key-pkey" to store Contours data
     */
    void SecDec::Contours(const string & key, const string & pkey) {
        if(IsZero) return;
        if(Minimizer==NULL) Minimizer = new HookeJeeves();

        if(key != "") {
            set_precision(NNDigits);
            ostringstream garfn;
            garfn << key << ".ci.gar";
            archive ar;
            ifstream in(garfn.str());
            in >> ar;
            in.close();
            auto c = ar.unarchive_ex("c");
            if(c!=19790923) throw Error("Contours: *.ci.gar error!");
            FT_N_XN = ex_to<lst>(ar.unarchive_ex("ftnxn"));
            reset_precision();
        }
        
        //change 2->1 from GiNaC 1.7.7
        if(FT_N_XN.nops()<1) return;
        if(Verbose > 0) cout << Color_HighLight << "  Contours @ " << now() << RESET << endl;
        
        exvector ftnxn_vec;
        //change 1->0 from GiNaC 1.7.7
        for(int i=0; i<FT_N_XN.nops(); i++) ftnxn_vec.push_back(FT_N_XN.op(i));
        
        auto pid = getpid();
        ostringstream cmd;
        cmd << "mkdir -p " << pid;
        if(!dir_exists(to_string(pid))) auto rc = system(cmd.str().c_str());
        
        ostringstream fsofn;
        if(key != "") {
            fsofn << key << "F.so";
        } else {
            fsofn << pid << "F.so";
        }
        void* module = dlopen(fsofn.str().c_str(), RTLD_NOW);
        if (module == nullptr) {
            module = dlopen(("./"+fsofn.str()).c_str(), RTLD_NOW);
            if (module == nullptr) {
                cerr << ErrColor << "Contours: could not open compiled module!" << RESET << endl;
                cout << "dlerror(): " << dlerror() << endl;
                exit(1);
            }
        }

        GiNaC_Parallel_NP["LAS"]=0;
        GiNaC_Parallel_RM["LAS"] = !Debug;
        auto res =
        GiNaC_Parallel(ftnxn_vec.size(), [&ftnxn_vec,this,module](int idx)->ex {
            // return lst{ ft_n, lst{lambda-i, lambda-max} }
            // with I*[lambda-i]*lambda, lambda < lambda-max
            // note that lambda sequence only matches to x sequence in F-term
            if(Verbose>5) cout << "     λ: " << idx << " / " << ftnxn_vec.size() << flush;
            auto ftnxn = ftnxn_vec[idx];
            int npara = -1;
            for(auto kv : Parameter) if(npara<kv.first) npara = kv.first;
            dREAL paras[npara+1];
            lst plRepl;
            for(auto kv : Parameter) {
                paras[kv.first] = (dREAL)ex2q(kv.second);
                plRepl.append(PL(kv.first)==kv.second);
            }
            
            auto ft = ftnxn.op(0);
            ft = ft.subs(plRepl);
            if(xSign(ft)>0 || CTLaMax<0) {
                if(Verbose>5) {
                    cout << "\r                                                    \r";
                    cout << "     λ: xSign>0 or CTLaMax<0, back to REAL mode!" << endl;
                }
                return lst{ ftnxn.op(1), 1979 }; // ft_id, las
            }
            
            int nvars = ex_to<numeric>(ftnxn.op(2)).to_int();
            ostringstream fname;
            
            if(CTMaxF) {
                MinimizeBase::FunctionType fp = NULL, fp2 = NULL;
                fname.clear();
                fname.str("");
                fname << "minF_"<<ftnxn.op(1);
                fp = (MinimizeBase::FunctionType)dlsym(module, fname.str().c_str());
                fname.clear();
                fname.str("");
                fname << "minFM_"<<ftnxn.op(1);
                fp2 = (MinimizeBase::FunctionType)dlsym(module, fname.str().c_str());
                if(fp==NULL || fp2==NULL) {
                    cerr << ErrColor << "Contours: fp==NULL/fp2==NULL" << RESET << endl;
                    cout << "dlerror(): " << dlerror() << endl;
                    exit(1);
                }
                
                dREAL maxf = Minimizer->FindMinimum(nvars, fp, paras);
                if(maxf>0) {
                    if(Verbose>5) {
                        cout << "\r                                                    \r";
                        cout << "     λ: F>0 Found, back to REAL mode!" << endl;
                    }
                    return lst{ ftnxn.op(1), 1979 }; // ft_id, las
                }
                
                maxf = Minimizer->FindMinimum(nvars, fp2, paras);
                if(maxf>0) {
                    if(Verbose>5) {
                        cout << "\r                                                    \r";
                        cout << "     λ: F<0 Found, back to REAL mode!" << endl;
                    }
                    return lst{ ftnxn.op(1), 1979 }; // ft_id, las
                }
            }
            
            dREAL nlas[nvars];
            dREAL max_df = -1, max_f;
            for(int i=0; i<nvars; i++) {
                MinimizeBase::FunctionType dfp = NULL;
                fname.clear();
                fname.str("");
                fname << "dirC_"<<ftnxn.op(1)<<"_"<<i;
                dfp = (MinimizeBase::FunctionType)dlsym(module, fname.str().c_str());
                if(dfp==NULL) {
                    cerr << ErrColor << "Contours: dfp==NULL" << RESET << endl;
                    cout << "dlerror(): " << dlerror() << endl;
                    exit(1);
                }
                
                dREAL maxdf = Minimizer->FindMinimum(nvars, dfp, paras);
                maxdf = -maxdf;
                nlas[i] = maxdf;
                if(max_df<maxdf) max_df = maxdf;
            }
            
            for(int i=0; i<nvars; i++) { // From FIESTA
                if(nlas[i] < 1) nlas[i] = 1;
                else nlas[i] = 1/nlas[i];
            }
            
            lst las;
            for(int i=0; i<nvars; i++) {
                las.append(q2ex(nlas[i]));
            }
            
            MinimizeBase::FunctionType fp;
            fname.clear();
            fname.str("");
            fname << "imgF_"<<ftnxn.op(1);
            fp = (MinimizeBase::FunctionType)dlsym(module, fname.str().c_str());
            if(fp==NULL) {
                cerr << ErrColor << "Contours: fp==NULL" << RESET << endl;
                cout << "dlerror(): " << dlerror() << endl;
                exit(1);
            }
            
            dREAL laBegin = 0, laEnd = CTLaMax, min;
            dREAL UB[nvars+1];
            for(int i=0; i<nvars+1; i++) UB[i] = 1;

            min = laEnd;
            while(true) {
                UB[nvars] = min;
                dREAL res = Minimizer->FindMinimum(nvars+1, fp, paras, nlas, UB, NULL, NULL, true, CTTryPTS, CTSavePTS);
                if(res < -1E-5) laEnd = min;
                else laBegin = min;
                
                if(laEnd < 1E-10) {
                    cout << "\r                                                    \r";
                    cout << "     λ: " << ErrColor << "too small lambda!" << RESET << endl;
                    break;
                }
                
                if(laEnd-laBegin <= 0.01*laEnd) break;
                min = (laBegin + laEnd) / 2.0;
            }
            min = laBegin;
        
            las.append(q2ex(min));
            if(Verbose>5) {
                cout << "\r                                                    \r";
                cout << "     λ: " << NN(las,3) << endl;
            }

            return lst{ ftnxn.op(1), las }; // ft_id, las
        
        }, "LAS");
        
        if(use_dlclose) dlclose(module);
        
        ostringstream garfn;
        if(key != "") {
            garfn << key;
            if(pkey != "") garfn << "-" << pkey;
            garfn << ".las.gar";
            lst gar_res;
            for(auto &item : res) gar_res.append(item);
            archive ar;
            ar.archive_ex(gar_res, "res");
            ar.archive_ex(19790923, "c");
            ofstream out(garfn.str());
            out << ar;
            out.close();
        } else {
            for(auto &item : res) LambdaMap[item.op(0)] = item.op(1);
        }
        
        cmd.clear();
        cmd.str("");
        cmd << "rm -rf " << pid;
        if(!Debug) auto rc = system(cmd.str().c_str());
    }

}
