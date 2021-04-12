/**
 * @file
 * @brief Functions for Numerical Integration
 */
 
#include "SD.h"
#include <math.h>
#include <cmath>

namespace HepLib::SD {

    /**
     * @brief ContinuousWRA, note that here we need to provide the specific Parameter
     * @param expr_in expression containing WRA
     * @param nc similar to NRCLog, the number of try
     * @return result with WRA removed
     */
    ex SecDec::ContinuousWRA(ex expr_in, int nc) {
        auto expr = expr_in;
        static symbol xwra("xwra");
        expr = xyz_pow_simplify(expr);
        exset pows_set;
        expr.find(sqrt(w), pows_set);
        expr.find(pow(w1,w2), pows_set);
        exmap pow_map;
        for(auto item : pows_set) {
            if(item.has(WRA(w))) {
                if(item.match(pow(w1,w2)) && !item.op(1).info(info_flags::integer)) {
                    pow_map[item] = exp(log(item.op(0)) * item.op(1));
                } else if(item.match(sqrt(w))) {
                    pow_map[item] = exp(log(item.op(0))/2);
                }
            }
        }

        expr = expr.subs(pow_map); 
        expr = exp_simplify(expr);;
        
        int log_id = 0;
        exset logs_set;
        expr.find(log(w), logs_set);
        exmap log_map;
        lst id_logz_lst;
        for(auto item : logs_set){
            if(item.has(WRA(w))) {
                log_id++;
                id_logz_lst.append(iWF(log_id,item.op(0)));
                log_map[item] = iWF(log_id);
            }
        }
        if(log_id<1) return expr.subs(WRA(w)==w);
        
        expr = expr.subs(log_map).subs(WRA(w)==w);
        
        exmap log_map2;
        for(auto id_logz : id_logz_lst) {
            auto id = id_logz.op(0);
            auto zz = id_logz.op(1);
            exset wra_set;
            zz.find(WRA(w), wra_set);
            if(wra_set.size()<1) return expr;
            if(wra_set.size()!=1) {
                cout << zz << endl;
                throw Error("ContinuousWRA: too many WRA.");
            }
            zz = zz.subs(WRA(w)==xwra);
            ex wra = (*(wra_set.begin())).op(0);

            ex ret = log(zz.subs(xwra==wra));
            if(nc<1) {
                log_map2[iWF(id)] = ret;
                break;
            }
            int total=0;
            int ReIm[nc][2];
            for(int k=0; k<=nc; k++) {
                ex zzk;
                if(k==0) zzk = NN(zz.subs(xwra==wra/(25*nc)));
                else zzk = NN(zz.subs(xwra==k*wra/nc));
                if(!is_a<numeric>(zzk)) throw Error("ContinuousWRA: zzk is not numeric: "+ex2str(zzk));
                auto nzzk = ex_to<numeric>(zzk);
                auto curR = real(nzzk);
                auto curI = imag(nzzk);
                if(is_zero(curI) && k==0 && curR<0) ReIm[total][1] = -1;
                else if(is_zero(curR) || is_zero(curR)) continue; 
                else ReIm[total][1] = curI>0 ? 1 : -1;
                ReIm[total][0] = curR>0 ? 1 : -1;
                total++;
            }
            
            int cutN = 0;
            for(int k=0; k<total-1; k++) {
                if(ReIm[k][0]*ReIm[k+1][0]<0 && ReIm[k][1]*ReIm[k+1][1]<0) throw Error("ContinuousWRA: 1<->3 or 2<->4 happened.");
                if(ReIm[k][0]<0 && ReIm[k+1][0]<0 && ReIm[k][1]*ReIm[k+1][1]<0) {
                    if(ReIm[k][1]>0) cutN++;
                    else cutN--;
                }
            }
            if(cutN!=0) ret += I * cutN * 2 * Pi;
            log_map2[iWF(id)] = ret;
        }
        expr = expr.subs(log_map2);
        if(expr.has(iWF(w))) throw Error("ContinuousWRA: iWF(w) still exists in final result.");
        return expr;
    }
    
    /**
     * @brief Contours, note that here we need to provide the specific Parameter
     * @param key the to retrive data from CIPrepares call
     * @param pkey secondary key, final key will "key-pkey" to store Contours data
     * @param kid only the kid-th will be evaluated and updated the original "key-pkey" data
     */
    void SecDec::Integrates(const string & key, const string & pkey, int kid) {
        if(IsZero) return;
        if(Integrator==NULL) Integrator = new HCubature();
                
        if(Verbose>0) cout << Color_HighLight << "  Integrates @ " << now() << RESET << endl;
        
        lst lstRE;
        auto pid = getpid();
        ostringstream fsofn, sofn, cmd;
        if(key == "") {
            sofn << pid << ".so";
            fsofn << pid << "F.so";
        } else {
            auto oDigits = Digits;
            Digits = NNDigits; // a fix to float overflow
            sofn << key << ".so";
            ostringstream garfn;
            garfn << key << ".ci.gar";
            archive ar;
            ifstream in(garfn.str());
            in >> ar;
            in.close();
            auto c = ar.unarchive_ex(GiNaC_archive_Symbols, "c");
            if(c!=19790923) throw Error("Integrates: *.ci.gar error!");
            auto gl = ar.unarchive_ex(GiNaC_archive_Symbols, "soLimit");
            auto epn = ar.unarchive_ex(GiNaC_archive_Symbols, "epN");
            auto epsn = ar.unarchive_ex(GiNaC_archive_Symbols, "epsN");
            soLimit = ex2int(gl);
            epN = ex2int(epn);
            epsN = ex2int(epsn);
            auto res = ar.unarchive_ex(GiNaC_archive_Symbols, "res");
            ciResult.clear();
            for(auto item : ex_to<lst>(res)) ciResult.push_back(ex_to<lst>(item));
            garfn.clear();
            garfn.str("");
            garfn << key;
            if(pkey != "") garfn << "-" << pkey;
            garfn << ".las.gar";
            if(file_exists(garfn.str().c_str())) {
                archive la_ar;
                ifstream la_in(garfn.str());
                la_in >> la_ar;
                la_in.close();
                auto la_c = la_ar.unarchive_ex(GiNaC_archive_Symbols, "c");
                auto la_res = la_ar.unarchive_ex(GiNaC_archive_Symbols, "res");
                if(la_c!=19790923) throw Error("Integrates: *.ci.gar error!");
                for(auto item : ex_to<lst>(la_res)) {
                    LambdaMap[item.op(0)] = item.op(1);
                }
            }
            
            if(kid>0) {
                garfn.clear();
                garfn.str("");
                garfn << key;
                if(pkey != "") garfn << "-" << pkey;
                garfn << ".res.gar";
                if(!file_exists(garfn.str().c_str())) {
                    throw Error("Integrates: File Not Found: " + garfn.str());
                }
                
                archive res_ar;
                ifstream res_in(garfn.str());
                res_in >> res_ar;
                res_in.close();
                auto res_c = res_ar.unarchive_ex(GiNaC_archive_Symbols, "c");
                auto relst = res_ar.unarchive_ex(GiNaC_archive_Symbols, "relst");
                if(res_c!=19790923) throw Error("*.res.gar error with kid!");
                lstRE = ex_to<lst>(relst);
            }
            Digits = oDigits;
        }
        
        void* main_module = dlopen(sofn.str().c_str(), RTLD_NOW);
        if(main_module == nullptr) {
            main_module = dlopen(("./"+sofn.str()).c_str(), RTLD_NOW);
            if(main_module == nullptr) {
                cout << "dlerror(): " << dlerror() << endl;
                throw Error("Integrates: could not open main module!");
            }
        }
        
        if(!debug && key == "") {
            if(file_exists(fsofn.str().c_str())) remove(fsofn.str().c_str());
            remove(sofn.str().c_str());
        }
        vector<void*> ex_modules;
        for(int n=1; true; n++) {
            ostringstream ex_sofn;
            if(key == "") {
                ex_sofn << pid << "X" << n << ".so";
            } else {
                ex_sofn << key << "X" << n << ".so";
            }
            if(file_exists(ex_sofn.str().c_str())) {
                void* module = dlopen(ex_sofn.str().c_str(), RTLD_NOW);
                if(module == nullptr) {
                    module = dlopen(("./"+ex_sofn.str()).c_str(), RTLD_NOW);
                    if(module == nullptr) {
                        cout << "dlerror(): " << dlerror() << endl;
                        throw Error("Integrates: could not open ex-module!");
                    }
                }
                ex_modules.push_back(module);
                if(!debug && key == "") remove(ex_sofn.str().c_str());
            } else break;
        }
        
        int npara = 0;
        lst plRepl;
        for(auto kv : Parameter) {
            plRepl.append(PL(kv.first)==kv.second);
            if(kv.first>npara) npara = kv.first;
        }
        plRepl.sort();
        plRepl.unique();
        
        int total = ciResult.size(), current = 0;
        qREAL stot = sqrtq(total*1.Q);

        ResultError = 0;
        //----------------------------------------------------------------
        for(auto &item : ciResult) {
            current++;
            if(kid>0 && current != kid) continue;
            if(Verbose>0) {
                cout << "\r                                           \r";
                cout << "  \\--Integrating [" <<current<<"/"<<total<< "] " << flush;
            } 
            
            unsigned int xsize = 0;
            ex co;
            vector<ex> xs, fxs;
            xsize = ex_to<numeric>(item.op(1)).to_int();
            co = item.op(2).subs(plRepl).subs(iEpsilon==iEpsilonN);
            if(co.has(WRA(w))) co = ContinuousWRA(co);
            
            if(xsize<1) { 
                // { expr, xs.size(), kvf.op(0), -1}
                ex exint = item.op(0).subs(plRepl).subs(iEpsilon==iEpsilonN);
                if(exint.has(WRA(w))) exint = ContinuousWRA(exint);
                ex res = VE(NN(exint),0);
                if(Verbose>5) {
                    cout << "XDim=" << xsize << endl;
                    cout << Color_HighLight << "     IRes = "<< HepLib::SD::VEResult(VESimplify(res)) << RESET << endl;
                }
                if(kid>0) {
                    lstRE.let_op(kid-1) = res*co;
                    break;
                }
                lstRE.append(res*co);
                continue;
            }
            
            //{ idx, xs.size(), kvf.op(0), ft_n }
            ex intid = item.op(0);
            ex ftid = item.op(3);
            
            if(co.is_zero()) continue;
            co = collect_ex(co, eps);
            if(co.is_zero()) continue;
            if(co.has(PL(w))) throw Error("Integrates: PL found @ " + ex2str(co));
            qREAL cmax = -1;
            int reim = 0;
            if(ReIm==3) reim = 3;
            
            for(int si=co.ldegree(eps); si<=co.degree(eps); si++) {
                auto tmp = co.coeff(eps, si);
                if(tmp.has(eps)) throw Error("Integrates: eps found @ " + ex2str(tmp));
                tmp = collect_ex(tmp, ep);
                for(int i=tmp.ldegree(ep); i<=tmp.degree(ep); i++) {
                    auto ccRes = NN(tmp.coeff(ep, i)).expand();
                    lst css;
                    css.append(ccRes);
                    if(is_a<add>(ccRes)) {
                        for(auto item : ccRes) css.append(item);
                    }

                    for(int ci=0; ci<css.nops(); ci++) {
                        auto nt = NN(css.op(ci).subs(epz==1).subs(log(vs)==1).subs(vs==1).subs(nReplacements).subs(lst{
                            epsID(w)==1, CV(w1,w2)==w2, ep==ex(1)/111, eps==ex(1)/1111
                        }));
                        if(nt.has(ep)) throw Error("Integrates: ep found @ nt = "+ex2str(nt));
                        
                        lst nt_lst;
                        if(!is_a<numeric>(nt)) {
                            auto cv_lst = collect_lst(nt,[](const ex &e)->bool{return Symbol::has(e);});
                            for(auto nti : cv_lst) {
                                auto nnt = nti.op(0);
                                if(!is_a<numeric>(nnt)) throw Error("Integrates: Not a number with nt = "+ex2str(nnt));
                                nt_lst.append(nnt);
                            }
                        } else nt_lst.append(nt);
                        
                        for(auto nnt : nt_lst) {
                            if(ReIm!=3 && reim!=3) {
                                if(ex_to<numeric>(nnt).imag()==0) {
                                    if(reim==2) reim = 3;
                                    else reim = 1;
                                } else if(ex_to<numeric>(nnt).real()==0) {
                                    if(reim==1) reim = 3;
                                    else reim = 2;
                                } else {
                                    reim = 3;
                                }
                            }
                            nnt = NN(abs(nnt)); // no PL here, then nReplacements
                            
                            qREAL qnt = ex2q(nnt);
                            if(qnt > cmax) cmax = qnt;
                        }
                    }
                }
            }
            if(cmax<=0) {
                while(true) {
                    auto co2 = subs(co,lst{exp(w1*log(w2)+w3)==pow(w2,w1)*exp(w3),exp(w1*log(w2))==pow(w2,w1)});
                    if(is_zero(co2-co)) break;
                    co = co2;
                }
                if(normal(co).is_zero()) continue;
                throw Error("Integrates: cmax<=0 with co = "+ex2str(co));
            }
            if(reim!=3 && ReIm!=3) {
                if(ReIm==2) {
                    if(reim==1) reim=2;
                    else if(reim==2) reim=1;
                }
            }
            
            if(Verbose>5) cout << "XDim=" << xsize << ", EpsAbs=" << (double)(EpsAbs/cmax/stot) << "/" << (double)cmax << endl;
            
            auto las = LambdaMap[ftid];
            bool hasF = (ftid>0);
            if(hasF && las.is_zero()) throw Error("Integrates: lambda with the key(ft_n=" + ex2str(ftid) + ") is NOT found!");
            
            if(hasF && !is_a<lst>(las)) {
                if(!is_zero(las-ex(1979))) { // the convention for xPositive or explict real mode
                    throw Error("Integrates: something is wrong with the F-term @ ft_n = "+ex2str(ftid) + ", las=" + ex2str(las));
                } else {
                    hasF = false;
                }
            }
            
            IntegratorBase::SD_Type fp = nullptr, fpQ = nullptr, fpMP = nullptr;
            IntegratorBase::FT_Type ftp = nullptr;
            int idx = ex_to<numeric>(intid).to_int();
            auto module = main_module;
            if(idx>=soLimit) module = ex_modules[idx/soLimit-1];
            ostringstream fname;
            if(hasF) fname << "C";
            fname << "SDD_" << idx;
            fp = (IntegratorBase::SD_Type)dlsym(module, fname.str().c_str());
            if(fp==NULL) {
                cout << "dlerror(): " << dlerror() << endl;
                throw Error("Integrates: fp==NULL");
            }
            
            fname.clear();
            fname.str("");
            if(hasF) fname << "C";
            fname << "SDQ_" << idx;
            fpQ = (IntegratorBase::SD_Type)dlsym(module, fname.str().c_str());
            if(fpQ==NULL) {
                cout << "dlerror(): " << dlerror() << endl;
                throw Error("Integrates: fpQ==NULL");
            }

            fname.clear();
            fname.str("");
            if(hasF) fname << "C";
            fname << "SDMP_" << idx;
            fpMP = (IntegratorBase::SD_Type)dlsym(module, fname.str().c_str());
                            
            if(is_a<lst>(las)) {
                fname.clear();
                fname.str("");
                fname << "FT_" << idx;
                ftp = (IntegratorBase::FT_Type)dlsym(module, fname.str().c_str());
                if(ftp==NULL) {
                    cout << "dlerror(): " << dlerror() << endl;
                    throw Error("Integrates: ftp==NULL.");
                }
            }
            
            qREAL lambda[las.nops()];
            qREAL paras[npara+1];
            for(auto kv : Parameter) paras[kv.first] = ex2q(kv.second);
            
            Integrator->ReIm = reim;
            Integrator->MPDigits = MPDigits;
            Integrator->Integrand = fp;
            Integrator->IntegrandQ = fpQ;
            Integrator->IntegrandMP = fpMP;
            Integrator->FT = ftp;
            Integrator->Parameter = paras;
            Integrator->Lambda = lambda;
            Integrator->XDim = xsize;
            
            if(hasF) {
                qREAL lamax = ex2q(las.op(las.nops()-1));
                if(lamax > IntLaMax) lamax = IntLaMax;
                
                if(TryPTS<10000) TryPTS = 10000;
                Integrator->RunMAX = -5;
                Integrator->RunPTS = TryPTS/5;
                Integrator->EpsAbs = EpsAbs/cmax/stot/2;
                Integrator->EpsRel = 0;
                
                if(MinPTS[xsize]>0) Integrator->MinPTS = MinPTS[xsize];
                else if(MinPTS[0]>0) Integrator->MinPTS = MinPTS[0];
                else Integrator->MinPTS = RunPTS/10;
                                
                int ctryR = 0, ctry = 0, ctryL = 0;
                int smin = -1;
                ex min_err, min_res;
                long long min_eval;
                qREAL log_lamax = log10q(lamax);
                qREAL log_lamin = log_lamax-1.Q;
                
                ostringstream las_fn;
                las_fn << key;
                if(pkey != "") las_fn << "-" << pkey;
                las_fn << "-" << current << ".las";
                if(use_las && file_exists(las_fn.str().c_str())) {
                    std::ifstream las_ifs;
                    las_ifs.open(las_fn.str(), ios::in);
                    if (!las_ifs) throw Error("Integrates: failed to open *.las file!");
                    for(int i=0; i<las.nops()-1; i++) {
                        dREAL la_tmp;
                        las_ifs >> la_tmp;
                        lambda[i] = la_tmp;
                    }
                    las_ifs.close();
                    auto res = Integrator->Integrate();
                    auto res_tmp = res.subs(VE(w1, w2)==w2);
                    auto err = real_part(res_tmp);
                    if(err < imag_part(res_tmp)) err = imag_part(res_tmp);
                    min_err = err;
                    min_res = res;
                } else {
                // ---------------------------------------
                while(true) {
                    smin = -1;
                    min_err = 0;
                    ex lastResErr = 0;
                    bool err_break = false;
                    for(int s=0; s<=LambdaSplit; s++) {
                        if(Verbose>10 && s==0) {
                            if(ctryR>0 || ctry>0 || ctryL>0)
                                cout << "     ------------------------------" << endl;
                        }
                        auto log_cla = (log_lamin + s * (log_lamax-log_lamin) / LambdaSplit);
                        auto cla = powq(10.Q, log_cla);
                        if(cla < 1E-10) throw Error("NIntegrate: too small lambda.");
                        for(int i=0; i<las.nops()-1; i++) lambda[i] = ex2q(las.op(i)) * cla;
     
                        auto res = Integrator->Integrate();
                        if(Verbose>10) {
                            cout << "\r                                                    \r";
                            if(res.has(NaN)) cout << "     λ=" << (double)cla << "/" << Integrator->NEval << ": " << NaN << endl;
                            else cout << "     λ=" << (double)cla << "/" << Integrator->NEval << ": " << VEResult2(VESimplify(res)) << endl;
                        }
                        
                        if(res.has(NaN) && s==0) continue;
                        else if(res.has(NaN)) break;
                        ex res_abs = NN(abs(res.subs(VE(w1,w2)==w1)));
                        if(lastResErr.is_zero()) lastResErr = res;
                        auto diff = VESimplify(lastResErr - res);
                        diff = diff.subs(VE(0,0)==0);
                        exset ves;
                        diff.find(VE(w0, w1), ves);
                        for(auto ve : ves) {
                            auto ve0 = abs(ve.op(0));
                            if(ve0>ve.op(1)) { 
                                if(numeric("1.E10")*ve0<res_abs) continue; // avoid fluctuation aroud 0
                                if(numeric("1.E10")*ve0<q2ex(EpsAbs)) continue; // avoid fluctuation aroud 0
                                err_break = true;
                                break;
                            }
                        }
                        if(err_break) {
                            if(Verbose>10) cout << Color_HighLight << "     Error Break ..." << RESET << endl;
                            break;
                        }
                        lastResErr = res;
                        
                        auto res_tmp = res.subs(VE(w1, w2)==w2);
                        auto err = real_part(res_tmp);
                        if(err < imag_part(res_tmp)) err = imag_part(res_tmp);
                        if(smin<0 || err < min_err) {
                            min_err = err;
                            min_res = res;
                            min_eval = Integrator->NEval;
                            smin = s;
                        } 
                        if(s>0 && min_err < q2ex(EpsAbs/cmax/stot)) { 
                            // s>0 make sure at least 2 λs compatiable
                            if(Verbose>5) {
                                cout << Color_HighLight << "     λ=" << (double)cla << "/" << min_eval << ": " << HepLib::SD::VEResult(VESimplify(min_res)) << RESET << endl;
                            }
                            
                            smin = -2;
                            if(kid>0) lstRE.let_op(kid-1) = co * min_res;
                            else lstRE.append(co * min_res);
                            break;
                        }
                        if(err > 100 * min_err) break; // s>0 make sure at least 2 λs compatiable
                    }
                    if(smin == -2) break;
                    if(smin == -1) throw Error("Integrates: smin = -1, too small lambda (<1E-10)!");
                    
                    if(smin <= 0) {
                        if((!err_break) && (ctryL >= CTryLeft || ctryR>0)) break;
                        log_lamax = log_lamin;
                        log_lamin -= 1.Q;
                        if(!err_break) ctryL++;
                    } else if(smin >= LambdaSplit) {
                        if(ctryR >= CTryRight || ctryL>0) break;
                        log_lamin = log_lamax;
                        log_lamax += log10q(CTryRightRatio);
                        ctryR++;
                    } else {
                        if(ctry >= CTry) break;
                        auto la1 = log_lamin + (smin-1) * (log_lamax-log_lamin) / LambdaSplit;
                        auto la2 = log_lamin + (smin+1) * (log_lamax-log_lamin) / LambdaSplit;
                        log_lamin = la1;
                        log_lamax = la2;
                        ctry++;
                    }
                }
                
                if(smin == -2) {
                    if(kid>0) break;
                    continue;
                }
                
                auto log_cla = (log_lamin + smin * (log_lamax-log_lamin) / LambdaSplit);
                auto cla = powq(10.Q, log_cla);
                if(Verbose>5) cout << Color_HighLight << "     Final λ = " << (double)cla << " / " << las.op(las.nops()-1) << RESET << endl;
                for(int i=0; i<las.nops()-1; i++) {
                    lambda[i] = ex2q(las.op(i)) * cla;
                }
                // ---------------------------------------
                }
                
                // ---------------------------------------
                // try HookeJeeves
                // ---------------------------------------
                if( use_ErrMin && (min_err > q2ex(1E5 * EpsAbs/cmax/stot)) ) {
                    ErrMin::lastResErr = min_res;
                    auto miner = new HookeJeeves();
                    ErrMin::miner = miner;
                    ErrMin::Integrator = Integrator;
                    dREAL oo[las.nops()-1], ip[las.nops()-1];
                    for(int i=0; i<las.nops()-1; i++) ip[i] = oo[i] = lambda[i];
                    ErrMin::lambda = oo;
                    ErrMin::err_max = 1E100;
                    auto oerrmin = ErrMin::err_min;
                    ErrMin::err_min = oerrmin < 0 ? -oerrmin * ex2q(min_err) : oerrmin/cmax;
                    ErrMin::RunRND = 0;
                    miner->Minimize(las.nops()-1, ErrMin::IntError, ip);
                    delete miner;
                    ErrMin::err_min = oerrmin;
                    for(int i=0; i<las.nops()-1; i++) {
                        lambda[i] = ErrMin::lambda[i];
                    }
                    
                    Integrator->Lambda = lambda; // Integrator->Lambda changed in ErrMin
                    if(Verbose>5) {
                        cout << Color_HighLight << "     Final λs: " << RESET;
                        for(int i=0; i<xsize; i++) {
                            char buffer[128];
                            quadmath_snprintf(buffer, sizeof buffer, "%.6QG", lambda[i]);
                            cout << buffer << " ";
                        }
                        cout << endl << "     ------------------------------" << endl;
                    }
                }
                // ---------------------------------------
                
                if(save_las) {
                    std::ofstream las_ofs;
                    las_ofs.open(las_fn.str(), ios::out);
                    if (las_ofs) {
                        for(int i=0; i<las.nops()-1; i++) {
                            dREAL la_tmp = lambda[i];
                            las_ofs << la_tmp << " ";
                        }
                        las_ofs << endl;
                        las_ofs.close();
                    }
                }
            }

            Integrator->RunMAX = RunMAX;
            Integrator->RunPTS = RunPTS;
            Integrator->EpsAbs = EpsAbs/cmax/stot;
            Integrator->EpsRel = 0;
            
            if(MinPTS[xsize]>0) Integrator->MinPTS = MinPTS[xsize];
            else if(MinPTS[0]>0) Integrator->MinPTS = MinPTS[0];
            else Integrator->MinPTS = RunPTS/10;
            
            auto res = Integrator->Integrate();
            if(Verbose>5) {
                cout << Color_HighLight << "     IRes = "<< HepLib::SD::VEResult(VESimplify(res)) << RESET << endl;
            }
            if(res.has(NaN)) {
                ResultError = NaN;
                if(kid>0) lstRE.let_op(kid-1) = NaN;
                else lstRE.append(NaN);
                break;
            } else {
                if(kid>0) {
                    lstRE.let_op(kid-1) = co * res;
                    break;
                } else lstRE.append(co * res);
            }
        }
        //----------------------------------------------------------------
                
        if(use_dlclose) {
            for(auto module : ex_modules) dlclose(module);
            dlclose(main_module);
        }
        if(total>0 && Verbose>1) cout << "@" << now(false) << endl;
        
        if(!ResultError.is_equal(NaN)) {
            ResultError = 0;
            for(auto item : lstRE) ResultError += item;
            ResultError = VESimplify(ResultError,epN,epsN);
        }

        if(key != "") {
            ostringstream garfn;
            garfn << key;
            if(pkey != "") garfn << "-" << pkey;
            garfn << ".res.gar";
            archive ar;
            ar.archive_ex(ResultError, "res");
            ar.archive_ex(lstRE, "relst");
            ar.archive_ex(19790923, "c");
            ofstream out(garfn.str());
            out << ar;
            out.close();
        }
    }
    
    /**
     * @brief Contours, note that here we need to provide the specific Parameter
     * @param key the to retrive data from CIPrepares call
     * @param pkey secondary key, final key will "key-pkey" to store Contours data
     * @param err only the item ( with error > err ) will be evaluated and updated the original "key-pkey" data
     */
    void SecDec::ReIntegrates(const string & key, const string & pkey, qREAL err) {
        if(IsZero) return;
        if(Integrator==NULL) Integrator = new HCubature();
                
        if(Verbose>0) cout << Color_HighLight << "  Integrates @ " << now() << RESET << endl;
        
        lst lstRE;
        auto pid = getpid();
        ostringstream fsofn, sofn, cmd;
        if(true) {
            auto oDigits = Digits;
            Digits = NNDigits; // a fix to float overflow
            sofn << key << ".so";
            ostringstream garfn;
            garfn << key << ".ci.gar";
            archive ar;
            ifstream in(garfn.str());
            in >> ar;
            in.close();
            auto c = ar.unarchive_ex(GiNaC_archive_Symbols, "c");
            if(c!=19790923) throw Error("Integrates: *.ci.gar error!");
            auto gl = ar.unarchive_ex(GiNaC_archive_Symbols, "soLimit");
            auto epn = ar.unarchive_ex(GiNaC_archive_Symbols, "epN");
            auto epsn = ar.unarchive_ex(GiNaC_archive_Symbols, "epsN");
            soLimit = ex2int(gl);
            epN = ex2int(epn);
            epsN = ex2int(epsn);
            auto res = ar.unarchive_ex(GiNaC_archive_Symbols, "res");
            ciResult.clear();
            for(auto item : ex_to<lst>(res)) ciResult.push_back(ex_to<lst>(item));
            garfn.clear();
            garfn.str("");
            garfn << key;
            if(pkey != "") garfn << "-" << pkey;
            garfn << ".las.gar";
            if(file_exists(garfn.str().c_str())) {
                archive la_ar;
                ifstream la_in(garfn.str());
                la_in >> la_ar;
                la_in.close();
                auto la_c = la_ar.unarchive_ex(GiNaC_archive_Symbols, "c");
                auto la_res = la_ar.unarchive_ex(GiNaC_archive_Symbols, "res");
                if(la_c!=19790923) throw Error("Integrates: *.ci.gar error!");
                for(auto item : ex_to<lst>(la_res)) {
                    LambdaMap[item.op(0)] = item.op(1);
                }
            }
            
            if(true) {
                garfn.clear();
                garfn.str("");
                garfn << key;
                if(pkey != "") garfn << "-" << pkey;
                garfn << ".res.gar";
                if(!file_exists(garfn.str().c_str())) {
                    throw Error("Integrates: File Not Found: " + garfn.str());
                }
                
                archive res_ar;
                ifstream res_in(garfn.str());
                res_in >> res_ar;
                res_in.close();
                auto res_c = res_ar.unarchive_ex(GiNaC_archive_Symbols, "c");
                auto relst = res_ar.unarchive_ex(GiNaC_archive_Symbols, "relst");
                if(res_c!=19790923) throw Error("*.res.gar error with ReIntegrates!");
                lstRE = ex_to<lst>(relst);
            }
            Digits = oDigits;
        }
        
        void* main_module = dlopen(sofn.str().c_str(), RTLD_NOW);
        if(main_module == nullptr) {
            main_module = dlopen(("./"+sofn.str()).c_str(), RTLD_NOW);
            if(main_module == nullptr) {
                cout << "dlerror(): " << dlerror() << endl;
                throw Error("Integrates: could not open main module!");
            }
        }
        
        vector<void*> ex_modules;
        for(int n=1; true; n++) {
            ostringstream ex_sofn;
            ex_sofn << key << "X" << n << ".so";
            if(file_exists(ex_sofn.str().c_str())) {
                void* module = dlopen(ex_sofn.str().c_str(), RTLD_NOW);
                if(module == nullptr) {
                    module = dlopen(("./"+ex_sofn.str()).c_str(), RTLD_NOW);
                    if(module == nullptr) {
                        cout << "dlerror(): " << dlerror() << endl;
                        throw Error("Integrates: could not open ex-module!");
                    }
                }
                ex_modules.push_back(module);
                if(!debug && key == "") remove(ex_sofn.str().c_str());
            } else break;
        }
        
        int npara = 0;
        lst plRepl;
        for(auto kv : Parameter) {
            plRepl.append(PL(kv.first)==kv.second);
            if(kv.first>npara) npara = kv.first;
        }
        plRepl.sort();
        plRepl.unique();
        
        int total = ciResult.size(), current = 0;
        qREAL stot = sqrtq(total*1.Q);

        ResultError = 0;
        //----------------------------------------------------------------
        for(auto &item : ciResult) {
            current++;
            auto cmerr = ex2q(VEMaxErr(lstRE.op(current-1)));
            if(cmerr < err) continue;
            if(Verbose>10) {
                char es[64];
                quadmath_snprintf(es, sizeof es, "%.10QG", cmerr);
                cout << "  \\--Current Err: " << es << endl;
            }
            if(Verbose>0) {
                cout << "\r                                           \r";
                cout << "  \\--Integrating [" <<current<<"/"<<total<< "] " << flush;
            }
            
            unsigned int xsize = 0;
            ex co;
            vector<ex> xs, fxs;
            xsize = ex_to<numeric>(item.op(1)).to_int();
            co = item.op(2).subs(plRepl).subs(iEpsilon==iEpsilonN);
            if(co.has(WRA(w))) co = ContinuousWRA(co);
            
            if(xsize<1) {
                // { expr, xs.size(), kvf.op(0), -1}
                ex exint = item.op(0).subs(plRepl).subs(iEpsilon==iEpsilonN);
                if(exint.has(WRA(w))) exint = ContinuousWRA(exint);
                ex res = VE(NN(exint),0);
                if(Verbose>5) {
                    cout << "XDim=" << xsize << endl;
                    cout << Color_HighLight << "     IRes = "<< HepLib::SD::VEResult(VESimplify(res)) << RESET << endl;
                }
                lstRE.let_op(current-1) = res*co;
                continue;
            }
            
            //{ idx, xs.size(), kvf.op(0), ft_n }
            ex intid = item.op(0);
            ex ftid = item.op(3);
            
            if(co.is_zero()) continue;
            co = collect_ex(co, eps);
            if(co.is_zero()) continue;
            if(co.has(PL(w))) throw Error("Integrates: PL found @ " + ex2str(co));
            qREAL cmax = -1;
            int reim = 0;
            if(ReIm==3) reim = 3;
            
            for(int si=co.ldegree(eps); si<=co.degree(eps); si++) {
                auto tmp = co.coeff(eps, si);
                if(tmp.has(eps)) throw Error("Integrates: eps found @ " + ex2str(tmp));
                tmp = collect_ex(tmp, ep);
                for(int i=tmp.ldegree(ep); i<=tmp.degree(ep); i++) {
                    auto ccRes = NN(tmp.coeff(ep, i)).expand();
                    lst css;
                    css.append(ccRes);
                    if(is_a<add>(ccRes)) {
                        for(auto item : ccRes) css.append(item);
                    }

                    for(int ci=0; ci<css.nops(); ci++) {
                        auto nt = NN(css.op(ci).subs(epz==1).subs(log(vs)==1).subs(vs==1).subs(nReplacements).subs(lst{
                            epsID(w)==1, CV(w1,w2)==w2, ep==ex(1)/111, eps==ex(1)/1111
                        }));
                        if(nt.has(ep)) throw Error("Integrates: ep found @ nt = "+ex2str(nt));
                        
                        lst nt_lst;
                        if(!is_a<numeric>(nt)) {
                            auto cv_lst = collect_lst(nt,[](const ex &e)->bool{return Symbol::has(e);});
                            for(auto nti : cv_lst) {
                                auto nnt = nti.op(0);
                                if(!is_a<numeric>(nnt)) throw Error("Integrates: Not a number with nt = "+ex2str(nnt));
                                nt_lst.append(nnt);
                            }
                        } else nt_lst.append(nt);
                        
                        for(auto nnt : nt_lst) {
                            if(ReIm!=3 && reim!=3) {
                                if(ex_to<numeric>(nnt).imag()==0) {
                                    if(reim==2) reim = 3;
                                    else reim = 1;
                                } else if(ex_to<numeric>(nnt).real()==0) {
                                    if(reim==1) reim = 3;
                                    else reim = 2;
                                } else {
                                    reim = 3;
                                }
                            }
                            nnt = NN(abs(nnt)); // no PL here, then nReplacements
                            
                            qREAL qnt = ex2q(nnt);
                            if(qnt > cmax) cmax = qnt;
                        }
                    }
                }
            }
            if(cmax<=0) {
                while(true) {
                    auto co2 = subs(co,lst{exp(w1*log(w2)+w3)==pow(w2,w1)*exp(w3),exp(w1*log(w2))==pow(w2,w1)});
                    if(is_zero(co2-co)) break;
                    co = co2;
                }
                if(normal(co).is_zero()) continue;
                throw Error("Integrates: cmax<=0 with co = "+ex2str(co));
            }
            if(reim!=3 && ReIm!=3) {
                if(ReIm==2) {
                    if(reim==1) reim=2;
                    else if(reim==2) reim=1;
                }
            }
            
            if(Verbose>5) cout << "XDim=" << xsize << ", EpsAbs=" << (double)(EpsAbs/cmax/stot) << "/" << (double)cmax << endl;
            
            auto las = LambdaMap[ftid];
            bool hasF = (ftid>0);
            if(hasF && las.is_zero()) throw Error("Integrates: lambda with the key(ft_n=" + ex2str(ftid) + ") is NOT found!");
            
            if(hasF && !is_a<lst>(las)) {
                if(!is_zero(las-ex(1979))) { // the convention for xPositive or explict real mode
                    throw Error("Integrates: something is wrong with the F-term @ ft_n = "+ex2str(ftid) + ", las=" + ex2str(las));
                } else {
                    hasF = false;
                }
            }
            
            IntegratorBase::SD_Type fp = nullptr, fpQ = nullptr, fpMP = nullptr;
            IntegratorBase::FT_Type ftp = nullptr;
            int idx = ex_to<numeric>(intid).to_int();
            auto module = main_module;
            if(idx>=soLimit) module = ex_modules[idx/soLimit-1];
            ostringstream fname;
            if(hasF) fname << "C";
            fname << "SDD_" << idx;
            fp = (IntegratorBase::SD_Type)dlsym(module, fname.str().c_str());
            if(fp==NULL) {
                cout << "dlerror(): " << dlerror() << endl;
                throw Error("Integrates: fp==NULL");
            }
            
            fname.clear();
            fname.str("");
            if(hasF) fname << "C";
            fname << "SDQ_" << idx;
            fpQ = (IntegratorBase::SD_Type)dlsym(module, fname.str().c_str());
            if(fpQ==NULL) {
                cout << "dlerror(): " << dlerror() << endl;
                throw Error("Integrates: fpQ==NULL");
            }

            fname.clear();
            fname.str("");
            if(hasF) fname << "C";
            fname << "SDMP_" << idx;
            fpMP = (IntegratorBase::SD_Type)dlsym(module, fname.str().c_str());
                            
            if(is_a<lst>(las)) {
                fname.clear();
                fname.str("");
                fname << "FT_" << idx;
                ftp = (IntegratorBase::FT_Type)dlsym(module, fname.str().c_str());
                if(ftp==NULL) {
                    cout << "dlerror(): " << dlerror() << endl;
                    throw Error("Integrates: ftp==NULL.");
                }
            }
            
            qREAL lambda[las.nops()];
            qREAL paras[npara+1];
            for(auto kv : Parameter) paras[kv.first] = ex2q(kv.second);
            
            Integrator->ReIm = reim;
            Integrator->MPDigits = MPDigits;
            Integrator->Integrand = fp;
            Integrator->IntegrandQ = fpQ;
            Integrator->IntegrandMP = fpMP;
            Integrator->FT = ftp;
            Integrator->Parameter = paras;
            Integrator->Lambda = lambda;
            Integrator->XDim = xsize;
            
            if(hasF) {
                qREAL lamax = ex2q(las.op(las.nops()-1));
                if(lamax > IntLaMax) lamax = IntLaMax;
                
                if(TryPTS<10000) TryPTS = 10000;
                Integrator->RunMAX = -5;
                Integrator->RunPTS = TryPTS/5;
                Integrator->EpsAbs = EpsAbs/cmax/stot/2;
                Integrator->EpsRel = 0;
                
                if(MinPTS[xsize]>0) Integrator->MinPTS = MinPTS[xsize];
                else if(MinPTS[0]>0) Integrator->MinPTS = MinPTS[0];
                else Integrator->MinPTS = RunPTS/10;
                                
                int ctryR = 0, ctry = 0, ctryL = 0;
                int smin = -1;
                ex min_err, min_res;
                long long min_eval;
                qREAL log_lamax = log10q(lamax);
                qREAL log_lamin = log_lamax-1.Q;
                
                ostringstream las_fn;
                las_fn << key;
                if(pkey != "") las_fn << "-" << pkey;
                las_fn << "-" << current << ".las";
                if(use_las && file_exists(las_fn.str().c_str())) {
                    std::ifstream las_ifs;
                    las_ifs.open(las_fn.str(), ios::in);
                    if (!las_ifs) throw Error("Integrates: failed to open *.las file!");
                    for(int i=0; i<las.nops()-1; i++) {
                        dREAL la_tmp;
                        las_ifs >> la_tmp;
                        lambda[i] = la_tmp;
                    }
                    las_ifs.close();
                    auto res = Integrator->Integrate();
                    auto res_tmp = res.subs(VE(w1, w2)==w2);
                    auto err = real_part(res_tmp);
                    if(err < imag_part(res_tmp)) err = imag_part(res_tmp);
                    min_err = err;
                    min_res = res;
                } else {
                // ---------------------------------------
                while(true) {
                    smin = -1;
                    min_err = 0;
                    ex lastResErr = 0;
                    bool err_break = false;
                    for(int s=0; s<=LambdaSplit; s++) {
                        if(Verbose>10 && s==0) {
                            if(ctryR>0 || ctry>0 || ctryL>0)
                                cout << "     ------------------------------" << endl;
                        }
                        auto log_cla = (log_lamin + s * (log_lamax-log_lamin) / LambdaSplit);
                        auto cla = powq(10.Q, log_cla);
                        if(cla < 1E-10) throw Error("NIntegrate: too small lambda.");
                        for(int i=0; i<las.nops()-1; i++) lambda[i] = ex2q(las.op(i)) * cla;
     
                        auto res = Integrator->Integrate();
                        if(Verbose>10) {
                            cout << "\r                                                    \r";
                            if(res.has(NaN)) cout << "     λ=" << (double)cla << "/" << Integrator->NEval << ": " << NaN << endl;
                            else cout << "     λ=" << (double)cla << "/" << Integrator->NEval << ": " << VEResult2(VESimplify(res)) << endl;
                        }
                        
                        if(res.has(NaN) && s==0) continue;
                        else if(res.has(NaN)) break;
                        ex res_abs = NN(abs(res.subs(VE(w1,w2)==w1)));
                        if(lastResErr.is_zero()) lastResErr = res;
                        auto diff = VESimplify(lastResErr - res);
                        diff = diff.subs(VE(0,0)==0);
                        exset ves;
                        diff.find(VE(w0, w1), ves);
                        for(auto ve : ves) {
                            auto ve0 = abs(ve.op(0));
                            if(ve0>ve.op(1)) {
                                if(numeric("1.E10")*ve0<res_abs) continue; // avoid fluctuation aroud 0
                                if(numeric("1.E10")*ve0<q2ex(EpsAbs)) continue; // avoid fluctuation aroud 0
                                err_break = true;
                                break;
                            }
                        }
                        if(err_break) {
                            if(Verbose>10) cout << Color_HighLight << "     Error Break ..." << RESET << endl;
                            break;
                        }
                        lastResErr = res;
                        
                        auto res_tmp = res.subs(VE(w1, w2)==w2);
                        auto err = real_part(res_tmp);
                        if(err < imag_part(res_tmp)) err = imag_part(res_tmp);
                        if(smin<0 || err < min_err) {
                            min_err = err;
                            min_res = res;
                            min_eval = Integrator->NEval;
                            smin = s;
                        }
                        if(s>0 && min_err < q2ex(EpsAbs/cmax/stot)) {
                            // s>0 make sure at least 2 λs compatiable
                            if(Verbose>5) {
                                cout << Color_HighLight << "     λ=" << (double)cla << "/" << min_eval << ": " << HepLib::SD::VEResult(VESimplify(min_res)) << RESET << endl;
                            }
                            
                            smin = -2;
                            lstRE.let_op(current-1) = co * min_res;
                            break;
                        }
                        if(err > 100 * min_err) break; // s>0 make sure at least 2 λs compatiable
                    }
                    if(smin == -2) break;
                    if(smin == -1) throw Error("Integrates: smin = -1, too small lambda (<1E-10)!");
                    
                    if(smin <= 0) {
                        if((!err_break) && (ctryL >= CTryLeft || ctryR>0)) break;
                        log_lamax = log_lamin;
                        log_lamin -= 1.Q;
                        if(!err_break) ctryL++;
                    } else if(smin >= LambdaSplit) {
                        if(ctryR >= CTryRight || ctryL>0) break;
                        log_lamin = log_lamax;
                        log_lamax += log10q(CTryRightRatio);
                        ctryR++;
                    } else {
                        if(ctry >= CTry) break;
                        auto la1 = log_lamin + (smin-1) * (log_lamax-log_lamin) / LambdaSplit;
                        auto la2 = log_lamin + (smin+1) * (log_lamax-log_lamin) / LambdaSplit;
                        log_lamin = la1;
                        log_lamax = la2;
                        ctry++;
                    }
                }
                
                if(smin == -2) continue;
                
                auto log_cla = (log_lamin + smin * (log_lamax-log_lamin) / LambdaSplit);
                auto cla = powq(10.Q, log_cla);
                if(Verbose>5) cout << Color_HighLight << "     Final λ = " << (double)cla << " / " << las.op(las.nops()-1) << RESET << endl;
                for(int i=0; i<las.nops()-1; i++) {
                    lambda[i] = ex2q(las.op(i)) * cla;
                }
                // ---------------------------------------
                }
                
                // ---------------------------------------
                // try HookeJeeves
                // ---------------------------------------
                if( use_ErrMin && (min_err > q2ex(1E5 * EpsAbs/cmax/stot)) ) {
                    ErrMin::lastResErr = min_res;
                    auto miner = new HookeJeeves();
                    ErrMin::miner = miner;
                    ErrMin::Integrator = Integrator;
                    dREAL oo[las.nops()-1], ip[las.nops()-1];
                    for(int i=0; i<las.nops()-1; i++) ip[i] = oo[i] = lambda[i];
                    ErrMin::lambda = oo;
                    ErrMin::err_max = 1E100;
                    auto oerrmin = ErrMin::err_min;
                    ErrMin::err_min = oerrmin < 0 ? -oerrmin * ex2q(min_err) : oerrmin/cmax;
                    ErrMin::RunRND = 0;
                    miner->Minimize(las.nops()-1, ErrMin::IntError, ip);
                    delete miner;
                    ErrMin::err_min = oerrmin;
                    for(int i=0; i<las.nops()-1; i++) {
                        lambda[i] = ErrMin::lambda[i];
                    }
                    
                    Integrator->Lambda = lambda; // Integrator->Lambda changed in ErrMin
                    if(Verbose>5) {
                        cout << Color_HighLight << "     Final λs: " << RESET;
                        for(int i=0; i<xsize; i++) {
                            char buffer[128];
                            quadmath_snprintf(buffer, sizeof buffer, "%.6QG", lambda[i]);
                            cout << buffer << " ";
                        }
                        cout << endl << "     ------------------------------" << endl;
                    }
                }
                // ---------------------------------------
                
                if(save_las) {
                    std::ofstream las_ofs;
                    las_ofs.open(las_fn.str(), ios::out);
                    if (las_ofs) {
                        for(int i=0; i<las.nops()-1; i++) {
                            dREAL la_tmp = lambda[i];
                            las_ofs << la_tmp << " ";
                        }
                        las_ofs << endl;
                        las_ofs.close();
                    }
                }
            }

            Integrator->RunMAX = RunMAX;
            Integrator->RunPTS = RunPTS;
            Integrator->EpsAbs = EpsAbs/cmax/stot;
            Integrator->EpsRel = 0;
            
            if(MinPTS[xsize]>0) Integrator->MinPTS = MinPTS[xsize];
            else if(MinPTS[0]>0) Integrator->MinPTS = xsize * MinPTS[0];
            else Integrator->MinPTS = RunPTS/10;
            
            auto res = Integrator->Integrate();
            if(Verbose>5) {
                cout << Color_HighLight << "     IRes = "<< HepLib::SD::VEResult(VESimplify(res)) << RESET << endl;
            }
            if(res.has(NaN)) {
                ResultError = NaN;
                lstRE.let_op(current-1) = NaN;
            } else {
                lstRE.let_op(current-1) = co * res;
            }
        }
        //----------------------------------------------------------------
                
        if(use_dlclose) {
            for(auto module : ex_modules) dlclose(module);
            dlclose(main_module);
        }
        if(total>0 && Verbose>1) cout << "@" << now(false) << endl;
        
        if(!ResultError.is_equal(NaN)) {
            ResultError = 0;
            for(auto item : lstRE) ResultError += item;
            ResultError = VESimplify(ResultError,epN,epsN);
        }

        if(key != "") {
            ostringstream garfn;
            garfn << key;
            if(pkey != "") garfn << "-" << pkey;
            garfn << ".res.gar";
            archive ar;
            ar.archive_ex(ResultError, "res");
            ar.archive_ex(lstRE, "relst");
            ar.archive_ex(19790923, "c");
            ofstream out(garfn.str());
            out << ar;
            out.close();
        }
    }
    
}
