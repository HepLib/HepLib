/**
 * @file
 * @brief Functions for Numerical Integration
 * @author F. Feng
 * @version 1.0.0
 * @date 2020-04-21
 */
 
#include "SD.h"
#include <math.h>
#include <cmath>

namespace HepLib::SD {
    
    // need Parameter
    void SecDec::Integrates(const string & key, const string & pkey, int kid) {
        if(IsZero) return;
        if(Integrator==NULL) Integrator = new HCubature();
                
        if(Verbose > 1) cout << now() << " - Integrates ..." << endl << flush;
        
        lst lstRE;
        auto pid = getpid();
        ostringstream fsofn, sofn, cmd;
        if(key == "") {
            sofn << pid << ".so";
            fsofn << pid << "F.so";
        } else {
            sofn << key << ".so";
            ostringstream garfn;
            garfn << key << ".ci.gar";
            archive ar;
            ifstream in(garfn.str());
            in >> ar;
            in.close();
            auto c = ar.unarchive_ex(GiNaC_archive_Symbols, "c");
            if(c!=19790923) throw runtime_error("*.ci.gar error!");
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
                if(la_c!=19790923) throw runtime_error("*.ci.gar error!");
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
                    cerr << Color_Error << "Integrates: File Not Found: " << garfn.str() << RESET << endl;
                    exit(1);
                }
                
                archive res_ar;
                ifstream res_in(garfn.str());
                res_in >> res_ar;
                res_in.close();
                auto res_c = res_ar.unarchive_ex(GiNaC_archive_Symbols, "c");
                auto relst = res_ar.unarchive_ex(GiNaC_archive_Symbols, "relst");
                if(res_c!=19790923) throw runtime_error("*.res.gar error with kid!");
                lstRE = ex_to<lst>(relst);
            }
        }
        
        void* main_module = dlopen(sofn.str().c_str(), RTLD_NOW);
        if(main_module == nullptr) {
            cerr << Color_Error << "Integrates: could not open main module!" << RESET << endl;
            cout << "dlerror(): " << dlerror() << endl;
            exit(1);
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
                    cerr << Color_Error << "Integrates: could not open ex-module!" << RESET << endl;
                    cout << "dlerror(): " << dlerror() << endl;
                    exit(1);
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

        for(auto &item : ciResult) {
            current++;
            if(kid>0 && current != kid) continue;
            if(Verbose > 1) {
                cout << "\r  \\--Evaluating [" <<current<<"/"<<total<< "] ... " << flush;
            }
            
            unsigned int xsize = 0;
            ex co, exint, exft;
            vector<ex> xs, fxs;
            xsize = ex_to<numeric>(item.op(1)).to_int();
            if(xsize<1) {
                Digits = 35;
                if(kid>0) {
                    lstRE.let_op(kid-1) = VE(item.op(0).subs(plRepl).evalf(),0) * item.op(2).subs(plRepl);
                    break;
                }
                ResultError +=  VE(item.op(0).subs(plRepl).evalf(),0) * item.op(2).subs(plRepl);
                lstRE.append(VE(item.op(0).subs(plRepl).evalf(),0) * item.op(2).subs(plRepl));
                continue;
            }
            co = item.op(2).subs(plRepl).subs(iEpsilon==0);
            
            if(co.is_zero()) continue;
            if(co.has(PL(w))) {
                cerr << Color_Error << "Integrates: PL found @ " << co << RESET << endl;
                exit(1);
            }
            qREAL cmax = -1;
            int reim = 0;
            if(ReIm==3) reim = 3;
            co = mma_collect(co, eps);
auto coo = co;
            if(co.is_zero()) {cout << "co=" << coo << endl; continue;}
            for(int si=co.ldegree(eps); si<=co.degree(eps); si++) {
                auto tmp = co.coeff(eps, si);
                if(tmp.has(eps)) {
                    cerr << Color_Error << "Integrates: eps found @ " << tmp << RESET << endl;
                    exit(1);
                }
                tmp = mma_collect(tmp, ep);
                for(int i=tmp.ldegree(ep); i<=tmp.degree(ep); i++) {
                    Digits = 50;
                    auto ccRes = tmp.coeff(ep, i).evalf().expand();
                    lst css;
                    css.append(ccRes);
                    if(is_a<add>(ccRes)) {
                        for(auto item : ccRes) css.append(item);
                    }

                    for(int ci=0; ci<css.nops(); ci++) {
                        auto nt = css.op(ci).subs(epz==1).subs(log(vs)==1).subs(vs==1).subs(nReplacements).subs(lst{
                            epsID(w)==1, CV(w1,w2)==w2, ep==ex(1)/111, eps==ex(1)/1111
                        }).evalf();
                        if(!is_a<numeric>(nt)) {
                            cerr << Color_Error << "Integrates: Not a number with nt = " << nt << RESET << endl;
                            exit(1);
                        }
                        if(nt.has(ep)) {
                            cerr << Color_Error << "Integrates: ep found @ nt = " << nt << RESET << endl;
                            exit(1);
                        }
                        
                        if(ReIm!=3 && reim!=3) {
                            if(ex_to<numeric>(nt).imag()==0) {
                                if(reim==2) reim = 3;
                                else reim = 1;
                            } else if(ex_to<numeric>(nt).real()==0) {
                                if(reim==1) reim = 3;
                                else reim = 2;
                            } else {
                                reim = 3;
                            }
                        }
                        nt = abs(nt).evalf(); // no PL here, then nReplacements
                        
                        qREAL qnt = CppFormat::ex2q(nt);
                        if(qnt > cmax) cmax = qnt;
                    }
                }
            }
            if(cmax<=0) {
                cerr << Color_Error << "Integrates: cmax<=0 with co = " << co << RESET <<endl;
                exit(1);
            }
            if(reim!=3 && ReIm!=3) {
                if(reim==1 && ReIm==2) reim=2;
                else if(reim==2 && ReIm==2) reim=1;
            }
            
            if(Verbose > 3) cout << "XDim=" << xsize << ", EpsAbs=" << (double)(EpsAbs/cmax/stot) << "/" << (double)cmax << endl;
            
            auto las = LambdaMap[item.op(3)];
            bool hasF = item.op(3)>0;
            if(hasF && las.is_zero()) {
                cerr << Color_Error << "Integrates: lambda with the key(ft_n=" << item.op(3) << ") is NOT found!" << RESET << endl;
                exit(1);
            }
            
            if(hasF && !is_a<lst>(las)) {
                if(!is_zero(las-ex(1979))) {
                    cerr << Color_Error << "Integrates: something is wrong with the F-term @ ft_n = " << item.op(3) << ", las=" << las << RESET << endl;
                    exit(1);
                } else {
                    hasF = false;
                }
            }
            
            IntegratorBase::SD_Type fp = nullptr, fpQ = nullptr, fpMP = nullptr;
            IntegratorBase::FT_Type ftp = nullptr;
            int idx = ex_to<numeric>(item.op(0)).to_int();
            auto module = main_module;
            if(idx>=GccLimit) module = ex_modules[idx/GccLimit-1];
            ostringstream fname;
            if(hasF) fname << "C";
            fname << "SDD_" << idx;
            fp = (IntegratorBase::SD_Type)dlsym(module, fname.str().c_str());
            if(fp==NULL) {
                cerr << Color_Error << "Integrates: fp==NULL" << RESET << endl;
                cout << "dlerror(): " << dlerror() << endl;
                exit(1);
            }
            
            fname.clear();
            fname.str("");
            if(hasF) fname << "C";
            fname << "SDQ_" << idx;
            fpQ = (IntegratorBase::SD_Type)dlsym(module, fname.str().c_str());
            if(fpQ==NULL) {
                cerr << Color_Error << "Integrates: fpQ==NULL" << RESET << endl;
                cout << "dlerror(): " << dlerror() << endl;
                exit(1);
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
                    cerr << Color_Error << "Integrates: ftp==NULL" << RESET << endl;
                    cout << "dlerror(): " << dlerror() << endl;
                    exit(1);
                }
            }
            
            qREAL lambda[las.nops()];
            qREAL paras[npara+1];
            for(auto kv : Parameter) paras[kv.first] = CppFormat::ex2q(kv.second);
            
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
                qREAL lamax = CppFormat::ex2q(las.op(las.nops()-1));
                if(lamax > LambdaMax) lamax = LambdaMax;
                
                if(TryPTS<10000) TryPTS = 10000;
                Integrator->RunMAX = -5;
                Integrator->RunPTS = TryPTS/5;
                Integrator->EpsAbs = EpsAbs/cmax/2/stot;
                Integrator->EpsRel = 0;
                
                if(MinPTS[xsize]>0) Integrator->MinPTS = MinPTS[xsize];
                else if(MinPTS[0]>0) Integrator->MinPTS = MinPTS[0];
                else Integrator->MinPTS = RunPTS/10;
                                
                int smin = -1;
                int ctryR = 0, ctry = 0, ctryL = 0;
                ex cerr;
                qREAL log_lamax = log10q(lamax);
                qREAL log_lamin = log_lamax-2.5Q;
                
                ostringstream las_fn;
                las_fn << key;
                if(pkey != "") las_fn << "-" << pkey;
                las_fn << "-" << current << ".las";
                if(use_las && file_exists(las_fn.str().c_str())) {
                    std::ifstream las_ifs;
                    las_ifs.open(las_fn.str(), ios::in);
                    if (!las_ifs) throw runtime_error("failed to open *.las file!");
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
                    ErrMin::lastResErr = res;
                    cerr = err;
                } else {
                // ---------------------------------------
                while(true) {
                    smin = -1;
                    ex emin = 0;
                    ex lastResErr = 0;
                    for(int s=0; s<=LambdaSplit; s++) {
                        if(Verbose>10 && s==0) {
                            if(ctryR>0 || ctry>0 || ctryL>0)
                                cout << "     ------------------------------" << endl;
                        }
                        auto log_cla = (log_lamin + s * (log_lamax-log_lamin) / LambdaSplit);
                        auto cla = powq(10.Q, log_cla);
                        if(cla < 1E-10) continue;
                        for(int i=0; i<las.nops()-1; i++) {
                            lambda[i] = CppFormat::ex2q(las.op(i)) * cla;
                        }
     
                        auto res = Integrator->Integrate();
                        if(Verbose>10) {
                            auto oDigits = Digits;
                            Digits = 3;
                            cout << "\r                                                    \r";
                            if(res.has(NaN)) cout << "     λ=" << (double)cla << "/" << Integrator->NEval << ": " << NaN << endl;
                            else cout << "     λ=" << (double)cla << "/" << Integrator->NEval << ": " << HepLib::SD::VEResult(VESimplify(res)) << endl;
                            Digits = oDigits;
                        }
                        
                        if(res.has(NaN) && s==0) continue;
                        else if(res.has(NaN)) break;
                        if(lastResErr.is_zero()) lastResErr = res;
                        auto diff = VESimplify(lastResErr - res);
                        diff = diff.subs(VE(0,0)==0);
                        exset ves;
                        diff.find(VE(w0, w1), ves);
                        bool err_break = false;
                        for(auto ve : ves) {
                            if(abs(ve.op(0)) > ve.op(1)) {
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
                        if(smin<0 || err < emin) {
                            ErrMin::lastResErr = res;
                            cerr = err;
                            smin = s;
                            emin = err;
                            if(emin < CppFormat::q2ex(EpsAbs/cmax/stot)) {
                                smin = -2;
                                if(kid>0) {
                                    lstRE.let_op(kid-1) = co * res;
                                    break;
                                }
                                ResultError += co * res;
                                lstRE.append(co * res);
                                if(Verbose>5) {
                                    cout << Color_HighLight;
                                    cout << "     λ=" << (double)cla << "/" << Integrator->NEval << ": " << HepLib::SD::VEResult(VESimplify(res)) << endl;
                                    cout << RESET;
                                }
                                break;
                            }
                        } else if(err > 1.E3 * emin) {
                            break;
                        }
                    }
                    if(smin == -2) break;
                    if(smin == -1) {
                        std::cerr << Color_Error << "Integrates: smin = -1, optimized lambda NOT found!" << RESET << endl;
                        exit(1);
                    }
                    
                    if(smin <= 0) {
                        if(ctryL >= CTryLeft || ctryR>0) break;
                        log_lamax = log_lamin;
                        log_lamin -= 1.5Q;
                        ctryL++;
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
                if(Verbose > 7) cout << Color_HighLight << "     Final λ = " << (double)cla << " / " << las.op(las.nops()-1) << RESET << endl;
                for(int i=0; i<las.nops()-1; i++) {
                    lambda[i] = CppFormat::ex2q(las.op(i)) * cla;
                }
                // ---------------------------------------
                }
                
                // ---------------------------------------
                // try HookeJeeves
                // ---------------------------------------
                if( use_ErrMin && (cerr > CppFormat::q2ex(1E5 * EpsAbs/cmax/stot)) ) {
                    auto miner = new HookeJeeves();
                    ErrMin::miner = miner;
                    ErrMin::Integrator = Integrator;
                    dREAL oo[las.nops()-1], ip[las.nops()-1];
                    for(int i=0; i<las.nops()-1; i++) ip[i] = oo[i] = lambda[i];
                    ErrMin::lambda = oo;
                    ErrMin::err_max = 1E100;
                    auto err_min = ErrMin::err_min;
                    ErrMin::err_min = err_min < 0 ? -err_min * ex_to<numeric>(cerr).to_double() : err_min/cmax;
                    ErrMin::RunRND = 0;
                    miner->Minimize(las.nops()-1, ErrMin::IntError, ip);
                    delete miner;
                    ErrMin::err_min = err_min;
                    for(int i=0; i<las.nops()-1; i++) {
                        lambda[i] = ErrMin::lambda[i];
                    }
                    
                    Integrator->Lambda = lambda; // Integrator->Lambda changed in ErrMin
                    if(Verbose > 7) {
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
                cout << Color_HighLight;
                cout << "     Res = "<< HepLib::SD::VEResult(VESimplify(res)) << endl;
                cout << RESET;
            }
            if(res.has(NaN)) {
                ResultError = NaN;
                if(kid>0) {
                    lstRE.let_op(kid-1) = NaN;
                } else {
                    lstRE.append(NaN);
                }
                break;
            } else {
                if(kid>0) {
                    lstRE.let_op(kid-1) = co * res;
                    break;
                } else {
                    ResultError += co * res;
                    lstRE.append(co * res);
                }
            }
        }
        
        if(use_dlclose) {
            dlclose(main_module);
            for(auto module : ex_modules) dlclose(module);
        }
        if(total>0 && Verbose > 1) cout << "@" << now(false) << endl;
        
        if(kid>0) {
            ResultError = 0;
            for(auto item : lstRE) ResultError += item;
        }
        ResultError = VESimplify(ResultError,epN,epsN);
        
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
