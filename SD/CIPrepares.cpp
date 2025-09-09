/**
 * @file
 * @brief Functions for Contour/Integration Preparation
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
        if(expResult.size()<1) IsZero = true;
        if(IsZero) return;
        int rc;
        
        if(Verbose > 0) cout << Color_HighLight << "  CIPrepares @ " << now() << RESET << endl;
        auto pid = getpid();
        
        GiNaC_Parallel_RM["FCI-F"] = false;
        auto resf =
        GiNaC_Parallel(expResult.size(), [this](int idx)->ex {
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
                exmap eps_map;
                ex epn = ex(1)/111;
                for(auto epi : eps_lst) {
                    eps_map[epi.op(0)] = epn;
                    epn = epn / 13;
                }
                auto tmp = ft.subs(nReplacement).subs(eps_map).subs(CV(w1,w2)==w2).expand();
                if(is_a<add>(tmp)) {
                    for(auto item : tmp) {
                        if(!is_a<numeric>(item.subs(x(w)==1))) {
                            throw Error("CIPrepares: (!is_a<numeric>(item.subs(x(w)==1)))");
                        }
                        if(item.subs(x(w)==1) < 0) {
                            need_contour_deformation = true;
                            break;
                        }
                    }
                } else {
                    if(!is_a<numeric>(tmp.subs(x(w)==1))) {
                        throw Error("CIPrepares: (!is_a<numeric>(tmp.subs(x(w)==1)))");
                    }
                    if(tmp.subs(x(w)==1) < 0) need_contour_deformation = true;
                }
                if(!need_contour_deformation) ft = 1; //note the difference with SDPrepare
            } else if(!ft.has(x(w))) {
                ft = 1;
            }
            
            ft = collect_common_factors(ft);
            return lst{ kv.op(0), kv.op(1), ft};
            
        }, "FCI-F");
        

    //============================================================================================================
        lst fts;
        for(auto item : resf) {
            if(item.op(2).has(x(w))) {
                fts.append(item.op(2));
            }
        }
        fts.sort();
        fts.unique();
        
        exvector ftnvec;
        map<ex,int,ex_is_less> ftnmap;
        int ft_n = 1;
        FT_N_XN.remove_all();
                
        for(auto item : fts) {
            ftnvec.push_back(lst{item, ft_n});
            ftnmap[item] = ft_n;
            FT_N_XN.append(lst{item, ft_n, get_xy_from(item).size()});
            ft_n++;
        }
        //ftnvec item: lst { ft, ft-id }
        
        exvector res_vec;
        exmap cf_int;
        
        for(auto &item : resf) {
            auto ii = ex_to<lst>(item);
            if(ii.op(2)==1) {
                ii.append(-1);
            } else {
                int ft_n = ftnmap[item.op(2)];
                if(ft_n==0) throw Error("CIPrepares: ft_n==0, " + ex2str(item.op(2)));
                ii.append(ft_n);
            }
            if(ii.op(0).is_zero() || ii.op(1).subs(FTX(w1,w2)==1).is_zero()) continue;
            if(IBF==0) res_vec.push_back(ii);
            else if(IBF==1) { // by coefficient and F-term
                ex key = ii;
                key[1] = 1;
                cf_int[key] += ii.op(1);
            } else { // by F-term
                auto cvs = collect_lst(ii.op(0), [](const ex & e)->bool { return has_symbol(e); });
                for(auto const & cv : cvs) {
                    ex key = ii;
                    key[0] = cv.op(1);
                    key[1] = 1;
                    cf_int[key] += cv.op(0)*ii.op(1);
                }
            }
        }
        
        if(IBF!=0) {
            for(auto kv : cf_int) {
                lst ii = ex_to<lst>(kv.first);
                ii[1] = kv.second;
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
        GiNaC_Parallel_RM["FCI-C"] = false;
        GiNaC_Parallel(ftnvec.size(), [&ftnvec,pid](int idx)->ex {
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
                plRepl.append(PL(i) == Symbol(pl.str()));
            }
            
            ex DFs[fxs.size()], DDFs[fxs.size()][fxs.size()];
            for(int i=0; i<fxs.size(); i++) {
                auto df = diff_ex(ft, fxs[i]);
                DFs[i] = collect_common_factors(df);
                ostringstream ilaos;
                ilaos << "ilas[" << i << "]";
                symbol ila(ilaos.str());
                for(int j=0; j<fxs.size(); j++) {
                    auto ddf = diff_ex(DFs[i], fxs[j]);
                    DDFs[i][j] = collect_common_factors(ddf);
                }
            }
            
            if(!dir_exists(to_string(pid))) auto rc = system(("mkdir -p "+to_string(pid)).c_str());
            ostringstream cppfn, sofn;
            cppfn << pid << "/F" << ft_n << ".cpp";
            sofn << pid << "/F" << ft_n << ".o";
            std::ofstream ofs;
            ofs.open(cppfn.str(), ios::out);
            if (!ofs) throw Error("failed to open *.cpp file! (1)");
            
            lst cxRepl, czRepl;
            for (int i=0; i<fxs.size(); i++) {
                ostringstream sx, sz;
                sx << "x[" << i << "]";
                cxRepl.append(fxs[i] == Symbol(sx.str()));
                sz << "z[" << i << "]";
                czRepl.append(fxs[i] == Symbol(sz.str()));
            }

            /*----------------------------------------------*/
            ofs << "#include \"NFunctions.h\"" << endl;
            /*----------------------------------------------*/
            auto cppL = CppFormat(ofs, "L");
            auto cppQ = CppFormat(ofs, "Q");
            auto cppMP = CppFormat(ofs, "MP");
            
            // FL_fid
            ofs << "dREAL FL_" << ft_n << "(const dREAL* x, const dREAL *pl) {" << endl;
            ofs << "dREAL yy = ";
            EvalL(ft.subs(plRepl).subs(cxRepl)).print(cppL);
            ofs << ";" << endl;
            ofs << "return yy;" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // FQ_fid
            ofs << "qREAL FQ_" << ft_n << "(const qREAL* x, const qREAL *pl) {" << endl;
            ofs << "qREAL yy = ";
            EvalQ(ft.subs(plRepl).subs(cxRepl)).print(cppQ);
            ofs << ";" << endl;
            ofs << "return yy;" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // FMP_fid
            ofs << "mpREAL FMP_" << ft_n << "(const mpREAL* x, const mpREAL *pl) {" << endl;
            ofs << "mpREAL yy = ";
            EvalMP(ft.subs(plRepl).subs(cxRepl)).print(cppMP);
            ofs << ";" << endl;
            ofs << "return yy;" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // D's FL_fid
            ofs << "dREAL FL_" << ft_n << "(const int i, const dREAL* x, const dREAL *pl) {" << endl;
            for(int i=0; i<fxs.size(); i++) {
                ofs << "if("<<i<<"==i) return ";
                EvalL(DFs[i].subs(plRepl).subs(cxRepl)).print(cppL);
                ofs << ";" << endl;
            }
            ofs << "return 0;" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // D's FQ_fid
            ofs << "qREAL FQ_" << ft_n << "(const int i, const qREAL* x, const qREAL *pl) {" << endl;
            for(int i=0; i<fxs.size(); i++) {
                ofs << "if("<<i<<"==i) return ";
                EvalQ(DFs[i].subs(plRepl).subs(cxRepl)).print(cppQ);
                ofs << ";" << endl;
            }
            ofs << "return 0;" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // D's FMP_fid
            ofs << "mpREAL FMP_" << ft_n << "(const int i, const mpREAL* x, const mpREAL *pl) {" << endl;
            for(int i=0; i<fxs.size(); i++) {
                ofs << "if("<<i<<"==i) return ";
                EvalMP(DFs[i].subs(plRepl).subs(cxRepl)).print(cppMP);
                ofs << ";" << endl;
            }
            ofs << "return 0;" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // DD's FL_fid
            ofs << "dREAL FL_" << ft_n << "(const int i, const int j, const dREAL* x, const dREAL *pl) {" << endl;
            for(int i=0; i<fxs.size(); i++) {
            for(int j=0; j<fxs.size(); j++) {
                ofs << "if("<<i<<"==i && "<<j<<"==j) return ";
                EvalL(DDFs[i][j].subs(plRepl).subs(cxRepl)).print(cppL);
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
                EvalQ(DDFs[i][j].subs(plRepl).subs(cxRepl)).print(cppQ);
                ofs << ";" << endl;
            }}
            ofs << "return 0;" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // DD's FMP_fid
            ofs << "mpREAL FMP_" << ft_n << "(const int i, const int j, const mpREAL* x, const mpREAL *pl) {" << endl;
            for(int i=0; i<fxs.size(); i++) {
            for(int j=0; j<fxs.size(); j++) {
                ofs << "if("<<i<<"==i && "<<j<<"==j) return ";
                EvalMP(DDFs[i][j].subs(plRepl).subs(cxRepl)).print(cppMP);
                ofs << ";" << endl;
            }}
            ofs << "return 0;" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // X2ZL_fid
            ofs << "void X2ZL_" << ft_n << "(const dREAL* x, dCOMPLEX* z, dCOMPLEX* r, dREAL* dff, const dREAL* pl, const dREAL* las) {" << endl;
            ofs << "X2Z("<<fxs.size()<<",FL_"<<ft_n<<",FL_"<<ft_n<<",x,z,r,dff,pl,las);" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // X2ZQ_fid
            ofs << "void X2ZQ_" << ft_n << "(const qREAL* x, qCOMPLEX* z, qCOMPLEX* r, qREAL* dff, const qREAL* pl, const qREAL* las) {" << endl;
            ofs << "X2Z("<<fxs.size()<<",FQ_"<<ft_n<<",FQ_"<<ft_n<<",x,z,r,dff,pl,las);" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // X2ZMP_fid
            ofs << "void X2ZMP_" << ft_n << "(const mpREAL* x, mpCOMPLEX* z, mpCOMPLEX* r, mpREAL* dff, const mpREAL* pl, const mpREAL* las) {" << endl;
            ofs << "X2Z("<<fxs.size()<<",FMP_"<<ft_n<<",FMP_"<<ft_n<<",x,z,r,dff,pl,las);" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // MatL_id
            ofs << "void MatL_"<<ft_n<<"(dCOMPLEX* mat, const dREAL* x, const dREAL* dff, const dREAL* pl, const dREAL* las) {" << endl;
            ofs << "Mat("<<fxs.size()<<",FL_"<<ft_n<<",mat,x,dff,pl,las);" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // MatQ_fid
            ofs << "void MatQ_"<<ft_n<<"(qCOMPLEX *mat, const qREAL* x, const qREAL* dff, const qREAL *pl, const qREAL *las) {" << endl;
            ofs << "Mat("<<fxs.size()<<",FQ_"<<ft_n<<",mat,x,dff,pl,las);" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // MatMP_fid
            ofs << "void MatMP_"<<ft_n<<"(mpCOMPLEX *mat, const mpREAL* x, const mpREAL* dff, const mpREAL *pl, const mpREAL *las) {" << endl;
            ofs << "Mat("<<fxs.size()<<",FMP_"<<ft_n<<",mat,x,dff,pl,las);" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // for Minimization of F(z)-image, x[xn-1] is the lambda
            ofs << "extern \"C\" " << endl;
            ofs << "dREAL imgF_"<<ft_n<<"(const int xn, const dREAL* x, const dREAL *pl, const dREAL *las_in) {" << endl;
            ofs << "dREAL las[xn-1];" <<endl;
            ofs << "for(int i=0; i<xn-1; i++) las[i] = las_in[i]*x[xn-1];" <<endl;
            ofs << "dCOMPLEX z[xn], r[xn];" << endl;
            ofs << "dREAL dff[xn+1];" << endl;
            ofs << "X2ZL_"<<ft_n<<"(x,z,r,dff,pl,las);" << endl;
            ofs << "dCOMPLEX zf = ";
            EvalL(ft.subs(plRepl).subs(czRepl)).print(cppL);
            ofs << ";" << endl;
            ofs << "return -zf.imag()/x[xn-1];" << endl; // find max image part, check with 0
            ofs << "}" << endl;
            ofs << endl;
            
            // for Minimization of F
            ofs << "extern \"C\" " << endl;
            ofs << "dREAL minF_"<<ft_n<<"(const int xn, const dREAL* x, const dREAL *pl, const dREAL *las_in) {" << endl;
            ofs << "return FL_"<<ft_n<<"(x,pl);" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // for Minimization of -F
            ofs << "extern \"C\" " << endl;
            ofs << "dREAL minFM_"<<ft_n<<"(const int xn, const dREAL* x, const dREAL *pl, const dREAL *las_in) {" << endl;
            ofs << "return 0.L-FL_"<<ft_n<<"(x,pl);" << endl;
            ofs << "}" << endl;
            ofs << endl;
            
            // for Minimization of DF-i
            for(int i=0; i<fxs.size(); i++) {
                ofs << "extern \"C\" " << endl;
                ofs << "dREAL dirC_"<<ft_n<<"_"<<i<<"(const int xn, const dREAL* x, const dREAL *pl, const dREAL *las) {" << endl;
                ofs << "int i = " << i << ";" << endl;
                ofs << "dREAL yy = x[i]*(1-x[i])*FL_"<<ft_n<<"(i, x, pl);" << endl;
                ofs << "return -fabs(yy);" << endl;
                ofs << "}" << endl;
                ofs << endl;
            }
            
            ostringstream cmd;
            cmd << cpp << " -fPIC -c " << INC_FLAGS << " -o " << sofn.str() << " " << cppfn.str();
            auto rc = system(cmd.str().c_str());
            
            if(!file_exists(sofn.str().c_str())) {
                cmd.clear();
                cmd.str("");
                cmd << "cp " << sofn.str() << " . ";
                rc = system(cmd.str().c_str());
            }
            
            return 0;
        
        }, "FCI-C");
        
        bool hasF = (ftnvec.size()>0);
        if(hasF) {
            ostringstream sofn, cmd;
            if(key != "") {
                sofn << key << "F.so";
            } else {
                sofn << pid << "F.so";
            }
            cmd << cpp << " " << LIB_FLAGS <<  " -Wl,-rpath,. -rdynamic -fPIC -shared -lHepLib -lquadmath -lmpfr -lgmp " << " -o " << sofn.str() << " " << pid << "/F*.o";
            cmd << " -lHepLib -lquadmath -lmpfr -lgmp";
            cmd << " 1> /dev/null 2> /dev/null";
            rc = system(cmd.str().c_str());
            
            cmd.clear();
            cmd.str("");
            cmd << "rm -rf " << pid;
            if(!Debug) rc = system(cmd.str().c_str());
        }


    //============================================================================================================
        // Compile the null.o
        if(true) {
            ostringstream cmd;
            if(!dir_exists(to_string(pid))) cmd << "mkdir -p " << pid;
            rc = system(cmd.str().c_str());
            cmd.clear();
            cmd.str("");
            cmd << "echo ''>" << pid << "/null.cpp;";
            cmd << cpp << " -fPIC -c -o " << pid << "/null.o " << pid << "/null.cpp";
            rc = system(cmd.str().c_str());
        }

        // Prepare Integrand
        GiNaC_Parallel_RM["FCI-I"] = false;
        auto res =
        GiNaC_Parallel(res_vec.size(), [&res_vec,pid](int idx)->ex {
            // return lst{ no-x-result, xn, x-indepent prefactor, ft_n }
            // or     lst{ id(SD(D|Q)_id in .so), xn, x-indepent prefactor, ft_n }
            
            static symbol xwra("xwra");
            auto kvf = res_vec[idx];
            auto expr = kvf.op(1);
            auto xs = get_xy_from(expr);
            auto ft_n = kvf.op(3);
            bool hasF = (ft_n>0);
            
            if(xs.size()<1) {
                ostringstream cmd;
                cmd << "cp " << pid << "/null.o " << pid << "/" << idx << ".o";
                auto rc = system(cmd.str().c_str());
                return lst{ expr.subs(lst{FTX(w1,w2)==1}), xs.size(), kvf.op(0), -1};
            }
            
            if(xs.size()<1) xs.push_back(x(0));
            
            auto ft = kvf.op(2);
            auto fxs = get_xy_from(ft);
            
            exset ftxset;
            expr.find(FTX(w1,w2), ftxset);
            lst ftxlst;
            for(auto it : ftxset) ftxlst.append(it);
            expr = collect_ex(expr, FTX(w1,w2));
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
                cxRepl.append(fxs[i] == Symbol(xs.str()));
                czRepl.append(fxs[i] == Symbol(zs.str()));
                czzRepl.append(fxs[i] == Symbol(zzs.str()));
            }
            int count = fxs.size();
            for(auto xi : xs) {
                auto xii = xi.subs(czRepl);
                if(is_zero(xii-xi)) {
                    ostringstream xs, zs;
                    xs << "x[" << count << "]";
                    cxRepl.append(xi == Symbol(xs.str()));
                    czRepl.append(xi == Symbol(xs.str()));
                    czzRepl.append(xi == Symbol(xs.str()));
                    count++;
                }
            }
            if(count!=xs.size()) throw Error("CIPrepares: (count!=xs.size())");
            auto pls = get_pl_from(expr);
            int npls = pls.size()>0 ? ex_to<numeric>(pls[pls.size()-1].subs(lst{PL(w)==w})).to_int() : -1;
            for(int i=0; i<npls+1; i++) {
                ostringstream pl;
                pl << "pl[" << i << "]";
                plRepl.append(PL(i) == Symbol(pl.str()));
            }
            
            if(!dir_exists(to_string(pid))) auto rc = system(("mkdir -p "+to_string(pid)).c_str());
            ostringstream cppfn;
            cppfn << pid << "/" << idx << ".cpp";
            std::ofstream ofs;
            ofs.open(cppfn.str(), ios::out);
            if (!ofs) throw Error("failed to open *.cpp file! (2)");

            /*----------------------------------------------*/
            ofs << "#include \"NFunctions.h\"" << endl;
            /*----------------------------------------------*/
            if(hasF) {
                ofs << "qREAL FQ_"<<ft_n<<"(const qREAL*, const qREAL*);" << endl; // for FT only
                ofs << "dREAL FL_"<<ft_n<<"(const int, const dREAL*, const dREAL*);" << endl;
                ofs << "qREAL FQ_"<<ft_n<<"(const int, const qREAL*, const qREAL*);" << endl;
                ofs << "qREAL FMP_"<<ft_n<<"(const int, const mpREAL*, const mpREAL*);" << endl;
                ofs << "dREAL FL_"<<ft_n<<"(const int, const int, const dREAL*, const dREAL*);" << endl;
                ofs << "qREAL FQ_"<<ft_n<<"(const int, const int, const qREAL*, const qREAL*);" << endl;
                ofs << "qREAL FMP_"<<ft_n<<"(const int, const int, const mpREAL*, const mpREAL*);" << endl;
                ofs << "void X2ZL_"<<ft_n<<"(const dREAL*, dCOMPLEX*, dCOMPLEX*, dREAL*, const dREAL*, const dREAL*);" << endl;
                ofs << "void X2ZQ_"<<ft_n<<"(const qREAL*, qCOMPLEX*, qCOMPLEX*, qREAL*, const qREAL*, const qREAL*);" << endl;
                ofs << "void X2ZMP_"<<ft_n<<"(const mpREAL*, mpCOMPLEX*, mpCOMPLEX*, mpREAL*, const mpREAL*, const mpREAL*);" << endl;
                ofs << "void MatL_"<<ft_n<<"(dCOMPLEX*, const dREAL*, const dREAL*, const dREAL*, const dREAL*);" << endl;
                ofs << "void MatQ_"<<ft_n<<"(qCOMPLEX*, const qREAL*, const qREAL*, const qREAL*, const qREAL*);" << endl;
                ofs << "void MatMP_"<<ft_n<<"(mpCOMPLEX*, const mpREAL*, const mpREAL*, const mpREAL*, const mpREAL*);" << endl;
                ofs << endl << endl;
            }

            /*----------------------------------------------*/
            // long double
            /*----------------------------------------------*/
            auto cppL = CppFormat(ofs, "L");
            // alwasy export non-complex function
            if(true) {
                ofs << "extern \"C\" " << endl;
                ofs << "int SDD_"<<idx<<"(const unsigned int xn, const dREAL x[], const unsigned int yn, dREAL y[], const dREAL pl[], const dREAL las[]) {" << endl;
                if(Debug) {
                    auto tmp = expr.subs(FTX(w1,w2)==1).subs(cxRepl).subs(plRepl);
                    ofs << "//for debug, intg: " << endl << "//" << tmp << endl;
                }
                
                auto intg = expr.subs(FTX(w1,w2)==1);
                intg = xyz_pow_simplify(intg);
                bool hasF2 = intg.has(iEpsilon) || intg.has(I);
                
                // WickRotation
                hasF2 = hasF2 || intg.has(WRA(w));
                if(intg.has(WRA(w))) {
                    exset wras;
                    find(intg, WRA(w), wras);
                    if(wras.size()!=1) throw Error("CIPrepares: Too many WRA(w).");
                    auto wra = (*(wras.begin())).op(0);
                    intg = intg.subs(WRA(w)==xwra);
                    ofs << "dREAL wra = "; EvalL(wra).print(cppL); ofs << ";" << endl;
                    
                    exset pows_set;
                    find(intg, pow(w1,w2), pows_set);
                    find(intg, sqrt(w), pows_set);
                    lst pow_subs;
                    for(auto item : pows_set) {
                        if(item.has(xwra)) {
                            if(item.match(pow(w1,w2)) && !item.op(1).info(info_flags::integer)) {
                                pow_subs.append(item == exp(item.op(1)*log(item.op(0))));
                            } else if(item.match(sqrt(w))) {
                                pow_subs.append(item == exp(log(item.op(0))/2));
                            }
                        }
                    }
                    if(pow_subs.nops()>0) intg = exp_simplify(intg);
                    
                    exset logs_set;
                    find(intg, log(w), logs_set);
                    lst logs;
                    for(auto item : logs_set) {
                        if(item.has(xwra)) logs.append(item.op(0));
                    }
                    
                    if(logs.nops()>0) {
                        ofs << "int nlog = "<<logs.nops()<<";" << endl;
                        ofs << "dCOMPLEX CLog[nlog], LogZ[nlog][NRCLog+1];" << endl;
                        cseParser cse;
                        lst clogs = ex_to<lst>(cse.Parse(logs));
                        
                        ofs << "for(int ti=0; ti<=NRCLog; ti++) {" << endl;
                        ofs << "dREAL xwra;" << endl;
                        ofs << "if(ti==0) xwra = wra/(25*NRCLog);" << endl;
                        ofs << "else xwra = ti*wra/NRCLog;" << endl;
                        ofs << "dCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                        for(auto kv : cse.os()) {
                            ofs <<cse.oc<< "["<<kv.first<<"] = ";
                            EvalL(kv.second.subs(cxRepl).subs(plRepl)).print(cppL);
                            ofs << ";" << endl;
                        }
                        
                        exmap log_subs;
                        for(int i=0; i<clogs.nops(); i++) {
                            ofs << "LogZ["<<i<<"][ti] = ";
                            EvalL(clogs.op(i).subs(cxRepl).subs(plRepl)).print(cppL);
                            ofs << ";" << endl;
                            log_subs[log(logs.op(i))] = Symbol("CLog["+to_string(i)+"]");
                        }
                        ofs << "}" << endl;
                        ofs << "for(int li=0; li<nlog; li++) CLog[li] = RCLog(LogZ[li],NRCLog);" << endl;
                        intg = intg.subs(log_subs);
                    }
                    ofs << "dREAL xwra = wra;" << endl; 
                }
                
                cseParser cse;
                intg = cse.Parse(intg);
                if(hasF2) ofs << "dCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                else ofs << "dREAL "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                for(auto kv : cse.os()) {
                    ofs <<cse.oc<< "["<<kv.first<<"] = ";
                    EvalL(kv.second.subs(cxRepl).subs(plRepl)).print(cppL);
                    ofs << ";" << endl;
                }
                
                ofs << "dCOMPLEX yy = ";
                EvalL(intg.subs(cxRepl).subs(plRepl)).print(cppL);
                ofs << ";" << endl;
                
                ofs << "y[0] = yy.real();" << endl;
                ofs << "y[1] = yy.imag();" << endl;
                ofs << "return 0;" << endl;
                ofs << "}" << endl;
                ofs << endl;
            }
            
            if(hasF) {
                ofs << "extern \"C\" " << endl;
                ofs << "int CSDD_"<<idx<<"(const unsigned int xn, const dREAL x[], const unsigned int yn, dREAL y[], const dREAL pl[], const dREAL las[]) {" << endl;
                ofs << "dREAL x0[xn];" << endl;
                ofs << "dCOMPLEX z[xn],zz[xn],r[xn];" << endl;
                ofs << "dREAL dff[xn+1];" << endl;
                ofs << "dCOMPLEX yy=0, ytmp, det;" << endl;
                ofs << "int ii, nfxs="<<fxs.size()<<";" << endl;
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
                    ofs  << "det = MatDet(mat, nfxs);" << endl;
                    
                    ex intg = kv.second;
                    
                    if(true) { // RCLog
                        exset pows_set;
                        find(intg, pow(w1,w2), pows_set);
                        find(intg, sqrt(w), pows_set);
                        lst pow_subs;
                        for(auto item : pows_set) {
                            if(!item.op(0).match(x(w))) {
                                if(item.match(pow(w1,w2)) && !item.op(1).info(info_flags::integer)) {
                                    pow_subs.append(item == exp(item.op(1)*log(item.op(0))));
                                } else if(item.match(sqrt(w))) {
                                    pow_subs.append(item == exp(log(item.op(0))/2));
                                }
                            }
                        }
                        if(pow_subs.nops()>0) intg = exp_simplify(intg);
                    
                        exset logs_set;
                        find(intg, log(w), logs_set);
                        if(logs_set.size()>0) {
                            lst logs;
                            for(auto item : logs_set) {
                                if(!item.op(0).match(x(w))) logs.append(item.op(0));
                            }
                            if(logs.nops()>0) {
                                ofs << "int nlog = "<<logs.nops()<<";" << endl;
                                ofs << "dCOMPLEX CLog[nlog], LogZ[nlog][NRCLog+1];" << endl;
                                cseParser cse;
                                lst clogs = ex_to<lst>(cse.Parse(logs));
                                
                                ofs << "for(int ti=0; ti<=NRCLog; ti++) {" << endl;
                                ofs << "if(ti==0) { for(int i=0; i<xn; i++) zz[i] = complex<dREAL>(z[i].real(), z[i].imag()/(25*NRCLog)); }" << endl;
                                ofs << "else { for(int i=0; i<xn; i++) zz[i] = complex<dREAL>(z[i].real(), ti*z[i].imag()/NRCLog); }" << endl;
                                ofs << "dCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                                for(auto kv : cse.os()) {
                                    ofs <<cse.oc<< "["<<kv.first<<"] = ";
                                    EvalL(kv.second.subs(czzRepl).subs(plRepl)).print(cppL);
                                    ofs << ";" << endl;
                                }
                                
                                exmap log_subs;
                                for(int i=0; i<clogs.nops(); i++) {
                                    ofs << "LogZ["<<i<<"][ti] = ";
                                    EvalL(clogs.op(i).subs(czzRepl).subs(plRepl)).print(cppL);
                                    ofs << ";" << endl;
                                    log_subs[log(logs.op(i))] = Symbol("CLog["+to_string(i)+"]");
                                }
                                ofs << "}" << endl;
                                ofs << "for(int li=0; li<nlog; li++) CLog[li] = RCLog(LogZ[li],NRCLog);" << endl;
                                intg = intg.subs(log_subs);
                            }
                        }
                    }
                    
                    cseParser cse;
                    intg = cse.Parse(intg);
                    ofs << "dCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                    for(auto kv : cse.os()) {
                        ofs <<cse.oc<< "["<<kv.first<<"] = ";
                        EvalL(kv.second.subs(czRepl).subs(plRepl)).print(cppL);
                        ofs << ";" << endl;
                    }
                    
                    ofs << "ytmp = ";
                    EvalL(intg.subs(czRepl).subs(plRepl)).print(cppL);
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
            auto cppQ = CppFormat(ofs, "Q");
            
            // always export non-complex function
            if(true) {
                ofs << "extern \"C\" " << endl;
                ofs << "int SDQ_"<<idx<<"(const unsigned int xn, const qREAL x[], const int unsigned yn, qREAL y[], const qREAL pl[], const qREAL las[]) {" << endl;
                
                auto intg = expr.subs(FTX(w1,w2)==1);
                intg = xyz_pow_simplify(intg);
                bool hasF2 = intg.has(iEpsilon) || intg.has(I);
                
                // WickRotation
                hasF2 = hasF2 || intg.has(WRA(w));
                if(intg.has(WRA(w))) {
                    exset wras;
                    find(intg, WRA(w), wras);
                    if(wras.size()!=1) throw Error("CIPrepares: Too many WRA(w).");
                    auto wra = (*(wras.begin())).op(0);
                    intg = intg.subs(WRA(w)==xwra);
                    ofs << "qREAL wra = "; EvalQ(wra).print(cppQ); ofs << ";" << endl;
                    
                    exset pows_set;
                    find(intg, pow(w1,w2), pows_set);
                    find(intg, sqrt(w), pows_set);
                    lst pow_subs;
                    for(auto item : pows_set) {
                        if(item.has(xwra)) {
                            if(item.match(pow(w1,w2)) && !item.op(1).info(info_flags::integer)) {
                                pow_subs.append(item == exp(item.op(1)*log(item.op(0))));
                            } else if(item.match(sqrt(w))) {
                                pow_subs.append(item == exp(log(item.op(0))/2));
                            }
                        }
                    }
                    if(pow_subs.nops()>0) intg = exp_simplify(intg);
                    
                    exset logs_set;
                    find(intg, log(w), logs_set);
                    lst logs;
                    for(auto item : logs_set) {
                        if(item.has(xwra)) logs.append(item.op(0));
                    }
                    
                    if(logs.nops()>0) {
                        ofs << "int nlog = "<<logs.nops()<<";" << endl;
                        ofs << "qCOMPLEX CLog[nlog], LogZ[nlog][NRCLog+1];" << endl;
                        cseParser cse;
                        lst clogs = ex_to<lst>(cse.Parse(logs));
                        
                        ofs << "for(int ti=0; ti<=NRCLog; ti++) {" << endl;
                        ofs << "qREAL xwra;" << endl;
                        ofs << "if(ti==0) xwra = wra/(25*NRCLog);" << endl;
                        ofs << "else xwra = ti*wra/NRCLog;" << endl;
                        ofs << "qCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                        for(auto kv : cse.os()) {
                            ofs <<cse.oc<< "["<<kv.first<<"] = ";
                            EvalQ(kv.second.subs(cxRepl).subs(plRepl)).print(cppQ);
                            ofs << ";" << endl;
                        }
                        
                        exmap log_subs;
                        for(int i=0; i<clogs.nops(); i++) {
                            ofs << "LogZ["<<i<<"][ti] = ";
                            EvalQ(clogs.op(i).subs(cxRepl).subs(plRepl)).print(cppQ);
                            ofs << ";" << endl;
                            log_subs[log(logs.op(i))] = Symbol("CLog["+to_string(i)+"]");
                        }
                        ofs << "}" << endl;
                        ofs << "for(int li=0; li<nlog; li++) CLog[li] = RCLog(LogZ[li],NRCLog);" << endl;
                        intg = intg.subs(log_subs);
                    } 
                    ofs << "qREAL xwra = wra;" << endl;
                }
                
                cseParser cse;
                intg = cse.Parse(intg);
                if(hasF2) ofs << "qCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                else ofs << "qREAL "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                for(auto kv : cse.os()) {
                    ofs <<cse.oc<< "["<<kv.first<<"] = ";
                    EvalQ(kv.second.subs(cxRepl).subs(plRepl)).print(cppQ);
                    ofs << ";" << endl;
                }
                
                ofs << "qCOMPLEX yy = ";
                EvalQ(intg.subs(cxRepl).subs(plRepl)).print(cppQ);
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
                ofs << "qREAL dff[xn+1];" << endl;
                ofs << "qCOMPLEX yy=0, ytmp, det;" << endl;
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
                    ofs  << "det = MatDet(mat, nfxs);" << endl;
                    
                    ex intg = kv.second;
                    
                    if(true) { // RCLog
                        exset pows_set;
                        find(intg, pow(w1,w2), pows_set);
                        find(intg, sqrt(w), pows_set);
                        lst pow_subs;
                        for(auto item : pows_set) {
                            if(!item.op(0).match(x(w))) {
                                if(item.match(pow(w1,w2)) && !item.op(1).info(info_flags::integer)) {
                                    pow_subs.append(item == exp(item.op(1)*log(item.op(0))));
                                } else if(item.match(sqrt(w))) {
                                    pow_subs.append(item == exp(log(item.op(0))/2));
                                }
                            }
                        }
                        if(pow_subs.nops()>0) intg = exp_simplify(intg);
                        
                        exset logs_set;
                        find(intg, log(w), logs_set);
                        if(logs_set.size()>0) {
                            lst logs;
                            for(auto item : logs_set) {
                                if(!item.op(0).match(x(w))) logs.append(item.op(0));
                            }
                            if(logs.nops()>0) {
                                ofs << "int nlog = "<<logs.nops()<<";" << endl;
                                ofs << "qCOMPLEX CLog[nlog], LogZ[nlog][NRCLog+1];" << endl;
                                cseParser cse;
                                lst clogs = ex_to<lst>(cse.Parse(logs));
                                
                                ofs << "for(int ti=0; ti<=NRCLog; ti++) {" << endl;
                                ofs << "if(ti==0) { for(int i=0; i<xn; i++) zz[i] = crealq(z[i]) + 1.Qi*cimagq(z[i])/(25.Q*NRCLog); }" << endl;
                                ofs << "else { for(int i=0; i<xn; i++) zz[i] = crealq(z[i]) + 1.Qi*ti*cimagq(z[i])/NRCLog; }" << endl;
                                ofs << "qCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                                for(auto kv : cse.os()) {
                                    ofs <<cse.oc<< "["<<kv.first<<"] = ";
                                    EvalQ(kv.second.subs(czzRepl).subs(plRepl)).print(cppQ);
                                    ofs << ";" << endl;
                                }
                                
                                exmap log_subs;
                                for(int i=0; i<clogs.nops(); i++) {
                                    ofs << "LogZ["<<i<<"][ti] = ";
                                    EvalQ(clogs.op(i).subs(czzRepl).subs(plRepl)).print(cppQ);
                                    ofs << ";" << endl;
                                    log_subs[log(logs.op(i))] = Symbol("CLog["+to_string(i)+"]");
                                }
                                ofs << "}" << endl;
                                ofs << "for(int li=0; li<nlog; li++) CLog[li] = RCLog(LogZ[li],NRCLog);" << endl;
                                intg = intg.subs(log_subs);
                            }
                        }
                    }
                    
                    cseParser cse;
                    intg = cse.Parse(intg);
                    ofs << "qCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                    for(auto kv : cse.os()) {
                        ofs <<cse.oc<< "["<<kv.first<<"] = ";
                        EvalQ(kv.second.subs(czRepl).subs(plRepl)).print(cppQ);
                        ofs << ";" << endl;
                    }
                    
                    ofs << "ytmp = ";
                    EvalQ(intg.subs(czRepl).subs(plRepl)).print(cppQ);
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
            auto cppMP =  CppFormat(ofs, "MP");
            
            // always export non-complex function
            if(true) {
                ofs << "extern \"C\" " << endl;
                ofs << "int SDMP_"<<idx<<"(const unsigned int xn, const mpREAL x[], const unsigned int yn, mpREAL y[], const mpREAL pl[], const mpREAL las[]) {" << endl;
                auto intg = expr.subs(FTX(w1,w2)==1);
                intg = xyz_pow_simplify(intg);
                bool hasF2 = intg.has(iEpsilon) || intg.has(I);
                
                // WickRotation
                hasF2 = hasF2 || intg.has(WRA(w));
                if(intg.has(WRA(w))) {
                    exset wras;
                    find(intg, WRA(w), wras);
                    if(wras.size()!=1) throw Error("CIPrepares: Too many WRA(w).");
                    auto wra = (*(wras.begin())).op(0);
                    intg = intg.subs(WRA(w)==xwra);
                    ofs << "mpREAL wra = "; EvalMP(wra).print(cppMP); ofs << ";" << endl;
                    
                    exset pows_set;
                    find(intg, pow(w1,w2), pows_set);
                    find(intg, sqrt(w), pows_set);
                    lst pow_subs;
                    for(auto item : pows_set) {
                        if(item.has(xwra)) {
                            if(item.match(pow(w1,w2)) && !item.op(1).info(info_flags::integer)) {
                                pow_subs.append(item == exp(item.op(1)*log(item.op(0))));
                            } else if(item.match(sqrt(w))) {
                                pow_subs.append(item == exp(log(item.op(0))/2));
                            }
                        }
                    }
                    if(pow_subs.nops()>0) intg = exp_simplify(intg);
                    
                    exset logs_set;
                    find(intg, log(w), logs_set);
                    lst logs;
                    for(auto item : logs_set) {
                        if(item.has(xwra)) logs.append(item.op(0));
                    }
                    
                    if(logs.nops()>0) {
                        ofs << "int nlog = "<<logs.nops()<<";" << endl;
                        ofs << "mpCOMPLEX CLog[nlog], LogZ[nlog][NRCLog+1];" << endl;
                        cseParser cse;
                        lst clogs = ex_to<lst>(cse.Parse(logs));
                        
                        ofs << "for(int ti=0; ti<=NRCLog; ti++) {" << endl;
                        ofs << "mpREAL xwra;" << endl;
                        ofs << "if(ti==0) xwra = wra/(25*NRCLog);" << endl;
                        ofs << "else xwra = ti*wra/NRCLog;" << endl; 
                        ofs << "mpCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                        for(auto kv : cse.os()) {
                            ofs <<cse.oc<< "["<<kv.first<<"] = ";
                            EvalMP(kv.second.subs(cxRepl).subs(plRepl)).print(cppMP);
                            ofs << ";" << endl;
                        }
                        
                        exmap log_subs;
                        for(int i=0; i<clogs.nops(); i++) {
                            ofs << "LogZ["<<i<<"][ti] = ";
                            EvalMP(clogs.op(i).subs(cxRepl).subs(plRepl)).print(cppMP);
                            ofs << ";" << endl;
                            log_subs[log(logs.op(i))] = Symbol("CLog["+to_string(i)+"]");
                        }
                        ofs << "}" << endl;
                        ofs << "for(int li=0; li<nlog; li++) CLog[li] = RCLog(LogZ[li],NRCLog);" << endl;
                        intg = intg.subs(log_subs);
                    } 
                    ofs << "mpREAL xwra = wra;" << endl;
                }
                
                cseParser cse;
                intg = cse.Parse(intg);
                if(hasF2) ofs << "mpCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                else ofs << "mpREAL "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                for(auto kv : cse.os()) {
                    ofs <<cse.oc<< "["<<kv.first<<"] = ";
                    EvalMP(kv.second.subs(cxRepl).subs(plRepl)).print(cppMP);
                    ofs << ";" << endl;
                }
                
                ofs << "mpCOMPLEX yy = ";
                EvalMP(intg.subs(cxRepl).subs(plRepl)).print(cppMP);
                ofs << ";" << endl;
                ofs << "y[0] = yy.real();" << endl;
                ofs << "y[1] = yy.imag();" << endl;
                ofs << "return 0;" << endl;
                ofs << "}" << endl;
                ofs << endl;
            }
            
            if(hasF) {
                ofs << "extern \"C\" " << endl;
                ofs << "int CSDMP_"<<idx<<"(const unsigned int xn, const mpREAL x[], const unsigned int yn, mpREAL y[], const mpREAL pl[], const mpREAL las[]) {" << endl;
                ofs << "mpREAL x0[xn];" << endl;
                ofs << "mpCOMPLEX z[xn],zz[xn],r[xn];" << endl;
                ofs << "mpREAL dff[xn+1];" << endl;
                ofs << "mpCOMPLEX yy=mpREAL(0), ytmp, det;" << endl;
                ofs << "int ii, nfxs="<<fxs.size()<<";" << endl;
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
                    ofs  << "det = MatDet(mat, nfxs);" << endl;
                    
                    ex intg = kv.second;
                    
                    if(true) { // RCLog
                        exset pows_set;
                        find(intg, pow(w1,w2), pows_set);
                        find(intg, sqrt(w), pows_set);
                        lst pow_subs;
                        for(auto item : pows_set) {
                            if(!item.op(0).match(x(w))) {
                                if(item.match(pow(w1,w2)) && !item.op(1).info(info_flags::integer)) {
                                    pow_subs.append(item == exp(item.op(1)*log(item.op(0))));
                                } else if(item.match(sqrt(w))) {
                                    pow_subs.append(item == exp(log(item.op(0))/2));
                                }
                            }
                        }
                        if(pow_subs.nops()>0) intg = exp_simplify(intg);
                        
                        exset logs_set;
                        find(intg, log(w), logs_set);
                        if(logs_set.size()>0) {
                            lst logs;
                            for(auto item : logs_set) {
                                if(!item.op(0).match(x(w))) logs.append(item.op(0));
                            }
                            if(logs.nops()>0) {
                                ofs << "int nlog = "<<logs.nops()<<";" << endl;
                                ofs << "mpCOMPLEX CLog[nlog], LogZ[nlog][NRCLog+1];" << endl;
                                cseParser cse;
                                lst clogs = ex_to<lst>(cse.Parse(logs));
                                
                                ofs << "for(int ti=0; ti<=NRCLog; ti++) {" << endl;
                                ofs << "if(ti==0) { for(int i=0; i<xn; i++) zz[i] = complex<mpREAL>(z[i].real(), z[i].imag()/(25*NRCLog)); }" << endl;
                                ofs << "else { for(int i=0; i<xn; i++) zz[i] = complex<mpREAL>(z[i].real(), ti*z[i].imag()/NRCLog); }" << endl;
                                ofs << "mpCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                                for(auto kv : cse.os()) {
                                    ofs <<cse.oc<< "["<<kv.first<<"] = ";
                                    EvalMP(kv.second.subs(czzRepl).subs(plRepl)).print(cppMP);
                                    ofs << ";" << endl;
                                }
                                
                                exmap log_subs;
                                for(int i=0; i<clogs.nops(); i++) {
                                    ofs << "LogZ["<<i<<"][ti] = ";
                                    EvalMP(clogs.op(i).subs(czzRepl).subs(plRepl)).print(cppMP);
                                    ofs << ";" << endl;
                                    log_subs[log(logs.op(i))] = Symbol("CLog["+to_string(i)+"]");
                                }
                                ofs << "}" << endl;
                                ofs << "for(int li=0; li<nlog; li++) CLog[li] = RCLog(LogZ[li],NRCLog);" << endl;
                                intg = intg.subs(log_subs);
                            }
                        }
                    }
                    
                    cseParser cse;
                    intg = cse.Parse(intg);
                    ofs << "mpCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                    for(auto kv : cse.os()) {
                        ofs <<cse.oc<< "["<<kv.first<<"] = ";
                        EvalMP(kv.second.subs(czRepl).subs(plRepl)).print(cppMP);
                        ofs << ";" << endl;
                    }
                    
                    ofs << "ytmp = ";
                    EvalMP(intg.subs(czRepl).subs(plRepl)).print(cppMP);
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
            
            ofs.close();
            ostringstream ofn, cmd;
            ofn << pid << "/" << idx << ".o";
            cmd << cpp << " -pipe -fPIC " << INC_FLAGS <<  " -c -o " << ofn.str() << " " << cppfn.str();
            auto rc = system(cmd.str().c_str());
            if(!Debug) remove(cppfn.str().c_str());
            return lst{ idx, xs.size(), kvf.op(0), ft_n };
        }, "FCI-I");
        
//        if(true) { // try other complilation method
//            ostringstream cmd;
//            cmd << "cd " << pid << ";find . -name '*.cpp' | xargs -n 100 -P " << CpuCores() << " " << cpp << " -pipe -fPIC " << INC_FLAGS <<  " -c";
//            auto rc = system(cmd.str().c_str());
//            //if(!Debug) remove(cppfn.str().c_str());
//        }
                
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
            ar.archive_ex(soLimit,"soLimit");
            ar.archive_ex(eps_lst,"eps_lst");
            ofstream out(garfn.str());
            out << ar;
            out.close();
            cmd.clear();
            cmd.str("");
            cmd << "rm -f " << key << "X*.so";
            rc = system(cmd.str().c_str());
        } else {
            fsofn << pid << "F.so";
            sofn << pid << ".so";
            for(auto &item : res) ciResult.push_back(ex_to<lst>(item));
            cmd.clear();
            cmd.str("");
            cmd << "rm -f " << pid << "X*.so";
            rc = system(cmd.str().c_str());
        }
        
        int res_size = res.size();
        if(soLimit<100) soLimit = 100;
        if(res_size>soLimit) {
            cmd.clear();
            cmd.str("");
            cmd << cpp << " " << LIB_FLAGS <<  " -Wl,-rpath,. -rdynamic -fPIC -shared -lHepLib -lquadmath -lmpfr -lgmp";
            if(hasF) cmd << " " << fsofn.str();
            cmd << " -o " << sofn.str() << " $(seq -f '" << pid << "/%g.o' 0 " << (soLimit-1) << ")";
            cmd << " -lHepLib -lquadmath -lmpfr -lgmp";
            cmd << " 1> /dev/null 2> /dev/null";
            rc = system(cmd.str().c_str());
            
            for(int n=1; true; n++) {
                int start = n*soLimit;
                int end = (n+1)*soLimit-1;
                if(end>res_size-1) end = res_size-1;
                sofn.clear();
                sofn.str("");
                if(key != "") sofn << key << "X" << n << ".so";
                else sofn << pid << "X" << n << ".so";
                cmd.clear();
                cmd.str("");
                cmd << cpp << " " << LIB_FLAGS <<  " -Wl,-rpath,. -rdynamic -fPIC -shared -lHepLib -lquadmath -lmpfr -lgmp";
                if(hasF) cmd << " " << fsofn.str();
                cmd << " -o " << sofn.str() << " $(seq -f '" << pid << "/%g.o' " << start << " " << end << ")";
                if(hasF) cmd << " " << fsofn.str();
                cmd << " -lHepLib -lquadmath -lmpfr -lgmp";
                cmd << " 1> /dev/null 2> /dev/null";
                rc = system(cmd.str().c_str());
                if(end>=res_size-1) break;
            }
        } else {
            cmd.clear();
            cmd.str("");
            cmd << cpp << " " << LIB_FLAGS <<  " -Wl,-rpath,. -rdynamic -fPIC -shared -lHepLib -lquadmath -lmpfr -lgmp";
            if(hasF) cmd << " " << fsofn.str();
            cmd << " -o " << sofn.str() << " " << pid << "/*.o";
            if(hasF) cmd << " " << fsofn.str();
            cmd << " -lHepLib -lquadmath -lmpfr -lgmp";
            cmd << " 1> /dev/null 2> /dev/null";
            rc = system(cmd.str().c_str());
        }
        cmd.clear();
        cmd.str("");
        if(Debug) cmd << "rm -rf " << key << "_debug;mv -f " << pid << " " << key << "_debug;rm -f " << key << "_debug/*.o";
        else cmd << "rm -rf " << pid;
        rc = system(cmd.str().c_str());

    }
    

}
