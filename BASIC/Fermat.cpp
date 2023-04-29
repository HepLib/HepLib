// Functions related to fermat program
#include "BASIC.h"
#include "cln/cln.h"

namespace HepLib {

    Fermat& Fermat::get() {
        static map<pid_t, Fermat> fermat_map;
        auto pid = getpid();
        if((fermat_map.find(pid)==fermat_map.end())) { // init section
            fermat_map[pid].Init();
            fermat_map[pid].vmax = 0;
        }
        return fermat_map[pid];
    }

    /**
     * @brief return the numerator and denominator after normalization
     * @param expr the input expression
     * @param dfactor true for factorize on the denominator
     * @return a list of { numer, denom }
     */
    ex numer_denom_fermat(const ex & expr, bool dfactor) {
        int _fermat_using_array = fermat_using_array;
        bool use_ncheck = false;
        
        Fermat &fermat = Fermat::get();
        int &v_max = fermat.vmax;
        
        auto expr_in = expr;
        exmap map_rat;
        expr_in = expr_in.to_rational(map_rat);
        
        lst rep_vs;
        if(true) {
            map<ex,long long,ex_is_less> s2c;
            for(const_preorder_iterator i = expr_in.preorder_begin(); i != expr_in.preorder_end(); ++i) {
                if(is_a<symbol>(*i)) s2c[*i]++;
            }
            exvector sv1, sv2, sv3; // sv2:fermat_weight, sv3:map_rat
            for(auto kv : s2c) {
                auto fw = fermat_weight.find(kv.first);
                if(fw!=fermat_weight.end()) sv2.push_back(lst{fw->second, fw->first});
                else if(map_rat.find(kv.first)!=map_rat.end()) sv3.push_back(lst{kv.second, kv.first});
                else sv1.push_back(lst{kv.second, kv.first});
            }
            sort_vec(sv1);
            sort_vec(sv2);
            sort_vec(sv3);
            for(auto sv : sv1) rep_vs.append(sv.op(1));
            for(auto sv : sv2) rep_vs.append(sv.op(1));
            for(auto sv : sv3) rep_vs.append(sv.op(1));
        }
                
        exmap v2f, f2v;
        exmap nn_map;
        auto nn_pi1 = cln::nextprobprime(3);
        auto nn_pi2 = cln::nextprobprime(3);
        int fvi = 0;
        for(auto vi : rep_vs) {
            auto name = "v" + to_string(fvi);
            Symbol s(name);
            v2f[vi] = s;
            f2v[s] = vi;
            fvi++;
            nn_pi1 = cln::nextprobprime(nn_pi2+1);
            nn_pi2 = cln::nextprobprime(nn_pi1+1);
            nn_map[s] = numeric(nn_pi1)/numeric(nn_pi2);
        }
        
        stringstream ss;
        if(fvi>111) {
            cout << rep_vs << endl;
            throw Error("Fermat: Too many variables.");
        }
        if(fvi>v_max) {
            for(int i=v_max; i<fvi; i++) ss << "&(J=v" << i << ");" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
            v_max = fvi;
        }
        
        ex nn_chk=0, num, den;
        lst item;
        if(!is_a<add>(expr_in)) item = lst{ expr_in };
        else for(auto ii : expr_in) item.append(ii);
        //sort_lst(item); // no need
        if(item.nops()>999999) _fermat_using_array = 0;
        if(_fermat_using_array) ss << "Array m[" << item.nops() << "];" << endl;
        else ss << "res:=0;" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        for(int i=0; i<item.nops(); i++) {
            ex tt = item.op(i).subs(v2f, subs_options::no_pattern);
            if(use_ncheck) nn_chk += tt.subs(nn_map, subs_options::no_pattern);
            if(_fermat_using_array) ss << "m[" << (i+1) << "]:=";
            else ss << "item:=";
            ss << tt << ";" << endl;
            if(!_fermat_using_array) ss << "res:=*res+*item;" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
        }
        if(_fermat_using_array) {
            if(_fermat_using_array==1) ss << "res:=Sumup([m]);" << endl;
            else if(_fermat_using_array==2) ss << "res:=Sigma<i=1,"<<item.nops()<<">(*m[i]);" << endl;
            else {
                ss << "n:=" << item.nops() << ";" << endl;
                ss << "while n>1 do n2:=n\\2; for i=1,n2 do m[i]:=*m[i]+*m[n+1-i] od; &_G; if (n|2)=0 then n:=n2 else n:=n2+1 fi od;" << endl;
                ss << "res:=*m[1];" << endl;
            }
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
        }
        
        static string bstr("[-begin-]"), estr("[-end-]");
        ss << "&(U=1);" << endl; // ugly printing, the whitespace matters
        ss << "!('" <<bstr<< "','{',Numer(^res),',',Denom(^res),'}','" <<estr<< "')" << endl;
        auto ostr = fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        // note the order,(normal_fermat will be called again in factor_form)
        ss << "&(U=0);" << endl; // disable ugly printing
        if(_fermat_using_array) ss << "@(res,[m]);" << endl;
        else ss << "@(res,item);" << endl;
        ss << "&_G;" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");

        // make sure last char is 0
        if(ostr[ostr.length()-1]!='0') throw Error("fermat_together: last char is NOT 0.");
        ostr = ostr.substr(0, ostr.length()-1);
        auto cpos = ostr.find(bstr);
        if(cpos==string::npos) throw Error(bstr+" NOT Found.");
        ostr = ostr.substr(cpos+bstr.length(),string::npos);
        cpos = ostr.find(estr);
        if(cpos==string::npos) throw Error(estr+" NOT Found.");
        ostr = ostr.substr(0,cpos);
        string_trim(ostr);

        symtab st;
        Parser fp(st);
        auto ret = fp.Read(ostr);
        ReShare(ret,expr);
        num = ret.op(0);
        if(dfactor) den = factor_form(ret.op(1));
        else den = ret.op(1);
        //fermat.Exit();
        
        if(use_ncheck) {
            auto nn_ret = subs(num/den,nn_map);
            if(nn_chk-nn_ret!=0) {
                cout << nn_chk << " : " << nn_ret << endl;
                throw Error("fermat_together: N Check Failed.");
            }
        }
        
        num = num.subs(f2v,subs_options::no_pattern).subs(map_rat,nopat);
        den = den.subs(f2v,subs_options::no_pattern).subs(map_rat,nopat);
        return lst{num, den};
        
    }
    
    ex numer_fermat(const ex & expr) {
        int _fermat_using_array = fermat_using_array;
        bool use_ncheck = false;

        Fermat &fermat = Fermat::get();
        int &v_max = fermat.vmax;
        
        auto expr_in = expr;
        exmap map_rat;
        expr_in = expr_in.to_polynomial(map_rat);
        
        lst rep_vs;
        if(true) {
            map<ex,long long,ex_is_less> s2c;
            for(const_preorder_iterator i = expr_in.preorder_begin(); i != expr_in.preorder_end(); ++i) {
                if(is_a<symbol>(*i)) s2c[*i]++;
            }
            exvector sv1, sv2, sv3; // sv2:fermat_weight, sv3:map_rat
            for(auto kv : s2c) {
                auto fw = fermat_weight.find(kv.first);
                if(fw!=fermat_weight.end()) sv2.push_back(lst{fw->second, fw->first});
                else if(map_rat.find(kv.first)!=map_rat.end()) sv3.push_back(lst{kv.second, kv.first});
                else sv1.push_back(lst{kv.second, kv.first});
            }
            sort_vec(sv1);
            sort_vec(sv2);
            sort_vec(sv3);
            for(auto sv : sv1) rep_vs.append(sv.op(1));
            for(auto sv : sv2) rep_vs.append(sv.op(1));
            for(auto sv : sv3) rep_vs.append(sv.op(1));
        }
                
        exmap v2f, f2v;
        exmap nn_map;
        auto nn_pi1 = cln::nextprobprime(3);
        auto nn_pi2 = cln::nextprobprime(3);
        int fvi = 0;
        for(auto vi : rep_vs) {
            auto name = "v" + to_string(fvi);
            Symbol s(name);
            v2f[vi] = s;
            f2v[s] = vi;
            fvi++;
            nn_pi1 = cln::nextprobprime(nn_pi2+1);
            nn_pi2 = cln::nextprobprime(nn_pi1+1);
            nn_map[s] = numeric(nn_pi1)/numeric(nn_pi2);
        }
        
        stringstream ss;
        if(fvi>111) {
            cout << rep_vs << endl;
            throw Error("Fermat: Too many variables.");
        }
        if(fvi>v_max) {
            for(int i=v_max; i<fvi; i++) ss << "&(J=v" << i << ");" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
            v_max = fvi;
        }
        
        ex nn_chk=0, nres;
        lst item;
        if(!is_a<add>(expr_in)) item = lst{ expr_in };
        else for(auto ii : expr_in) item.append(ii);
        //sort_lst(item); // no need
        if(item.nops()>999999) _fermat_using_array = 0;
        if(_fermat_using_array) ss << "Array m[" << item.nops() << "];" << endl;
        else ss << "res:=0;" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        for(int i=0; i<item.nops(); i++) {
            ex tt = item.op(i).subs(v2f, subs_options::no_pattern);
            if(use_ncheck) nn_chk += tt.subs(nn_map, subs_options::no_pattern);
            if(_fermat_using_array) ss << "m[" << (i+1) << "]:=";
            else ss << "item:=";
            ss << tt << ";" << endl;
            if(!_fermat_using_array) ss << "res:=*res+*item;" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
        }
        if(_fermat_using_array) {
            if(_fermat_using_array==1) ss << "res:=Sumup([m]);" << endl;
            else if(_fermat_using_array==2) ss << "res:=Sigma<i=1,"<<item.nops()<<">(*m[i]);" << endl;
            else {
                ss << "n:=" << item.nops() << ";" << endl;
                ss << "while n>1 do n2:=n\\2; for i=1,n2 do m[i]:=*m[i]+*m[n+1-i] od; &_G; if (n|2)=0 then n:=n2 else n:=n2+1 fi od;" << endl;
                ss << "res:=*m[1];" << endl;
            }
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
        }
        
        static string bstr("[-begin-]"), estr("[-end-]");
        ss << "&(U=1);" << endl; // ugly printing, the whitespace matters
        ss << "!('" <<bstr<< "',res,'" <<estr<< "')" << endl;
        auto ostr = fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        // note the order,(normal_fermat will be called again in factor_form)
        ss << "&(U=0);" << endl; // disable ugly printing
        if(_fermat_using_array) ss << "@(res,[m]);" << endl;
        else ss << "@(res,item);" << endl;
        ss << "&_G;" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");

        // make sure last char is 0
        if(ostr[ostr.length()-1]!='0') throw Error("fermat_together: last char is NOT 0.");
        ostr = ostr.substr(0, ostr.length()-1);
        auto cpos = ostr.find(bstr);
        if(cpos==string::npos) throw Error(bstr+" NOT Found.");
        ostr = ostr.substr(cpos+bstr.length(),string::npos);
        cpos = ostr.find(estr);
        if(cpos==string::npos) throw Error(estr+" NOT Found.");
        ostr = ostr.substr(0,cpos);
        string_trim(ostr);

        symtab st;
        Parser fp(st);
        auto ret = fp.Read(ostr);
        //fermat.Exit();
        
        if(use_ncheck) {
            auto nn_ret = subs(ret,nn_map);
            if(nn_chk-nn_ret!=0) {
                cout << nn_chk << " : " << nn_ret << endl;
                throw Error("fermat_together: N Check Failed.");
            }
        }
        
        ret = ret.subs(f2v,subs_options::no_pattern).subs(map_rat,subs_options::no_pattern);
        return ret;
        
    }

    /**
     * @brief return the numerator and denominator after normalization
     * @param expr the input expression
     * @return fermat evaluated expression
     */
    ex fermat_eval(const ex & expr) {
        Fermat &fermat = Fermat::get();
        int &v_max = fermat.vmax;
        
        auto expr_in = expr;
        exmap map_rat;
        expr_in = expr_in.to_rational(map_rat);
        
        lst rep_vs;
        if(true) {
            map<ex,long long,ex_is_less> s2c;
            for(const_preorder_iterator i = expr_in.preorder_begin(); i != expr_in.preorder_end(); ++i) {
                if(is_a<symbol>(*i)) s2c[*i]++;
            }
            exvector sv1, sv2, sv3; // sv2:fermat_weight, sv3:map_rat
            for(auto kv : s2c) {
                auto fw = fermat_weight.find(kv.first);
                if(fw!=fermat_weight.end()) sv2.push_back(lst{fw->second, fw->first});
                else if(map_rat.find(kv.first)!=map_rat.end()) sv3.push_back(lst{kv.second, kv.first});
                else sv1.push_back(lst{kv.second, kv.first});
            }
            sort_vec(sv1);
            sort_vec(sv2);
            sort_vec(sv3);
            for(auto sv : sv1) rep_vs.append(sv.op(1));
            for(auto sv : sv2) rep_vs.append(sv.op(1));
            for(auto sv : sv3) rep_vs.append(sv.op(1));
        }
                
        exmap v2f, f2v;
        int fvi = 0;
        for(auto vi : rep_vs) {
            auto name = "v" + to_string(fvi);
            Symbol s(name);
            v2f[vi] = s;
            f2v[s] = vi;
            fvi++;
        }
        expr_in = expr_in.subs(v2f);
        
        stringstream ss;
        if(fvi>111) {
            cout << rep_vs << endl;
            throw Error("Fermat: Too many variables.");
        }
        if(fvi>v_max) {
            for(int i=v_max; i<fvi; i++) ss << "&(J=v" << i << ");" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
            v_max = fvi;
        }
        
        ex res;
        ss << "res:=" << expr_in << ";" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        static string bstr("[-begin-]"), estr("[-end-]");
        ss << "&(U=1);" << endl; // ugly printing, the whitespace matters
        ss << "!('" <<bstr<< "',res,'" <<estr<< "')" << endl;
        auto ostr = fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        // note the order,(normal_fermat will be called again in factor_form)
        ss << "&(U=0);" << endl; // disable ugly printing
        ss << "@(res);" << endl;
        ss << "&_G;" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");

        // make sure last char is 0
        if(ostr[ostr.length()-1]!='0') throw Error("fermat_together: last char is NOT 0.");
        ostr = ostr.substr(0, ostr.length()-1);
        auto cpos = ostr.find(bstr);
        if(cpos==string::npos) throw Error(bstr+" NOT Found.");
        ostr = ostr.substr(cpos+bstr.length(),string::npos);
        cpos = ostr.find(estr);
        if(cpos==string::npos) throw Error(estr+" NOT Found.");
        ostr = ostr.substr(0,cpos);
        string_trim(ostr);

        symtab st;
        Parser fp(st);
        res = fp.Read(ostr);
        res = res.subs(f2v,nopat);
        return res.subs(map_rat,nopat);
    }
    
    matrix fermat_Redrowech(const matrix & mat_in) {
        Fermat &fermat = Fermat::get();
        int &v_max = fermat.vmax;
        
        exmap map_rat;
        int nrow = mat_in.rows();
        int ncol = mat_in.cols();
        matrix mat(nrow, ncol);
        for(int r=0; r<nrow; r++) for(int c=0; c<ncol; c++) mat(r,c) = mat_in(r,c).to_rational(map_rat);
        
        lst rep_vs;
        if(true) {
            map<ex,long long,ex_is_less> s2c;
            ex expr_in = mat;
            for(const_preorder_iterator i = expr_in.preorder_begin(); i != expr_in.preorder_end(); ++i) {
                if(is_a<symbol>(*i)) s2c[*i]++;
            }
            exvector sv1, sv2, sv3; // sv2:fermat_weight, sv3:map_rat
            for(auto kv : s2c) {
                auto fw = fermat_weight.find(kv.first);
                if(fw!=fermat_weight.end()) sv2.push_back(lst{fw->second, fw->first});
                else if(map_rat.find(kv.first)!=map_rat.end()) sv3.push_back(lst{kv.second, kv.first});
                else sv1.push_back(lst{kv.second, kv.first});
            }
            sort_vec(sv1);
            sort_vec(sv2);
            sort_vec(sv3);
            for(auto sv : sv1) rep_vs.append(sv.op(1));
            for(auto sv : sv2) rep_vs.append(sv.op(1));
            for(auto sv : sv3) rep_vs.append(sv.op(1));
        }
        
        exmap v2f;
        symtab st;
        int fvi = 0;
        for(auto vi : rep_vs) {
            auto name = "v" + to_string(fvi);
            v2f[vi] = Symbol(name);
            st[name] = vi;
            fvi++;
        }
        
        stringstream ss;
        if(fvi>111) {
            cout << rep_vs << endl;
            throw Error("Fermat: Too many variables.");
        }
        if(fvi>v_max) {
            for(int i=v_max; i<fvi; i++) ss << "&(J=v" << i << ");" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
            v_max = fvi;
        }
        
        ss << "Array m[" << nrow << "," << ncol << "];" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        ss << "[m]:=[(";
        for(int c=0; c<ncol; c++) {
            for(int r=0; r<nrow; r++) {
                ss << mat(r,c).subs(iEpsilon==0,nopat).subs(v2f,nopat) << ",";
            }
        }
        ss << ")];" << endl;
        ss << "Redrowech([m]);" << endl;
        auto tmp = ss.str();
        string_replace_all(tmp,",)]",")]");
        fermat.Execute(tmp);
        ss.clear();
        ss.str("");

        ss << "&(U=1);" << endl; // ugly printing, the whitespace matters
        ss << "![m" << endl;
        auto ostr = fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        //fermat.Exit();
        
        // note the order, before exfactor (normal_fermat will be called again here)
        ss << "&(U=0);" << endl; // disable ugly printing
        ss << "@([m]);" << endl;
        ss << "&_G;" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");

        // make sure last char is 0
        if(ostr[ostr.length()-1]!='0') throw Error("Direc::Export, last char is NOT 0.");
        ostr = ostr.substr(0, ostr.length()-1);
        string_trim(ostr);
        
        ostr.erase(0, ostr.find(":=")+2);
        string_replace_all(ostr, "[", "{");
        string_replace_all(ostr, "]", "}");
        Parser fp(st);
        matrix mr(nrow, ncol);
        auto res = fp.Read(ostr);
        for(int r=0; r<nrow; r++) {
            auto cur = res.op(r);
            for(int c=0; c<ncol; c++) mr(r,c) = cur.op(c).subs(map_rat);
        }
        return mr;
    }
    
    ex fermat_Det(const matrix & mat_in) {
        Fermat &fermat = Fermat::get();
        int &v_max = fermat.vmax;
        
        exmap map_rat;
        int nrow = mat_in.rows();
        int ncol = mat_in.cols();
        matrix mat(nrow, ncol);
        for(int r=0; r<nrow; r++) for(int c=0; c<ncol; c++) mat(r,c) = mat_in(r,c).to_rational(map_rat);
        
        lst rep_vs;
        if(true) {
            map<ex,long long,ex_is_less> s2c;
            ex expr_in = mat;
            for(const_preorder_iterator i = expr_in.preorder_begin(); i != expr_in.preorder_end(); ++i) {
                if(is_a<symbol>(*i)) s2c[*i]++;
            }
            exvector sv1, sv2, sv3; // sv2:fermat_weight, sv3:map_rat
            for(auto kv : s2c) {
                auto fw = fermat_weight.find(kv.first);
                if(fw!=fermat_weight.end()) sv2.push_back(lst{fw->second, fw->first});
                else if(map_rat.find(kv.first)!=map_rat.end()) sv3.push_back(lst{kv.second, kv.first});
                else sv1.push_back(lst{kv.second, kv.first});
            }
            sort_vec(sv1);
            sort_vec(sv2);
            sort_vec(sv3);
            for(auto sv : sv1) rep_vs.append(sv.op(1));
            for(auto sv : sv2) rep_vs.append(sv.op(1));
            for(auto sv : sv3) rep_vs.append(sv.op(1));
        }
        
        exmap v2f;
        symtab st;
        int fvi = 0;
        for(auto vi : rep_vs) {
            auto name = "v" + to_string(fvi);
            v2f[vi] = Symbol(name);
            st[name] = vi;
            fvi++;
        }
        
        stringstream ss;
        if(fvi>111) {
            cout << rep_vs << endl;
            throw Error("Fermat: Too many variables.");
        }
        if(fvi>v_max) {
            for(int i=v_max; i<fvi; i++) ss << "&(J=v" << i << ");" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
            v_max = fvi;
        }
        
        ss << "Array m[" << nrow << "," << ncol << "];" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        ss << "[m]:=[(";
        for(int c=0; c<ncol; c++) {
            for(int r=0; r<nrow; r++) {
                ss << mat(r,c).subs(iEpsilon==0,nopat).subs(v2f,nopat) << ",";
            }
        }
        ss << ")];" << endl;
        ss << "res:=Det[m];" << endl;
        auto tmp = ss.str();
        string_replace_all(tmp,",)]",")]");
        fermat.Execute(tmp);
        ss.clear();
        ss.str("");

        static string bstr("[-begin-]"), estr("[-end-]");
        ss << "&(U=1);" << endl; // ugly printing, the whitespace matters
        ss << "!('" <<bstr<< "',res,'" <<estr<< "')" << endl;
        auto ostr = fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        // note the order,(normal_fermat will be called again in factor_form)
        ss << "&(U=0);" << endl; // disable ugly printing
        ss << "@(res,[m]);" << endl;
        ss << "&_G;" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");

        // make sure last char is 0
        if(ostr[ostr.length()-1]!='0') throw Error("fermat_together: last char is NOT 0.");
        ostr = ostr.substr(0, ostr.length()-1);
        auto cpos = ostr.find(bstr);
        if(cpos==string::npos) throw Error(bstr+" NOT Found.");
        ostr = ostr.substr(cpos+bstr.length(),string::npos);
        cpos = ostr.find(estr);
        if(cpos==string::npos) throw Error(estr+" NOT Found.");
        ostr = ostr.substr(0,cpos);
        string_trim(ostr);
        
        Parser fp(st);
        auto res = fp.Read(ostr);
        res = res.subs(map_rat);
        return res;
    }
        
    matrix fermat_inv(const matrix & mat_in) {
        Fermat &fermat = Fermat::get();
        int &v_max = fermat.vmax;
        
        exmap map_rat;
        int nrow = mat_in.rows();
        int ncol = mat_in.cols();
        matrix mat(nrow, ncol);
        for(int r=0; r<nrow; r++) for(int c=0; c<ncol; c++) mat(r,c) = mat_in(r,c).to_rational(map_rat);
        
        lst rep_vs;
        if(true) {
            map<ex,long long,ex_is_less> s2c;
            ex expr_in = mat;
            for(const_preorder_iterator i = expr_in.preorder_begin(); i != expr_in.preorder_end(); ++i) {
                if(is_a<symbol>(*i)) s2c[*i]++;
            }
            exvector sv1, sv2, sv3; // sv2:fermat_weight, sv3:map_rat
            for(auto kv : s2c) {
                auto fw = fermat_weight.find(kv.first);
                if(fw!=fermat_weight.end()) sv2.push_back(lst{fw->second, fw->first});
                else if(map_rat.find(kv.first)!=map_rat.end()) sv3.push_back(lst{kv.second, kv.first});
                else sv1.push_back(lst{kv.second, kv.first});
            }
            sort_vec(sv1);
            sort_vec(sv2);
            sort_vec(sv3);
            for(auto sv : sv1) rep_vs.append(sv.op(1));
            for(auto sv : sv2) rep_vs.append(sv.op(1));
            for(auto sv : sv3) rep_vs.append(sv.op(1));
        }
        
        exmap v2f;
        symtab st;
        int fvi = 0;
        for(auto vi : rep_vs) {
            auto name = "v" + to_string(fvi);
            v2f[vi] = Symbol(name);
            st[name] = vi;
            fvi++;
        }
        
        stringstream ss;
        if(fvi>111) {
            cout << rep_vs << endl;
            throw Error("Fermat: Too many variables.");
        }
        if(fvi>v_max) {
            for(int i=v_max; i<fvi; i++) ss << "&(J=v" << i << ");" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
            v_max = fvi;
        }
        
        ss << "Array m[" << nrow << "," << ncol << "];" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        ss << "[m]:=[(";
        for(int c=0; c<ncol; c++) {
            for(int r=0; r<nrow; r++) {
                ss << mat(r,c).subs(iEpsilon==0,nopat).subs(v2f,nopat) << ",";
            }
        }
        ss << ")];" << endl;
        ss << "[m]:=1/[m];" << endl;
        auto tmp = ss.str();
        string_replace_all(tmp,",)]",")]");
        fermat.Execute(tmp);
        ss.clear();
        ss.str("");

        ss << "&(U=1);" << endl; // ugly printing, the whitespace matters
        ss << "![m" << endl;
        auto ostr = fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        //fermat.Exit();
        
        // note the order, before exfactor (normal_fermat will be called again here)
        ss << "&(U=0);" << endl; // disable ugly printing
        ss << "@([m]);" << endl;
        ss << "&_G;" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");

        // make sure last char is 0
        if(ostr[ostr.length()-1]!='0') throw Error("Direc::Export, last char is NOT 0.");
        ostr = ostr.substr(0, ostr.length()-1);
        string_trim(ostr);
        
        ostr.erase(0, ostr.find(":=")+2);
        string_replace_all(ostr, "[", "{");
        string_replace_all(ostr, "]", "}");
        Parser fp(st);
        matrix mr(nrow, ncol);
        auto res = fp.Read(ostr);
        for(int r=0; r<nrow; r++) {
            auto cur = res.op(r);
            for(int c=0; c<ncol; c++) mr(r,c) = cur.op(c).subs(map_rat);
        }
        return mr;
    }
    
    matrix fermat_mul(const matrix & m1, const matrix & m2) {
        Fermat &fermat = Fermat::get();
        int &v_max = fermat.vmax;
        
        exmap map_rat;
        int nr1 = m1.rows();
        int nc1 = m1.cols();
        matrix mat1(nr1, nc1);
        for(int r=0; r<nr1; r++) for(int c=0; c<nc1; c++) mat1(r,c) = m1(r,c).to_rational(map_rat);
        int nr2 = m2.rows();
        int nc2 = m2.cols();
        matrix mat2(nr2, nc2);
        for(int r=0; r<nr2; r++) for(int c=0; c<nc2; c++) mat2(r,c) = m2(r,c).to_rational(map_rat);
        
        lst rep_vs;
        if(true) {
            map<ex,long long,ex_is_less> s2c;
            ex expr_in = mat1;
            for(const_preorder_iterator i = expr_in.preorder_begin(); i != expr_in.preorder_end(); ++i) {
                if(is_a<symbol>(*i)) s2c[*i]++;
            }
            expr_in = mat2;
            for(const_preorder_iterator i = expr_in.preorder_begin(); i != expr_in.preorder_end(); ++i) {
                if(is_a<symbol>(*i)) s2c[*i]++;
            }
            exvector sv1, sv2, sv3; // sv2:fermat_weight, sv3:map_rat
            for(auto kv : s2c) {
                auto fw = fermat_weight.find(kv.first);
                if(fw!=fermat_weight.end()) sv2.push_back(lst{fw->second, fw->first});
                else if(map_rat.find(kv.first)!=map_rat.end()) sv3.push_back(lst{kv.second, kv.first});
                else sv1.push_back(lst{kv.second, kv.first});
            }
            sort_vec(sv1);
            sort_vec(sv2);
            sort_vec(sv3);
            for(auto sv : sv1) rep_vs.append(sv.op(1));
            for(auto sv : sv2) rep_vs.append(sv.op(1));
            for(auto sv : sv3) rep_vs.append(sv.op(1));
        }
        
        exmap v2f;
        symtab st;
        int fvi = 0;
        for(auto vi : rep_vs) {
            auto name = "v" + to_string(fvi);
            v2f[vi] = Symbol(name);
            st[name] = vi;
            fvi++;
        }
        
        stringstream ss;
        if(fvi>111) {
            cout << rep_vs << endl;
            throw Error("Fermat: Too many variables.");
        }
        if(fvi>v_max) {
            for(int i=v_max; i<fvi; i++) ss << "&(J=v" << i << ");" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
            v_max = fvi;
        }
        
        ss << "Array m1[" << nr1 << "," << nc1 << "];" << endl;
        ss << "Array m2[" << nr2 << "," << nc2 << "];" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        ss << "[m1]:=[(";
        for(int c=0; c<nc1; c++) {
            for(int r=0; r<nr1; r++) {
                ss << mat1(r,c).subs(iEpsilon==0,nopat).subs(v2f,nopat) << ",";
            }
        }
        ss << ")];" << endl;
        
        ss << "[m2]:=[(";
        for(int c=0; c<nc2; c++) {
            for(int r=0; r<nr2; r++) {
                ss << mat2(r,c).subs(iEpsilon==0,nopat).subs(v2f,nopat) << ",";
            }
        }
        ss << ")];" << endl;
        
        ss << "[m]:=[m1]*[m2];" << endl;
        auto tmp = ss.str();
        string_replace_all(tmp,",)]",")]");
        fermat.Execute(tmp);
        ss.clear();
        ss.str("");

        ss << "&(U=1);" << endl; // ugly printing, the whitespace matters
        ss << "![m" << endl;
        auto ostr = fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        //fermat.Exit();
        
        // note the order, before exfactor (normal_fermat will be called again here)
        ss << "&(U=0);" << endl; // disable ugly printing
        ss << "@([m1],[m2],[m]);" << endl;
        ss << "&_G;" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");

        // make sure last char is 0
        if(ostr[ostr.length()-1]!='0') throw Error("Direc::Export, last char is NOT 0.");
        ostr = ostr.substr(0, ostr.length()-1);
        string_trim(ostr);
        
        ostr.erase(0, ostr.find(":=")+2);
        string_replace_all(ostr, "[", "{");
        string_replace_all(ostr, "]", "}");
        Parser fp(st);
        matrix mr(nr1, nc2);
        auto res = fp.Read(ostr);
        for(int r=0; r<nr1; r++) {
            auto cur = res.op(r);
            for(int c=0; c<nc2; c++) mr(r,c) = cur.op(c).subs(map_rat);
        }
        return mr;
    }
    
    matrix fermat_pow(const matrix & mat_in, int n) {
        Fermat &fermat = Fermat::get();
        int &v_max = fermat.vmax;
        
        exmap map_rat;
        int nrow = mat_in.rows();
        int ncol = mat_in.cols();
        matrix mat(nrow, ncol);
        for(int r=0; r<nrow; r++) for(int c=0; c<ncol; c++) mat(r,c) = mat_in(r,c).to_rational(map_rat);
        
        lst rep_vs;
        if(true) {
            map<ex,long long,ex_is_less> s2c;
            ex expr_in = mat;
            for(const_preorder_iterator i = expr_in.preorder_begin(); i != expr_in.preorder_end(); ++i) {
                if(is_a<symbol>(*i)) s2c[*i]++;
            }
            exvector sv1, sv2, sv3; // sv2:fermat_weight, sv3:map_rat
            for(auto kv : s2c) {
                auto fw = fermat_weight.find(kv.first);
                if(fw!=fermat_weight.end()) sv2.push_back(lst{fw->second, fw->first});
                else if(map_rat.find(kv.first)!=map_rat.end()) sv3.push_back(lst{kv.second, kv.first});
                else sv1.push_back(lst{kv.second, kv.first});
            }
            sort_vec(sv1);
            sort_vec(sv2);
            sort_vec(sv3);
            for(auto sv : sv1) rep_vs.append(sv.op(1));
            for(auto sv : sv2) rep_vs.append(sv.op(1));
            for(auto sv : sv3) rep_vs.append(sv.op(1));
        }
        
        exmap v2f;
        symtab st;
        int fvi = 0;
        for(auto vi : rep_vs) {
            auto name = "v" + to_string(fvi);
            v2f[vi] = Symbol(name);
            st[name] = vi;
            fvi++;
        }
        
        stringstream ss;
        if(fvi>111) {
            cout << rep_vs << endl;
            throw Error("Fermat: Too many variables.");
        }
        if(fvi>v_max) {
            for(int i=v_max; i<fvi; i++) ss << "&(J=v" << i << ");" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
            v_max = fvi;
        }
        
        ss << "Array m[" << nrow << "," << ncol << "];" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        ss << "[m]:=[(";
        for(int c=0; c<ncol; c++) {
            for(int r=0; r<nrow; r++) {
                ss << mat(r,c).subs(iEpsilon==0,nopat).subs(v2f,nopat) << ",";
            }
        }
        ss << ")];" << endl;
        ss << "[m]:=[m]^(" << n << ");" << endl;
        auto tmp = ss.str();
        string_replace_all(tmp,",)]",")]");
        fermat.Execute(tmp);
        ss.clear();
        ss.str("");

        ss << "&(U=1);" << endl; // ugly printing, the whitespace matters
        ss << "![m" << endl;
        auto ostr = fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        //fermat.Exit();
        
        // note the order, before exfactor (normal_fermat will be called again here)
        ss << "&(U=0);" << endl; // disable ugly printing
        ss << "@([m]);" << endl;
        ss << "&_G;" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");

        // make sure last char is 0
        if(ostr[ostr.length()-1]!='0') throw Error("Direc::Export, last char is NOT 0.");
        ostr = ostr.substr(0, ostr.length()-1);
        string_trim(ostr);
        
        ostr.erase(0, ostr.find(":=")+2);
        string_replace_all(ostr, "[", "{");
        string_replace_all(ostr, "]", "}");
        Parser fp(st);
        matrix mr(nrow, ncol);
        auto res = fp.Read(ostr);
        for(int r=0; r<nrow; r++) {
            auto cur = res.op(r);
            for(int c=0; c<ncol; c++) mr(r,c) = cur.op(c).subs(map_rat);
        }
        return mr;
    }
    
    namespace {
        exmap mat_map_rat;
        symtab mat_st;
    }
    
    void fermat_mat(const matrix & mat_in, const string & _name) {
        string name = _name;
        string_replace_all(name, "[", "");
        string_replace_all(name, "]", "");
        
        Fermat &fermat = Fermat::get();
        int &v_max = fermat.vmax;
        
        int nrow = mat_in.rows();
        int ncol = mat_in.cols();
        matrix mat(nrow, ncol);
        for(int r=0; r<nrow; r++) for(int c=0; c<ncol; c++) mat(r,c) = mat_in(r,c).to_rational(mat_map_rat);
        
        lst rep_vs;
        if(true) {
            map<ex,long long,ex_is_less> s2c;
            ex expr_in = mat;
            for(const_preorder_iterator i = expr_in.preorder_begin(); i != expr_in.preorder_end(); ++i) {
                if(is_a<symbol>(*i)) s2c[*i]++;
            }
            exvector sv1, sv2, sv3; // sv2:fermat_weight, sv3:map_rat
            for(auto kv : s2c) {
                auto fw = fermat_weight.find(kv.first);
                if(fw!=fermat_weight.end()) sv2.push_back(lst{fw->second, fw->first});
                else if(mat_map_rat.find(kv.first)!=mat_map_rat.end()) sv3.push_back(lst{kv.second, kv.first});
                else sv1.push_back(lst{kv.second, kv.first});
            }
            sort_vec(sv1);
            sort_vec(sv2);
            sort_vec(sv3);
            for(auto sv : sv1) rep_vs.append(sv.op(1));
            for(auto sv : sv2) rep_vs.append(sv.op(1));
            for(auto sv : sv3) rep_vs.append(sv.op(1));
        }
        
        exmap v2f;
        int fvi = 0;
        for(auto vi : rep_vs) {
            auto name = "v" + to_string(fvi);
            v2f[vi] = Symbol(name);
            mat_st[name] = vi;
            fvi++;
        }
        
        stringstream ss;
        if(fvi>111) {
            cout << rep_vs << endl;
            throw Error("Fermat: Too many variables.");
        }
        if(fvi>v_max) {
            for(int i=v_max; i<fvi; i++) ss << "&(J=v" << i << ");" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
            v_max = fvi;
        }
        
        ss << "Array " << name << "[" << nrow << "," << ncol << "];" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        ss << "["<<name<<"]:=[(";
        for(int c=0; c<ncol; c++) {
            for(int r=0; r<nrow; r++) {
                ss << mat(r,c).subs(iEpsilon==0,nopat).subs(v2f,nopat) << ",";
            }
        }
        ss << ")];" << endl;
        //ss << "Redrowech([m]);" << endl;
        auto tmp = ss.str();
        string_replace_all(tmp,",)]",")]");
        fermat.Execute(tmp);
        ss.clear();
        ss.str("");
    }
    
    matrix fermat_mat(const string & _name) {
        string name = _name;
        string_replace_all(name, "[", "");
        string_replace_all(name, "]", "");
        
        Fermat &fermat = Fermat::get();
        int &v_max = fermat.vmax;
        
        stringstream ss;
        ss << "&(U=1);" << endl; // ugly printing, the whitespace matters
        ss << "![" << name << endl;
        auto ostr = fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        //fermat.Exit();
        
        // note the order, before exfactor (normal_fermat will be called again here)
        ss << "&(U=0);" << endl; // disable ugly printing
        ss << "&_G;" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");

        // make sure last char is 0
        if(ostr[ostr.length()-1]!='0') throw Error("Direc::Export, last char is NOT 0.");
        ostr = ostr.substr(0, ostr.length()-1);
        string_trim(ostr);
        ostr.erase(0, ostr.find(":=")+2);
        string_replace_all(ostr, "[", "{");
        string_replace_all(ostr, "]", "}");
        Parser fp(mat_st);
        auto res = fp.Read(ostr);
        int nrow = res.nops();
        int ncol = res.op(0).nops();
        matrix mr(nrow, ncol);
        for(int r=0; r<nrow; r++) {
            auto cur = res.op(r);
            for(int c=0; c<ncol; c++) mr(r,c) = cur.op(c).subs(mat_map_rat);
        }
        return mr;
    }
    
    void fermat_eval(const string & fcmd) {
        Fermat &fermat = Fermat::get();
        int &v_max = fermat.vmax;
        
        stringstream ss;
        ss << fcmd << endl; // ugly printing, the whitespace matters
        auto ostr = fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        //fermat.Exit();
    }
    
}
