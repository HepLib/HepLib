/**
 * @file
 * @brief IBP with KIRA
 */
 
#include "IBP.h"
#include <cmath>

namespace HepLib {

    /**
     * @brief export integral to KIRA form
     */
    string KIRA::Fout(const ex & expr) {
        ex f = expr;
        if(is_a<lst>(f)) f = F(f);
        else if(expr.match(F(w1,w2))) f = F(f.op(1));
        string fstr = ex2str(f);
        string_replace_all(fstr,"{","");
        string_replace_all(fstr,"}","");
        string_replace_all(fstr,"(","[");
        string_replace_all(fstr,")","]");
        return fstr;
    }
    
    /**
     * @brief import integral from KIRA
     */
    ex KIRA::Fin(const string & expr) {
        string fstr = expr;
        string_replace_all(fstr,"[","("+to_string(ProblemNumber)+",{");
        string_replace_all(fstr,"]","})");
        return str2ex(fstr);
    }

    /**
     * @brief Export input data for KIRA
     */
    void KIRA::Export() {

        if(Integral.nops()<1) return;
        
        if(WorkingDir.length()<1) WorkingDir = to_string(getpid())+"IBP";
        string job_dir = WorkingDir + "/" + to_string(ProblemNumber);
        system(("rm -rf "+job_dir).c_str());
        if(!dir_exists(job_dir)) system(("mkdir -p "+job_dir).c_str());
        
        string config_dir = job_dir + "/config";
        system(("rm -rf "+config_dir).c_str());
        if(!dir_exists(config_dir)) system(("mkdir -p "+config_dir).c_str());
        
        ostringstream oss;
        
        // integralfamilies.yaml
        oss << "integralfamilies:" << endl;
        oss << "  - name: F" << endl;
        oss << "    loop_momenta: [";
        for(int i=0; i<Internal.nops(); i++) {
            oss << Internal.op(i);
            if(i+1<Internal.nops()) oss << ", ";
            else oss << "]" << endl;
        }
        
        vector<int> symbolic_ibp;
        for(auto kv : Shift) {
            if(kv.second.is_zero()) continue;
            symbolic_ibp.push_back(kv.first);
        }
        if(symbolic_ibp.size()>0) {
            oss << "    symbolic_ibp: [";
            for(int i=0; i<symbolic_ibp.size(); i++) {
                if(i>0) oss << ",";
                oss << symbolic_ibp[i];
            }
            oss << "]" << endl;
        }
        
        //oss << "    top_level_sectors: [127]" << endl;
        oss << "    propagators:" << endl;
        for(int i=0; i<Propagator.nops(); i++) {
            oss << "     - [ \"" << Propagator.op(i) << "\", 0]" << endl;
        }
        if(Cut.nops()>0) {
            oss << "    cut_propagators: [";
            for(auto i=0; i<Cut.nops(); i++) {
                oss << Cut.op(i);
                if(i+1<Cut.nops()) oss << ", ";
                else oss << "]" << endl;
            }
        }
        ofstream if_out(config_dir+"/integralfamilies.yaml");
        if_out << oss.str() << endl;
        if_out.close();
        oss.str("");
        oss.clear();
        
        
        // kinematics.yaml
        oss << "kinematics:" << endl;
        oss << "  incoming_momenta: [";
        for(int i=0; i<External.nops(); i++) {
            oss << External.op(i);
            if(i+1<External.nops()) oss << ", ";
            else oss << "]" << endl;
        }
        //oss << "  outgoing_momenta: []" << endl;
        //oss << "  momentum_conservation: [p0,-p1-p2]" << endl;
        
        auto vars = gather_symbols(lst{Propagator, Replacement});
        exset vset_all;
        for(auto vi : vars) vset_all.insert(vi);
        exset vset_mom;
        for(auto vi : Internal) vset_mom.insert(vi);
        for(auto vi : External) vset_mom.insert(vi);
        exset vset;
        for(auto vi : vset_all) {
            if(vset_mom.find(vi) == vset_mom.end()) vset.insert(vi);
        }
        if(vset.size()>0) {
            oss << "  kinematic_invariants:" << endl;
            for(auto vi : vset) oss << "   - [ " << vi << ", 0 ]" << endl;
            for(auto i : symbolic_ibp) oss << "   - [ b" << (i-1) << ", 0 ]" << endl;
        }
        
        oss << "  scalarproduct_rules:" << endl;
        for(int i=0; i<External.nops(); i++) {
            auto pi = External.op(i);
            for(int j=i; j<External.nops(); j++) {
                auto pj = External.op(j);
                oss << "    - [ [" << pi << "," << pj << "], ";
                oss << subs(pi*pj, Replacement) << "]" << endl;
            }
        }
        
        //oss << "  symbol_to_replace_by_one: s" << endl;

        ofstream km_out(config_dir+"/kinematics.yaml");
        km_out << oss.str() << endl;
        km_out.close();
        oss.str("");
        oss.clear();
        
        // handle rmax, smax, dmax
        int _rmax = -1;
        int _smax = -1;
        if(_rmax < 0 || _smax < 0) {
            int rrmax = 0;
            int ssmax = 0;
            for(auto integral : Integral) {
                int rr = 0;
                int ss = 0;
                for(auto item : integral) {
                    if(item>0) rr += ex2int(item);
                    else ss -= ex2int(item);
                }
                if(rrmax<rr) rrmax = rr;
                if(ssmax<ss) ssmax = ss;
            }
            if(rrmax > _rmax) _rmax = rrmax;
            if(ssmax > _smax) _smax = ssmax;
        }
        _rmax += ra;
        _smax += sa;
        
        // job.yaml
        oss << "jobs:" << endl;
        oss << "  - reduce_sectors:" << endl;
        oss << "      reduce: " << endl;
        oss << "        - {r: " << _rmax << ", s: " << _smax << "}" << endl;
        oss << "      select_integrals: " << endl;
        oss << "        select_mandatory_recursively: " << endl;
        oss << "          - {r: " << _rmax << ", s: " << _smax;
        if(dmax>0)  oss << ", d: " << dmax;
        oss << "}" << endl;
        oss << "      run_initiate: true" << endl;
        oss << "      run_triangular: true" << endl;
        oss << "      run_back_substitution: true" << endl;
        if(PIntegral.nops()>0) {
            ostringstream oss2;
            int nn = PIntegral.nops();
            for(int i=0; i<nn; i++) oss2 << Fout(PIntegral.op(i)) << endl;
            ofstream pref_out(job_dir+"/preferred");
            pref_out << oss2.str() << endl;
            pref_out.close();
            oss << "      preferred_masters: preferred" << endl;
        }
        oss << "  - kira2file:" << endl;
        oss << "      target:" << endl;
        oss << "        - [F,integrals]" << endl;
        ofstream job_out(job_dir+"/job.yaml");
        job_out << oss.str() << endl;
        job_out.close();
        oss.str("");
        oss.clear();
        
        // integrals
        for(auto integral : Integral) oss << Fout(integral) << endl;
        ofstream intg_out(job_dir+"/integrals");
        intg_out << oss.str() << endl;
        intg_out.close();
        oss.str("");
        oss.clear();
    }
    
    /**
     * @brief invoke kira program for reduction
     */
    void KIRA::Run() {
        string job_dir = WorkingDir + "/" + to_string(ProblemNumber);
        ostringstream cmd;
        cmd << "cd " << job_dir << " && kira " << KArgs << " --silent job.yaml >/dev/null 2>&1";
        system(cmd.str().c_str());
    }

    /**
     * @brief import kira result
     */
    void KIRA::Import() {
        string job_dir = WorkingDir + "/" + to_string(ProblemNumber);
        ostringstream fn;
        fn << job_dir << "/results/F/kira_integrals.kira";
        auto strvec = file2strvec(fn.str());
        
        exmap bMAP;
        for(auto kv : Shift) {
            Symbol bi("b"+to_string(kv.first-1));
            bMAP[bi] = kv.second;
        }
        
        ex exL=0, exR=0;
        map<ex,int,ex_is_less> flags;
        lst exRs;
        for(auto intg : Integral) flags[F(ProblemNumber,intg)] = 1;
        Rules.remove_all();
        for(auto line : strvec) {
            if(line.size()==0) {
                if(!is_zero(exL)) {
                    exR = exR.subs(bMAP);
                    Rules.append(exL==exR);
                    flags[exL] = 0;
                    exRs.append(exR);
                }
                exL = exR = 0;
            } else if(is_zero(exL)) {
                exL -= Fin(line);
                if(!exL.match(F(w1,w2))) {
                    cout << line << endl;
                    throw Error("KIRA::Import error found.");
                }
            } else {
                exR += Fin(line);
            }
        }
        if(!is_zero(exL)) {
            exR = exR.subs(bMAP);
            Rules.append(exL==exR);
            flags[exL] = 0;
            exRs.append(exR);
        }
        MIntegral.remove_all();
        for(auto kv : flags) {
            if(kv.second!=0) MIntegral.append(kv.first);
        }
        exset miset;
        find(exRs,F(w1,w2),miset);
        for(auto mi : miset) MIntegral.append(mi);
        MIntegral.sort();
        MIntegral.unique();

    }

}
