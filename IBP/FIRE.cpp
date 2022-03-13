/**
 * @file
 * @brief IBP with FIRE
 */
 
#include "IBP.h"
#include <cmath>

namespace HepLib {

    // FROM FIRE6.m
    namespace {
        unsigned long long HardCodedPrimes[128] = {
            2017llu, 18446744073709551557llu, 18446744073709551533llu,
            18446744073709551521llu, 18446744073709551437llu, 18446744073709551427llu,
            18446744073709551359llu, 18446744073709551337llu, 18446744073709551293llu,
            18446744073709551263llu, 18446744073709551253llu, 18446744073709551191llu,
            18446744073709551163llu, 18446744073709551113llu, 18446744073709550873llu,
            18446744073709550791llu, 18446744073709550773llu, 18446744073709550771llu,
            18446744073709550719llu, 18446744073709550717llu, 18446744073709550681llu,
            18446744073709550671llu, 18446744073709550593llu, 18446744073709550591llu,
            18446744073709550539llu, 18446744073709550537llu, 18446744073709550381llu,
            18446744073709550341llu, 18446744073709550293llu, 18446744073709550237llu,
            18446744073709550147llu, 18446744073709550141llu, 18446744073709550129llu,
            18446744073709550111llu, 18446744073709550099llu, 18446744073709550047llu,
            18446744073709550033llu, 18446744073709550009llu, 18446744073709549951llu,
            18446744073709549861llu, 18446744073709549817llu, 18446744073709549811llu,
            18446744073709549777llu, 18446744073709549757llu, 18446744073709549733llu,
            18446744073709549667llu, 18446744073709549621llu, 18446744073709549613llu,
            18446744073709549583llu, 18446744073709549571llu, 18446744073709549519llu,
            18446744073709549483llu, 18446744073709549441llu, 18446744073709549363llu,
            18446744073709549331llu, 18446744073709549327llu, 18446744073709549307llu,
            18446744073709549237llu, 18446744073709549153llu, 18446744073709549123llu,
            18446744073709549067llu, 18446744073709549061llu, 18446744073709549019llu,
            18446744073709548983llu, 18446744073709548899llu, 18446744073709548887llu,
            18446744073709548859llu, 18446744073709548847llu, 18446744073709548809llu,
            18446744073709548703llu, 18446744073709548599llu, 18446744073709548587llu,
            18446744073709548557llu, 18446744073709548511llu, 18446744073709548503llu,
            18446744073709548497llu, 18446744073709548481llu, 18446744073709548397llu,
            18446744073709548391llu, 18446744073709548379llu, 18446744073709548353llu,
            18446744073709548349llu, 18446744073709548287llu, 18446744073709548271llu,
            18446744073709548239llu, 18446744073709548193llu, 18446744073709548119llu,
            18446744073709548073llu, 18446744073709548053llu, 18446744073709547821llu,
            18446744073709547797llu, 18446744073709547777llu, 18446744073709547731llu,
            18446744073709547707llu, 18446744073709547669llu, 18446744073709547657llu,
            18446744073709547537llu, 18446744073709547521llu, 18446744073709547489llu,
            18446744073709547473llu, 18446744073709547471llu, 18446744073709547371llu,
            18446744073709547357llu, 18446744073709547317llu, 18446744073709547303llu,
            18446744073709547117llu, 18446744073709547087llu, 18446744073709547003llu,
            18446744073709546897llu, 18446744073709546879llu, 18446744073709546873llu,
            18446744073709546841llu, 18446744073709546739llu, 18446744073709546729llu,
            18446744073709546657llu, 18446744073709546643llu, 18446744073709546601llu,
            18446744073709546561llu, 18446744073709546541llu, 18446744073709546493llu,
            18446744073709546429llu, 18446744073709546409llu, 18446744073709546391llu,
            18446744073709546363llu, 18446744073709546337llu, 18446744073709546333llu,
            18446744073709546289llu, 18446744073709546271llu};
    }
        
    void FIRE::RRTables(const string & filename, int pnum) {
        if(pnum<1) return;
        exvector _tables(pnum);
        unsigned int nmin = 0;
        for(int i=0; i<pnum; i++) {
            string fn = filename;
            string_replace_all(fn, ".tables", "-"+to_string(i+1)+".tables");
            ifstream ifs(fn);
            string ostr((istreambuf_iterator<char>(ifs)), (istreambuf_iterator<char>()));
            ifs.close();
            string_replace_all(ostr, "\"", "");
            Parser tp;
            _tables[i] = tp.Read(ostr);
            auto nn = _tables[i].op(1).nops();
            if(i==0 || nmin>nn) nmin = nn;
        }
        exvector tables;
        vector<numeric> primes;
        ex convs = 0;
        for(int i=0; i<_tables.size(); i++) {
            auto nn = _tables[i].op(1).nops();
            if(nmin<nn) {
                cout << "RRTables: current table will be dropped." << endl;
                continue;
            }
            if(!is_a<lst>(convs)) convs = _tables[i].op(1); // second part is equal
            else if(!convs.is_equal(_tables[i].op(1))) {
                throw Error("RRT: convs are different.");
            }
            tables.push_back(_tables[i].op(0));
            primes.push_back(HardCodedPrimes[i+1]);
        }
        
        vector<pair<ex,vector<pair<ex,exvector>>>> int_mi_cs_vec;
        int ci = 0;
        for(int i=0; i<tables.size(); i++) {
            if(i==0) {
                for(auto ti : tables[0]) {
                    pair<ex,vector<pair<ex,exvector>>> int_mi_cs;
                    int_mi_cs.first = ti.op(0);
                    for(auto tii : ti.op(1)) {
                        pair<ex,exvector> mi_cs;
                        mi_cs.first = tii.op(0);
                        mi_cs.second.push_back(tii.op(1));
                        int_mi_cs.second.push_back(mi_cs);
                    }
                    int_mi_cs_vec.push_back(int_mi_cs);
                }
            } else {
                for(int i2=0; i2<int_mi_cs_vec.size(); i2++) {
                    auto ti = tables[i].op(i2);
                    auto & int_mi_cs = int_mi_cs_vec[i2];
                    if(int_mi_cs.first != ti.op(0)) throw Error("E1");
                    for(int j=0; j<int_mi_cs.second.size(); j++) {
                        auto & mi_cs = int_mi_cs.second[j];
                        if(mi_cs.first != ti.op(1).op(j).op(0)) throw Error("E2");
                        mi_cs.second.push_back(ti.op(1).op(j).op(1));
                    }
                }
            }
        }
        
        std::ofstream ofs;
        ofs.open(filename, ios::out);
        ofs << "{" << endl;
        ofs << "    {" << endl;
        for(int i=0; i<int_mi_cs_vec.size(); i++) {
            auto imc = int_mi_cs_vec[i];
            ofs << "        {" << imc.first << "," << endl;
            ofs << "            {" << endl;
            for(int j=0; j<imc.second.size(); j++) {
                auto item = imc.second[j];
                ofs << "                {" << item.first << ", ";
                vector<numeric> nvec;
                for(auto e : item.second) nvec.push_back(ex_to<numeric>(e));
                ofs << "\"" << RationalReconstruct(nvec,primes) << "\"}";
                if(j+1<imc.second.size()) ofs << ",";
                ofs << endl;
            }
            ofs << "            }" << endl;
            ofs << "        }";
            if(i+1<int_mi_cs_vec.size()) ofs << ",";
            ofs << endl;
        }
        ofs << "    }," << endl << "    {" << endl;
        for(int i=0; i<convs.nops(); i++) {
            auto item = convs.op(i);
            ofs << "        {" << item.op(0) << ", " << item.op(1) << "}";
            if(i+1<convs.nops()) ofs << ",";
            ofs << endl;
        }
        ofs << "    }" << endl;
        ofs << "}" << endl;
        ofs.close();
    }
    
    void FIRE::ThieleTables(const string & filename, int si, int ei) {
        int pnum = 40;
        map<int,ex> _tables;
        int nd;
        unsigned int nmin = 0;
        for(int i=si; i<ei; i++) {
            string fn = filename;
            string_replace_all(fn, ".tables", "-"+to_string(i)+".tables");
            ifstream ifs(fn);
            string ostr((istreambuf_iterator<char>(ifs)), (istreambuf_iterator<char>()));
            ifs.close();
            string_replace_all(ostr, "\"", "");
            Parser tp;
            _tables[i] = tp.Read(ostr);
            auto nn = _tables[i].op(1).nops();
            if(i==si || nmin>nn) nmin = nn;
        }
        exvector tables;
        exvector keys;
        ex convs = 0;
        for(int i=si; i<ei; i++) {
            auto nn = _tables[i].op(1).nops();
            if(nmin<nn) {
                cout << "ThieleTables: current table will be dropped." << endl;
                continue;
            }
            if(!is_a<lst>(convs)) convs = _tables[i].op(1); // second part is equal
            else if(!convs.is_equal(_tables[i].op(1))) {
                throw Error("RRT: convs are different.");
            }
            tables.push_back(_tables[i].op(0));
            keys.push_back(i);
        }
        
        vector<pair<ex,vector<pair<ex,exvector>>>> int_mi_cs_vec;
        int ci = 0;
        for(int i=0; i<tables.size(); i++) {
            if(i==0) {
                for(auto ti : tables[0]) {
                    pair<ex,vector<pair<ex,exvector>>> int_mi_cs;
                    int_mi_cs.first = ti.op(0);
                    for(auto tii : ti.op(1)) {
                        pair<ex,exvector> mi_cs;
                        mi_cs.first = tii.op(0);
                        mi_cs.second.push_back(tii.op(1));
                        int_mi_cs.second.push_back(mi_cs);
                    }
                    int_mi_cs_vec.push_back(int_mi_cs);
                }
            } else {
                for(int i2=0; i2<int_mi_cs_vec.size(); i2++) {
                    auto ti = tables[i].op(i2);
                    auto & int_mi_cs = int_mi_cs_vec[i2];
                    if(int_mi_cs.first != ti.op(0)) throw Error("E1");
                    for(int j=0; j<int_mi_cs.second.size(); j++) {
                        auto & mi_cs = int_mi_cs.second[j];
                        if(mi_cs.first != ti.op(1).op(j).op(0)) throw Error("E2");
                        mi_cs.second.push_back(ti.op(1).op(j).op(1));
                    }
                }
            }
        }
        
        std::ofstream ofs;
        ofs.open(filename, ios::out);
        ofs << "{" << endl;
        ofs << "    {" << endl;
        for(int i=0; i<int_mi_cs_vec.size(); i++) {
            auto imc = int_mi_cs_vec[i];
            ofs << "        {" << imc.first << "," << endl;
            ofs << "            {" << endl;
            for(int j=0; j<imc.second.size(); j++) {
                auto item = imc.second[j];
                ofs << "                {" << item.first << ", ";
                exvector vals(item.second.begin(), item.second.end());
                auto res = Thiele(keys, vals, d);
                ofs << "\"" << res << "\"}";
                if(j+1<imc.second.size()) ofs << ",";
                ofs << endl;
            }
            ofs << "            }" << endl;
            ofs << "        }";
            if(i+1<int_mi_cs_vec.size()) ofs << ",";
            ofs << endl;
        }
        ofs << "    }," << endl << "    {" << endl;
        for(int i=0; i<convs.nops(); i++) {
            auto item = convs.op(i);
            ofs << "        {" << item.op(0) << ", " << item.op(1) << "}";
            if(i+1<convs.nops()) ofs << ",";
            ofs << endl;
        }
        ofs << "    }" << endl;
        ofs << "}" << endl;
        ofs.close();
    }

    /**
     * @brief Export start config intgral etc. files
     */
    void FIRE::Export() {
        if(WorkingDir=="") WorkingDir = to_string(getpid());
        // check already exported, used while restarting from check point
        if(file_exists(WorkingDir+"/"+to_string(ProblemNumber)+".intg")) return;
    
        int pdim = Propagators.nops();
        if(Integrals.nops()<1) return;
        int pn = 0; // to avoid unsigned short overflow in FIRE
        lst InExternal;
        for(auto ii : Internal) InExternal.append(ii);
        for(auto ii : External) InExternal.append(ii);
        
        if(ISP.nops()<1) {
            for(auto it : Internal) {
                for(auto ii : InExternal) ISP.append(it*ii);
            }
            ISP.sort();
            ISP.unique();
        }
        
        if(ISP.nops() > pdim) {
            cout << "ISP = " << ISP << endl;
            cout << "Propagators = " << Propagators << endl;
            throw Error("FIRE::Export: #(ISP) > #(Propagators).");
        }
        
        lst sp2s, s2sp, ss;
        int _pic=0;
        for(auto item : ISP) {
            _pic++;
            Symbol si("P"+to_string(_pic));
            ss.append(si);
            sp2s.append(item==si);
            s2sp.append(si==item);
        }
        
        lst eqns;
        for(int i=0; i<ISP.nops(); i++) { // note NOT pdim
            auto eq = expand(Propagators.op(i)).subs(iEpsilon==0); // drop iEpsilon
            eq = eq.subs(sp2s, algbr);
            eq = eq.subs(Replacements, algbr);
            if(eq.has(iWF(w))) throw Error("FIRE::Export, iWF used in eq.");
            eqns.append(eq == iWF(i));
        }
        auto s2p = lsolve(eqns, ss);
        if(s2p.nops() != ISP.nops()) {
            cout << ISP << endl << s2p << endl << eqns << endl;
            throw Error("FIRE::Export: lsolve failed.");
        }
        
        if(DSP.nops()<1) {
            for(auto p1 : Internal)
            for(auto p2 : InExternal)
            DSP.append(lst{p1,p2});
        }

        vector<exmap> ibps;
        exvector IBPvec;
        lst ns0;
        for(int i=0; i<pdim; i++) ns0.append(0);
        for(auto sp : DSP) {
            symbol ss;
            auto ilp = sp.op(0);
            auto iep = sp.op(1);
            lst dp_lst;
            for(int i=0; i<pdim; i++) {
                dp_lst.append(Propagators.op(i).subs(ilp==ss).diff(ss).subs(ss==ilp));
            } 
            
            exmap nc_map;
            for(int i=0; i<pdim; i++) { // diff on each propagator
                auto ns = ns0;
                ns.let_op(i) = ns.op(i)+1; // note the covention
                auto tmp = dp_lst.op(i) * iep;
                tmp = expand(tmp);
                tmp = tmp.subs(Replacements, algbr);
                tmp = tmp.subs(sp2s, algbr);
                tmp = tmp.subs(s2p, algbr);
                tmp = tmp.subs(Replacements, algbr);
                
                tmp = ex(0) - (Shift[i+1]+a(i+1))*tmp; // note Shift here

                for(int j=0; j<pdim; j++) {
                    auto cj = tmp.coeff(iWF(j));
                    if(is_zero(cj)) continue;
                    lst cns = ns;
                    cns.let_op(j) = cns.op(j)-1; // note the covention
                    nc_map[cns] = nc_map[cns] + cj;
                }
                tmp = tmp.subs(iWF(w)==0); // constant term
                if(!is_zero(tmp)) nc_map[ns] = nc_map[ns] + tmp;
            }
            
            auto cilp = iep.coeff(ilp);
            if(!is_zero(cilp)) nc_map[ns0] = nc_map[ns0] + d*cilp;
            bool ok = false;
            for(auto nc : nc_map) {
                if(!is_zero(nc.second)) {
                    ok = true;
                    IBPvec.push_back(nc.second);
                }
            }
            if(ok) ibps.push_back(nc_map);
        }

        auto Variables = gather_symbols(IBPvec);
        for(auto e : External) {
            if(Variables.has(e)) {
                cout << "IBPvec: " << IBPvec << endl;
                cout << "Replacements: " << Replacements << endl;
                cout << "Variables: " << Variables << endl;
                throw Error("FIRE: Replacement NOT work!");
            }
        }
        
        ostringstream start;

        start << "ExampleDimension[" << pn << "]=" << pdim << endl << endl;
        start << "ProblemNumber=" << pn << endl << endl;
        
        // .start - SBasisL
        if(Version==5) {
            PermutationsR(2, pdim, [pdim,pn,&start](const int *ns) {
                start << "SBasisL[" << pn << ",{";
                for(int i=0; i<pdim; i++) start << (ns[i]<1 ? -1 : 1) << (i<pdim-1 ? "," : "");
                start << "}]=0" << endl << endl;
            });
        }
        
        // .start - SBasis0L, SBasis0D && SBasis0C
        start << "SBasis0L[" << pn << "]=" << ibps.size() << endl << endl;
        ostringstream oss;
        for(int i=0; i<ibps.size(); i++) {
            start << "SBasis0D[" << pn << "," << (i+1) << "]=";
            lst items;
            for(auto kv : ibps[i]) {
                items.append(kv.first);
                if(Version==5) {
                    oss << "SBasis0C[" << pn << "," << (i+1) << "," << kv.first << "]=" << 
                    collect_common_factors(kv.second.normal()) << endl << endl;
                } else {
                    lst olst;
                    auto cv_lst = collect_lst(kv.second, a(w));
                    for(auto item : cv_lst) {
                        auto cc = item.op(0);
                        auto cv = item.op(1);
                        cc = collect_common_factors(cc.normal());
                        if(is_zero(cc)) continue;
                        if(is_zero(cv-1)) cv=0;
                        else cv = cv.subs(a(w)==w);
                        olst.append(lst{cc, cv});
                    }
                    oss << "SBasis0C[" << pn << "," << (i+1) << "," << kv.first << "]=" << olst << endl << endl;
                }
            }
            start << items << endl << endl;
        }
        start << oss.str();
        
        // .start - SBasisS
        start << "SBasisS[" << pn << "]={{{";
        // 1,2,3,...
        for(int i=0; i<pdim; i++) start << (i+1) << (i<pdim-1 ? "," : "");
        start << "},{";
        // 1,1,1,...
        for(int i=0; i<pdim; i++) start << 1 << (i<pdim-1 ? "," : "");
        start << "},{";
        // 0,0,0,...
        for(int i=0; i<pdim; i++) start << 0 << (i<pdim-1 ? "," : "");
        start << "}}}" << endl << endl;
        
        // .start - SBasisR
        lst Rlst;
        Rlst.append(lst{});
        for(int i=0; i<pdim; i++) {
            let_op_append(Rlst, 0, -1);
        }
        for(auto lpi : Internal) {
            vector<int> ns_vec;
            lst ns0;
            for(int i=0; i<pdim; i++) ns0.append(1);
            for(int i=0; i<pdim; i++) {
                if(Propagators.op(i).has(lpi)) ns0.let_op(i) = -1;
                else ns_vec.push_back(i);
            }
            long long tot = std::pow(2LL,ns_vec.size());
            for(long long n=0; n<tot; n++) {
                int cn = n;
                lst ns1 = ns0;
                for(int j=0; j<ns_vec.size(); j++) {
                    if((cn%2)==1) ns1.let_op(ns_vec[j]) = -1;
                    cn /= 2;
                }
                Rlst.append(ns1);
            } 
        }
        
        // Lee Zero Sector
        if(true) {
            exset sectors;
            IsAlwaysZero = true;
            lst ns0;
            for(int i=0; i<pdim; i++) ns0.append(1);
            long long tot = std::pow(2LL,pdim);
            
            for(long long n=0; n<tot; n++) {
                int cn = n;
                lst sector = ns0;
                for(int j=0; j<pdim; j++) {
                    if((cn%2)==1) sector.let_op(j) = 0;
                    cn /= 2;
                }
                sectors.insert(sector);
            }
            
            while(!sectors.empty()) {
                auto first = *(sectors.begin());
                sectors.erase(first);
                
                auto sector = ex_to<lst>(first);
                if(IsZero(sector)) { // from LiteRed: all subsector is zero
                    lst ns1 = sector;
                    for(int j=0; j<pdim; j++) {
                        if(sector.op(j).is_zero()) ns1.let_op(j) = -1;
                    }
                    Rlst.append(ns1);
                    exset subsets;
                    for(auto item : sectors) {
                        bool ok = true;
                        for(int j=0; j<pdim; j++) {
                            if(sector.op(j)<item.op(j)) { ok = false; break; }
                        }
                        if(ok) subsets.insert(item);
                    }
                    for(auto item : subsets) {
                        sectors.erase(item);
                        lst ns1 = ex_to<lst>(item);
                        for(int j=0; j<pdim; j++) {
                            if(item.op(j).is_zero()) ns1.let_op(j) = -1;
                        }
                        Rlst.append(ns1);
                    }
                } else { // from LiteRed: all supersector is non-zero
                    if(IsAlwaysZero) IsAlwaysZero = false;
                    exset supsets;
                    for(auto item : sectors) {
                        bool ok = true;
                        for(int j=0; j<pdim; j++) {
                            if(sector.op(j)>item.op(j)) { ok = false; break; }
                        }
                        if(ok) supsets.insert(item);
                    }
                    for(auto item : supsets) sectors.erase(item);
                }
            }
            
            if(IsAlwaysZero) {
                Rules.remove_all();
                for(auto ii : Integrals) Rules.append(F(ProblemNumber, ii)==0);
                return;
            }
        }
        
        // handle Cut Propagators
        if(Cuts.nops()>0) {
            for(auto cx : Cuts) {
                int ci = ex_to<numeric>(cx-1).to_int(); // start from 1 in Cuts
                lst ns0;
                for(int i=0; i<pdim; i++) ns0.append(1);
                ns0.let_op(ci) = -1;
                for(int n=0; n<std::pow(2,pdim-1); n++) {
                    int cn = n;
                    lst ns1 = ns0;
                    for(int j=0; j<pdim; j++) {
                        if(ci==j) continue;
                        if((cn%2)==1) ns1.let_op(j) = -1;
                        cn /= 2;
                    }
                    Rlst.append(ns1);
                } 
            }
        }
        
        Rlst.sort();
        Rlst.unique();
        //sort_lst(Rlst);
        
        for(auto iR : Rlst) {
            start << "SBasisR[" << pn << ",{";
            for(int i=0; i<pdim; i++) start << iR.op(i) << (i<pdim-1 ? "," : "");
            start << "}]=True" << endl << endl;
        }
        
        // .start - Others
        start << "SBasisRL[" << pn << "]=0" << endl << endl;
        start << "HPI[" << pn << "]={}" << endl << endl;
        
        string sss = start.str();
        string_replace_all(sss, "=", " = ");
        string_replace_all(sss, ",", ", ");
        
        string spn = to_string(ProblemNumber);
        if(!dir_exists(WorkingDir)) system(("mkdir -p " + WorkingDir).c_str());
        
        if(WorkingDir.length()<1) WorkingDir = to_string(getpid());
        if(!dir_exists(WorkingDir)) system(("mkdir -p "+WorkingDir).c_str());
        ofstream start_out(WorkingDir+"/"+spn+".start");
        start_out << sss << endl;
        start_out.close();
        
        // .config
        int ct = 1; // 1-poly, 1-poly+prime, 2-poly+prime+mpq, since fire-config provided, we only generate poly
        for(int ci=0; ci<ct; ci++) { // 0-poly, 1-prime, 2-mpq
            ostringstream config;
            if(Version>5) config << "#compressor none" << endl;
            if(Version==5) config << "#bucket 20" << endl;
            config << "#threads " << Threads << endl;
            if(fThreads>0) config << "#fthreads " << fThreads << endl;
            else if(ci==2) config << "#fthreads 1" << endl;
            else config << "#fthreads " << Threads << endl;
            if(Version>5) { // FIRE6 or lator
                if(sThreads>0) config << "#sthreads " << sThreads << endl;
                if(lThreads>1) config << "#lThreads " << lThreads << endl;
                if(pos_pref!=1) config << "#pos_pref "<< pos_pref << endl;
                if(ci==1) config << "#prime 1" << endl;
            }
            config << "#variables ";
            bool first = true;
            exvector ev_sort;
            for(auto v : Variables) {
                auto fw = fermat_weight.find(v);
                if(fw!=fermat_weight.end()) ev_sort.push_back(lst{ fw->second, v });
                else ev_sort.push_back(lst{ 0, v });
            }
            sort_vec(ev_sort);
            for(auto nv : ev_sort) { 
                const symbol & s = ex_to<symbol>(nv.op(1));
                if(!islower(s.get_name()[0])) {
                    cout << "Replacements: " << Replacements << endl;
                    cout << "IBPvec: " << IBPvec << endl;
                    throw Error("FIRE: Fermat requires a name must begin with a lower case letter: "+s.get_name());
                }
                config << (first ? "" : ",") << s; 
                first=false; 
            }
            config << endl;
            config << "#database db" << ProblemNumber << endl;
            config << "#start" << endl;
            config << "#problem " << pn << " " << ProblemNumber << ".start" << endl;
            if(PIntegrals.nops()>0) {
                ostringstream oss;
                oss << "{";
                int nn = PIntegrals.nops();
                for(int i=0; i<nn; i++) {
                    if(PIntegrals.op(i).nops()!=pdim) throw Error("FIRE::Export@1, Index dimension NOT match Propagators.");
                    oss << "{" << pn << "," << PIntegrals.op(i) << (i<nn-1 ? "}," : "}");
                }
                oss << "}";
                ofstream pref_out(WorkingDir+"/"+spn+".pref");
                pref_out << oss.str() << endl;
                pref_out.close();
                config << "#preferred " << ProblemNumber << ".pref" << endl;
            }
            config << "#integrals " << ProblemNumber << ".intg" << endl;
            config << "#output " << ProblemNumber << ".tables" << endl;
            
            string cpos = "";
            if(ci==1) cpos = "p";
            else if(ci==2) cpos = "q";
            ofstream config_out(WorkingDir+"/"+spn+cpos+".config");
            config_out << config.str() << endl;
            config_out.close();
        }
        
        // *.intg
        ostringstream intg;
        intg << "{";
        for(int i=0; i<Integrals.nops(); i++) {
            if(Integrals[i].nops()!=pdim) throw Error("FIRE::Export@2, Index dimension NOT match Propagators.");
            intg << "{" << pn << "," << Integrals[i] << (i<Integrals.nops()-1 ? "}," : "}");
        }
        intg << "}" << endl;
        
        ofstream intg_out(WorkingDir+"/"+spn+".intg");
        intg_out << intg.str() << endl;
        intg_out.close();
    }
    
    /**
     * @brief Run FIRE reduction 
     */
    void FIRE::Run() {
        if(IsAlwaysZero) return;
        if(file_exists(WorkingDir + "/" + to_string(ProblemNumber) + ".tables")) return;
        int tried = 0;
        while(tried<3) {
            tried++;
            ostringstream cmd;
            cmd << "cd " << WorkingDir << " && $(which FIRE" << Version << ")";
            if(Version>5) cmd << " -silent -parallel";
            cmd << " -c " << ProblemNumber << " >/dev/null";
            system(cmd.str().c_str());
            system(("rm -rf "+WorkingDir+"/db"+to_string(ProblemNumber)).c_str());
            if(file_exists(WorkingDir + "/" + to_string(ProblemNumber) + ".tables")) break;
            sleep(3);
        }
        if(!file_exists(WorkingDir + "/" + to_string(ProblemNumber) + ".tables")) {
            cout << "Propagators: " << Propagators << endl;
            throw Error("FIRE::Run failed!");
        }    
    }
    
    /**
     * @brief Import tables  
     */
    void FIRE::Import() {
        if(IsAlwaysZero) return;
        string spn = to_string(ProblemNumber);
        ifstream ifs(WorkingDir+"/"+spn+".tables");
        string ostr((istreambuf_iterator<char>(ifs)), (istreambuf_iterator<char>()));
        ifs.close();
        
        string_replace_all(ostr, "\"", "");
        
        Parser tp;
        auto tp_lst = tp.Read(ostr);
        exmap id2F;
        for(auto item : tp_lst.op(1)) {
            if(!is_zero(item.op(1).op(0))) throw Error("FIRE::Reduce: pn is NOT 0.");
            id2F[item.op(0)] = F(ProblemNumber, item.op(1).op(1));
        }
        
        for(auto item : tp_lst.op(0)) {
            ex left = item.op(0).subs(id2F);
            ex right = 0;
            for(auto it : item.op(1)) {
                right += it.op(0).subs(id2F) * it.op(1);
            }
            if(left.is_equal(right)) MIntegrals.append(left);
            else Rules.append(left==right);
        }
        MIntegrals.sort();
        MIntegrals.unique();
        sort_lst(MIntegrals);
     
        // handle Cuts not equal 1, using #preferred
        if(reCut && Cuts.nops()>0) {
            lst ois, vis;
            auto mis = MIntegrals;
            MIntegrals.remove_all();
            int nrun = -1;
            while(true) {
                nrun++;
                if(nrun>100) break;
                
                lst iis, pis;
                bool allOK = true;
                for(auto item : mis) {
                    lst mi = ex_to<lst>(item.op(1));
                    bool isOK = true;
                    for(auto cx : Cuts) {
                        if(!is_zero(mi.op(ex_to<numeric>(cx).to_int()-1)-1)) {
                            isOK = false;
                            allOK = false;
                            break;
                        }
                    }
                    if(nrun==0 && isOK) MIntegrals.append(item);
                    else if(!isOK) {
                        if(nrun==0) {
                            ois.append(item);
                            vis.append(item);
                        }
                        iis.append(mi);
                        auto pi = mi;
                        vector<int> ipos, ineg;
                        for(int i=0; i<mi.nops(); i++) {
                            bool isCut = false;
                            for(auto cx : Cuts) {
                                if(is_zero(cx-1-i)) {
                                    isCut = true;
                                    break;
                                }
                            }
                            if(isCut) pi.let_op(i)=1;
                            else if(mi.op(i)<=0) ineg.push_back(i);
                            else ipos.push_back(i);
                        }
                        
                        lst tis;
                        int max = 2;
                        ex total = pow(numeric(max), ineg.size());
                        for(numeric in=0; in<total; in++) {
                            auto cin = in;
                            auto pi2 = pi;
                            for(int i=0; i<ineg.size(); i++) {
                                int re = mod(cin,max).to_int();
                                pi2.let_op(ineg[i]) = ex(0)-re;
                                cin = (cin-re)/max;
                            }
                            tis.append(pi2);
                        }
                        
                        int max2 = 2*max+1;
                        total = pow(numeric(max2), ipos.size());
                        for(auto pi : tis) {
                            for(numeric in=0; in<total; in++) {
                                auto cin = in;
                                auto pi2 = pi;
                                for(int i=0; i<ipos.size(); i++) {
                                    int re = mod(cin,max2).to_int();
                                    pi2.let_op(ipos[i]) = re-max;
                                    cin = (cin-re)/max2;
                                }
                                pis.append(pi2);
                            }
                        }
                    }
                }
                if(allOK) break;
                
                // Reduce again
                pis.sort();
                pis.unique();
                iis.sort();
                iis.unique();
                if(pis.nops()>0) {
                    FIRE fire;
                    fire.Propagators = Propagators;
                    fire.Internal = Internal;
                    fire.External = External;
                    fire.Replacements = Replacements;
                    fire.ProblemNumber = ProblemNumber;
                    fire.ISP = ISP;
                    fire.DSP = DSP;
                    fire.Cuts = Cuts;
                    fire.reCut = false;
                    fire.Shift = Shift;
                    fire.Integrals = iis;
                    fire.PIntegrals = pis;
                    fire.WorkingDir = WorkingDir + "_C"+to_string(ProblemNumber);
                    fire.Reduce();
                    system(("rm -rf " + fire.WorkingDir).c_str());
                    
                    auto vis_chk = vis;
                    vis_chk = ex_to<lst>(subs(vis_chk, fire.Rules));
                    if(vis_chk.is_equal(vis)) break;
                    vis = vis_chk;
                    exset fs;
                    find(vis,F(w1,w2),fs);
                    mis.remove_all();
                    for(auto item : fs) mis.append(item);
                } else break;
            }
            
            exmap c2m;
            for(int i=0; i<ois.nops(); i++) c2m[ois.op(i)] = vis.op(i);
            for(auto item : mis) MIntegrals.append(item);
            MIntegrals.sort();
            MIntegrals.unique();
            sort_lst(MIntegrals);
            
            auto rules = Rules;
            Rules.remove_all();
            lst rIntegrals;
            for(auto r : rules) {
                auto ri = r.op(0);
                bool isMI = false;
                for(auto item : MIntegrals) {
                    if(item.is_equal(ri)) {
                        isMI = true;
                        rIntegrals.append(ri);
                        break;
                    }
                }
                if(isMI) continue;
                auto rv = r.op(1);
                rIntegrals.append(rv);
                Rules.append(ri==rv.subs(c2m));
            }
            
            rIntegrals.sort();
            rIntegrals.unique();
            exset fs;
            find(rIntegrals,F(w1,w2),fs);
            MIntegrals.remove_all();
            for(auto fi : fs) MIntegrals.append(fi);
            sort_lst(MIntegrals);
        }
    }
    
}
