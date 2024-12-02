/**
 * @file
 * @brief IBP with FIRE
 */
 
#include "IBP.h"
#include <cmath>

#include "exFlint.h"

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
    
    inline lst syms(const exvector & ev) {
        exset ss;
        for(auto e : ev) {
            for(const_preorder_iterator i=e.preorder_begin(); i!=e.preorder_end(); ++i)
                if(is_a<symbol>(*i)) ss.insert(*i);
        }
        lst ls;
        for(auto item : ss) ls.append(item);
        return ls;
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
        for(int i=si; i<=ei; i++) {
            cout << "\r                               \r" << flush;
            cout << "Reading tables: " << ei << "|" << i << flush;
            string fn = filename;
            string_replace_all(fn, ".tables", "-"+to_string(i)+".tables");
            ifstream ifs(fn);
            string ostr((istreambuf_iterator<char>(ifs)), (istreambuf_iterator<char>()));
            ifs.close();
            for(auto & c : ostr) if(c=='\"') c = ' ';
            Parser tp;
            _tables[i] = tp.Read(ostr);
            auto nn = _tables[i].op(1).nops();
            if(i==si || nmin>nn) nmin = nn;
        }
        cout << endl;
        
        exvector tables;
        exvector keys;
        ex convs = 0;
        for(int i=si; i<=ei; i++) {
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
            cout << "\r                               \r" << flush;
            cout << "Thiele: " << int_mi_cs_vec.size() << "|" << i+1 << flush;
            auto imc = int_mi_cs_vec[i];
            ofs << "        {" << imc.first << "," << endl;
            ofs << "            {" << endl;
            
            auto nnn = imc.second.size();
            vector<fmpz_mpoly_q_struct*> res_vec(nnn);
            vector<lst> xs_vec(nnn);
            #pragma omp parallel for schedule(dynamic, 1)
            for(int j=0; j<nnn; j++) {
                auto item = imc.second[j];
                xs_vec[j] = syms(item.second);
                fmpz_mpoly_ctx_t ctx;
                fmpz_mpoly_ctx_init(ctx, xs_vec[j].nops(), ORD_LEX);
                res_vec[j] = (fmpz_mpoly_q_struct*) flint_malloc(sizeof(fmpz_mpoly_q_struct));
                fmpz_mpoly_q_init(res_vec[j], ctx);
                int n = keys.size();
                vector<fmpz_mpoly_q_struct*> key_vec(n);
                vector<fmpz_mpoly_q_struct*> val_vec(n);
                vector<fmpz_mpoly_q_struct*> coeff_vec(n);
                for(int k=0; k<n; k++) {
                    key_vec[j] = (fmpz_mpoly_q_struct*) flint_malloc(sizeof(fmpz_mpoly_q_struct));
                    fmpz_mpoly_q_init(key_vec[k], ctx);
                    _to_(xs_vec[j], key_vec[k], ctx, item.first.op(k));
                    val_vec[j] = (fmpz_mpoly_q_struct*) flint_malloc(sizeof(fmpz_mpoly_q_struct));
                    fmpz_mpoly_q_init(val_vec[k], ctx);
                    _to_(xs_vec[j], val_vec[j], ctx, item.second[k]);
                    coeff_vec[j] = (fmpz_mpoly_q_struct*) flint_malloc(sizeof(fmpz_mpoly_q_struct));
                    fmpz_mpoly_q_init(coeff_vec[k], ctx);
                    fmpz_mpoly_q_zero(coeff_vec[k], ctx);
                }
                
                fmpz_mpoly_q_t t1, t2;
                fmpz_mpoly_q_init(t1, ctx);
                fmpz_mpoly_q_init(t2, ctx);
                for(int i=0; i<n; i++) {
                    fmpz_mpoly_q_set(t1, val_vec[i], ctx);
                    for(int j=0; j<i; j++) {
                        if(fmpz_mpoly_q_equal(t1, coeff_vec[j], ctx)) {
                            n = i-1;
                            goto out;
                        }
                        fmpz_mpoly_q_sub(t2, key_vec[i], key_vec[j], ctx);
                        fmpz_mpoly_q_sub(t1, t1, coeff_vec[j], ctx);
                        fmpz_mpoly_q_div(t1, t2, t1, ctx);
                    }
                    fmpz_mpoly_q_set(coeff_vec[i], t1, ctx);
                }
                throw Error("Thiele unstable!");
                out: ;
                fmpz_mpoly_q_set(t1, coeff_vec[n], ctx);
                fmpz_mpoly_q_t dx;
                fmpz_mpoly_q_init(dx, ctx);
                _to_(xs_vec[j], dx, ctx, d);
                for(int i=n-1; i>=0; i--) {
                    fmpz_mpoly_q_sub(t2, dx, key_vec[i], ctx);
                    fmpz_mpoly_q_div(t2, t2, t1, ctx);
                    fmpz_mpoly_q_add(t1, coeff_vec[i], t2, ctx);
                }
                fmpz_mpoly_q_clear(dx, ctx);
                fmpz_mpoly_q_set(res_vec[i], t1, ctx);
                for(int k=0; k<n; k++) {
                    fmpz_mpoly_q_clear(key_vec[k], ctx);
                    fmpz_mpoly_q_clear(val_vec[k], ctx);
                    fmpz_mpoly_q_clear(coeff_vec[k], ctx);
                }
                fmpz_mpoly_q_clear(t1, ctx);
                fmpz_mpoly_q_clear(t2, ctx);
                fmpz_mpoly_ctx_clear(ctx);
            }
            for(int j=0; j<imc.second.size(); j++) {
                fmpz_mpoly_ctx_t ctx;
                fmpz_mpoly_ctx_init(ctx, xs_vec[j].nops(), ORD_LEX);
                auto item = imc.second[j];
                ex res;
                res = _to_(xs_vec[j], res_vec[j], ctx);
                fmpz_mpoly_q_clear(res_vec[j], ctx);
                ofs << "                {" << item.first << ", ";
                ofs << "\"" << res << "\"}";
                if(j+1<imc.second.size()) ofs << ",";
                ofs << endl;
            }
            ofs << "            }" << endl;
            ofs << "        }";
            if(i+1<int_mi_cs_vec.size()) ofs << ",";
            ofs << endl;
        }
        cout << endl;
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
        if(WorkingDir=="") WorkingDir = to_string(getpid())+"IBP";
        int pdim = Propagator.nops();
        string spn = to_string(ProblemNumber);
        int pn = 0; // to avoid unsigned short overflow in FIRE
        
        if(!file_exists(WorkingDir+"/"+spn+".start") || !file_exists(WorkingDir+"/"+spn+".config")) {
            
            if(Integral.nops()<1) return;
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
                cout << "Propagator = " << Propagator << endl;
                throw Error("FIRE::Export: #(ISP) > #(Propagator).");
            }
            
            lst sp2s, s2sp, ss;
            int _pic=0;
            for(auto item : ISP) {
                _pic++;
                Symbol si("P"+to_string(_pic));
                ss.append(si);
                sp2s.append(w*item==w*si);
                sp2s.append(item==si);
                s2sp.append(si==item);
            }
            
            if(!has_w(Replacement)) {
                lst repl;
                for(auto item : Replacement) {
                    repl.append(item);
                    repl.append(w*item.op(0) == w*item.op(1));
                }
                Replacement = repl;
            }
            
            lst eqns;
            for(int i=0; i<ISP.nops(); i++) { // note NOT pdim
                auto eq = expand(Propagator.op(i)).subs(iEpsilon==0); // drop iEpsilon
                eq = eq.subs(sp2s);
                eq = eq.subs(Replacement);
                if(eq.has(iWF(w))) throw Error("FIRE::Export, iWF used in eq.");
                eqns.append(eq == iWF(i));
            }
            auto s2p = lsolve(eqns, ss);
            if(s2p.nops() != ISP.nops()) {
                cout << "ISP=" << ISP << endl;
                cout << "Propagator=" << Propagator << endl;
                cout << "eqns=" << eqns << endl;
                throw Error("FIRE::Export: lsolve failed.");
            }
            
            // FIRE version
            if(true)
            if(DSP.nops()<1) {
                for(auto p1 : Internal)
                for(auto p2 : InExternal)
                DSP.append(lst{p1,p2});
            }
            
            // Lee version
            if(false)
            if(DSP.nops()<1) {
                for(int i=0; i<Internal.nops(); i++)
                    DSP.append(lst{Internal.op(i), Internal.op(i==Internal.nops()-1 ? 0 : i+1)});
                for(auto p2 : InExternal) DSP.append(lst{Internal.op(0), p2});
                DSP.append(lst{Internal, Internal});
                allIBP = true;
            }

            vector<exmap> ibps;
            exvector IBPvec;
            lst ns0;
            for(int i=0; i<pdim; i++) ns0.append(0);
            
            for(auto sp : DSP) {
                auto ilp_lst = sp.op(0);
                auto iep_lst = sp.op(1);
                if(!is_a<lst>(ilp_lst)) {
                    ilp_lst = lst{ ilp_lst };
                    iep_lst = lst{ iep_lst };
                }
                
                exmap nc_map;
                for(int ilst=0; ilst<ilp_lst.nops(); ilst++) {
                    symbol ss;
                    auto ilp = ilp_lst.op(ilst);
                    auto iep = iep_lst.op(ilst);
                    lst dp_lst;
                    for(int i=0; i<pdim; i++) {
                        dp_lst.append(Propagator.op(i).subs(ilp==ss).diff(ss).subs(ss==ilp));
                    }
                    
                    for(int i=0; i<pdim; i++) { // diff on each propagator
                        auto ns = ns0;
                        ns.let_op(i) = ns.op(i)+1; // note the covention
                        auto tmp = dp_lst.op(i) * iep;
                        tmp = expand(tmp);
                        tmp = tmp.subs(Replacement);
                        tmp = tmp.subs(sp2s);
                        tmp = tmp.subs(s2p);
                        tmp = tmp.subs(Replacement);
                        
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
                }
                
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
                    cout << "Replacement: " << Replacement << endl;
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
                    if(Propagator.op(i).has(lpi)) ns0.let_op(i) = -1;
                    else ns_vec.push_back(i);
                }
                size_t tot = std::pow(2LL,ns_vec.size());
                for(size_t n=0; n<tot; n++) {
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
                size_t tot = std::pow(2LL,pdim);
                
                for(size_t n=0; n<tot; n++) {
                    int cn = n;
                    lst sector = ns0;
                    for(int j=0; j<pdim; j++) {
                        if((cn%2)==1) sector.let_op(j) = 0;
                        cn /= 2;
                    }
                    if(SECTOR.nops()==pdim) {
                        for(int j=0; j<pdim; j++) if(sector.op(j)>SECTOR.op(j)) goto add_sector_done;
                    }
                    sectors.insert(sector);
                    add_sector_done: ;
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
                    for(auto ii : Integral) Rules.append(F(ProblemNumber, ii)==0);
                    return;
                }
            }
            
            // handle Cut Propagator
            if(Cut.nops()>0) {
                for(auto cx : Cut) {
                    int ci = ex_to<numeric>(cx-1).to_int(); // start from 1 in Cut
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
            
            if(!dir_exists(WorkingDir)) auto rc = system(("mkdir -p " + WorkingDir).c_str());
            
            ofstream start_out(WorkingDir+"/"+spn+".start");
            start_out << sss << endl;
            start_out.close();
        
            // .config
            {
                ostringstream config;
                
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
                    if(false && !islower(s.get_name()[0])) {
                        cout << "Replacement: " << Replacement << endl;
                        cout << "IBPvec: " << IBPvec << endl;
                        throw Error("FIRE: Fermat requires a name must begin with a lower case letter: "+s.get_name());
                    }
                    config << (first ? "" : ",") << s;
                    auto itr = NVariables.find(nv.op(1));
                    if(itr!=NVariables.end()) config << "->" << itr->second;
                    first=false; 
                }
                config << endl;
                
                if(PosPref!=1) config << "#pos_pref "<< PosPref << endl;
                config << "#database db" << ProblemNumber << endl;
                config << "##prime 111" << endl; // ## for comment
                if(allIBP) config << "#allIBP" << endl;
                config << "#start" << endl;
                
                if(SECTOR.nops()==pdim) {
                    int pmax = pdim;
                    for(int i=pdim-1; i>=0; i--) {
                        if(SECTOR.op(i)!=0) {
                            pmax = i+1;
                            break;
                        }
                    }
                    int pmin = -1;
                    for(int i=0; i<pdim; i++) {
                        if(SECTOR.op(i)!=0) {
                            pmin = i+1;
                            break;
                        }
                    }
                    config << "#problem " << pn << " ";
                    if(pmin>1) config << "|" << pmin << "," << pmax << "|";
                    else if(pmax!=pdim) config << "|" << pmax << "|";
                    config << ProblemNumber << ".start" << endl;
                } else config << "#problem " << pn << " " << ProblemNumber << ".start" << endl;
                
                config << "#output " << ProblemNumber << ".tables" << endl;
                
                if(PIntegral.nops()>0) {
                    ostringstream oss;
                    oss << "{";
                    int nn = PIntegral.nops();
                    for(int i=0; i<nn; i++) {
                        if(PIntegral.op(i).nops()!=pdim) throw Error("FIRE::Export@1, Index dimension NOT match Propagator.");
                        oss << "{" << pn << "," << PIntegral.op(i) << (i<nn-1 ? "}," : "}");
                    }
                    oss << "}";
                    ofstream pref_out(WorkingDir+"/"+spn+".pref");
                    pref_out << oss.str() << endl;
                    pref_out.close();
                    config << "#preferred " << ProblemNumber << ".pref" << endl;
                }
                
                config << "#integrals " << ProblemNumber << ".intg" << endl;
                
                ofstream config_out(WorkingDir+"/"+spn+".config");
                config_out << config.str() << endl;
                config_out.close();
            }

        }
        
        // *.intg
        if(!file_exists(WorkingDir+"/"+spn+".intg")) {
            ostringstream intg;
            intg << "{";
            for(int i=0; i<Integral.nops(); i++) {
                if(Integral[i].nops()!=pdim) throw Error("FIRE::Export@2, Index dimension NOT match Propagator.");
                intg << "{" << pn << "," << Integral[i] << (i<Integral.nops()-1 ? "}," : "}");
            }
            intg << "}" << endl;
            
            ofstream intg_out(WorkingDir+"/"+spn+".intg");
            intg_out << intg.str() << endl;
            intg_out.close();
        }
        
        // export .gar here
        if(!file_exists(WorkingDir+"/"+spn+".gar")) {
            garWrite(TO(), WorkingDir+"/"+spn+".gar");
        }
    }
    
    /**
     * @brief Run FIRE reduction 
     */
    void FIRE::Run() {
        if(IsAlwaysZero) return;
        if(Execute=="") Execute = InstallPrefix + "/FIRE/M/FIRE";
        if(file_exists(WorkingDir + "/" + to_string(ProblemNumber) + ".tables")) return;
        if(Integral.nops()<1) return;
        int tried = 0;
        while(tried<3) {
            tried++;
            ostringstream cmd;
            cmd << "cd " << WorkingDir << " && " << Execute;
            if(Version>5) {
                if(opt=="") cmd << " -silent -t " << T1 << "," << T2 << " -lt " << LT1 << "," << LT2;
                else cmd << " " << opt;
            }
            cmd << " -c " << ProblemNumber << " 1>/dev/null 2>&1";
            auto rc = system(cmd.str().c_str());
            rc = system(("rm -rf "+WorkingDir+"/db"+to_string(ProblemNumber)).c_str());
            if(file_exists(WorkingDir + "/" + to_string(ProblemNumber) + ".tables")) break;
            sleep(3);
        }
        if(!file_exists(WorkingDir + "/" + to_string(ProblemNumber) + ".tables")) {
            cout << "Propagator: " << Propagator << endl;
            throw Error("FIRE::Run failed!");
        }    
    }
    
    /**
     * @brief Import tables  
     */
    void FIRE::Import() {
        if(IsAlwaysZero) return;
        if(Integral.nops()<1) return;
        string spn = to_string(ProblemNumber);
        ifstream ifs(WorkingDir+"/"+spn+".tables");
        string ostr((istreambuf_iterator<char>(ifs)), (istreambuf_iterator<char>()));
        ifs.close();
        
        //string_replace_all(ostr, "\"", "");
        for(auto & c : ostr) if(c=='"') c = ' '; 
        
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
            if(left.is_equal(right)) MIntegral.append(left);
            else Rules.append(left==right);
        }
        MIntegral.sort();
        MIntegral.unique();
        sort_lst(MIntegral);
     
        // handle Cut not equal 1, using #preferred
        if(reCut && Cut.nops()>0) {
            lst ois, vis;
            auto mis = MIntegral;
            MIntegral.remove_all();
            int nrun = -1;
            while(true) {
                nrun++;
                if(nrun>100) break;
                
                lst iis, pis;
                bool allOK = true;
                for(auto item : mis) {
                    lst mi = ex_to<lst>(item.op(1));
                    bool isOK = true;
                    for(auto cx : Cut) {
                        if(!is_zero(mi.op(ex_to<numeric>(cx).to_int()-1)-1)) {
                            isOK = false;
                            allOK = false;
                            break;
                        }
                    }
                    if(nrun==0 && isOK) MIntegral.append(item);
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
                            for(auto cx : Cut) {
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
                    fire.Propagator = Propagator;
                    fire.Internal = Internal;
                    fire.External = External;
                    fire.Replacement = Replacement;
                    fire.ProblemNumber = ProblemNumber;
                    fire.ISP = ISP;
                    fire.DSP = DSP;
                    fire.Cut = Cut;
                    fire.reCut = false;
                    fire.Shift = Shift;
                    fire.Integral = iis;
                    fire.PIntegral = pis;
                    fire.WorkingDir = WorkingDir + "_C"+to_string(ProblemNumber);
                    fire.Reduce();
                    auto rc = system(("rm -rf " + fire.WorkingDir).c_str());
                    
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
            for(auto item : mis) MIntegral.append(item);
            MIntegral.sort();
            MIntegral.unique();
            sort_lst(MIntegral);
            
            auto rules = Rules;
            Rules.remove_all();
            lst rIntegral;
            for(auto r : rules) {
                auto ri = r.op(0);
                bool isMI = false;
                for(auto item : MIntegral) {
                    if(item.is_equal(ri)) {
                        isMI = true;
                        rIntegral.append(ri);
                        break;
                    }
                }
                if(isMI) continue;
                auto rv = r.op(1);
                rIntegral.append(rv);
                Rules.append(ri==rv.subs(c2m));
            }
            
            rIntegral.sort();
            rIntegral.unique();
            exset fs;
            find(rIntegral,F(w1,w2),fs);
            MIntegral.remove_all();
            for(auto fi : fs) MIntegral.append(fi);
            sort_lst(MIntegral);
        }
    }
    
}
