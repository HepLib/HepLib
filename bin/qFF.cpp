#include "SD.h"

// Computation for quark Fragmentation Function

using namespace HepLib;

auto ep = SD::ep;
auto eps = SD::eps;
auto iEpsilon = SD::iEpsilon;
auto vs = SD::vs;
auto vz = SD::vz;

lst para_sym = lst { ep, eps, vs, vz, iEpsilon };

const char* cm_path = "cm";
const char* SD_path = "SD";
int verb = 0;
int epN = 0;
bool useq = false;
ex nL = 3;
ex nH = 1;

void ExportNull(const char* fn) {
    ofstream ofs;
    ofs.open(fn, ios::out);
    if (!ofs) throw runtime_error("failed to open final null file!");
    ofs << "Result is Null." << endl;
    ofs.close();
}

void Prepare(int idx) {
    
    symbol p("p"), n("n");
    symbol q1("q1"), q2("q2");
    symbol k1("k1"), k2("k2"), k3("k3");
    symbol K1("K1"), K2("K2"), K3("K3");
    
    FeynmanParameter fp;
    fp.LoopMomenta = lst{ q1 };
    fp.tLoopMomenta = lst{ K1 };
    ex S = 1; // phase space symmetry factor, note the ONE Quark
    
    ex kp = 1;
    ex m = 1;
    ex m2=m*m;
    ex zz = PL(0);
    
    ex k1p = z(1)*(1-zz)*kp;
    ex k1m = (pow(K1,2)+m2)/(2*k1p); // k1^2 = m^2
    ex k2p = z(2)*(1-zz)*kp;
    ex k2m = (pow(K2,2))/(2*k2p);
    ex k3p = z(3)*(1-zz)*kp;
    ex k3m = (pow(K3,2))/(2*k3p);
    
    if(fp.tLoopMomenta.nops()<1) zz = 1;
    ex pp = zz/2*kp;
    ex pm = m2/(2*pp);
    
    symtab table;
    table["ep"] = ep;
    table["eps"] = eps;
    table["iEpsilon"] = iEpsilon;
    
    table["p"] = p;
    table["n"] = n;
    table["k1"] = k1;
    table["K1"] = K1;
    table["k2"] = k2;
    table["K2"] = K2;
    table["k3"] = k3;
    table["K3"] = K3;
    table["q1"] = q1;
    table["q2"] = q2;
    table["m"] = m;
    table["zz"] = zz;
    table["nL"] = nL;
    table["nH"] = nH;
    table["ApartNull"] = 1;
    
    table["kp"] = kp;
    table["pp"] = pp;
    table["pm"] = pm;
    parser reader(table);

    ostringstream ifn;
    ifn << cm_path << "/" << idx << ".txt";
    ifstream in(ifn.str());
    
    auto cms = reader(in);
    
    fp.Prefactor = pow(2*Pi, (2*ep-4)*fp.LoopMomenta.nops()) * cms.op(2);
    fp.Propagators = ex_to<lst>(cms.op(0));
    fp.Exponents = ex_to<lst>(cms.op(1));
        
    ex NCS = pow(zz,(1-2*ep))/(8*3*kp*Pi);
    ex M = 2*m;
    ex nts = fp.tLoopMomenta.nops();
    ex psFactor;
    if(nts>0) psFactor = NCS*M/S*(-((pow(2,(2-nts-(3-2*ep)*nts))*pow(Pi,(1-(3-2*ep)*nts)))/(kp*(-1+zz))));
    else psFactor = NCS*M*4*Pi/kp;
    fp.Prefactor = fp.Prefactor * psFactor;
    
    fp.lReplacements[p*p] = m2;
    fp.lReplacements[n*n] = 0;
    fp.lReplacements[p*n] = pp;
    fp.lReplacements[k1*k1] = m2;
    fp.lReplacements[k2*k2] = 0;
    fp.lReplacements[k3*k3] = 0;
    fp.lReplacements[n*k1] = k1p;
    fp.lReplacements[n*k2] = k2p;
    fp.lReplacements[n*k3] = k3p;
    fp.lReplacements[p*k1] = k1p*pm+k1m*pp;
    fp.lReplacements[p*k2] = k2p*pm+k2m*pp;
    fp.lReplacements[p*k3] = k3p*pm+k3m*pp;
    fp.lReplacements[k1*k2] = k1p*k2m+k2p*k1m-K1*K2;
    fp.lReplacements[k1*k3] = k1p*k3m+k3p*k1m-K1*K3;
    fp.lReplacements[k3*k2] = k3p*k2m+k2p*k3m-K3*K2;
    
    fp.tReplacements[p*p] = m2;
    fp.tReplacements[n*n] = 0;
    fp.tReplacements[p*n] = pp;
    
    if(nts>0) fp.nReplacements[z(wild())] = ex(1)/nts;
    
    fp.nReplacements[ep] = ex(1)/11;
    fp.nReplacements[eps] = ex(1)/111;
    
    SD work;
    work.epN = epN;
    work.Verbose = verb;
    work.ParallelSymbols = para_sym;
    work.ParameterLB[0] = 0;
    work.ParameterUB[0] = 1;
    
    char *CFLAGS = getenv("SD_CFLAGS");
    work.CFLAGS = CFLAGS;
    work.CheckF1 = true;
    
    work.Initialize(fp);
    
    if(work.IsZero) {
        ostringstream ifn;
        ifn << SD_path << "/" << idx << ".null";
        ExportNull(ifn.str().c_str());
        return;
    }
    
    if(work.SecDec==NULL) work.SecDec = new SecDecG();
    if(work.Minimizer==NULL) work.Minimizer = new MinUit();
    
    int xn = work.Deltas[0].nops();
    exmap z2x;
    
    lst zs;
    ex zFactor = 1;
    for(int i=1; i<=nts; i++) {
        z2x[z(i)] = x(xn+i-1);
        zs.append(x(xn+i-1));
        zFactor /= x(xn+i-1);
    }
    if(nts>0) work.Deltas.push_back(zs);
    
    if(nts>0)
    for(auto &kv : work.FunExp) {
        
        kv.first.let_op(0) = kv.first.op(0) * zFactor;
        kv.first = lstHelper::subs(kv.first, z2x);
        
        auto tmp = collect_common_factors(kv.first.op(0));
        if(tmp.has(x(wild())) && is_a<mul>(tmp)) {
            ex rem = 1;
            for(auto item : tmp) {
                if(!item.has(x(wild()))) {
                    rem *= item;
                } else if(item.match(pow(wild(0), wild(1)))) {
                    kv.first.append(item.op(0));
                    kv.second.append(item.op(1));
                } else {
                    kv.first.append(item);
                    kv.second.append(1);
                }
            }
            kv.first.let_op(0) = rem;
        } else if(tmp.has(x(wild())) && tmp.match(pow(wild(0), wild(1)))) {
            kv.first.let_op(0) = 1;
            kv.first.append(tmp.op(0));
            kv.second.append(tmp.op(1));
        }
        
        if(false)
        for(int i=1; i<=nts; i++) {
            kv.first.append(x(xn+i-1));
            kv.second.append(eps);
        }
        
    }
    
    work.Normalizes();
    work.XReOrders();
    work.Normalizes();
    work.Scalelesses();
    work.RemoveDeltas();
    
    work.SDPrepares();
    work.EpsEpExpands();
    ostringstream ikey;
    ikey << SD_path << "/" << idx;
    work.CIPrepares(ikey.str().c_str());
    
    delete work.SecDec;
    delete work.Minimizer;
    
    if(work.IsZero) {
        ostringstream ifn;
        ifn << SD_path << "/" << idx << ".null";
        ExportNull(ifn.str().c_str());
        return;
    }
}

void Contour(int idx, numeric zz) {
    SD work;
    char *CFLAGS = getenv("SD_CFLAGS");
    work.CFLAGS = CFLAGS;
    work.ParallelSymbols = para_sym;
    work.Verbose = verb;
    work.Parameter[0] = zz;
    work.ParallelProcess = 0;
    
    ostringstream ikey;
    ikey << SD_path << "/" << idx;
    
    stringstream skey;
    skey << setprecision(2) << zz.to_double();
    const char *pkey = skey.str().c_str();
    
    if(work.Minimizer==NULL) work.Minimizer = new MinUit();
    work.Contours(ikey.str().c_str(), pkey);
    
    delete work.Minimizer;
}

ex Integrate(int idx, numeric zz, int ii = -1) {
    SD work;
    char *CFLAGS = getenv("SD_CFLAGS");
    work.CFLAGS = CFLAGS;
    work.ParallelSymbols = para_sym;
    work.Verbose = verb;
    work.Parameter[0] = zz;
    work.epN = epN;
    
    work.use_ErrMin = false;
    ErrMin::err_min = 1E-3;
    ErrMin::MaxRND = 20;
    ErrMin::hjRHO = 0.5;
    
    work.EpsAbs = 5E-6;
    work.RunPTS = 1000000;
    work.RunMAX = 10;
    work.LambdaSplit = 5;
    work.TryPTS = 1000000;
    work.CTryLeft = 3;
    work.CTryRight = 3;
    work.CTry = 1;
    work.CTryRightRatio = 5;
    work.ReIm = 1;
    
    ostringstream ikey;
    ikey << SD_path << "/" << idx;
    
    stringstream skey;
    skey << setprecision(2) << zz.to_double();
    const char *pkey = skey.str().c_str();
    
    if(work.Integrator==NULL) {
        auto intor = new HCubature();
        intor->use_last = true;
        work.Integrator = intor;
        intor->UseQ = useq;
    }
    work.Integrates(ikey.str().c_str(), pkey, ii);
    
    delete work.Integrator;
    return work.ResultError;
}

int main(int argc, char** argv) {

    //SD::debug = true;
    char* dc = getenv("SD_DLCLOSE");
    if(dc!=NULL && !strcmp(dc, "0")) SD::use_dlclose = false;
    
    const char* arg_a = NULL;
    const char* arg_n = NULL;
    const char* arg_i = NULL;
    const char* arg_z = "3/5";
    const char* arg_v = "1";
    const char* arg_e = "1E-2";
    bool is_o = false;
    const char* arg_o = "";
    
    for (int opt; (opt = getopt(argc, argv, "a:n:i:z:o:v:qe:")) != -1;) {
        switch (opt) {
            case 'a': arg_a = optarg; break;
            case 'n': arg_n = optarg; break;
            case 'i': arg_i = optarg; break;
            case 'z': arg_z = optarg; break;
            case 'o': arg_o = optarg; is_o = true; break;
            case 'v': arg_v = optarg; break;
            case 'e': arg_e = optarg; break;
            case 'q': useq = true; break;
            default:
                cout << "supported options: -a A -n N -i I -z Z -o O -v V -q -e E" << endl;
                cout << "A: can be p(prepare), c(contour), i(integrate), r(result), ar(analyse result)." << endl;
                cout << "N: only apply A on prefix N." << endl;
                cout << "I: only apply A on part I in prefix N." << endl;
                cout << "Z: z pararmeter in Fragmentation Function." << endl;
                cout << "O: output filename." << endl;
                cout << "V: verbose level." << endl;
                cout << "E: print when error > E." << endl;
                cout << "-q: use quadruple precision." << endl;
                exit(1);
        }
    }
    argc -= optind;
    argv += optind;
    
    ostringstream cmd;
    cmd << "mkdir -p " << SD_path;
    system(cmd.str().c_str());
    cmd.clear();
    cmd.str("");
    cmd << "ls " << cm_path << "|wc -l";
    auto owc = RunOS(cmd.str().c_str());
    int nmi;
    stringstream ss;
    ss << owc;
    ss >> nmi;
    
    int in = 0;
    if(arg_n != NULL) in = stoi(arg_n);
    int ii = -1;
    if(arg_n != NULL && arg_i != NULL) ii = stoi(arg_i);
    verb = stoi(arg_v);
    
    if(arg_a == NULL || !strcmp(arg_a, "p")) {
        vector<int> fid;
        for(int i=1; i<=nmi; i++) {
            if(in>0 && i!=in) continue;
            try {
                cout << "Preparing " << i << "/" << nmi << endl;
                Prepare(i);
            } catch(...) {
                fid.push_back(i);
            }
            cout << endl;
        }
        
        if(fid.size()>0) {
            cout << "Failed IDs: {";
            for(auto item : fid) cout << item << ", ";
            cout << "}" << endl;
            return 0;
        }
    }
    
    numeric zz(arg_z);
    numeric ee(arg_e);
    stringstream skey;
    skey << setprecision(2) << zz.to_double();
    const char *pkey = skey.str().c_str();

    if(arg_a == NULL || !strcmp(arg_a, "c") || !strcmp(arg_a, "ci")) {
        for(int i=1; i<=nmi; i++) {
            if(in>0 && i!=in) continue;
            stringstream ss;
            ss << SD_path << "/" << i << ".ci.gar";
            if(!file_exists(ss.str().c_str())) continue;
            cout << "Contouring @ z=" << arg_z << " - " << i << "/" << nmi << endl;
            Contour(i, zz);
            cout << endl;
        }
    }
    
    if(arg_a == NULL || !strcmp(arg_a, "i") || !strcmp(arg_a, "ci")) {
        ex fRes = 0;
        vector<int> nid;
        for(int i=1; i<=nmi; i++) {
            if(in>0 && i!=in) continue;
            stringstream ss;
            ss << SD_path << "/" << i << ".ci.gar";
            if(!file_exists(ss.str().c_str())) continue;
            cout << "Integrating @ z=" << arg_z << " - " << i << "/" << nmi << endl;
            auto res = Integrate(i, zz, ii);
            fRes += res;
            //if(!res.is_zero()) cout << VEResult(res) << endl;
            cout << endl;
            if(res.has(SD::NaN)) nid.push_back(i);
            
            auto err = res.subs(lst{VE(wild(1),wild(2))==wild(2)});
            if(abs(err.coeff(ep,0))>ee) {
                cout << "n = " << i << endl;
                cout << res.subs(I*wild()==0, subs_options::algebraic) << endl << endl;
            }
        }
        
        if(nid.size()>0) {
            cout << "NaN IDs: {";
            for(auto item : nid) cout << item << ", ";
            cout << "}" << endl;
            return 0;
        }
    }
    
    if(arg_a == NULL || !strcmp(arg_a, "r")) {
        ex fRes = 0;
        vector<int> nid;
        for(int i=1; i<=nmi; i++) {
            if(in>0 && i!=in) continue;
            stringstream ss;
            ss << SD_path << "/" << i << ".ci.gar";
            if(!file_exists(ss.str().c_str())) {
                ostringstream ifn;
                ifn << SD_path << "/" << i << ".null";
                if(!file_exists(ifn.str().c_str())) {
                    cout << RED << "File NOT Found: " << ifn.str() << RESET << endl;
                    assert(false);
                }
                continue;
            }
            
            ss.clear();
            ss.str("");
            ss << SD_path << "/" << i << "-" << pkey << ".res.gar";
            if(!file_exists(ss.str().c_str())) {
                cout << RED << "File NOT Found: " << ss.str() << RESET << endl;
                assert(false);
            }
            auto res = garResult(ss.str().c_str(), para_sym);
            fRes += res;
            if(res.has(SD::NaN)) nid.push_back(i);
            
            auto err = res.subs(lst{VE(wild(1),wild(2))==wild(2)});
            if(abs(err.coeff(ep,0))>ee) {
                cout << "n = " << i << endl;
                cout << res.subs(I*wild()==0, subs_options::algebraic) << endl << endl;
            }
        }
        
        if(nid.size()>0) {
            cout << "NaN IDs: {";
            for(auto item : nid) cout << item << ", ";
            cout << "}" << endl;
            return 0;
        }
        
        fRes = VESimplify(fRes, 0);
        cout << endl << "Final Result @ z=" << arg_z << endl;
        cout << "------" << endl;
        cout << fRes.subs(I*wild()==0, subs_options::algebraic) << endl;
        cout << "------" << endl;
        cout << VEResult(fRes.subs(I*wild()==0, subs_options::algebraic)) << endl << endl << endl;
        
        if(is_o) {
            ofstream ofs;
            ofs.open(arg_o, ios::out);
            if (!ofs) throw runtime_error("failed to open final out file!");
            ofs << fRes << endl;
            ofs.close();
        }
    }
    
    if(arg_a != NULL && !strcmp(arg_a, "ar")) {
        if(in<1) {
            cout << "n should be setted with -a ar option!" << endl;
            return 1;
        }
        stringstream ss;
        ss << SD_path << "/" << in << "-" << pkey << ".res.gar";
        archive ar;
        ifstream in(ss.str().c_str());
        in >> ar;
        in.close();
        auto c = ar.unarchive_ex(para_sym, "c");
        if(c!=19790923) {
            cout << "gar file: " << ss.str() << endl;
            cout << "c=" << c << ", different from 19790923!" << endl;
            assert(false);
        }
        auto relst = ex_to<lst>(ar.unarchive_ex(para_sym, "relst"));
        int max_index;
        ex max_re = -1;
        for(int i=0; i<relst.nops(); i++) {
            auto tmp = relst.op(i).subs(lst{VE(wild(1),wild(2))==wild(2)});
            tmp = tmp.subs(VE(wild(1), wild(2))==wild(2));
            tmp = tmp.expand().coeff(ep, epN).evalf();
            if(abs(tmp)>max_re) {
                max_re = abs(tmp);
                max_index = i;
            }
        }
        cout << "Max Index: " << max_index+1 << " / " << relst.nops() << endl;
        cout << VESimplify(relst.op(max_index).expand().coeff(ep, epN)) << endl;
    }

    return 0;
}
