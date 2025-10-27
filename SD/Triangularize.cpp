/**
 * @file
 * @brief Triangularize Function
 */

#include "SD.h"
extern "C" {
    #include <libqhull_r/qhull_ra.h>
}

namespace HepLib::SD {
    
    namespace {
        vector<vector<int>> QDelaunay(const ex & pts, int nshift=1) {
            vector<vector<int>> ret;
            int npts = pts.nops();
            int dim = pts.op(0).nops()-nshift;

            qhT qh_qh;
            qhT *qh= &qh_qh;
            
            coordT* cpts = (coordT*)malloc(sizeof(coordT) * npts * dim);
            for(int r=0; r<npts; r++) {
                for(int c=0; c<dim; c++) {
                    cpts[r*dim + c] = ex_to<numeric>(pts.op(r).op(c)).to_double();
                }
            }

            char opts[32];
            sprintf(opts, "qhull d Qbb QJ Qt i");
            FILE* dev_null = fopen("/dev/null", "w");
            int curlong, totlong;
            boolT ismalloc = false;
            
            QHULL_LIB_CHECK
            qh_zero(qh, dev_null);

            int exit_code = qh_new_qhull(qh, dim, npts, cpts, ismalloc, opts, NULL, dev_null);
            if(exit_code!=0) throw Error("QDelaunay: qhull exit code "+to_string(exit_code));
            fclose(dev_null);
            free(cpts);
                    
            facetT *facet;
            vertexT *vertex, **vertexp;
            FORALLfacets {
                if(!facet->upperdelaunay) {
                    vector<int> lvec;
                    FOREACHvertex_(facet->vertices) {
                        lvec.push_back(qh_pointid(qh, vertex->point));
                    }
                    ret.push_back(lvec);
                }
            }

            qh_freeqhull(qh, !qh_ALL);
            qh_memfreeshort(qh, &curlong, &totlong);
            if (curlong || totlong) {
                cout << "QDelaunay: Non-freed " << totlong << " bytes of long memory(" << curlong << " pieces)" << endl;
            }

            if(ret.size()<=0) {
                cerr << ErrColor << "QDelaunay: (ret.size()<=0)" << RESET << endl;
                throw Error("QDelaunay failed.");
            }

            return ret;
        }
    }

    /**
     * @brief to Triangularize the domain  with each xi from 0 to +infinity
     * @param fs_in is the list containing the liner equations w.r.t. xs_in
     * @param xs_in is a list of x's in Delta function
     * @param subs to replace the parameters with numerical values
     * @return vector of tranformation matrixes. for each matrix M, one has X = M.Y, with Y the new variables
     */
    vector<matrix> Triangularize(const lst & fs_in, const ex & xs_in, const lst & nsubs) {
        int nx = xs_in.nops();
        lst fs = fs_in;
        if(fs.nops()<1) {
            vector<matrix> mats;
            matrix mat(nx,nx);
            for(int i=0; i<nx; i++) mat(i,i) = 1;
            mats.push_back(mat);
            return mats;
        }
        lst xs;
        ex xsum = 0;
        for(int i=0; i<nx; i++) {
            symbol xi("x"+to_string(i));
            xs.append(xi);
            xsum += xi;
        }
        exmap x2x;
        for(int i=0; i<nx; i++) {
            fs.append(xs_in.op(i));
            x2x[xs_in.op(i)] = xs.op(i);
        }
        for(int i=0; i<fs.nops(); i++) {
            lst tls = {fs.op(i), -fs.op(i)};
            sort_lst(tls);
            fs.let_op(i) = tls.op(0);
        }
        fs.sort();
        fs.unique();
        fs = ex_to<lst>(subs(fs,x2x));
        int n = fs.nops();
        lst pts;
        int nmax = 1;
        int nx1 = nx-1;
        Combinations(fs.nops(), nx1, [&](const int is[]){
            lst eqs = lst{ xsum==nmax };
            for(int i=0; i<nx1; i++) eqs.append(fs.op(is[i])==0);
            auto sol = lsolve(eqs, xs);
            if(sol.nops()>0) {
                auto ps = subs(xs,ex_to<lst>(sol));
                bool ok = true;
                for(auto pi : ps) {
                    ex npi = pi.subs(nsubs);
                    if(!is_a<numeric>(npi)) {
                        cout << "Found non-numeric(after nsubs): " << pi.subs(nsubs) << endl;
                        cout << eqs << endl;
                        cout << sol << endl;
                        cout << endl << endl;
                    }
                    if(npi<0 || npi>nmax) { ok = false; break; }
                }
                if(ok) pts.append(ps);
            }
        });
        pts.sort();
        pts.unique();
  
        exset pts_set;
        pts_set.insert(pts);
        for(auto const & fi : fs) { // divide pts by the planes, i.e., equations
            auto cset = pts_set;
            pts_set.clear();
            for(auto const & si : cset) {
                lst lstn, lst0, lstp;
                for(auto const & item : si) {
                    exmap xmap;
                    for(int i=0; i<nx; i++) xmap[xs.op(i)] = item.op(i);
                    ex c = exnormal(fi.subs(xmap).subs(nsubs)); // nsubs
                    if(!is_a<numeric>(c)) {
                        cout << "c(after nsbus)=" << c << endl;
                        throw Error("TriSector: c is NOT a number.");
                    }
                    if(c.is_zero()) lst0.append(item);
                    else if(c>0) lstp.append(item);
                    else lstn.append(item);
                }
                if(lstn.nops()==0 || lstp.nops()==0) {
                    pts_set.insert(si);
                    continue;
                }
                for(auto const & item : lst0) {
                    lstp.append(item);
                    lstn.append(item);
                }
                pts_set.insert(lstp);
                pts_set.insert(lstn);
            }
        }

        vector<matrix> mats;
        for(auto const & pts : pts_set) { // redefine pts here
            vector<vector<int>> res;
            if(pts.nops()>nx) res = QDelaunay(pts.subs(nsubs)); // nsubs
            else {
                vector<int> n_vec(nx);
                for(int i=0; i<nx; i++) n_vec[i] = i;
                res.push_back(n_vec);
            }
            for(auto item : res) {
                exvector c_vec(nx);
                exmap c2i;
                for(int i=0; i<nx; i++) {
                    symbol ci("c"+to_string(i));
                    c_vec[i] = ci;
                    c2i[ci] = 1;
                }
                
                matrix mat(nx,nx);
                for(int i=0; i<nx; i++) {
                    exvector pt_vec(nx);
                    for(int j=0; j<nx; j++) pt_vec[j] = pts.op(item[(i+j)%nx]);
                    lst eqs;
                    for(int j=1; j<nx; j++) { // no eq0, pt_vec[0] used later
                        ex eq = 0;
                        for(int k=0; k<nx; k++) eq += c_vec[k] * pt_vec[j][k];
                        eqs.append(eq==0);
                    }
                    lst c_lst = vec2lst(c_vec);
                    ex sol = lsolve(eqs, c_lst);
                    sol = subs(c_lst, ex_to<lst>(sol)).subs(c2i); // one ci remains, just set it to 1
                    ex sgn = 0;
                    for(int j=0; j<nx; j++) sgn += sol.op(j) * pt_vec[0].op(j);
                    sgn = sgn.subs(nsubs);
                    if(!is_a<numeric>(sgn)) {
                        cout << "sgn(after nsubs) = " << sgn << endl;
                        throw Error("TriSector: sgn is NOT a number!");
                    }
                    if(sgn<0) sgn = -1;
                    else sgn = 1;
                    for(int j=0; j<nx; j++) mat(i,j) = sgn * sol.op(j);
                }
                mats.push_back(mat.inverse());
            }
        }

        return mats;
    }
    
    
    /**
     * @brief to Triangularize the domain  with each xi from 0 to +infinity, and update FunExp directly
     * @param FunExp the same as FunExp in SecDec class, will be updated
     * @param fs is the list containing the liner equations w.r.t. xs_in
     * @param xs is a list of x's in Delta function
     * @param subs to replace the parameters with numerical values
     * @return vector of tranformation matrixes. for each matrix M, one has X = M.Y, with Y the new variables
     */
    void Triangularize(exvector & FunExp, const lst & fs, const ex & xs, const lst & nsubs) {
        auto mats = Triangularize(fs, xs, nsubs);
        ex xsum = 0;
        for(auto xi : xs) xsum += xi;
        auto fun_exp_lst = FunExp;
        FunExp.clear();
        for(auto fun_exp : fun_exp_lst) {
            for(auto const & mat : mats) {
                lst funs = ex_to<lst>(fun_exp.op(0));
                lst exps = ex_to<lst>(fun_exp.op(1));
                lst dlts = ex_to<lst>(fun_exp.op(2));
                exmap x2y;
                for(int r=0; r<xs.nops(); r++) {
                    ex xr = 0;
                    for(int c=0; c<xs.nops(); c++) xr += mat(r,c)*y(c);
                    x2y[xs.op(r)] = xr;
                }
                funs = ex_to<lst>(subs(funs,x2y));
                ex det = mat.determinant();
                ex ndet = det.subs(nsubs);
                if(!is_a<numeric>(ndet)) {
                    cout << "ndet(after nsbus)=" << ndet << endl;
                    throw Error("Triangularize: ndet is NOT a number.");
                }
                if(ndet<0) det = -det;
                ex ysum = xsum.subs(x2y);
                exmap y2x;
                for(int i=0; i<xs.nops(); i++) {
                    ex c = ysum.coeff(y(i));
                    det /= c;
                    y2x[y(i)] = xs.op(i)/c;
                }
                funs = ex_to<lst>(subs(funs,y2x));
                funs.append(det);
                exps.append(1);
                FunExp.push_back(lst{funs, exps, dlts});
            }
        }
    }
    

}
