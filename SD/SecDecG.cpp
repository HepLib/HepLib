#include "SD.h"

extern "C" {
    #include <libqhull/qhull_a.h>
}

namespace HepLib::SD {

vector<vector<int>> SecDecG::RunQHull(const matrix &pts) {

    vector<vector<int>> ret;
    
    int npts = pts.rows();
    int dim = pts.cols();
    ex imax = -1;
    for(int r=0; r<npts; r++) {
        for(int c=0; c<dim; c++) {
            if(imax<abs(pts(r,c))) imax = abs(pts(r,c));
        }
    }
    
    coordT cpts[npts * dim];
    for(int r=0; r<npts; r++) {
        for(int c=0; c<dim; c++) {
            if(imax<1000) cpts[r*dim + c] = ex_to<numeric>(pts(r,c)).to_int();
            else cpts[r*dim + c] = ex_to<numeric>(pts(r,c)).to_double();
        }
    }
    
    char opts[32];
    sprintf(opts, "qhull Fv");
    FILE* dev_null = fopen("/dev/null", "w");
    int curlong, totlong;
    int exit_code = qh_new_qhull(dim, npts, cpts, 0, opts, NULL, dev_null);
    if(exit_code) {
        qh_freeqhull(!qh_ALL);
        qh_memfreeshort(&curlong, &totlong);
        char opts[32];
        sprintf(opts, "qhull QbB Fv");
        exit_code = qh_new_qhull(dim, npts, cpts, 0, opts, NULL, dev_null);
        if(exit_code) {
            cout << Color_Error << "qhull return code : " << exit_code << RESET << endl;
            cout << "input for qhull Fv:" << endl;
            cout << dim << endl;
            cout << npts << endl;
            for(int r=0; r<npts; r++) {
                for(int c=0; c<dim; c++) cout << pts(r,c) << " ";
                cout << endl;
            }
            fclose(dev_null);
            qh_freeqhull(!qh_ALL);
            qh_memfreeshort(&curlong, &totlong);
            throw "Qhull Error!";
        }
    }
    fclose(dev_null);
    
    facetT *facet;
    vertexT *vertex, **vertexp;
    FORALLfacets {
        vector<int> lvec;
        FOREACHvertex_(facet->vertices) {
            lvec.push_back(qh_pointid(vertex->point));
        }
        ret.push_back(lvec);
    }

    qh_freeqhull(!qh_ALL);
    qh_memfreeshort(&curlong, &totlong);
    if (curlong || totlong) {
        cout << "Qhull: Non-freed " << totlong << " bytes of long memory(" << curlong << " pieces)" << endl;
    }
    //qh_freeqhull(qh_ALL);
    
    if(ret.size()<=0) {
        cerr << Color_Error << "RunQHull: (ret.size()<=0)" << RESET << endl;
        exit(1);
    }
    return ret;
}

vector<matrix> SecDecG::ZeroFaces(const matrix &pts) {
    if(pts.rows()<=0) {
        cerr << Color_Error << "ZeroFaces: (pts.rows()<=0)" << RESET << endl;
        exit(1);
    }
    auto zri = MatHelper::zero_row_index(pts);
    if(zri.size() <= 0) {
        cerr << Color_Error << "ZeroFaces: (zri.size() <= 0)" << RESET << endl;
        exit(1);
    }
    int zpos = zri[0];
    auto QH = RunQHull(pts);
    vector<matrix> ret;
    for(auto vi : QH) {
        bool iszero = find(vi.begin(), vi.end(), zpos) != vi.end();
        if(iszero) {
            matrix zmat(vi.size(), pts.cols());
            for(int r=0; r<zmat.rows(); r++) {
                for(int c=0; c<zmat.cols(); c++) {
                    zmat(r, c) = pts(vi[r], c);
                }
            }
            ret.push_back(zmat);
        }
    }
    
    return ret;
}

matrix SecDecG::NormalVectors(const vector<matrix> &zfs) {
    int ncols = zfs[0].cols();
    int nrows = zfs.size();
    matrix ret(nrows, ncols);
    for(int ii=0; ii<nrows; ii++) {
    
        matrix dir;
        if(ii==0) {
            dir = MatHelper::sub(MatHelper::remove_zero_rows(zfs[1]), 0, 1, 0, ncols);
        } else {
            dir = MatHelper::sub(MatHelper::remove_zero_rows(zfs[0]), 0, 1, 0, ncols);
        }
        
        matrix tmat = MatHelper::remove_zero_rows(zfs[ii]);
        if(tmat.rows() >= zfs[ii].rows()) {
            cerr << Color_Error << "NormalVectors: (tmat.rows() >= zfs[ii].rows())" << RESET << endl;
            exit(1);
        }
        
        for(int ti=0; ti<tmat.rows(); ti++) {
            int pos = -1;
            for(int tc=0; tc<tmat.cols(); tc++) {
                if(tmat(ti, tc) !=0 ) {
                    pos = tc;
                    break;
                }
            }
            if(pos == -1) continue;
            
            for(int tj=ti+1; tj<tmat.rows(); tj++) {
                if(tmat(tj,pos)!=0) {
                    matrix matj = MatHelper::sub(tmat, tj,1, 0,tmat.cols());
                    for(int tc=0; tc<tmat.cols(); tc++) {
                        tmat(tj, tc)=tmat(ti,pos)*matj(0,tc)-matj(0,pos)*tmat(ti,tc);
                    }
                }
            }
        }
        
        tmat = MatHelper::remove_zero_rows(tmat);
        if(tmat.rows()!=tmat.cols()-1) {
            cerr << Color_Error << "NormalVectors: (tmat.rows()!=tmat.cols()-1)" << RESET << endl;
            exit(1);
        }
        matrix tmat2(tmat.cols(), tmat.cols());
        for(int r=0; r<tmat.rows(); r++) {
            for(int c=0; c<tmat.cols(); c++) {
                tmat2(r,c) = tmat(r,c);
            }
        }
        
        for(int c=0; c<tmat.cols(); c++) {
            tmat2(tmat.cols()-1,c) = 1;
        }
        tmat = tmat2;
        
        auto det = tmat.determinant();
        matrix vec(tmat.cols(), 1);
        vec(tmat.cols()-1, 0) = det;
        vec = tmat.inverse().mul(vec);
        det = 0;
        for(int r=0; r<vec.rows(); r++) {
            det = gcd(det, vec(r,0));
        }
        
        for(int r=0; r<vec.rows(); r++) {
            vec(r,0) = vec(r,0)/det;
        }
        
        if(dir.mul(vec)(0,0) < 0) vec=vec.mul_scalar(-1);
        
        for(int r=0; r<vec.rows(); r++) {
            ret(ii,r) = vec(r,0);
        }
    }
    return ret;
}

matrix SecDecG::DualCone(const matrix &pts) {
    auto zfs = ZeroFaces(pts);
    if(zfs.size()<1 || zfs.size()<pts.cols()) return matrix(0,0);
    auto nvs = NormalVectors(zfs);
    ex sum_vec[nvs.rows()];
    ex sum_lcm = 1;
    for(int r=0; r<nvs.rows(); r++) {
        sum_vec[r] = 0;
        for(int c=0; c<nvs.cols(); c++) sum_vec[r] += nvs(r,c);
        sum_lcm = lcm(sum_lcm, sum_vec[r]);
    }
    
    matrix ret(nvs.rows(), nvs.cols());
    for(int r=0; r<nvs.rows(); r++) {
        for(int c=0; c<nvs.cols(); c++) ret(r, c) = sum_lcm/sum_vec[r]*nvs(r,c);
    }
    return ret;
}

vector<vector<int>> SecDecG::QHull(const matrix &dc, int dim) {
    if(dim >= dc.cols()) return RunQHull(dc);
    
    matrix mdc(dc.rows()-1, dc.cols());
    for(int r=0; r<mdc.rows(); r++) {
        for(int c=0; c<mdc.cols(); c++) mdc(r, c) = dc(r+1, c) - dc(0, c);
    }
    
    int n = mdc.cols();
    mdc = mdc.transpose();
    for(int i=0; i<n; i++) {
        if(MatHelper::is_zero_row(mdc, i)) continue;
        int pos = -1;
        for(int c=0; c<mdc.cols(); c++) {
            if(mdc(i,c)!=0) {
                pos = c;
                break;
            }
        }
        
        for(int j=i+1; j<n; j++) {
            if(mdc(j,pos)!=0) {
                matrix matj = MatHelper::sub(mdc, j,1, 0,mdc.cols());
                for(int c=0; c<mdc.cols(); c++) {
                    mdc(j, c)=matj(0,c) - mdc(i,c) * matj(0, pos)/mdc(i, pos);
                }
            }
        }
    }
    
    mdc = MatHelper::remove_zero_rows(mdc);
    mdc = mdc.transpose();

    matrix tmp(mdc.rows()+1, mdc.cols());
    for(int r=0; r<mdc.rows(); r++) {
        for(int c=0; c<mdc.cols(); c++) {
            tmp(r+1,c) = mdc(r,c);
        }
    }
    for(int c=0; c<mdc.cols(); c++) {
        tmp(0,c) = 0;
    }
    mdc = tmp;

    ex all_lcm = 1;
    for(int r=0; r<mdc.rows(); r++) {
        for(int c=0; c<mdc.cols(); c++) {
            all_lcm = lcm(all_lcm, normal(mdc(r,c)).denom());
        }
    }
    mdc = mdc.mul_scalar(all_lcm);
    return RunQHull(mdc);
}

vector<matrix> SecDecG::Simplexify(const matrix &dc, int dim) {
    if(dc.rows()-dim<=0) {
        cerr << Color_Error << "Simplexify: (dc.rows()-dim<=0)" << RESET << endl;
        exit(1);
    }
    vector<matrix> ret;
    if(dc.rows() == dim+1) {
        ret.push_back(dc);
        return ret;
    }
    
    auto tmp = QHull(dc, dim);
    vector<matrix> vmat;
    for(auto vi : tmp) {
        bool found = find(vi.begin(), vi.end(), 0) != vi.end();
        if(found) continue;
        matrix mat(vi.size(), dc.cols());
        for(int r=0; r<mat.rows(); r++) {
            for(int c=0; c<mat.cols(); c++) {
                mat(r,c) = dc(vi[r], c);
            }
        }
        vmat.push_back(mat);
    }
    
    for(auto mat : vmat) {
        auto tmp = SimplexifyR(mat, dim-1);
        for(auto it : tmp) ret.push_back(it);
    }
    
    for(int i=0; i<ret.size(); i++) {
        auto tmp = ret[i];
        matrix mat(tmp.rows()+1, tmp.cols());
        for(int r=1; r<mat.rows(); r++) {
            for(int c=0; c<mat.cols(); c++) {
                mat(r,c) = tmp(r-1, c);
            }
        }
        for(int c=0; c<mat.cols(); c++) {
            mat(0,c) = dc(0, c);
        }
        ret[i] = mat;
    }
    vmat.clear();
    vmat.shrink_to_fit();
    return ret;
}

vector<matrix> SecDecG::SimplexifyR(const matrix &dc, int dim) {
    int min = -1;
    int i = 0;
    auto tmat = dc;
    vector<matrix> ret;

    while (i<tmat.rows() && (min<0 || min>2)) {
        auto tmp_ret = Simplexify(tmat, dim);
        if(tmp_ret.size()<min || min<0) {
            min = tmp_ret.size();
            ret = tmp_ret;
        }
        i++;
        matrix tmp(tmat.rows(), tmat.cols());
        for(int r=1; r<tmat.rows(); r++) {
            for(int c=0; c<tmat.cols(); c++) {
                tmp(r,c) = tmat(r-1, c);
            }
        }
        for(int c=0; c<tmat.cols(); c++) {
            tmp(0,c) = tmat(tmat.rows()-1, c);
        }
        tmat = tmp;
    }
    return ret;
}

vector<matrix> SecDecG::SimplexCones(matrix pts) {
    auto ds = DualCone(pts);
    if(ds.rows() == 0) return vector<matrix>();
    
    if(ds.rows() < ds.cols()) {
        cerr << Color_Error << "SimplexCones: (ds.rows() < ds.cols())" << RESET << endl;
        exit(1);
    }
    if(ds.rank()<ds.cols()) return vector<matrix>();
    return SimplexifyR(ds, ds.cols()-1);
}

// return a replacement/transformation, using x(-1) as key for determinant
vector<exmap> SecDecG::x2y(const ex &xpol) {
    if(xpol.has(y(w))) {
        cerr << Color_Error << "SecDecG::x2y, y(w) found @ " << xpol << RESET << endl;
        exit(1);
    }
    auto xs = get_xy_from(xpol);
    int nx = xs.size();
    auto pol = mma_collect(xpol, x(w), true);
    pol = pol.subs(coCF(w)==1);
    int np = is_a<add>(pol) ? pol.nops() : 1;
    
    if(nx<2 || np<2) {
        vector<exmap> vmap;
        exmap nmap;
        nmap[x(-1)] = 1;
        int i=0;
        for(auto xi : xs) nmap[xi] = y(i++);
        vmap.push_back(nmap);
        return vmap;
    }
    
    vector<ex> vpols;
    for(auto item : pol) vpols.push_back(item);
    sort(vpols.begin(), vpols.end(), [&](const auto &ain, const auto &bin){
        //auto a=ain; auto b=bin; // < < <
        auto a=bin; auto b=ain; // > > >
        int tai=0, tbi=0;
        int ab = 0;
        for(auto xi : xs) {
            tai += a.degree(xi);
            tbi += b.degree(xi);
            if(ab==0 && tai!=tbi) ab = (tai<tbi) ? 1 : -1;
        }
        if(tai!=tbi) return (tai<tbi);
        if(ab!=0) return (ab>0);
        return ex_is_less()(a,b);
    });
        
    matrix deg_mat(np, nx);
    for(int r=0; r<np; r++) {
        auto tmp = vpols[r];
        for(int c=0; c<nx; c++) {
            deg_mat(r, c) = tmp.degree(xs[c]);
        }
    }

    vector<matrix> vmat;
    for(int r=0; r<deg_mat.rows(); r++) {
        matrix tmp(deg_mat.rows()+xs.size(), deg_mat.cols());
        for(int rr=0; rr<deg_mat.rows(); rr++) {
            for(int c=0; c<deg_mat.cols(); c++) tmp(rr,c) = deg_mat(rr,c) - deg_mat(r, c);
        }
        
        for(int rr=0; rr<tmp.cols(); rr++) {
            tmp(rr+deg_mat.rows(),rr) = 1;
        }
        auto sc = SimplexCones(tmp);
        for(auto isc : sc) vmat.push_back(isc);
    }

    vector<map<ex,ex,ex_is_less>> ret;
    for(int vi=0; vi<vmat.size(); vi++) {
        matrix &tmp = vmat[vi];
        for(int r=0; r<tmp.rows(); r++) {
            ex row_gcd = 0;
            for(int c=0; c<tmp.cols(); c++) row_gcd = gcd(row_gcd, tmp(r, c));
            for(int c=0; c<tmp.cols(); c++) tmp(r, c) = tmp(r, c)/row_gcd;
        }
        if(tmp.determinant()<0) {
            ex row0[tmp.cols()];
            for(int c=0; c<tmp.cols(); c++) {
                row0[c] = tmp(0, c);
                tmp(0, c) = tmp(1, c);
            }
            for(int c=0; c<tmp.cols(); c++) tmp(1,c) = row0[c];
        }
        
        tmp = tmp.transpose();
        matrix Dxy(nx, nx);
        exmap transmap;
        for(int n=0; n<nx; n++) {
            ex tt = 1;
            for(int m=0; m<nx; m++) {
                tt = tt*pow(y(m), tmp(n,m));
            }
            transmap[xs[n]] = tt;
            for(int m=0; m<nx; m++) Dxy(n, m) = mma_diff(tt, y(m), 1, false);
        }
        transmap[x(-1)] = Dxy.determinant();
        ret.push_back(transmap);
    }
    
    vmat.clear();
    vmat.shrink_to_fit();
    return ret;
    
}


}
