#include "SD.h"
#include <cmath>

namespace HepLib {

bool SD::KillPowersWithDelta(lst fe, int kpi) {
    if(fe.op(0).op(fe.op(0).nops()-1)==WF(1) && fe.op(1).op(fe.op(1).nops()-1).is_zero()) {
        FunExp.push_back(fe);
        return false;
    }
    ex ft = fe.op(0).op(1);
    ft = Factor(ft);
    lst fts;
    if(is_a<mul>(ft)) {
        for(auto item : ft) fts.append(item);
    } else {
        fts.append(ft);
    }
    ex eqn;
    
    //-----------------------------------------------------
    bool ok2 = true;
    for(auto ftitem : fts) {
        auto xs = get_x_from(ftitem);
        for(int i=0; i<xs.size(); i++) {
            for(int j=i+1; j<xs.size(); j++) {
                bool has_xij = false;
                for(auto delta : ex_to<lst>(fe.op(2))) {
                    if(delta.has(xs[i]) && delta.has(xs[j])) {
                        has_xij = true;
                        break;
                    }
                }
                if(!has_xij) continue;
                
                symbol xi("xi"), xj("xj");
                auto ftij = ftitem.subs(lst{xs[i]==xi, xs[j]==xj});
                auto xs2 = get_x_from(ftij);
                auto ftt = Factor(ftij.subs(x(wild())==0));
                lst fts2;
                if(is_a<mul>(ftt)) {
                    for(auto item : ftt) fts2.append(item);
                } else {
                    fts2.append(ftt);
                }
                
                for(auto item : fts2) {
                    if(!item.has(xi) || !item.has(xj)) continue;
                    if(item.has(xi) && item.has(xj)) {
                        if(item.match(pow(wild(1),wild(2)))) eqn = item.op(0).expand();
                        else eqn = item;
                        if(eqn.degree(xi)==1 && eqn.degree(xj)==1) {
                            ex ci = eqn.coeff(xi);
                            ex cj = eqn.coeff(xj);
                            if((ci*xi+cj*xj-eqn).is_zero() && is_a<numeric>(ci*cj) && (ci*cj)<0) {
                                eqn = eqn.subs(lst{xi==xs[i], xj==xs[j]});
                                ok2 = false;
                                goto OK2;
                            }
                        }
                    }
                }
            }
        }
    }

    OK2:
    if(!ok2) {
        auto xij = get_x_from(eqn);
        ex xi = xij[0];
        ex xj = xij[1];
        ex ci = eqn.coeff(xi);
        ex cj = eqn.coeff(xj);
        
        // handle eqn==ci xi - cj xj
        if((ci*xi+cj*xj-eqn).is_zero() && is_a<numeric>(ci * cj) && (ci*cj)<0) {
            if(Verbose>10) cout << "  \\--" << WHITE << "KillPowers-δ ["<<kpi<<"]: "  << eqn << RESET << endl;
            ci = abs(ci);
            cj = abs(cj);
            symbol yi,yj;
            // Part I: ci xi-cj xj>0, i.e., xi>cj/ci xj
            auto f1 = fe.op(0);
            auto e1 = fe.op(1);
            ex c1 = cj/ci;
            for(int i=0; i<f1.nops(); i++) {
                f1.let_op(i) = f1.op(i).subs(lst{xi==c1*yj/(1+c1)+yi,xj==yj/(1+c1)}).subs(lst{yi==xi,yj==xj});
            }
            f1.let_op(0) = f1.op(0)/(1+c1); // Jaccobi
            FunExp.push_back(lst{f1,e1,fe.op(2)});
            // Part II: ci xi-cj xj<0, i.e., i.e., xj>ci/cj xi
            auto f2 = fe.op(0);
            auto e2 = fe.op(1);
            ex c2 = ci/cj;
            for(int i=0; i<f2.nops(); i++) {
                f2.let_op(i) = f2.op(i).subs(lst{xj==c2*yi/(1+c2)+yj,xi==yi/(1+c2)}).subs(lst{yi==xi,yj==xj});
            }
            f2.let_op(0) = f2.op(0)/(1+c2); // Jaccobi
            FunExp.push_back(lst{f2,e2,fe.op(2)});
        } else {
            assert(false);
        }
        return true;
    }
            
    // when ok2 is true
    let_op_append(fe, 0, WF(1));
    let_op_append(fe, 1, 0);
    FunExp.push_back(fe);
    return false;
}

bool SD::KillPowersWithoutDelta(lst fe, int kpi, int bits) {
    if(fe.op(0).op(fe.op(0).nops()-1)==WF(1) && fe.op(1).op(fe.op(1).nops()-1).is_zero()) {
        FunExp.push_back(fe);
        return false;
    }
    ex ft = fe.op(0).op(1);
    ft = Factor(ft);
    lst fts;
    if(is_a<mul>(ft)) {
        for(auto item : ft) fts.append(item);
    } else {
        fts.append(ft);
    }
    ex eqn;
    
    //-----------------------------------------------------
    bool ok2 = true;
    if(((bits/2)%2)==1) {
    for(auto ftitem : fts) {
        auto xs = get_x_from(ftitem);
        for(int i=0; i<xs.size(); i++) {
            for(int j=i+1; j<xs.size(); j++) {
                symbol xi("xi"), xj("xj");
                auto ftij = ftitem.subs(lst{xs[i]==xi, xs[j]==xj});
                auto xs2 = get_x_from(ftij);
                for(int nn=0; nn<std::pow(2,xs2.size()); nn++) {
                    int tnn = nn;
                    lst xsubs;
                    for(int ni=0; ni<xs2.size(); ni++) {
                        if(tnn%2==1) xsubs.append(xs2[ni]==1);
                        else xsubs.append(xs2[ni]==0);
                        tnn /= 2;
                    }
                    auto ftt = Factor(ftij.subs(xsubs));
                    lst fts2;
                    if(is_a<mul>(ftt)) {
                        for(auto item : ftt) fts2.append(item);
                    } else {
                        fts2.append(ftt);
                    }
                    
                    int NN = 100;
                    for(auto item : fts2) {
                        if(!item.has(xi) || !item.has(xj)) continue;
                        if(item.match(pow(wild(1),wild(2))) && (item.has(xi) || item.has(xj))) {
                            eqn = item.op(0);
                            auto t1 = eqn.subs(lst{xi==1/ex(11), xj==1/ex(19)});
                            if(t1.is_zero()) t1 = eqn.subs(lst{xi==1/ex(3), xj==1/ex(23)});
                            if(t1.is_zero()) t1 = eqn.subs(lst{xi==1/ex(13), xj==1/ex(37)});
                            
                            bool ook = true;
                            for(int ni=0; ni<=NN; ni++) {
                                for(int nj=0; nj<=NN; nj++) {
                                    auto t2 = eqn.subs(lst{xi==ni/ex(NN), xj==nj/ex(NN)});
                                    if(t1*t2 < 0) {
                                        ook = false;
                                        break;
                                    }
                                }
                                if(!ook) break;
                            }
                            if(ook) continue;
                            
                            if(eqn.degree(xi)>1 || eqn.degree(xj)>1 || eqn.coeff(xi).has(xj)) {
                                cout << MAGENTA << "Warning: Not handled with eqn1=" << eqn << RESET << endl;
                                continue;
                            }

                            if(!ook) {
                                eqn = eqn.subs(lst{xi==xs[i], xj==xs[j]});
                                ok2 = false;
                                goto OK2;
                            }
                        } else if(item.degree(xi)==1 && item.degree(xj)==1) {
                            ex ci = item.coeff(xi);
                            ex cj = item.coeff(xj);
                            if((ci*xi+cj*xj-item).is_zero() && is_a<numeric>(ci*cj) && (ci*cj)<0) {
                                eqn = item.subs(lst{xi==xs[i], xj==xs[j]});
                                ok2 = false;
                                goto OK2;
                            }
                        }
                    }
                }
            }
        }
    }}

    OK2:
    if(!ok2) {
        auto xij = get_x_from(eqn);
        ex xi = xij[0];
        ex xj = xij[1];
        
        if(false)
        if((eqn-(xi+xj-1)).is_zero() || (eqn+(xi+xj-1)).is_zero()) {
            symbol xx;
            auto f1 = fe.op(0);
            auto e1 = fe.op(1);
            for(int i=0; i<f1.nops(); i++) f1.let_op(i) = f1.op(i).subs(xi==1-xx).subs(xx==xi);
            FunExp.push_back(lst{f1,e1});
            return true;
        }
        
        ex ci = eqn.coeff(xi);
        ex cj = eqn.coeff(xj);
        
        // handle eqn==ci xi - cj xj
        if((ci*xi+cj*xj-eqn).is_zero() && is_a<numeric>(ci * cj) && (ci*cj)<0) {
            if(Verbose>10) cout << "  \\--" << WHITE << "KillPowers ["<<kpi<<"]: "  << eqn << RESET << endl;
            ci = abs(ci);
            cj = abs(cj);
            if(ci==cj) {
                symbol xx;
                // Part I: xi<xj
                auto f1 = ex_to<lst>(fe.op(0));
                auto e1 = ex_to<lst>(fe.op(1));
                for(int i=0; i<f1.nops(); i++) f1.let_op(i) = f1.op(i).subs(xi==xx*xj).subs(xx==xi);
                f1.append(xj); // Jaccobi
                e1.append(1);
                FunExp.push_back(lst{f1,e1});
                // Part II: xj<xi
                auto f2 = ex_to<lst>(fe.op(0));
                auto e2 = ex_to<lst>(fe.op(1));
                for(int i=0; i<f2.nops(); i++) f2.let_op(i) = f2.op(i).subs(xj==xx*xi).subs(xx==xj);
                f2.append(xi); // Jaccobi
                e2.append(1);
                FunExp.push_back(lst{f2,e2});
            } else {
                // we set c1 > c2 always
                ex c1 = ci;
                ex x1 = xi;
                ex c2 = cj;
                ex x2 = xj;
                if(c1 < c2) {
                    c1 = cj;
                    x1 = xj;
                    c2 = ci;
                    x2 = xi;
                }
                symbol xx;
                // Part I: 0<x2<1, c2/c1<x1<1
                auto f1 = ex_to<lst>(fe.op(0));
                auto e1 = ex_to<lst>(fe.op(1));
                for(int i=0; i<f1.nops(); i++) f1.let_op(i) = f1.op(i).subs(x1==xx*(c1-c2)/c1+c2/c1).subs(xx==x1);
                f1.let_op(0) = f1.op(0)*(c1-c2)/c1; // Jaccobi
                FunExp.push_back(lst{f1,e1});
                // Part II: x1>c2/c1 x2, i.e., x2<c1/c2 x1, 0<x1<c2/c1
                auto f2 = ex_to<lst>(fe.op(0));
                auto e2 = ex_to<lst>(fe.op(1));
                for(int i=0; i<f2.nops(); i++) f2.let_op(i) = f2.op(i).subs(x1==xx*c2/c1).subs(xx==x1);
                f2.let_op(0) = f2.op(0)*c2/c1;
                // now x2<x1, 0<x1<1
                for(int i=0; i<f2.nops(); i++) f2.let_op(i) = f2.op(i).subs(x2==xx*x1).subs(xx==x2);
                f2.append(x1); // Jaccobi
                e2.append(1);
                FunExp.push_back(lst{f2,e2});
                // Part III: x1<c2/c1 x2
                auto f3 = ex_to<lst>(fe.op(0));
                auto e3 = ex_to<lst>(fe.op(1));
                for(int i=0; i<f3.nops(); i++) f3.let_op(i) = f3.op(i).subs(x1==xx*x2*c2/c1).subs(xx==x1);
                f3.let_op(0) = f3.op(0)*c2/c1; // Jaccobi
                f3.append(x2);
                e3.append(1);
                FunExp.push_back(lst{f3,e3});
            }
        } else if( (eqn-(xi+xj-1)).is_zero() || (eqn+(xi+xj-1)).is_zero() ) {
            if(Verbose>10) cout << "  \\--" << WHITE << "KillPowers ["<<kpi<<"]: "  << eqn << RESET << endl;
            symbol xx, yy, zz;
            // Part I: xi+xj-1>0
            auto f1 = ex_to<lst>(fe.op(0));
            auto e1 = ex_to<lst>(fe.op(1));
            for(int i=0; i<f1.nops(); i++) f1.let_op(i) = f1.op(i).subs(xj==xx+1-xi).subs(xx==xj);
            // now 0<xi<1, 0<xj<xi
            for(int i=0; i<f1.nops(); i++) f1.let_op(i) = f1.op(i).subs(xj==xx*xi).subs(xx==xj);
            f1.append(xi); // Jaccobi
            e1.append(1);
            FunExp.push_back(lst{f1,e1});
            // Part IIa: 1-xi-xj>0, (a): xj<xi
            auto f2 = ex_to<lst>(fe.op(0));
            auto e2 = ex_to<lst>(fe.op(1));
            for(int i=0; i<f2.nops(); i++) f2.let_op(i) = f2.op(i).subs(lst{xi==(1+zz)*yy/2,xj==(1-zz)*yy/2}).subs(lst{yy==xi, zz==xj});
            f2.append(xi/2); // Jaccobi
            e2.append(1);
            FunExp.push_back(lst{f2,e2});
            // Part IIb: 1-xi-xj>0, (a): xj>xi
            auto f3 = ex_to<lst>(fe.op(0));
            auto e3 = ex_to<lst>(fe.op(1));
            for(int i=0; i<f3.nops(); i++) f3.let_op(i) = f3.op(i).subs(lst{xj==(1+zz)*yy/2,xi==(1-zz)*yy/2}).subs(lst{yy==xi, zz==xj});
            f3.append(xi/2); // Jaccobi
            e3.append(1);
            FunExp.push_back(lst{f3,e3});
        } else {
            let_op_append(fe, 0, WF(1));
            let_op_append(fe, 1, 0);
            FunExp.push_back(fe);
            cout << MAGENTA << "Warning: Not handled with eqn2=" << eqn << RESET << endl;
        }
        return true;
    }
    
    //-----------------------------------------------------
    bool ok1 = true;
    if((bits%2)==1) {
    for(auto ftitem : fts) {
        auto xs = get_x_from(ftitem);
        for(int i=0; i<xs.size(); i++) {
            symbol xi("xi");
            auto fti = ftitem.subs(lst{xs[i]==xi});
            auto xs2 = get_x_from(fti);
            for(int nn=0; nn<std::pow(2,xs2.size()); nn++) {
                int tnn = nn;
                lst xsubs;
                for(int ni=0; ni<xs2.size(); ni++) {
                    if(tnn%2==1) xsubs.append(xs2[ni]==1);
                    else xsubs.append(xs2[ni]==0);
                    tnn /= 2;
                }
                auto ftt = Factor(fti.subs(xsubs));
                lst fts2;
                if(is_a<mul>(ftt)) {
                    for(auto item : ftt) fts2.append(item);
                } else {
                    fts2.append(ftt);
                }
                
                int NN = 100;
                for(auto item : fts2) {
                    if(item.match(pow(wild(1),wild(2))) && item.has(xi)) {
                        eqn = item.op(0);
                        auto t1 = eqn.subs(lst{xi==1/ex(11)});
                        if(t1.is_zero()) t1 = eqn.subs(lst{xi==1/ex(3)});
                        if(t1.is_zero()) t1 = eqn.subs(lst{xi==1/ex(13)});
                        
                        bool ook = true;
                        for(int ni=0; ni<=NN; ni++) {
                            auto t2 = eqn.subs(lst{xi==ni/ex(NN)});
                            if(t1*t2 < 0) {
                                ook = false;
                                break;
                            }
                        }
                        if(ook) continue;
                        
                        if(eqn.degree(xi)>1) {
                            cout << MAGENTA << "Warning: Not handled with eqn3=" << eqn << RESET << endl;
                            continue;
                        }
                        
                        if(!ook) {
                            eqn = eqn.subs(lst{xi==xs[i]});
                            ok1 = false;
                            goto OK1;
                        }
                    }
                }
            }
        }
    }}

    OK1:
    if(!ok1) {
        auto xij = get_x_from(eqn);
        ex xi = xij[0];
        ex ci = eqn.coeff(xi);
        ex c0 = eqn.subs(lst{xi==0});
        // handle eqn==ci xi - c0
        if((ci*xi+c0-eqn).is_zero() && is_a<numeric>(ci*c0) && (ci*c0)<0 && abs(c0)<abs(ci)) {
            if(Verbose>10) cout << "  \\--" << WHITE << "KillPowers ["<<kpi<<"]: "  << eqn << RESET << endl;
            ci = abs(ci);
            c0 = abs(c0);
            ex cc = c0/ci;
            symbol xx;
            // Part I: xi<cc
            auto f1 = fe.op(0);
            auto e1 = fe.op(1);
            for(int i=0; i<f1.nops(); i++) f1.let_op(i) = f1.op(i).subs(xi==xx*cc).subs(xx==xi);
            f1.let_op(0) = f1.op(0)*cc; // Jaccobi
            FunExp.push_back(lst{f1,e1});
            // Part II: xi>cc
            auto f2 = fe.op(0);
            auto e2 = fe.op(1);
            for(int i=0; i<f2.nops(); i++) f2.let_op(i) = f2.op(i).subs(xi==(1-cc)*xx+cc).subs(xx==xi);
            f2.let_op(0) = f2.op(0)*(1-cc); // Jaccobi
            FunExp.push_back(lst{f2,e2});
        } else {
            cout << MAGENTA << "Warning: Not handled with eqn4=" << eqn << RESET << endl;
            FunExp.push_back(fe);
        }
        return true;
    }
    
    // when ok2 && ok1 is true
    let_op_append(fe, 0, WF(1));
    let_op_append(fe, 1, 0);
    FunExp.push_back(fe);
    return false;
}

void SD::KillPowers(int bits) {
    int kpi = 0;
    while(kpi<30) {
        kpi++;
        if(kpi>10 && (kpi % 0)==0) {
            cout << "Warning: kip>10, (kpi=" << kpi << ") maybe a dead loop!" << endl;
        }
        
        vector<lst> funexp;
        for(auto fe : FunExp) {
            funexp.push_back(fe);
        }
        FunExp.clear();
        FunExp.shrink_to_fit();

        bool repeat = false;
        for(auto &fe : funexp) {
            bool ret;
            if(fe.nops()>2) ret = KillPowersWithDelta(fe, kpi);
            else ret = KillPowersWithoutDelta(fe, kpi, bits);
            if(ret) repeat = true;
        }
        if(!repeat) break;
    }
        
    for(auto &fe : FunExp) {
        auto nn = fe.op(0).nops()-1;
        if(fe.op(0).op(nn)==WF(1) && fe.op(1).op(nn)==0) {
            let_op_remove_last(fe, 0);
            let_op_remove_last(fe, 1);
        }
    }
    Normalizes();
}


}
