/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "BASIC.h"

namespace HepLib {
    
    long long node_number(const ex & expr, int level) {
        if(expr.nops()<1) return level+1;
        long long tot = 0;
        for(auto item : expr) tot += node_number(item,level+1)+level;
        return tot;
    }
    
    bool ex_less(const ex &a, const ex &b) {
        static map<ex,bool,ex_is_less> cache;
        ex key = lst{a,b};
        if(using_cache && cache.find(key)!=cache.end()) return cache[key];
        if(a.is_equal(b)) return using_cache ? (cache[key]=false) : false;
        
        // numeric
        if(is_a<numeric>(a) && is_a<numeric>(b)) {
            auto ab = a-b;
            auto abr = real_part(ab);
            if(!is_zero(abr)) return using_cache ? (cache[key]=(abr<0)) : (abr<0);
            else return using_cache ? (cache[key]=(imag_part(ab)<0)) : (imag_part(ab)<0);
        }
        if(is_a<numeric>(a)) return using_cache ? (cache[key]=true) : true;
        if(is_a<numeric>(b)) return using_cache ? (cache[key]=false) : false;
        
        // symbol
        if(is_a<symbol>(a) && is_a<symbol>(b)) return using_cache ? (cache[key]=(ex2str(a) < ex2str(b))) : (ex2str(a) < ex2str(b));
        if(is_a<symbol>(a)) return using_cache ? (cache[key]=true) : true;
        if(is_a<symbol>(b)) return using_cache ? (cache[key]=false) : false;
        
        // matrix
        if(is_a<matrix>(a) && is_a<matrix>(b)) {
            auto ma = ex_to<matrix>(a);
            auto mb = ex_to<matrix>(b);
            if(ma.cols() != mb.cols()) return using_cache ? (cache[key]=(ma.cols() < mb.cols())) : (ma.cols() < mb.cols());
            if(ma.rows() != mb.rows()) return using_cache ? (cache[key]=(ma.rows() < mb.rows())) : (cache[key]=(ma.rows() < mb.rows()));
            for(int c=0; c<ma.cols(); c++) {
            for(int r=0; r<ma.rows(); r++) {
                if(ma(r,c).is_equal(mb(r,c))) continue;
                return using_cache ? (cache[key]=ex_less(ma(r,c),mb(r,c))) : ex_less(ma(r,c),mb(r,c));
            }}
            return using_cache ? (cache[key]=false) : false;
        }
        if(is_a<matrix>(b)) return using_cache ? (cache[key]=true) : true;
        if(is_a<matrix>(a)) return using_cache ? (cache[key]=false) : false;
        
        // lst
        if(is_a<lst>(a) && is_a<lst>(b)) {
            auto na = a.nops();
            auto nb = b.nops();
            if(na!=nb) return using_cache ? (cache[key]=(na<nb)) : (na<nb);
            for(int i=0; i<na; i++) {
                if(a.op(i).is_equal(b.op(i))) continue;
                return using_cache ? (cache[key]=ex_less(a.op(i), b.op(i))) : ex_less(a.op(i), b.op(i));
            }
            return using_cache ? (cache[key]=false) : false;
        }
        if(is_a<lst>(b)) return using_cache ? (cache[key]=true) : true;
        if(is_a<lst>(a)) return using_cache ? (cache[key]=false) : false;
        
        // atomic
        auto an = a.nops();
        auto bn = b.nops();
        if(an==0 && bn==0) {
            string na = a.return_type_tinfo().tinfo->name();
            string nb = b.return_type_tinfo().tinfo->name();
            auto nc = na.compare(nb);
            if(nc<0) return using_cache ? (cache[key]=true) : true;
            else if(nc>0) return using_cache ? (cache[key]=false) : false;
            else return using_cache ? (cache[key]=(ex2str(a) < ex2str(b))) : (ex2str(a) < ex2str(b));
        }
        
        // power
        if(is_a<GiNaC::power>(a) || is_a<GiNaC::power>(b)) {
            ex ae=a, be=b;
            ex an=1, bn=1;
            if(is_a<GiNaC::power>(a)) {
                ae = a.op(0);
                an = a.op(1);
            }
            if(is_a<GiNaC::power>(b)) {
                be = b.op(0);
                bn = b.op(1);
            }
            if(!ae.is_equal(be)) return using_cache ? (cache[key]=ex_less(ae,be)) : ex_less(ae,be);
            if(an.info(info_flags::real) && bn.info(info_flags::real)) return using_cache ? (cache[key]=(an<bn)) : (an<bn);
            if(!is_zero(an-bn)) return using_cache ? (cache[key]=ex_less(an,bn)) : ex_less(an,bn);
            return using_cache ? (cache[key]=false) : false;
        }
        
        // function
        if(is_a<GiNaC::function>(a) && is_a<GiNaC::function>(b)) {
            string na = ex_to<GiNaC::function>(a).get_name();
            string nb = ex_to<GiNaC::function>(b).get_name();
            auto nc = na.compare(nb);
            if(nc<0) return using_cache ? (cache[key]=true) : true;
            else if(nc>0) return using_cache ? (cache[key]=false) : false;
            if(an!=bn) return using_cache ? (cache[key]=(an < bn)) : (an < bn);
            for(int i=0; i<an; i++) {
                if(a.op(i).is_equal(b.op(i))) continue;
                return using_cache ? (cache[key]=ex_less(a.op(i), b.op(i))) : ex_less(a.op(i), b.op(i));
            }
            return using_cache ? (cache[key]=false) : false;
        }
        
        // add
        if(is_a<add>(a) && is_a<add>(b)) {
            auto as = add2lst(a);
            auto bs = add2lst(b);
            auto na = as.nops();
            auto nb = bs.nops();
            sort_lst(as,false);
            sort_lst(bs,false);
            int nn = ((na>nb) ? nb : na);
            for(int i=0; i<nn; i++) {
                if(as.op(i).is_equal(bs.op(i))) continue;
                return using_cache ? (cache[key]=ex_less(as.op(i), bs.op(i))) : ex_less(as.op(i), bs.op(i));
            }
            if(na!=nb) return using_cache ? (cache[key]=(na<nb)) : (na<nb);
            return using_cache ? (cache[key]=false) : false;
        }
        if(is_a<add>(a)) return using_cache ? (cache[key]=false) : false;
        if(is_a<add>(b)) return using_cache ? (cache[key]=true) : true;
        
        // mul
        if(is_a<mul>(a) && is_a<mul>(b)) {
            auto as = mul2lst(a);
            auto bs = mul2lst(b);
            auto na = as.nops();
            auto nb = bs.nops();
            sort_lst(as,false);
            sort_lst(bs,false);
            int nn = ((na>nb) ? nb : na);
            for(int i=0; i<nn; i++) {
                if(as.op(i).is_equal(bs.op(i))) continue;
                return using_cache ? (cache[key]=ex_less(as.op(i), bs.op(i))) : ex_less(as.op(i), bs.op(i));
            }
            if(na!=nb) return using_cache ? (cache[key]=(na<nb)) : (na<nb);
            return using_cache ? (cache[key]=false) : false;
        }
        if(is_a<mul>(a)) return using_cache ? (cache[key]=false) : false;
        if(is_a<mul>(b)) return using_cache ? (cache[key]=true) : true;
        
        // type
        string tna = a.return_type_tinfo().tinfo->name();
        string tnb = b.return_type_tinfo().tinfo->name();
        auto tnc = tna.compare(tnb);
        if(tnc<0) return using_cache ? (cache[key]=true) : true;
        else if(tnc>0) return using_cache ? (cache[key]=false) : false;
        
        // node_number
        auto nna = node_number(a);
        auto nnb = node_number(b);
        if(nna!=nnb) return using_cache ? (cache[key]=(nna < nnb)) : (nna < nnb);

        // all others
        if(an!=bn) return (an<bn);
        for(int i=0; i<an; i++) {
            if(a.op(i).is_equal(b.op(i))) continue;
            return using_cache ? (cache[key]=ex_less(a.op(i), b.op(i))) : ex_less(a.op(i), b.op(i));
        }
        
        return using_cache ? (cache[key]=(ex2str(a) < ex2str(b))) : (ex2str(a) < ex2str(b));
    }
    
    bool ex_less_cache(const ex &a, const ex &b, map<ex,bool,ex_is_less> &cache) {
        ex key = lst{a,b};
        if(cache.find(key)!=cache.end()) return cache[key];
        if(a.is_equal(b)) return cache[key]=false;
        
        // numeric
        if(is_a<numeric>(a) && is_a<numeric>(b)) {
            auto ab = a-b;
            auto abr = real_part(ab);
            if(!is_zero(abr)) return cache[key]=(abr<0);
            else return cache[key]=(imag_part(ab)<0);
        }
        if(is_a<numeric>(a)) return cache[key]=true;
        if(is_a<numeric>(b)) return cache[key]=false;
        
        // symbol
        if(is_a<symbol>(a) && is_a<symbol>(b)) return cache[key]=(ex2str(a) < ex2str(b));
        if(is_a<symbol>(a)) return cache[key]=true;
        if(is_a<symbol>(b)) return cache[key]=false;
        
        // matrix
        if(is_a<matrix>(a) && is_a<matrix>(b)) {
            auto ma = ex_to<matrix>(a);
            auto mb = ex_to<matrix>(b);
            if(ma.cols() != mb.cols()) return cache[key]=(ma.cols() < mb.cols());
            if(ma.rows() != mb.rows()) return cache[key]=(ma.rows() < mb.rows());
            for(int c=0; c<ma.cols(); c++) {
            for(int r=0; r<ma.rows(); r++) {
                if(ma(r,c).is_equal(mb(r,c))) continue;
                return cache[key]=ex_less_cache(ma(r,c),mb(r,c),cache);
            }}
            return cache[key]=false;
        }
        if(is_a<matrix>(b)) return cache[key]=true;
        if(is_a<matrix>(a)) return cache[key]=false;
        
        // lst
        if(is_a<lst>(a) && is_a<lst>(b)) {
            auto na = a.nops();
            auto nb = b.nops();
            if(na!=nb) return cache[key]=(na<nb);
            for(int i=0; i<na; i++) {
                if(a.op(i).is_equal(b.op(i))) continue;
                return cache[key]=ex_less_cache(a.op(i), b.op(i),cache);
            }
            return cache[key]=false;
        }
        if(is_a<lst>(b)) return cache[key]=true;
        if(is_a<lst>(a)) return cache[key]=false;
        
        // atomic
        auto an = a.nops();
        auto bn = b.nops();
        if(an==0 && bn==0) {
            string na = a.return_type_tinfo().tinfo->name();
            string nb = b.return_type_tinfo().tinfo->name();
            auto nc = na.compare(nb);
            if(nc<0) return cache[key]=true;
            else if(nc>0) return cache[key]=false;
            else return cache[key]=(ex2str(a) < ex2str(b));
        }
        
        // power
        if(is_a<GiNaC::power>(a) || is_a<GiNaC::power>(b)) {
            ex ae=a, be=b;
            ex an=1, bn=1;
            if(is_a<GiNaC::power>(a)) {
                ae = a.op(0);
                an = a.op(1);
            }
            if(is_a<GiNaC::power>(b)) {
                be = b.op(0);
                bn = b.op(1);
            }
            if(!ae.is_equal(be)) return cache[key]=ex_less_cache(ae,be,cache);
            if(an.info(info_flags::real) && bn.info(info_flags::real)) return cache[key]=(an<bn);
            if(!is_zero(an-bn)) return cache[key]=ex_less_cache(an,bn,cache);
            return cache[key]=false;
        }
        
        // function
        if(is_a<GiNaC::function>(a) && is_a<GiNaC::function>(b)) {
            string na = ex_to<GiNaC::function>(a).get_name();
            string nb = ex_to<GiNaC::function>(b).get_name();
            auto nc = na.compare(nb);
            if(nc<0) return cache[key]=true;
            else if(nc>0) return cache[key]=false;
            if(an!=bn) return cache[key]=(an < bn);
            for(int i=0; i<an; i++) {
                if(a.op(i).is_equal(b.op(i))) continue;
                return cache[key]=ex_less_cache(a.op(i), b.op(i), cache);
            }
            return cache[key]=false;
        }
        
        // add
        if(is_a<add>(a) && is_a<add>(b)) {
            auto as = add2lst(a);
            auto bs = add2lst(b);
            auto na = as.nops();
            auto nb = bs.nops();
            sort_lst(as,false);
            sort_lst(bs,false);
            int nn = ((na>nb) ? nb : na);
            for(int i=0; i<nn; i++) {
                if(as.op(i).is_equal(bs.op(i))) continue;
                return cache[key]=ex_less_cache(as.op(i), bs.op(i), cache);
            }
            if(na!=nb) return cache[key]=(na<nb);
            return cache[key]=false;
        }
        if(is_a<add>(a)) return cache[key]=false;
        if(is_a<add>(b)) return cache[key]=true;
        
        // mul
        if(is_a<mul>(a) && is_a<mul>(b)) {
            auto as = mul2lst(a);
            auto bs = mul2lst(b);
            auto na = as.nops();
            auto nb = bs.nops();
            sort_lst(as,false);
            sort_lst(bs,false);
            int nn = ((na>nb) ? nb : na);
            for(int i=0; i<nn; i++) {
                if(as.op(i).is_equal(bs.op(i))) continue;
                return cache[key]=ex_less_cache(as.op(i), bs.op(i), cache);
            }
            if(na!=nb) return cache[key]=(na<nb);
            return cache[key]=false;
        }
        if(is_a<mul>(a)) return cache[key]=false;
        if(is_a<mul>(b)) return cache[key]=true;
        
        // type
        string tna = a.return_type_tinfo().tinfo->name();
        string tnb = b.return_type_tinfo().tinfo->name();
        auto tnc = tna.compare(tnb);
        if(tnc<0) return cache[key]=true;
        else if(tnc>0) return cache[key]=false;
        
        // node_number
        auto nna = node_number(a);
        auto nnb = node_number(b);
        if(nna!=nnb) return cache[key]=(nna < nnb);

        // all others
        if(an!=bn) return (an<bn);
        for(int i=0; i<an; i++) {
            if(a.op(i).is_equal(b.op(i))) continue;
            return cache[key]=ex_less_cache(a.op(i), b.op(i), cache);
        }
        
        return cache[key]=(ex2str(a) < ex2str(b));
    }
     
     /**
      * @brief sort the list in less order, or the reverse
      * @param ilst input lst, will be updated after call
      * @param less true for less order
      */
     void sort_lst(lst & ilst, bool less) {
        auto ivec = lst2vec(ilst);
        sort_vec(ivec,less);
        for(auto i=0; i<ivec.size(); i++) ilst.let_op(i) = ivec[i];
     }
     
     void sort_lst(lst & ilst, map<ex,bool,ex_is_less> &cache, bool less) {
        auto ivec = lst2vec(ilst);
        sort_vec(ivec, cache, less);
        for(auto i=0; i<ivec.size(); i++) ilst.let_op(i) = ivec[i];
     }
     
     /**
      * @brief sort the list in less order, or the reverse
      * @param ilst input lst, will be updated after call
      * @param ki the sort key is at .op(n)
      * @param less true for less order
      */
     void sort_lst_by(lst & ilst, int ki, bool less) {
        auto ivec = lst2vec(ilst);
        std::sort(ivec.begin(), ivec.end(), [ki](const auto &as, const auto &bs){
            return ex_less(as.op(ki),bs.op(ki));
        });
        auto n = ivec.size();
        if(less) for(auto i=0; i<n; i++) ilst.let_op(i) = ivec[i];
        else for(auto i=0; i<n; i++) ilst.let_op(i) = ivec[n-1-i];
     }
     
     /**
      * @brief sort the list in less order, or the reverse
      * @param ivec input exvector, will be updated after call
      * @param less true for less order
      */
     void sort_vec(exvector & ivec, bool less) {
        std::sort(ivec.begin(), ivec.end(), [less](const auto &a, const auto &b){
            if(less) return ex_less(a,b);
            else return ex_less(b,a);
        });
     }
     
     void sort_vec(exvector & ivec, map<ex,bool,ex_is_less> &cache, bool less) {
        std::sort(ivec.begin(), ivec.end(), [&cache,less](const auto &a, const auto &b){
            if(less) return ex_less_cache(a,b,cache);
            else return ex_less_cache(b,a,cache);
        });
     }
     
     /**
      * @brief sort the list in less order, or the reverse
      * @param ivec input exvector, will be updated after call
      * @param ki the sort key is at .op(n)
      * @param less true for less order
      */
     void sort_vec_by(exvector & ivec, int ki, bool less) {
        std::sort(ivec.begin(), ivec.end(), [ki,less](const auto &as, const auto &bs){
            if(less) return ex_less(as.op(ki),bs.op(ki));
            else return ex_less(bs.op(ki),as.op(ki));
        });
     }
    
}
