/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "BASIC.h"

namespace HepLib {
        
    bool ex_less(const ex &a, const ex &b) {
        static ex_is_less less;
        if(a.is_equal(b)) return false;
        
        // numeric
        if(is_a<numeric>(a) && is_a<numeric>(b)) {
            auto ab = a-b;
            auto abr = real_part(ab);
            if(!is_zero(abr)) return abr<0;
            else return imag_part(ab)<0;
        }
        if(is_a<numeric>(a)) return true;
        if(is_a<numeric>(b)) return false;
        
        // symbol/Symbol
        if(is_a<symbol>(a) && is_a<symbol>(b)) {
            string sa = ex_to<symbol>(a).get_name();
            string sb = ex_to<symbol>(b).get_name();
            return sa < sb;
        }
        if(is_a<symbol>(a)) return true;
        if(is_a<symbol>(b)) return false;
        
        // lst
        if(is_a<lst>(a) && is_a<lst>(b)) {
            auto na = a.nops();
            auto nb = b.nops();
            if(na!=nb) return (na<nb);
            for(int i=0; i<na; i++) {
                if(a.op(i).is_equal(b.op(i))) continue;
                return ex_less(a.op(i), b.op(i));
            }
            return false;
        }
        if(is_a<lst>(b)) return true;
        if(is_a<lst>(a)) return false;
        
        return less(a,b);
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
     
     /**
      * @brief sort the list in less order, or the reverse
      * @param ilst input lst, will be updated after call
      * @param ki the sort key is at .op(n)
      * @param less true for less order
      */
     void sort_lst_by(lst & ilst, int ki, bool less) {
        auto ivec = lst2vec(ilst);
        sort_vec_by(ivec,ki,less);
        for(auto i=0; i<ivec.size(); i++) ilst.let_op(i) = ivec[i];
     }
    
    
}
