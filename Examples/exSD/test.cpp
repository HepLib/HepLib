#include "SD.h"
#include "mpreal.h"

using namespace HepLib;

int main(int argc, char** argv) {

    //SD::debug = true;

    auto ep = SD::ep;
    
    if(false) {
        archive a;
        ifstream in("test.gar");
        in >> a;
        in.close();
        
        lst sys_lst = {ep};

        auto tmp = a.unarchive_ex(sys_lst, "tmp");
        
        cout << tmp << endl;
        
        cout << tmp.subs(lst{ x(wild()) == 0 }) << endl;
        
        
        ex oo("ep(z(1))", lst{ep});
        cout << oo << endl;
    }
    
    if(false) {
        ex zz = tgamma(2+ep);
        cout << is_a<numeric>(zz.evalf()) << endl;
        cout << zz.evalf() << endl;
    }
    
    if(false) {
        XIntegrand xint;
        xint.Functions = lst{x(3),(-(x(3)*x(4)) - x(2)*(2*x(3) + x(4)))*(x(1)*x(3) + 4*(x(3)*x(4) + x(2)*(3*x(3) + x(4)))),x(2) + x(3),x(4), x(3)*x(4) + x(2)*(2*x(3) + x(4))};
        xint.Exponents = lst{1,-1 - 2*ep,2*ep,1,-1 + 3*ep};
        xint.Deltas.push_back(lst{x(1), x(2), x(3), x(4)});
        SD work;
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        //work.ParallelProcess = 0;
        work.epN = 0;
        work.CheckEnd = true;
        work.Evaluate(xint);
        cout << work.VEResult() << endl;
    }
    
    if(false) {
        XIntegrand xint;
        xint.Functions = lst{ 5, x(1)-ex(1)/2};
        xint.Exponents = lst{ 1, -1+ep };
        SD work;
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        //work.ParallelProcess = 0;
        work.Verbose = 0;
        work.epN = 2;
        work.CheckEnd = true;
        
        if(false) {
            auto cuba = new CUBA();
            cuba->VERBOSE = 1;
            cuba->Method = CUBA::VEGAS;
            cuba->VEGAS_NSTART = 1000;
            cuba->VEGAS_NINCREASE = 1000;
            cuba->VEGAS_NBATCH = 1000;
            work.Integrator = cuba;
            work.TryPTS = 1000000;
            work.RunPTS = 1000000;
            work.RunMAX = 10;
        }
        
        work.Evaluate(xint, "test");
        cout << work.VEResult() << endl;
    }
    
    if(false) {
        symtab table;
        table["ep"] = SD::ep;
        parser reader(table);
    
        ex res = x(1)*x(2)*x(3);
        
        auto w = wild();
                    auto w1 = wild(1);
                    auto w2 = wild(2);
        
        exset oo;
        find(res, x(w1)*x(w2), oo);
        cout << oo << endl;
        
    }
    
    if(false) {
        symbol s, ss;
        ex ex1 = VF(s+4);
        
        ex oo = ex1.diff(ss);
        cout << oo << endl;
        
        oo = series(ex1,ss,3);
        cout << oo << endl;
        
    }
    
    if(false) {
        lst oo = lst{ lst {0,1,2}, 3,4,5 };
        cout << oo << endl;
        oo.let_op(0).let_op(2) = x(1);
        cout << oo << endl;
    
    }
    
    if(false) {
        ex tt = 5+4*I;
        CppFormat cc(cout, "M");
        tt.print(cc);
        
        cout << endl;
        cout << tt << endl;
        cout << ("M" == "M") << endl;
    }
    
    if(false) { // without delta
        XIntegrand xint;
        xint.Functions = lst{ 1, pow(x(1)-3*x(2)+x(3),2)};
        xint.Exponents = lst{ 1, -1+ep };
        SD work;
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        work.Verbose = 100;
        work.epN = 1;
        work.CheckEnd = true;
        work.Evaluate(xint);
        cout << work.VEResult() << endl;
    }
    
    if(false) { // with delta
        XIntegrand xint;
        xint.Functions = lst{ 1, pow(x(1)-2*x(2),2)};
        xint.Exponents = lst{ 1, -1+ep };
        xint.Deltas.push_back(lst{x(1),x(2)});
        SD work;
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        work.Verbose = 100;
        work.epN = 3;
        work.RunPTS = 10000;
        work.CheckEnd = true;
        work.Evaluate(xint);
        cout << work.VEResult() << endl;
    }
    
    if(false) { // with delta
        XIntegrand xint;
        xint.Functions = lst{ 1, pow(x(1)-x(2)+x(3),2)};
        xint.Exponents = lst{ 1, -1+ep };
        xint.Deltas.push_back(lst{x(1),x(2),x(3)});
        SD work;
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        work.Verbose = 100;
        work.epN = 1;
        work.RunPTS = 10000;
        work.CheckEnd = true;
        work.Evaluate(xint);
        cout << work.VEResult() << endl;
    }
    
    if(false) {
        int digits = 64;
        cout.precision(digits);
        mpfr::mpreal::set_default_prec(mpfr::digits2bits(digits));
        mpfr::mpreal outs[10];
        
        cout << RED << "Not Wanted:" << RESET << endl;
        #pragma omp parallel for num_threads(omp_get_num_procs()) schedule(dynamic, 1)
        for(int i=0; i<10; i++) {
            //mpfr::mpreal::set_default_prec(mpfr::digits2bits(digits));
            outs[i] = mpfr::mpreal("1")/mpfr::mpreal("3");
            //cout << omp_get_thread_num() << endl;
        }
        for(int i=0; i<10; i++) {
            cout << outs[i] << endl;
        }
        
        cout << RED << "Wanted:" << RESET << endl;
        #pragma omp parallel for num_threads(omp_get_num_procs()) schedule(dynamic, 1)
        for(int i=0; i<10; i++) {
            mpfr::mpreal::set_default_prec(mpfr::digits2bits(digits));
            outs[i] = mpfr::mpreal("1")/mpfr::mpreal("3");
            //cout << omp_get_thread_num() << endl;
        }
        
        for(int i=0; i<10; i++) {
            cout << outs[i] << endl;
        }
    }
    
    if(false) { // without delta
        XIntegrand xint;
        xint.Functions = lst{ 1, pow(x(1)+x(2)-1,2), x(2)};
        xint.Exponents = lst{ 1, -1+ep, 2+ep};
        SD work;
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        work.Verbose = 100;
        work.epN = 3;
        work.CheckEnd = true;
        work.Evaluate(xint);
        cout << work.VEResult() << endl;
    }
    
    if(false) {
        cout << CV(1,2) << ": " << evalf(CV(1,2).hold()) << endl;
    }
    
    if(true) {
        #define Power pow
        symbol p1("p1"), p2("p2"), p3("p3"), p4("p4"), p5("p5"), l2("l2"), k2("k2");
        auto tt = SD::vs;
        ex oo = (x(1)*x(4) + x(2)*x(4) + x(3)*x(4) + x(1)*x(5) + x(2)*x(5) + x(3)*x(5) + x(4)*x(5) + x(1)*x(6) + x(2)*x(6) + x(3)*x(6) + x(5)*x(6) + x(1)*x(7) + x(2)*x(7) + x(3)*x(7) + x(5)*x(7))*(x(1)*x(3)*x(4) + x(1)*x(3)*x(5) + x(3)*x(4)*x(5) + x(1)*x(3)*x(6) + x(1)*x(4)*x(6) + x(2)*x(4)*x(6) + x(3)*x(4)*x(6) + x(1)*x(5)*x(6) + x(4)*x(5)*x(6) + x(1)*x(3)*x(7) + tt*x(2)*x(5)*x(7));

        cout << SD::PExpand(oo) << endl;
        
//        auto vecs = SecDecG::RunQHull(oo);
//
//
//        cout << endl << vecs.size() << endl;
//        for(auto ii : vecs) {
//            for( auto iii : ii) cout << iii << ", ";
//            cout << endl;
//        }
        
    }
    
    
    
    return 0;
}
