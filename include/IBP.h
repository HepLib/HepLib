#pragma once

#include "ExGiNaC.h"

#include <dlfcn.h>

#include <string>
#include <fcntl.h>
#include <signal.h>
#include <sys/syscall.h>
#include <sstream>
#include <ios>

namespace HepLib {

matrix RowReduce(matrix);
void RunFermat(char const * ifn);

DECLARE_FUNCTION_1P(P)
DECLARE_FUNCTION_1P(F)
DECLARE_FUNCTION_1P(L)
DECLARE_FUNCTION_1P(n)


class FormatCache {
public:
    map<ex,lst,ex_is_less> F2fmt;
    map<ex,lst,ex_is_less> I2fmt;
    map<lst,ex,ex_is_less> fmt2F;
    lst Cuts;
};

class IBP {
public:
    static symbol const ep;
    static symbol const eps;
    static symbol const d;
    
    static lst formatF(ex f, FormatCache &cache);
    static lst formatI(ex ibp, FormatCache &cache);
    static bool less(lst ls1, lst ls2);
    
    lst Cuts;
    lst Variables;
    int Verbose = 0;
    lst FSolution;
    
    void Prepare(lst loop, lst ext, lst prop, lst repl);
    void Generate(vector<lst> seeds);
    void Reduce();
    
private:
    lst preIBPs;
    lst IBPs;
    static ex collectF(ex);
};


}
