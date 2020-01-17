#include "SD.h"
#include "mpreal.h"

using namespace HepLib;

int main(int argc, char** argv) {

    //SD::debug = true;

    auto ep = SD::ep;
    
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
        cout << "check with:" << endl;
        cout << "ep^3*(-5.318503740159345 +- 6.09714E-6)+ep^2*(-2.6222491120986264 +- 5.60406E-6)+ep*(-1.23104906018664843647241070465427486519 +- 3.4071E-19)+(-0.5 +- 0)" << endl << endl;
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
        cout << "check with:" << endl;
        cout << "(-0.4999999912658149 +- 6.7716E-6)+ep*(-1.0000008702897987 +- 1.26966E-5)+ep^(-1)*(0 +- 1.96161464E-16)" << endl << endl;
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
        cout << "check with:" << endl;
        cout << "I*(-5.12124E-11 +- 2.94233E-9)+((0.24889315261249362 +- 1.17935E-5)+I*(3.50568E-10 +- 2.52899E-8))*ep^2+(6.4343E-9 +- 2.82284E-6)+ep^3*(I*(-7.5311E-10 +- 3.90023E-8)+(-9.854044005614545 +- 2.00735E-5))+ep*(I*(3.44408E-12 +- 1.12738E-8)+(-2.3550642800278156 +- 9.1018E-6))+ep^(-1)*(-0.5 +- 4.321953751E-47)" << endl << endl;
    }
    
    if(false) { // without delta
        XIntegrand xint;
        xint.Functions = lst{ 1, x(1)-x(2)};
        xint.Exponents = lst{ 1, -1+ep};
        SD work;
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        work.Verbose = 100;
        work.epN = 3;
        work.CheckEnd = true;
        work.Evaluate(xint);
        cout << work.VEResult() << endl;
        cout << "check with:" << endl;
        cout << "ep*(I*(-3.141592643497266 +- 4.4013E-6)+(4.934802200544682 +- 1.74201E-6))+ep^3*(I*(2.026119323775967 +- 9.9703E-6)+(0.8760900995073708 +- 1.20121E-5))+ep^2*((-4.934802184622438 +- 8.5754E-6)+I*(-2.026120013987933 +- 4.8383E-6))+I*(3.141592654 +- 1.146412548E-38)+(5.294242378E-41 +- 5.79961065E-39)" << endl << endl;
    }
    
    if(false) { // without delta
        XIntegrand xint;
        xint.Functions = lst{ 1, x(1)-2*x(2)};
        xint.Exponents = lst{ 1, -1+ep};
        SD work;
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        work.Verbose = 100;
        work.epN = 3;
        work.CheckEnd = true;
        work.Evaluate(xint);
        cout << work.VEResult() << endl;
        cout << "check with:" << endl;
        cout << "(-0.6931471805054186 +- 5.25342E-8)+I*(1.5707963267980893 +- 1.86634E-9)+((2.9203215085122483 +- 1.26295E-6)+I*(0.6067897683760062 +- 2.20715E-6))*ep+ep^2*((0.44471972382655117 +- 5.21156E-6)+I*(-2.435951434158803 +- 3.75704E-6))+((-1.298226402539505 +- 1.03668E-5)+I*(-0.9716668223978889 +- 1.17999E-5))*ep^3" << endl << endl;
    }
    
    if(false) { // without delta
        XIntegrand xint;
        xint.Functions = lst{ 1, pow(x(1)-2*x(2),5)};
        xint.Exponents = lst{ 1, -1+ep};
        SD work;
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        work.Verbose = 100;
        work.epN = 3;
        work.CheckEnd = true;
        work.Evaluate(xint);
        cout << work.VEResult() << endl;
        cout << "check with:" << endl;
        cout << "ep^2*(I*(-0.27735894763144603 +- 3.356E-6)+(0.23800256594336502 +- 3.4734E-7))+(0.07812500057646682 +- 2.843E-7)+ep*((0.2098138728313142 +- 1.06695E-6)+I*(-0.11453723397370566 +- 8.94135E-7))+ep^3*((0.3100030883895008 +- 2.80026E-7)+I*(-0.28367076724864015 +- 2.74675E-7))+I*(4.43946E-14 +- 9.0829E-11)" << endl << endl;
    }
    
    if(false) { // with delta
        XIntegrand xint;
        xint.Functions = lst{ 1, x(1)-2*x(2)};
        xint.Exponents = lst{ 1, -1+ep};
        SD work;
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        work.Verbose = 100;
        work.epN = 3;
        xint.Deltas.push_back(lst{x(1),x(2)});
        work.Evaluate(xint);
        cout << work.VEResult() << endl;
        cout << "check with:" << endl;
        cout << "(I*(-1.1358713887357565 +- 1.13014E-5)+(-0.9609533480725917 +- 9.073E-6))*ep^3+(-0.2310490602 +- 0)+I*(1.047197551 +- 0)+ep^2*((1.121680044071696 +- 1.29481E-6)+I*(-1.4710063303868242 +- 5.8686E-6))+ep*((1.564858564528526199027898078924914360594 +- 2.713E-18)+I*(0.725862030101200710166896299412537982445 +- 1.22963E-17))" << endl << endl;
    }
    
    if(false) { // without delta
        XIntegrand xint;
        xint.Functions = lst{ 1, 1-2*x(2)};
        xint.Exponents = lst{ 1, -1+ep};
        SD work;
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        work.Verbose = 100;
        work.epN = 3;
        work.Evaluate(xint);
        cout << work.VEResult() << endl;
        cout << "check with:" << endl;
        cout << "I*(1.5707963267945655 +- 1.84498E-6)+(2.888894917E-34 +- 1.228299845E-06)+ep*((2.4674011002723395 +- 1.4424E-10)+I*(1.06383E-17 +- 8.39734E-10))+ep^3*((-2.0293560632083842 +- 5.01876E-9)+I*(4.75046E-17 +- 7.7206E-9))+ep^2*((1.67107E-17 +- 1.62754E-9)+I*(-2.583856390024985 +- 3.4349E-9))" << endl << endl;
    }
    
    if(false) { // without delta
        XIntegrand xint;
        xint.Functions = lst{ 1, 3-2*x(2)};
        xint.Exponents = lst{ 1, -1+ep};
        SD work;
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        work.Verbose = 100;
        work.epN = 3;
        work.Evaluate(xint);
        cout << work.VEResult() << endl;
        cout << "check with:" << endl;
        cout << "(0.5493061443340549 +- 1.67324E-13)+ep^3*(I*(3.95703E-20 +- 1.39949E-11)+(0.030348454041803575 +- 2.55734E-12))+ep^2*((0.1104974133453256 +- 2.22172E-12)+I*(-4.58428E-20 +- 1.73448E-11))+I*(-8.9945E-21 +- 8.66814E-13)+ep*(I*(3.0894E-20 +- 8.6804E-12)+(0.3017372402031455 +- 1.0645E-12))" << endl << endl;
    }
    
    if(false) { // without delta
        XIntegrand xint;
        xint.Functions = lst{ 1, 1-2*x(2)};
        xint.Exponents = lst{ 1, -1+ep};
        SD work;
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        work.Verbose = 100;
        work.epN = 3;
        work.Evaluate(xint);
        cout << work.VEResult() << endl;
        cout << "check with:" << endl;
        cout << "ep^3*(I*(4.75046E-17 +- 7.7206E-9)+(-2.0293560632083842 +- 5.01876E-9))+(2.888894917E-34 +- 1.228299845E-06)+ep^2*(I*(-2.583856390024985 +- 3.4349E-9)+(1.67107E-17 +- 1.62754E-9))+(I*(1.06383E-17 +- 8.39734E-10)+(2.4674011002723395 +- 1.4424E-10))*ep+I*(1.5707963267945655 +- 1.84498E-6)" << endl << endl;
    }
    
    if(false) { // without delta
        XIntegrand xint;
        xint.Functions = lst{ 1, pow(1-2*x(2),5)};
        xint.Exponents = lst{ 1, -1+ep};
        SD work;
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        work.Verbose = 100;
        work.epN = 3;
        work.Evaluate(xint);
        cout << work.VEResult() << endl;
        cout << "check with:" << endl;
        cout << "(-2.996264821E-36 +- 2.146501408E-34)+((-0.7710628438 +- 4.869955915E-31)+I*(0.03237178235 +- 1.664824146E-33))*ep^3+((-0.6168502751 +- 5.296862617E-32)+I*(-0.4908738521 +- 7.354290341E-35))*ep^2+I*(0 +- 1.605965375E-38)+(I*(-0.3926990817 +- 1.589979186E-36)+(6.483586029E-37 +- 3.795654524E-33))*ep" << endl << endl;
    }
    
    
    if(false) { // without delta
        XIntegrand xint;
        xint.Functions = lst{ 1, pow(1-2*x(2),1)};
        xint.Exponents = lst{ 1, 5*(-1+ep)};
        SD work;
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        work.Verbose = 100;
        work.epN = 3;
        work.Evaluate(xint);
        cout << work.VEResult() << endl;
        cout << "check with:" << endl;
        cout << "(1.39333E-28 +- 1.19694E-11)+(I*(-1.9634954084936207 +- 3.47274E-10)+(-1.16233E-20 +- 3.1221E-10))*ep+ep^2*((-15.421256876702122 +- 4.07687E-9)+I*(-2.454369260617026 +- 4.48597E-9))+I*(-1.47993E-21 +- 1.30842E-11)+ep^3*(I*(77.6775506125095 +- 3.56385E-8)+(-19.276571095877653 +- 3.91E-8))" << endl << endl;
    }
    
    if(false) { // without delta
        XIntegrand xint;
        xint.Functions = lst{ 1, pow(1-2*x(2),4)};
        xint.Exponents = lst{ 1, (-1+ep)};
        SD work;
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        work.Verbose = 100;
        work.epN = 3;
        work.Evaluate(xint);
        cout << work.VEResult() << endl;
        cout << "check with:" << endl;
        cout << "ep^3*((-0.7901234568 +- 4.365732098E-32)+I*(1.412853787E-41 +- 6.79405539E-34))+ep^2*((-0.5925925926 +- 4.971404441E-33)+I*(0 +- 3.56600453E-35))+(-0.3333333333 +- 2.394812699E-35)+ep*((-0.4444444444 +- 3.72948807E-34)+I*(1.808452847E-39 +- 1.103138152E-36))+I*(0 +- 1.145611015E-38)" << endl << endl;
    }
    
    if(true) { // without delta
        XIntegrand xint;
        xint.Functions = lst{ 1, pow(1-2*x(2),1)};
        xint.Exponents = lst{ 1, 4*(-1+ep)};
        SD work;
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        work.Verbose = 100;
        work.epN = 3;
        work.Evaluate(xint);
        cout << work.VEResult() << endl;
        cout << "check with:" << endl;
        cout << "((12.56687994219322 +- 1.40373E-11)+I*(2.792526803190927 +- 1.47529E-11))*ep^2+ep*(I*(2.094395102393195 +- 1.82328E-6)+(-0.4444444444444465 +- 5.2678E-7))+I*(1.32119E-30 +- 7.1761E-8)+ep^3*(I*(-51.39890058294511 +- 9.77645E-11)+(16.75583992292429 +- 1.13734E-10))+(-0.33333333333333326 +- 2.06924E-8)" << endl << endl;
    }
    
    return 0;
}
