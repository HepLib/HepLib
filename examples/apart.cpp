#include "HepLib.h"

using namespace HepLib;

int main() {

    Symbol a("a"), b("b"), c("c"), x("x"), y("y"), z("z");
    Vector p1("p1"), p2("p2"), p3("p3"), q1("q1"), q2("q2");
    Symbol m("m"), s("s");
    
    { // a general case
        ex expr = x*(x+y+a)*(x+y+b);
        expr = 1/expr;
        ex res = Apart(expr, lst{x,y,z});
        res = collect_ex(res, ApartIR(w1,w2));
        cout << "input: " << expr << endl;
        cout << "aparted result:" << res << endl;
        cout << "result in normal expression: " << ApartIR2ex(res) << endl;
        cout << "check : 0 = " << normal(ApartIR2ex(res)-expr) << endl;
        cout << endl;
    }
    
    ex p4 = -(p1+p2+p3);
        
    { // item with iEpsilon
        ex pro1 = a*x+iEpsilon;
        ex pro2 = b*y+iEpsilon;
        ex pro3 = c*z+iEpsilon;
        ex pro4 = x+y+z+iEpsilon;
        ex expr = 1/(pro1*pro2*pro3*pro4);
        
        auto vars = lst{x,y,z};
        
        ex res = Apart(expr, vars);
        ex chk = ApartIR2ex(res);
        
        cout << "input: " << expr << endl;
        cout << "aparted result:" << chk << endl;
        cout << "check : 0 = " << normal(subs(expr-chk,iEpsilon==0)) << endl;
        cout << endl;
    }
        
    { // vars in SP form, 13 propagators
        ex pro1 = SP(p1+p2+p3+q1+q2)-m*m;
        ex pro2 = SP(p3+q1+q2)-m*m;
        ex pro3 = SP(p2+q1+q2)-s;
        ex pro4 = SP(p3+p2+q1+q2)-m*m;
        ex pro5 = SP(p1+p2+q1)-m*m+iEpsilon;
        ex pro6 = SP(p3+q2)-m*m;
        ex pro7 = SP(p4+q2)-s;
        ex pro8 = SP(p1+p2+q1+q2)-s;
        ex pro9 = SP(p4+q2+q1)-s;
        ex pro10 = SP(p2+q2+q1)-s;
        ex pro11 = SP(p4+p2+q2-q1)-m;
        ex pro12 = SP(p2+p2+q2-q1)-s*m;
        ex pro13 = SP(p1+p2+q2+q1)-m;
        
        ex expr = pow(pro1,-1) * pow(pro2,-2) * pow(pro3,-1) * pow(pro4,-1) * pow(pro5,-2) * pow(pro6,-2) * pow(pro7,-4) * pow(pro8,-1) * pow(pro9,-2) * pow(pro10,-2) * pow(pro11,-3) * pow(pro12,-2) * pow(pro13,-2);
        
        auto vars = lst{SP(q1), SP(q2,q2), SP(q1,q2), SP(q1,p1), SP(q1,p2), SP(q1,p3), SP(q1,p4), SP(q2,p1), SP(q2,p2), SP(q2,p3), SP(q2,p4)};
        
        exmap sign;
        sign[SP(q1)]=-1;
        sign[SP(q2)]=-1;
        
        ex res = Apart(expr, vars, sign);
        ex chk = ApartIR2ex(res);
        
        cout << "we omit the output due to large experssions." << endl;
        cout << "numerical check : 0 = " << normal(subs(expr-chk,lst{
            iEpsilon==0,SP(q1)==ex(1)/14,SP(q2)==ex(1)/21,
            SP(q1,p1)==ex(1)/41,SP(q1,p2)==12,SP(q1,p3)==ex(1)/41,SP(q1,p4)==ex(1)/15,
            SP(q2,p1)==ex(1)/50,SP(q2,p2)==ex(1)/15,SP(q2,p3)==ex(1)/42,SP(q2,p4)==71,
            SP(p1)==1,SP(p2)==ex(1)/35,SP(p3)==ex(1)/41,SP(p4)==ex(1)/51,
            SP(p1,p2)==ex(1)/9,SP(p1,p3)==1,SP(p1,p4)==ex(1)/9,SP(p2,p3)==ex(1)/2,
            SP(p2,p4)==11,SP(p3,p4)==15,s==ex(1)/2,m==ex(1)/11
        })) << endl;
        cout << endl;
    }
    
    
    { // double box, vars inferred from loop/external mommena
        ex pro1 = SP(q1)-m*m;
        ex pro2 = SP(q2)-m*m;
        ex pro3 = SP(q1+p1)-m*m;
        ex pro4 = SP(q2+p4)-m*m;
        ex pro5 = SP(q1-p3)-m*m;
        ex pro6 = SP(q1-p2)-m*m;
        ex pro7 = SP(q1+q2+p1+p4)-m*m;
        ex pro8 = SP(q1+q2+p1+p4+p2)-m*m;
        ex pro9 = SP(q1+q2+p1+p4+p2+p3)-m*m;
        
        letSP(p1)=m*m;
        letSP(p2)=m*m;
        letSP(p3)=m*m;
        letSP(p1,p2)=0;
        letSP(p1,p3)=s;
        letSP(p2,p3)=0;
        
        ex expr = pow(pro1,-1) * pow(pro2,-1) * pow(pro3,-1) * pow(pro4,-1) * pow(pro5,-1) * pow(pro6,-1) * pow(pro7,-1) * pow(pro8,-1) * pow(pro9,-1);
        ex num = SP(q1+q2+p1+p4+p3+p4, q1+q2+p1+p4) * SP(q1+q2+p1+p4+p3+p4, p4+q1+q2+p1+p4);
                
        ex res = Apart(num * expr, lst{q1, q2}, lst{p1,p2,p3});
        ex chk = ApartIR2ex(res);
        ex nchk = subs(num*expr-chk,SP_map); // SP_map is defined in HepLib
                
        cout << "we omit the output due to large experssions." << endl;
        cout << "numerical check : 0 = " << normal(subs(nchk,lst{
            iEpsilon==0,SP(q1)==ex(1)/14,SP(q2)==ex(1)/21,
            SP(q1,p1)==ex(1)/41,SP(q1,p2)==12,SP(q1,p3)==ex(1)/41,SP(q1,p4)==ex(1)/15,
            SP(q2,p1)==ex(1)/50,SP(q2,p2)==ex(1)/15,SP(q2,p3)==ex(1)/42,SP(q2,p4)==71,
            s==ex(1)/2,m==ex(1)/11
        })) << endl;
        cout << endl;
        
    }
                
    return 0;
}
