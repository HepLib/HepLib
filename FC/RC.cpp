/**
 * @file
 * @brief Functions for Renormalization Constants
 */
 
#include "FC.h"

namespace HepLib::FC {

    //-----------------------------------------------------------
    // RC namespace
    //-----------------------------------------------------------
    namespace RC {
        
        /**
         * @brief the Z2 renomalization constant
         * @param name a string to identify the particle
         * @param m the mass of the particle
         * @param loop the loop order
         * @return Z2 for particle of name at loop order
         */
        ex Z2(string name, ex m, int loop) {
            return Z2(Symbol(name), m, loop);
        }
        
        /**
         * @brief the Z2 renomalization constant
         * @param n a Symbol to identify the particle
         * @param m the mass of the particle
         * @param loop the loop order
         * @return Z2 for particle of n at loop order
         */
        ex Z2(Symbol n, ex m, int loop) {
            if(loop==0) return 1;
            ex rr = 1;
            ex Tf = ex(1)/2;
            ex as_pi = as/Pi;
            ex CA = NF;
            ex CF = NA/(2*NF);
            ex lmu = 2*log(mu/m);
            if(n==Symbol("g")) {
                if(loop>0) rr += (nH*Tf*(as_pi)*((-2*lmu)/ex(3) - (2*pow(ep, -1))/ex(3) - (ep*pow(lmu, 2))/ex(3) - (pow(ep, 2)*pow(lmu, 3))/ex(9) - (ep*pow(Pi, 2))/ex(18) - (lmu*pow(ep, 2)*pow(Pi, 2))/ex(18) + (2*pow(ep, 2)*zeta(3))/ex(9)))/ex(2);
                if(loop>1) rr += (nH*Tf*(CF*(-15/ex(4) - lmu - pow(ep, -1)/ex(2)) + nL*Tf*((-4*pow(ep, -2))/ex(9) - (4*lmu*pow(ep, -1))/ex(9) - (2*pow(lmu, 2))/ex(9) - pow(Pi, 2)/ex(27)) + nH*Tf*((4*lmu*pow(ep, -1))/ex(9) + (2*pow(lmu, 2))/ex(3) + pow(Pi, 2)/ex(27)) + CA*(13/ex(48) - (5*lmu)/ex(4) + (35*pow(ep, -2))/(36) - (5*pow(ep, -1))/ex(8) + (13*lmu*pow(ep, -1))/ex(18) + pow(lmu, 2)/ex(9) + (13*pow(Pi, 2))/ex(216)))*pow(as_pi, 2))/ex(4);
            } else if(n==Symbol("Q")) {
                if(loop>0) rr += (CF*(as_pi)*(-2 - 4*ep - (3*lmu)/ex(2) - 2*ep*lmu - (3*pow(ep, -1))/ex(2) - 8*pow(ep, 2) - 4*lmu*pow(ep, 2) - (3*ep*pow(lmu, 2))/ex(4) - pow(ep, 2)*pow(lmu, 2) - (pow(ep, 2)*pow(lmu, 3))/ex(4) - (ep*pow(Pi, 2))/ex(8) - (pow(ep, 2)*pow(Pi, 2))/ex(6) - (lmu*pow(ep, 2)*pow(Pi, 2))/ex(8) + (pow(ep, 2)*zeta(3))/ex(2)))/ex(2);
                if(loop>1) rr += (CF*pow(as_pi, 2)*(nH*Tf*(947/ex(72) + (11*lmu)/ex(6) + pow(ep, -1)/ex(4) + lmu*pow(ep, -1) + (3*pow(lmu, 2))/ex(2) - (5*pow(Pi, 2))/ex(4)) + nL*Tf*(113/ex(24) + (19*lmu)/ex(6) - pow(ep, -2)/ex(2) + (11*pow(ep, -1))/ex(12) + pow(lmu, 2)/ex(2) + pow(Pi, 2)/ex(3)) + CF*(433/ex(32) + (51*lmu)/ex(8) + (9*pow(ep, -2))/ex(8) + (51*pow(ep, -1))/ex(16) + (9*lmu*pow(ep, -1))/ex(4) + (9*pow(lmu, 2))/ex(4) - (49*pow(Pi, 2))/ex(16) + 4*log(2)*pow(Pi, 2) - 6*zeta(3)) + CA*(-1705/ex(96) - (215*lmu)/ex(24) + (11*pow(ep, -2))/ex(8) - (127*pow(ep, -1))/ex(48) - (11*pow(lmu, 2))/ex(8) + (5*pow(Pi, 2))/ex(4) - 2*log(2)*pow(Pi, 2) + 3*zeta(3))))/ex(4);
            } else if(n==Symbol("q")) {
                if(loop>1) rr += (CF*nH*Tf*(-5/ex(24) + lmu/ex(2) + pow(ep, -1)/ex(4))*pow(as_pi, 2))/ex(4);
            }
            
            return rr;
        }
        
        /**
         * @brief the Zm renomalization constant
         * @param m the mass of the particle
         * @param loop the loop order
         * @return Zm at loop order
         */
        ex Zm(ex m, int loop) {
            if(loop==0) return 1;
            ex rr = 1;
            ex Tf = ex(1)/2;
            ex as_pi = as/Pi;
            ex CA = NF;
            ex CF = NA/(2*NF);
            ex lmu = 2*log(mu/m);
            if(loop>0) rr += (CF*(as_pi)*(-2 - 4*ep - (3*lmu)/ex(2) - 2*ep*lmu - (3*pow(ep, -1))/ex(2) - 8*pow(ep, 2) - 4*lmu*pow(ep, 2) - (3*ep*pow(lmu, 2))/ex(4) - pow(ep, 2)*pow(lmu, 2) - (pow(ep, 2)*pow(lmu, 3))/ex(4) - (ep*pow(Pi, 2))/ex(8) - (pow(ep, 2)*pow(Pi, 2))/ex(6) - (lmu*pow(ep, 2)*pow(Pi, 2))/ex(8) + (pow(ep, 2)*zeta(3))/ex(2)))/ex(2);
            if(loop>1) rr += (CF*pow(as_pi, 2)*(nH*Tf*(143/ex(24) + (13*lmu)/ex(6) - pow(ep, -2)/ex(2) + (5*pow(ep, -1))/ex(12) + pow(lmu, 2)/ex(2) - (2*pow(Pi, 2))/ex(3)) + nL*Tf*(71/ex(24) + (13*lmu)/ex(6) - pow(ep, -2)/ex(2) + (5*pow(ep, -1))/ex(12) + pow(lmu, 2)/ex(2) + pow(Pi, 2)/ex(3)) + CF*(199/ex(32) + (45*lmu)/ex(8) + (9*pow(ep, -2))/ex(8) + (45*pow(ep, -1))/ex(16) + (9*lmu*pow(ep, -1))/ex(4) + (9*pow(lmu, 2))/ex(4) - (17*pow(Pi, 2))/ex(16) + 2*log(2)*pow(Pi, 2) - 3*zeta(3)) + CA*(-1111/ex(96) - (185*lmu)/ex(24) + (11*pow(ep, -2))/ex(8) - (97*pow(ep, -1))/ex(48) - (11*pow(lmu, 2))/ex(8) + pow(Pi, 2)/ex(3) - log(2)*pow(Pi, 2) + (3*zeta(3))/2)))/ex(4);
            return rr;
        }
        
        /**
         * @brief the bare strong cupling constant up to loop order, MSbar schema
         * @param loop the loop order
         * @return the bare strong cupling constant
         */
        ex asBare(int loop) {
            auto as0 = asLO();
            ex as_pi = as/Pi;
            ex CA = NF;
            ex CF = NA/(2*NF);
            ex Tf = ex(1)/2;
            ex rr = 1;
            if(loop>0) rr += -(((11*CA)/ex(3) - (4*(nH + nL)*Tf)/ex(3))*(as_pi)*pow(ep, -1))/ex(4);
            if(loop>1) rr += ((-(((-20*CA*(nH + nL)*Tf)/ex(3) - 4*CF*(nH + nL)*Tf + (34*pow(CA, 2))/ex(3))*pow(ep, -1))/ex(8) + (pow(ep, -2)*pow((11*CA)/ex(3) - (4*(nH + nL)*Tf)/ex(3), 2))/ex(4))*pow(as_pi, 2))/ex(4);
            return as0 * rr;
        }
        
        /**
         * @brief the bare strong cupling constant at tree level, MSbar schema
         * @return the bare strong cupling constant at LO
         */
        ex asLO() {
            return as * pow(mu, 2*ep) * exp(ep*Euler)*pow(4*Pi, -ep);
        }
    }

}

