#include "WSKernel.h"

namespace HepLib {
  
    const char * WSKernel::Error::what() const throw () {
        return msg.c_str();
    }
    
    WSKernel::Error::Error(WSLINK lp) {
        const char *err_msg = WSErrorMessage(lp);
        msg = err_msg;
        if(err_msg!=NULL) WSReleaseErrorMessage(lp, err_msg);
    };

    WSKernel::WSKernel(const string & open_str) {
        int error;
        ep = WSInitialize((WSParametersPointer)0);
        if(ep == (WSENV)0) throw Error(lp);
        
        string ostr = open_str;
        if(ostr.length()<1) {
            #ifdef __APPLE__
            ostr = "-linklaunch -linkname '/Applications/Mathematica.app/Contents/MacOS/WolframKernel -wstp'";
            #else
            ostr = "-linklaunch -linkname 'math -wstp'";
            #endif
        }
        
        lp = WSOpenString(ep, ostr.c_str(), &error);
        if(lp == (WSLINK)0 || error != WSEOK) throw Error(lp);
        Evaluate("Off[LinkConnect::linkc];");
    }

    WSKernel::~WSKernel() {
        WSClose(lp);
        WSDeinitialize(ep);
    }

    string WSKernel::Evaluate(const string & expr, const string & OutputForm) {
        WSPutFunction(lp, "EvaluatePacket", 1);
        WSPutFunction(lp, "ToString", 1);
        WSPutFunction(lp, OutputForm.c_str(), 1);
        WSPutFunction(lp, "ToExpression", 1);
        WSPutString(lp, expr.c_str());
        WSEndPacket(lp);
        WSFlush(lp);

        int pkt;
        while( (pkt=WSNextPacket(lp)) && (pkt!= RETURNPKT) ) WSNewPacket(lp);

        const char *mma_out;
        if(!WSGetString(lp, &mma_out)) throw Error(lp);
        
        string ret_string = mma_out;
        if(!WSEndPacket(lp)) throw Error(lp);
        WSReleaseString(lp, mma_out);
        WSFlush(lp);
        return ret_string;
    }

}
