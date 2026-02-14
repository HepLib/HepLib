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
            ostr = "-linkmode launch -linkname '/Applications/Mathematica.app/Contents/MacOS/WolframKernel -wstp'";
            #else
            ostr = "-linkmode launch -linkname 'math -wstp'";
            #endif
        }
        
        lp = WSOpenString(ep, ostr.c_str(), &error);
        //lp = WSOpenString(ep, "-linkcreate -linkprotocol IntraProcess", &error);
        //lp = WSOpenArgcArgv(ep, argc, argv, &error);
        if(lp == (WSLINK)0 || error != WSEOK) throw Error(lp);
    }

    WSKernel::~WSKernel() {
        WSClose(lp);
        WSDeinitialize(ep);
    }

    string WSKernel::Evaluate(const string & expr) {
        WSPutFunction(lp, "EvaluatePacket", 1);
        WSPutFunction(lp, "ToString", 1);
        WSPutFunction(lp, "InputForm", 1);
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
