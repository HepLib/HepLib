#include "SD.h"
#include "mpreal.h"

using namespace HepLib;

void let_op_append(ex & ex_in, int index, ex const item) {
    auto tmp = ex_to<lst>(ex_in.op(index));
    tmp.append(item);
    ex_in.let_op(index) = tmp;
}
void let_op_append(lst & ex_in, int index, ex const item) {
    auto tmp = ex_to<lst>(ex_in.op(index));
    tmp.append(item);
    ex_in.let_op(index) = tmp;
}

void let_op_prepend(ex & ex_in, int index, ex const item) {
    auto tmp = ex_to<lst>(ex_in.op(index));
    tmp.prepend(item);
    ex_in.let_op(index) = tmp;
}
void let_op_prepend(lst & ex_in, int index, ex const item) {
    auto tmp = ex_to<lst>(ex_in.op(index));
    tmp.prepend(item);
    ex_in.let_op(index) = tmp;
}

void let_op_remove_last(ex & ex_in, int index) {
    auto tmp = ex_to<lst>(ex_in.op(index));
    tmp.remove_last();
    ex_in.let_op(index) = tmp;
}
void let_op_remove_last(lst & ex_in, int index) {
    auto tmp = ex_to<lst>(ex_in.op(index));
    tmp.remove_last();
    ex_in.let_op(index) = tmp;
}

void let_op_remove_first(ex & ex_in, int index) {
    auto tmp = ex_to<lst>(ex_in.op(index));
    tmp.remove_first();
    ex_in.let_op(index) = tmp;
}
void let_op_remove_first(lst & ex_in, int index) {
    auto tmp = ex_to<lst>(ex_in.op(index));
    tmp.remove_first();
    ex_in.let_op(index) = tmp;
}

int main(int argc, char** argv) {

    //SD::debug = true;

    auto ep = SD::ep;
    
    if(true) {
        
        lst olst = lst{ lst{x(1),2}, lst{3,x(4)}};
        cout << olst << endl;
        
        exset xset;
        ex tmp = olst;
        find(olst, x(wild()), xset);
        
        cout << xset << endl;
    }
    
    
    
    return 0;
}
