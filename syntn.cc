#ifdef DEBUG

#include "syntn.h"

/*-------------------------------------------------------------------*/
// prepare for realizability checking
/*-------------------------------------------------------------------*/
void
synt::prepare()
{
    // first add finals
    _tseq.push_back(_dfa.get_finals());
    _wseq.push_back(_dfa.get_finals());
    //remove input and output variables
    //_rest_cube = bdd_exist(_dfa._label_cube, _in_cube & _out_cube);
    _result = synt_result::UNKNOWN;
}

synt::~synt()
{
}

/*-------------------------------------------------------------------*/
// compute pre image
/*-------------------------------------------------------------------*/
bdd
synt::pre_image(bdd & curr)
{
    bdd next = bdd_replace(curr, _dfa._curr_to_next_pairs);
    // has label_cube
	bdd pre = bdd_relprod(_dfa._trans, next, _dfa._next_cube);
    // only allow reachable states
    pre = pre & _dfa._reach;
    return pre;
}

/*-------------------------------------------------------------------*/
// fixed point computation for checking realizability
/*-------------------------------------------------------------------*/
void 
synt::is_realizable()
{
    unsigned curr = 0;
    while(true)
    {
        // t_{i + 1}(Z, Y) = t_{i}(Z, Y) | (!w_{i}(Z) &  ∀ X. (tr(X, Y, (w_{i})(Z))))
        // pre_image (X, Y, Z)
        bdd pre_w = pre_image(_wseq[curr]);
        /*
        cout << "pre image of curr" << endl;
        bdd_print_set(cout, _dfa._state_vars.get_dict(), pre_w);
        cout << endl;
        */
        pre_w = bdd_forall(pre_w, _in_cube);
        /*
        cout << "bdd_foall pre_w " << endl;
        bdd_print_set(cout, _dfa._state_vars.get_dict(), pre_w);
        cout << endl;
        */
        // now pre_w has only Z and Y
        bdd tp = _tseq[curr] | ((!_wseq[curr]) & pre_w);
        /*
        cout << "bdd_foall tp " << endl;
        bdd_print_set(cout, _dfa._state_vars.get_dict(), tp);
        cout << endl;
        */
        _tseq.push_back(tp);
        // w_{i+1}(Z) = ∃Y.ti+1(Z, Y)
        bdd wp = bdd_exist(tp, _out_cube);
        /*
        cout << "bdd_foall wp " << endl;
        bdd_print_set(cout, _dfa._state_vars.get_dict(), wp);
        cout << endl;
        */
        _wseq.push_back(wp);
        ++ curr;
        if(_wseq[curr] == _wseq[curr - 1])
        {
            break;
        }
    }
    bdd overlap = _wseq[curr] & _dfa._init;
    if(overlap != bddfalse)
    {
        _result = synt_result::YES;
        cout << "REALIZABLE" << endl;
    }else
    {
        _result = synt_result::NO;
        cout << "UNREALIZABLE" << endl;
    }
}
/*-------------------------------------------------------------------*/
// synthesize winning strategy
/*-------------------------------------------------------------------*/        
void 
synt::synthesize()
{
    
}

#endif
