#pragma once

#ifdef DEBUG
#include <vector>
#include <string>
#include <algorithm>

#include <spot/twaalgos/hoa.hh>

#include <spot/twa/twagraph.hh>

#include <spot/misc/bddlt.hh>

#include "dfwan.hh"
#include "dfwavar.h"

using namespace std;
using namespace spot;

enum class synt_result { YES, NO, UNKNOWN };

/*
 * (X, Y, Z)
 *   X is the set of input variables (_in_cube)
 *   Y is the set of output variables (_out_cube)
 *   Z is the set of state variables
 * */
class syntn
{
    private:
        void prepare();
        bdd pre_image(bdd & w);
    public:
        unsigned _copies;
        // (Z, Y) is the sequence of winning states and outputs
        vector<bdd> _tseq;
        // (Z) is the sequence of winning states
        vector<bdd> _wseq;
        // input
        bdd _in_cube;
        // output
        bdd _out_cube;
        // rest of propositions
        //bdd _rest_cube;
        bdd_dict_ptr _dict;
        
        synt_result _result;
        dfwan_ptr _dfa;
        
        syntn(dfwan_ptr dfa, bdd& input_cube, bdd& output_cube)
        : _dfa(dfa), _in_cube(input_cube), _out_cube(output_cube)
        {
            prepare();
        }
        
        ~syntn();
        
        void is_realizable();
        
        void synthesize();
};
#endif
