#pragma once

#include <vector>
#include <string>
#include <algorithm>

#include <spot/twaalgos/hoa.hh>

#include <spot/twa/twagraph.hh>

#include <spot/misc/bddlt.hh>

#include "dfwa.hh"
#include "dfwavar.hh"
#include "strategy.hh"

using namespace std;
using namespace spot;

enum class synt_result { YES, NO, UNKNOWN };

/*
 * (X, Y, Z)
 *   X is the set of input variables (_in_cube)
 *   Y is the set of output variables (_out_cube)
 *   Z is the set of state variables
 * */
const char* const REALIZABLE = "REALIZABLE";
const char* const UNREALIZABLE = "UNREALIZABLE";

class synt
{
    private:
        void prepare();
        bdd pre_image(bdd& w);
        
        // output transducer
        bdd _mealy_init;
        bdd _mealy_trans;
        bdd _mealy_out;
        bool _env_first;

        bool eval(bdd& dd, bdd& assign);
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
        dfwa_ptr _dfa;
        
        synt(dfwa_ptr dfa, bdd& input_cube, bdd& output_cube)
        : _dfa(dfa), _in_cube(input_cube), _out_cube(output_cube)
        {
            prepare();
        }
        
        ~synt();
        
        void is_realizable();
        
        void synthesize();

        void env_play_first();
};
