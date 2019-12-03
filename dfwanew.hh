#pragma once

#include <vector>
#include <string>
#include <set>
#include <iostream>

#include <spot/twa/twagraph.hh>
#include <spot/twa/bddprint.hh>
#include <spot/misc/bddlt.hh>

#include "dfwavar.hh"
#include "spotutil.hh"

// this class only uses one copy of the state variables
// see the paper titled "Symbolic LTLf Synthesis" about this
// state needs to start from 1, since every state needs at least one bit

class dfwa_new
{
    public:
        twa_graph_ptr _aut;
        dfwa_var _state_vars;
        bdd _curr_cube;
        bdd _label_cube;
        bdd_dict_ptr _dict;
        //bdd _reach;
        
        //transitions
        bddPair* _curr_to_bdd_pairs;
        bdd _init;
        bdd _finals;
        
        // probably not needed
        vector<bdd> _bitseq;

        // product
        dfwa_new(bdd_dict_ptr dict, bdd& label_cube);
        // from explicit graph
        dfwa_new(twa_graph_ptr aut, bdd& label_cube, set<unsigned>& finals, const char* name = "s");
        virtual ~dfwa_new();
        
        //bdd next_image(bdd& curr);
        
        bdd pre_image(bdd& curr);
        
        bdd get_init()
        {
            return _init;
        }
        
        bdd_dict_ptr get_dict()
        {
        	return _dict;
        }

        unsigned get_seq_size()
        {
            return _bitseq.size();
        }
        
        bdd get_bit_seq(unsigned id)
        {
            return _bitseq[id]; 
        }
        
        bdd get_finals()
        {
            return _finals;
        }
        
        bdd back_explore();

        //bdd explore();
        
        //dfwa& operator &(dfwa& other);
        void output(ostream& os);
        
};
typedef dfwa_new& dfwa_new_ref;
typedef dfwa_new* dfwa_new_ptr;
//ostream& operator<<(ostream& os, const dfwa& dfa);

void intersect_dfwan(dfwa_new_ref result, dfwa_new_ref op1, dfwa_new_ref op2);

dfwa_new_ref
product_dfwa_new_and(dfwa_new_ref op1, dfwa_new_ref op2);
dfwa_new_ref
product_dfwa_new_or(dfwa_new_ref op1, dfwa_new_ref op2);
dfwa_new_ref
product_dfwa_new_minus(dfwa_new_ref op1, dfwa_new_ref op2);
