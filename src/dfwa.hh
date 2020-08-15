/*-------------------------------------------------------------------*/
// Yong Li (liyong@ios.ac.cn)
/*-------------------------------------------------------------------*/

#pragma once

#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <queue>
#include <functional>

#include <spot/twa/twagraph.hh>
#include <spot/twa/bddprint.hh>
#include <spot/misc/bddlt.hh>

#include "dfwavar.hh"
#include "debug.hh"
#include "spotutil.hh"


// needs to input as pointers
class dfwa
{
    public:
        twa_graph_ptr _aut;
        dfwa_var _state_vars;
        bdd _curr_cube;
        bdd _next_cube;
        bdd _label_cube;
        bddPair* _curr_to_next_pairs;
        bddPair* _next_to_curr_pairs;
        
        bdd _init;
        bdd _trans;
        bdd _finals;
        bdd _reach;
        // label_cube should not be a pointer
        dfwa(bdd_dict_ptr dict, bdd label_cube);
        dfwa(twa_graph_ptr aut, bdd label_cube, set<unsigned>& finals, const char* name = "s");
        ~dfwa();
        
        bdd next_image(bdd curr);
        
        bdd pre_image(bdd curr);
        
        bdd get_init();
        
        bdd get_trans();
        
        bdd get_finals();
        
        bdd explore();
        bdd back_explore();
        
        bool is_empty();

        bdd_dict_ptr get_dict()
        {
        	return _state_vars.get_dict();
        }

        //dfwa& operator &(dfwa& other);
        void output(ostream& os);

        void output_dfwa(ostream& os);

        void make_complete();

};

//ostream& operator<<(ostream& os, const dfwa& dfa);
typedef dfwa& dfwa_ptr;

//typedef dfwa* dfwa_pointer;


void
intersect_dfwa(dfwa_ptr result, dfwa_ptr op1, dfwa_ptr op2);
dfwa_ptr
product_dfwa_and(dfwa_ptr op1, dfwa_ptr op2);
dfwa_ptr
product_dfwa_or(dfwa_ptr op1, dfwa_ptr op2);
dfwa_ptr
product_dfwa_minus(dfwa_ptr op1, dfwa_ptr op2);
