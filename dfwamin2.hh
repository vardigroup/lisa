#pragma once

#ifndef __ORIGINAL__
#define __ORIGINAL__ 1 // set original mode
#endif

#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <queue>
#include <unordered_map>
//#include <tuple>


#include "common.hh"

#include "dfwa.hh"
#include "dfwavar.hh"
#include "spotutil.hh"

//using namespace CUDD;
using namespace std;

//static unsigned nodes_num = 0;

/**
 * The input DFA has the ordering of A, (X, X')
 * The ordering of BDD variables in the DFA minimization is A, (X, X'), A' and then K (block variables)
 * the strategy for changing the order of a BDD data structure in BuDDy is as follows:
 * 		renaming A to A' in T(A, (X, X')) before doing DFA minimization
 *
 * CUDD is good at renaming and building BDDs, and BuDDy is good at state space exploration
 * */

class dfwa_min_bdd
{
    public:
    	dfwa_min_bdd(dfwa_ptr aut);
        // cube bdd
        buddy_bdd _curr_cube;
        buddy_bdd _next_cube;
        buddy_bdd _label_cube;
        
        // dfa (I(X) T(X, X', A), F(X))
        buddy_bdd _init;
        buddy_bdd _trans;
        buddy_bdd _finals;
        
        // first position after state variables 
        int _max_states_var = -1; //_last_pos_states;
        // first position after label variables
        unsigned _last_pos_labels;
        
        unsigned _num_min_states;

#ifdef __ORIGINAL__
        // state to anonymous
        bddPair* _state_to_anony_pair = nullptr;
        bddPair* _curr_to_next_state_pair = nullptr;
        bddPair* _next_to_curr_state_pair = nullptr;
        unordered_map<int, int> _map;
        //unordered_map<int, int> _anony_to_label;

#endif
        // state variables
        vector<vector<buddy_bdd>> _state_vars;
        // block variables
        vector<vector<buddy_bdd>> _block_vars;
        bddPair* _curr_to_next_block_pair = nullptr;
        bddPair* _next_to_curr_block_pair = nullptr;

        // new to label
        bddPair* _label_to_anony_pair = nullptr;
        bddPair* _anony_to_label_pair = nullptr;
        // label variables
        vector<buddy_bdd> _label_vars;
        
        dfwa_ptr _aut;

        ~dfwa_min_bdd();
        
        void minimize();
        
        dfwa_ptr move_dfwa();

        //dfwa& operator &(dfwa& other);
        void output(ostream& os);

        unsigned get_num_min_states()
        {
        	return _num_min_states;
        }
    private:
    
        bool is_state_variable(const unsigned var_index)
        {
            return var_index < _max_states_var;
        }

        // prepare for the minimization
        void prepare();
        void prepare_variables();

        // bdd_replace function for cudd
        buddy_bdd move_to(buddy_bdd dd);
        buddy_bdd move_to(buddy_bdd dd, unordered_map<int, buddy_bdd>& computed_table);

        // helper functions for minimization
        void new_block_variables(unsigned num = 1);
        buddy_bdd new_block_number(int block_number, int copy);
        buddy_bdd initialize_partition();
        pair<buddy_bdd, unsigned> refine_partition(buddy_bdd_ptr sig);
        buddy_bdd refine_partition_rec(buddy_bdd_ptr dd, unsigned & block_number, unordered_map<int, buddy_bdd> &computed_table);
        void compute_block_number_rec(buddy_bdd_ptr sig, unordered_map<cudd_node_ptr, int>& table, int & block_number);
        void reduce_block_variables(buddy_bdd_ptr partition);
        
        // print functions
        void print(buddy_bdd bdd);
        void generate_all_bits(vector<buddy_bdd>& vars, unsigned index, buddy_bdd_ptr dd, buddy_bdd temp);

        // state space exploration
        buddy_bdd next_image(buddy_bdd_ptr curr);
        buddy_bdd pre_image(buddy_bdd_ptr curr);
        buddy_bdd forward_explore();
        buddy_bdd backward_explore();
};
