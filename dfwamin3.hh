#pragma once

#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <queue>
#include <unordered_map>
//#include <tuple>


#include <bddx.h>
#include <sylvan.h>
#include <sylvan_obj.hpp>

#include "dfwa.hh"
#include "dfwavar.hh"
#include "spotutil.hh"

//using namespace CUDD;
using namespace std;
using namespace sylvan;

// this file will use sylvan

// CUDD and BuDDy both use BDD, we need new names for CUDD::BDD

typedef bdd buddy_bdd;

typedef bdd& buddy_bdd_ptr;

typedef Bdd sylvan_bdd;
typedef BddSet sylvan_bdd_set;
typedef Bdd& sylvan_bdd_ptr;

typedef sylvan::BDD sylvan_BDD;


//static unsigned nodes_num = 0;


/**
 * The ordering of BDD variables in the DFA: (X, X'), A and then K (block variables)
 * the strategy for moving a BDD data structure from BuDDy to CUDD is as follows:
 * 		1. create variables A', X, X', A, K in CUDD, where A', X, X' are used in migration
 * 		2. renaming A' to A before doing DFA minimization
 * */

typedef std::tuple<sylvan_BDD, sylvan_BDD, sylvan_BDD> key_triple_sylvan;
struct key_hash_sylvan: public std::unary_function<key_triple_sylvan, std::size_t>
{
	std::size_t operator()(const key_triple_sylvan& k) const
	{
		return (long)get<0>(k) ^ (long)get<1>(k) ^ (long)get<2>(k);
	}
};

struct key_equal_sylvan: public std::function<bool(const key_triple_sylvan& k1, const key_triple_sylvan& k2)>
{
	bool operator()(const key_triple_sylvan& k1, const key_triple_sylvan& k2) const
	{
		return get<0>(k1) == get<0>(k2)
			&& get<1>(k1) == get<1>(k2)
			&& get<2>(k1) == get<2>(k2);
	}
};

typedef unordered_map<const key_triple_sylvan, sylvan_bdd, key_hash_sylvan, key_equal_sylvan> map_t_sylvan;

class dfwa_min_sylvan
{
    public:
		dfwa_min_sylvan(dfwa_ptr aut)
    	: _aut (aut)
    	{
    		prepare();
    	}

		uint32_t _num_vars;
        // cube bdd
		sylvan_bdd_set _curr_cube;
		sylvan_bdd_set _next_cube;
		sylvan_bdd_set _label_cube;
		sylvan_bdd_set _curr_label_cube;
		sylvan_bdd_set _next_label_cube;
        
        // dfa (I(X) T(X, X', A), F(X))
        sylvan_bdd _init;
        sylvan_bdd _trans;
        sylvan_bdd _finals;
        
        // first position after state variables 
        uint32_t _num_states_vars; //_last_pos_states;
        // first position after label variables
        unsigned _last_pos_labels;
        
        unsigned _num_min_states;

        // state variables
        vector<vector<uint32_t>> _state_vars;
        // block variables
        vector<vector<uint32_t>> _block_vars;
        // label variables
        vector<uint32_t> _label_vars;
        
        // utility for mapping to and from BuDDy bdds
        vector<int> _sylvan_to_buddy_vars; 
        unordered_map<int, int> _buddy_to_sylvan_vars;

        dfwa_ptr _aut;

        // CUDD manager

        ~dfwa_min_sylvan();
        
        void minimize();
        
        dfwa_ptr move_dfwa();

        //dfwa& operator &(dfwa& other);
        void output(ostream& os);

        unsigned get_num_min_states()
        {
        	return _num_min_states;
        }
    private:
    
        bool is_state_variable(const uint32_t var_index)
        {
            return var_index < _num_states_vars;
        }

        // prepare for the minimization
        void prepare();
        void prepare_variables();        
        
        // move BDD to a different BDD manager
        sylvan_bdd move_to_sylvan(buddy_bdd dd);
        sylvan_bdd move_to_cudd(buddy_bdd dd, unordered_map<int, sylvan_bdd>& computed_table);

        // move BDD from CUDD to BuDDy
        buddy_bdd move_to_buddy(sylvan_bdd dd, unordered_map<int, int>& cudd_to_buddy_map);
        buddy_bdd move_to_buddy(sylvan_bdd_ptr dd, unordered_map<int, int>& cudd_to_buddy_map
        		, unordered_map<sylvan_BDD, buddy_bdd>& computed_table);

        // bdd_replace function for cudd
        //sylvan_bdd cudd_permute(sylvan_bdd_ptr dd, vector<sylvan_bdd>& from, vector<sylvan_bdd>& to);
        
        // helper functions for minimization
        void new_block_variables(unsigned num = 1);
        sylvan_bdd new_block_number(int block_number, int copy);
        sylvan_bdd initialize_partition();
        pair<sylvan_bdd, unsigned> refine_partition(sylvan_bdd_ptr sig);
        sylvan_bdd refine_partition_rec(sylvan_bdd_ptr dd, unsigned & block_number, unordered_map<sylvan_BDD, sylvan_bdd> &computed_table);
        void compute_block_number_rec(sylvan_bdd_ptr sig, unordered_map<sylvan_BDD, int>& table, int & block_number);
        void reduce_block_variables(sylvan_bdd_ptr partition);
        
        sylvan_bdd make_quotient_trans(sylvan_bdd_ptr trans, sylvan_bdd_ptr curr_partition, sylvan_bdd_ptr next_partition);
        sylvan_bdd make_quotient_trans(sylvan_bdd trans, sylvan_bdd curr_partition, sylvan_bdd next_partition, vector<set<uint32_t>>& var_sets, map_t_sylvan& computed_table);

        // print functions
        void print(sylvan_bdd bdd);
        void generate_all_bits(vector<sylvan_bdd>& vars, unsigned index, sylvan_bdd_ptr dd, sylvan_bdd temp);

        // state space exploration
        sylvan_bdd next_image(sylvan_bdd_ptr curr);
        sylvan_bdd pre_image(sylvan_bdd_ptr curr);
        sylvan_bdd forward_explore();
        sylvan_bdd backward_explore();
};
