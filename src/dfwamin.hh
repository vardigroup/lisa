/*-------------------------------------------------------------------*/
// Yong Li (liyong@ios.ac.cn)
/*-------------------------------------------------------------------*/

#pragma once

#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <queue>
#include <unordered_map>
//#include <tuple>


#include <bddx.h>
#include <cudd.h>
#include <cuddObj.hh>

#include "dfwa.hh"
#include "dfwavar.hh"
#include "spotutil.hh"

//using namespace CUDD;
using namespace std;

// CUDD and BuDDy both use BDD, we need new names for CUDD::BDD
typedef CUDD::BDD cudd_bdd;

typedef CUDD::BDD& cudd_bdd_ptr;

typedef CUDD::Cudd* cudd_ptr;

typedef CUDD::Cudd cudd;

typedef DdNode* cudd_node_ptr;

typedef bdd buddy_bdd;

typedef bdd& buddy_bdd_ptr;

static unsigned nodes_num = 0;

typedef std::tuple<cudd_node_ptr, cudd_node_ptr, cudd_node_ptr> key_triple;
struct key_hash: public std::unary_function<key_triple, std::size_t>
{
	std::size_t operator()(const key_triple& k) const
	{
		return (long)get<0>(k) ^ (long)get<1>(k) ^ (long)get<2>(k);
	}
};

struct key_equal: public std::function<bool(const key_triple& k1, const key_triple& k2)>
{
	bool operator()(const key_triple& k1, const key_triple& k2) const
	{
		return get<0>(k1) == get<0>(k2)
			&& get<1>(k1) == get<1>(k2)
			&& get<2>(k1) == get<2>(k2);
	}
};

typedef unordered_map<const key_triple, cudd_bdd, key_hash, key_equal> map_t;

/**
 * The ordering of BDD variables in the DFA: (X, X'), A and then K (block variables)
 * the strategy for moving a BDD data structure from BuDDy to CUDD is as follows:
 * 		1. create variables A', X, X', A, K in CUDD, where A', X, X' are used in migration
 * 		2. renaming A' to A before doing DFA minimization
 * */

class dfwa_min
{
    public:
    	dfwa_min(cudd_ptr manager, dfwa_ptr aut)
    	:_manager(manager), _aut (aut)
    	{
    		prepare();
    	}
        // cube bdd
        cudd_bdd _curr_cube;
        cudd_bdd _next_cube;
        cudd_bdd _label_cube; // seems to be redundant
        
        // dfa (I(X) T(X, X', A), F(X))
        cudd_bdd _init;
        cudd_bdd _trans;
        cudd_bdd _finals;
        
        // first position after state variables 
        unsigned _num_states_vars; //_last_pos_states;
        // first position after label variables
        unsigned _last_pos_labels;
        
        unsigned _num_min_states;

        // state variables
        vector<vector<cudd_bdd>> _state_vars;
        // block variables
        vector<vector<cudd_bdd>> _block_vars;
        // label variables
        vector<cudd_bdd> _label_vars;
        
        // utility for mapping to and from BuDDy bdds
        vector<int> _cudd_to_buddy_vars; 
        unordered_map<int, int> _buddy_to_cudd_vars;

        dfwa_ptr _aut;

        // CUDD manager
        cudd_ptr _manager = nullptr;

        ~dfwa_min();
        
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
            return var_index < _num_states_vars;
        }

        // prepare for the minimization
        void prepare();
        void prepare_variables();        
        
        // move BDD to a different BDD manager
        cudd_bdd move_to_cudd(buddy_bdd dd);
        cudd_bdd move_to_cudd(buddy_bdd dd, unordered_map<int, cudd_bdd>& computed_table);

        // move BDD from CUDD to BuDDy
        buddy_bdd move_to_buddy(cudd_bdd dd, unordered_map<int, int>& cudd_to_buddy_map);
        buddy_bdd move_to_buddy(cudd_bdd_ptr dd, unordered_map<int, int>& cudd_to_buddy_map
        		, unordered_map<cudd_node_ptr, buddy_bdd>& computed_table);

        // bdd_replace function for cudd
        cudd_bdd cudd_permute(cudd_bdd_ptr dd, vector<cudd_bdd>& from, vector<cudd_bdd>& to);
        
        // helper functions for minimization
        void new_block_variables(unsigned num = 1);
        cudd_bdd new_block_number(int block_number, int copy);
        cudd_bdd initialize_partition();
        pair<cudd_bdd, unsigned> refine_partition(cudd_bdd_ptr sig);
        cudd_bdd refine_partition_rec(cudd_bdd_ptr dd, unsigned & block_number, unordered_map<cudd_node_ptr, cudd_bdd> &computed_table);
        void compute_block_number_rec(cudd_bdd_ptr sig, unordered_map<cudd_node_ptr, int>& table, int & block_number);
        void reduce_block_variables(cudd_bdd_ptr partition);
        
        cudd_bdd make_quotient_trans(cudd_bdd_ptr trans, cudd_bdd_ptr curr_partition, cudd_bdd_ptr next_partition);
        cudd_bdd make_quotient_trans(cudd_bdd trans, cudd_bdd curr_partition, cudd_bdd next_partition, vector<set<unsigned>>& var_sets, map_t& computed_table);

        // print functions
        void print(cudd_bdd bdd);
        void generate_all_bits(vector<cudd_bdd>& vars, unsigned index, cudd_bdd_ptr dd, cudd_bdd temp);

        // state space exploration
        cudd_bdd next_image(cudd_bdd_ptr curr);
        cudd_bdd pre_image(cudd_bdd_ptr curr);
        cudd_bdd forward_explore();
        cudd_bdd backward_explore();
};
