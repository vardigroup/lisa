#pragma once

#include <bddx.h>
#include <unordered_map>

#include "spotutil.hh"

/**
 * here we construct a strategy represented by a bdd such that
 *  state s needs to give output o iff s /\ o = 1 in the bdd
 *  */

class strategy
{
private:

	bdd& _states_and_outputs;
	//buddy_bdd _output_cube;
	int _max_state_var = -1;
	bdd _anony_cube;
	bddPair* _label_to_anony_pair;
	bddPair* _anony_to_label_pair;

	bool is_state_variable(int index)
	{
		return index < _max_state_var;
	}

	//buddy_bdd synthesize_strategy_rec(buddy_bdd_ptr s_and_o, int index);
	bdd synthesize_strategy_rec(bdd& s_and_o, unordered_map<int, bdd>& computed_table);

public:

	strategy(bdd& s_and_o, bdd output_cube, bdd state_cube);
	~strategy();

	bdd synthesize_strategy();

};

