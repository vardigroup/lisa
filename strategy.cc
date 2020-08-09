/*-------------------------------------------------------------------*/
// Yong Li (liyong@ios.ac.cn)
/*-------------------------------------------------------------------*/

#include "strategy.hh"

// NOTE THAT the ordering of the variables is A, X, X'
// first traverse the possible values of
// if the ordering is X, A
// then get bdd_satone when first get to the label variables

strategy::strategy(bdd& s_and_o, bdd output_cube, bdd state_cube)
	: _states_and_outputs(s_and_o)
{
	vector<int> list_label_vars;
	get_list_var_indices(list_label_vars, output_cube);
	vector<int> list_state_vars;
	get_list_var_indices(list_state_vars, state_cube);
	// largest state variable occuring in _states_and_outputs
	_max_state_var = list_state_vars.back() + 1;
	_label_to_anony_pair = bdd_newpair();
	_anony_to_label_pair = bdd_newpair();
	// compute pairs
	_anony_cube = bddtrue;
	for(unsigned i = 0; i < list_label_vars.size(); i ++)
	{
		int anony_index = _max_state_var + i;
		bdd_setpair(_label_to_anony_pair, list_label_vars[i], anony_index);
		bdd_setpair(_anony_to_label_pair, anony_index, list_label_vars[i]);
		_anony_cube = _anony_cube & bdd_ithvar(anony_index);
	}
}

strategy::~strategy()
{
	if(_label_to_anony_pair != nullptr)
	{
		bdd_freepair(_label_to_anony_pair);
		_label_to_anony_pair = nullptr;
	}
	if(_anony_to_label_pair != nullptr)
	{
		bdd_freepair(_anony_to_label_pair);
		_anony_to_label_pair = nullptr;
	}
}

bdd
strategy::synthesize_strategy()
{
	bdd new_state_output = bdd_replace(_states_and_outputs, _label_to_anony_pair);
	unordered_map<int, bdd> computed_table;
	bdd func_output = synthesize_strategy_rec(new_state_output, computed_table);
	func_output = bdd_replace(func_output, _anony_to_label_pair);
	return func_output;
}

/*
buddy_bdd
strategy::synthesize_strategy_rec(buddy_bdd_ptr s_and_o, int index)
{
	return bddtrue;
}
*/
// if the ordering is X, A
// then get bdd_satone when first get to the label variables
bdd
strategy::synthesize_strategy_rec(bdd& s_and_o
		, unordered_map<int, bdd>& computed_table)
{
	if(s_and_o == bddfalse)
	{
		return bddfalse;
	}
	// check whether it has been computed
	int dd_id = s_and_o.id();
	unordered_map<int, bdd>::const_iterator it = computed_table.find(dd_id);
	if(it != computed_table.end())
	{
		return it->second;
	}
	// now check if it is true

	if(s_and_o == bddtrue )
	{
		// choose one sat
		bdd set = bdd_satoneset(s_and_o, _anony_cube, bddtrue);
		computed_table[dd_id] = set;
		return set;
	}
	bdd result;
	int top_var = bdd_var(s_and_o);
	if(is_state_variable(top_var))
	{
		// now do it recursively
		bdd low_branch = bdd_low(s_and_o);
		bdd low = synthesize_strategy_rec(low_branch, computed_table);
		bdd high_branch = bdd_high(s_and_o);
		bdd high = synthesize_strategy_rec(high_branch, computed_table);
		bdd var_dd = bdd_ithvar(top_var);
		result = bdd_ite(var_dd, high, low);
	}else
	{
		// reach non state variables
		result = bdd_satoneset(s_and_o, _anony_cube, bddtrue);
	}
	computed_table[dd_id] = result;
	return result;
}
