/*-------------------------------------------------------------------*/
// Yong Li (liyong@ios.ac.cn)
/*-------------------------------------------------------------------*/

#include "dfwamin2.hh"

dfwa_min_bdd::dfwa_min_bdd(dfwa_ptr aut)
    	: _aut (aut)
{
	//bdd_gbc();
	prepare();
}

/*-------------------------------------------------------------------*/
// prepare state, label and block variables for DFA minimization
/*-------------------------------------------------------------------*/
void
dfwa_min_bdd::prepare_variables()
{
	_label_to_anony_pair = bdd_newpair();
	_anony_to_label_pair = bdd_newpair();
	_curr_to_next_block_pair = bdd_newpair();
	_next_to_curr_block_pair = bdd_newpair();
    // 0. create duplicate label variables for dfa minimization
	vector<int> state_vars;
	get_list_var_indices(state_vars, _aut._curr_cube & _aut._next_cube);
    // A' duplicate variables for labels
	//unsigned num_state_var = _aut._state_vars.get_bdd_vars(0);
#ifdef __ORIGINAL__
	_curr_to_next_state_pair = bdd_newpair();
	_next_to_curr_state_pair = bdd_newpair();
	_state_to_anony_pair = bdd_newpair();

	vector<buddy_bdd> curr_state_vars = _aut._state_vars.get_bdd_vars(0);
	vector<buddy_bdd> next_state_vars = _aut._state_vars.get_bdd_vars(1);
    //get_list_var_indices(next_state_vars, _aut._next_cube );
    _curr_cube = bddtrue;
    _next_cube = bddtrue;
    // start of the state index
    int buddy_index = 0; // state_vars.back() + 1
    for(unsigned i = 0; i < curr_state_vars.size(); i ++)
    {
    	int curr_index = bdd_var(curr_state_vars[i]);
    	bdd_setpair(_state_to_anony_pair, curr_index, buddy_index);
    	_map[curr_index] = buddy_index;
    	//_map[buddy_index] = buddy_index;
    	_curr_cube = _curr_cube & bdd_ithvar(buddy_index);
    	++ buddy_index;
    	int next_index = bdd_var(next_state_vars[i]);
    	bdd_setpair(_state_to_anony_pair, next_index, buddy_index);
    	_map[next_index] = buddy_index;
    	//_map[buddy_index] = buddy_index;
    	_next_cube = _next_cube & bdd_ithvar(buddy_index);
    	bdd_setpair(_curr_to_next_state_pair, buddy_index - 1, buddy_index);
    	bdd_setpair(_next_to_curr_state_pair, buddy_index, buddy_index - 1);
    	++ buddy_index;
    }
    cout << "Number of state variables is " << buddy_index << endl;
	// the first position for label variables
	_max_states_var = buddy_index + 1;
#else
	_curr_cube = _aut._curr_cube;
	_next_cube = _aut._next_cube;
	_max_states_var = state_vars.back() + 1;
#endif
    vector<int> label_vars;
    get_list_var_indices(label_vars, _aut._label_cube);
    // create anonymous variables for labels
    _label_cube = bddtrue;
    for(unsigned i = 0; i < label_vars.size(); i ++)
    {
    	int buddy_index = _max_states_var + i;
    	buddy_bdd var_dd = bdd_ithvar(buddy_index);
    	_label_cube = _label_cube & var_dd;
    	_map[label_vars[i]] = buddy_index;
    	//_map[label_vars[i]] = label_vars[i];
    	bdd_setpair(_label_to_anony_pair, label_vars[i], buddy_index);
    	bdd_setpair(_anony_to_label_pair, buddy_index, label_vars[i]);
    	//_anony_to_label[buddy_index] = label_vars[i];
    }
    _last_pos_labels = _max_states_var + label_vars.size();
    
    // 3. create block vectors
    // K block variables
    vector<buddy_bdd> curr_block;
    _block_vars.push_back(curr_block);
    vector<buddy_bdd> next_block;
    _block_vars.push_back(next_block);
    int num_block_var = 2 * _aut._state_vars.get_var_num(0);
    for(int i = 0; i < num_block_var; i +=2)
    {
    	int buddy_index = _last_pos_labels + i;
        //cout << "block variable " << var_dd << endl;
    	buddy_bdd var_dd = bdd_ithvar(buddy_index);
    	_block_vars[0].push_back(var_dd);
    	var_dd = bdd_ithvar(buddy_index + 1);
    	_block_vars[1].push_back(var_dd);
    	bdd_setpair(_curr_to_next_block_pair, buddy_index, buddy_index + 1);
    	bdd_setpair(_next_to_curr_block_pair, buddy_index + 1, buddy_index);
    }
}

/*-------------------------------------------------------------------*/
// prepare CUDD DFA minimization
/*-------------------------------------------------------------------*/
void
dfwa_min_bdd::prepare()
{
    // create variables of dfwa_min_bdd::
	cout << "minimized before prepare " << _aut._label_cube << endl;
    prepare_variables();
    // copy the initial state
    clock_t c_start = clock();
#ifdef __ORIGINAL__
    _init = bdd_replace(_aut._init, _state_to_anony_pair);
#else
    _init = _aut._init;
#endif
    cout << "minimized init: " << _init << endl;
    int count = bdd_nodecount(_aut._trans);
    cout << "The number of nodes in transition after reach is " << bdd_nodecount(_aut._trans) << endl;
#ifdef __ORIGINAL__
    //bddPair* result_pair = bdd_mergepairs(_state_to_anony_pair, _label_to_anony_pair);
    //_trans = bdd_replace(_aut._trans, result_pair);
    _trans = move_to(_aut._trans);
    //bdd_freepair(result_pair);

#else
    _trans = bdd_replace(_aut._trans, _label_to_anony_pair);
#endif
    cout << "The number of nodes in transition after renaming is " << bdd_nodecount(_trans) << endl;
    //_trans = bdd_replace(_trans, _state_to_anony_pair);
    // copy the final states
    //cout << "minimized trans: " << _trans << endl;
#ifdef __ORIGINAL__
    _finals = bdd_replace(_aut._finals, _state_to_anony_pair);
#else
    _finals = _aut._finals;
#endif
    //cout << "minimized finals " << _finals << endl;
    //cout << "minimized labels " << _aut._label_cube << endl;
    clock_t c_end = clock();
        cout << "Finished renaming BDD in "
        	 << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << "ms ..."
			 << "node for trans" << count << endl;

}

buddy_bdd
dfwa_min_bdd::move_to(buddy_bdd dd)
{
	unordered_map<int, buddy_bdd> computed_table;
	return move_to(dd, computed_table);
}

buddy_bdd
dfwa_min_bdd::move_to(buddy_bdd dd, unordered_map<int, buddy_bdd>& computed_table)
{
	int dd_id = dd.id();
	unordered_map<int, buddy_bdd>::const_iterator it = computed_table.find(dd_id);
	if(it != computed_table.end())
	{
		return it->second;
	}
	if(dd == bddtrue || dd == bddfalse)
	{
		//cout << "true: node id = " << dd_id << endl;
		return dd;
	}else
	{
		//nodes_num ++;
		//cout << "visited " << nodes_num << endl;
		buddy_bdd high = move_to(bdd_high(dd), computed_table);
		buddy_bdd low  = move_to(bdd_low(dd), computed_table);
		int buddy_index = bdd_var(dd);
		int anony_var_index = _map[buddy_index];
		buddy_bdd var_dd = bdd_ithvar(anony_var_index);
		buddy_bdd result = bdd_ite(var_dd, high, low);
		computed_table[dd_id] = result;
		return result;
	}
}

dfwa_min_bdd::~dfwa_min_bdd()
{
	if (_label_to_anony_pair != nullptr)
	{
		bdd_freepair(_label_to_anony_pair);
		_label_to_anony_pair = nullptr;
	}
	if (_anony_to_label_pair != nullptr)
	{
		bdd_freepair(_anony_to_label_pair);
		_anony_to_label_pair = nullptr;
	}
	if (_curr_to_next_block_pair != nullptr)
	{
		bdd_freepair(_curr_to_next_block_pair);
		_curr_to_next_block_pair = nullptr;
	}
	if (_next_to_curr_block_pair != nullptr)
	{
		bdd_freepair(_next_to_curr_block_pair);
		_next_to_curr_block_pair = nullptr;
	}
#ifdef __ORIGINAL__
	if (_state_to_anony_pair != nullptr)
	{
		bdd_freepair(_state_to_anony_pair);
		_state_to_anony_pair = nullptr;
	}

	if(_curr_to_next_state_pair != nullptr)
	{
		bdd_freepair(_curr_to_next_state_pair);
		_curr_to_next_state_pair = nullptr;
	}
	if(_next_to_curr_state_pair != nullptr)
	{
		bdd_freepair(_next_to_curr_state_pair);
		_next_to_curr_state_pair = nullptr;
	}
#endif
}

void
dfwa_min_bdd::output(ostream& os)
{
    os << "dfa: " << endl;
    os << "init: " << endl;
    os << _init << endl;
    
    os << "trans: " << endl;
    os << _trans << endl;
    
    os << "finals: " << endl;
    os << _finals << endl;
}

// -------------- Symbolic DFA Minimization Algorithms ----------------

/*-------------------------------------------------------------------*/
// create num block variables
/*-------------------------------------------------------------------*/

buddy_bdd
dfwa_min_bdd::new_block_number(int block_number, int copy)
{
    // now check whether block_number can be represented with current
    // number of block variables
    //cout << "block_number = " << block_number << endl;
    buddy_bdd dd = bddtrue;
    int bit = 1;
    for (buddy_bdd& bit_var : _block_vars[copy])
    {
        buddy_bdd bit_var_not = !bit_var;
        dd = dd & ((block_number & bit) != 0 ? bit_var : bit_var_not);
        bit <<= 1;
    }
    //cout << "returned block dd = " << dd << endl;
    
    return dd;
}

/*-------------------------------------------------------------------*/
// traverse BDD structure to compute refined partition
// 1. increase the number of block variables on-the-fly: will need exponential
//    more times of BDD traversing
// 2. set to the number of state variables and then remove unessential variables:
//    much less times of BDD traversing
/*-------------------------------------------------------------------*/ 
buddy_bdd
dfwa_min_bdd::refine_partition_rec(buddy_bdd_ptr dd, unsigned & block_number
, unordered_map<int, buddy_bdd> &computed_table)
{
    int dd_id = dd.id();
    cout << "Compute new partition from sigf(s, a, k) in refine_partition_rec  " << endl;
    // check whether the node has been visited before
    unordered_map<int, buddy_bdd>::const_iterator it = computed_table.find(dd_id);
    if(it != computed_table.end())
    {
        return it->second;
    }
    // false means cannot reach a block
    if(dd == bddfalse)
    {
    	return dd;
    }

    buddy_bdd result;
    const int top_index = bdd_var(dd);
    // check whether it is a state variable
    if(! is_state_variable(top_index) || dd == bddtrue)
    {
    	// This is a node which representing a block
        cout << "block  " << block_number << endl;
        result = new_block_number(block_number, 0);
        ++ block_number;

        // debug info
        vector<buddy_bdd> label_block(_label_vars);
        label_block.insert(label_block.end(), _block_vars[0].begin(), _block_vars[0].end());
        //generate_all_bits(label_block, 0, dd, _manager->bddOne());
    }else
    {
        // Traverse the BDD further.
        cout << "State variables identified  " << endl;
        buddy_bdd dd_var = bdd_ithvar(top_index);
        buddy_bdd low_cofactor = bdd_low(dd);
        buddy_bdd low  = refine_partition_rec(low_cofactor, block_number, computed_table);
        buddy_bdd high_cofactor = bdd_high(dd);
        buddy_bdd high = refine_partition_rec(high_cofactor, block_number, computed_table);
        result = bdd_ite(dd_var, high, low);
    }
    computed_table[dd_id] = result;
    return result;
}

pair<buddy_bdd, unsigned>
dfwa_min_bdd::refine_partition(buddy_bdd_ptr sig)
{
	unsigned int block_number = 0;
    unordered_map<int, buddy_bdd> computed_table;
    buddy_bdd partition = refine_partition_rec(sig, block_number, computed_table);
    return make_pair<>(partition, block_number);
}

/*-------------------------------------------------------------------*/
// permute variables for from and to
// can be improved by outgoing transitions
/*-------------------------------------------------------------------*/ 

/*-------------------------------------------------------------------*/
// compute next step image of curr
// FIXED (note the returned image contains propositions)
/*-------------------------------------------------------------------*/
buddy_bdd
dfwa_min_bdd::next_image(buddy_bdd_ptr curr)
{
	buddy_bdd next = bdd_relprod(curr, _trans, _curr_cube & _label_cube);
	// next to current
#ifdef __ORIGINAL__
    next = bdd_replace(next, _next_to_curr_state_pair);
#else
    next = bdd_replace(next, _aut._next_to_curr_pairs);
#endif
    return next;
}

/*-------------------------------------------------------------------*/
// compute previous image of curr
/*-------------------------------------------------------------------*/
buddy_bdd
dfwa_min_bdd::pre_image(buddy_bdd_ptr curr)
{
#ifdef __ORIGINAL__
	buddy_bdd next = bdd_replace(curr, _curr_to_next_state_pair);
#else
	buddy_bdd next = bdd_replace(curr, _aut._curr_to_next_pairs);
#endif
    // and exist
	buddy_bdd pre = bdd_relprod(next, _trans, _next_cube  & _label_cube);
    return pre;
}
/*-------------------------------------------------------------------*/
// compute reachable state space
/*-------------------------------------------------------------------*/
buddy_bdd
dfwa_min_bdd::forward_explore()
{
	clock_t c_start = clock();
    buddy_bdd s = _init;
    buddy_bdd sp = bddfalse;
    unsigned count = 1;
    while(sp != s)
    {
    	cout << "Iteration number = " << count << endl;
        // record the states reached within last step
        sp = s;
        // compute image of next step
#ifdef DEBUG

        cout << "reachable states in product: " << endl;
        bdd_print_set(cout, _state_vars.get_dict(), sp);
        cout << endl;
#endif
        s = sp | next_image(sp);
        ++ count;
    }
    clock_t c_end = clock();
    cout << "Finished forward exploration in " << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms\n";
    return s;
}

buddy_bdd
dfwa_min_bdd::backward_explore()
{
	clock_t c_start = clock();
	buddy_bdd s = _finals;
	buddy_bdd sp = bddfalse;
    unsigned count = 1;
    while(sp != s)
    {
    	cout << "Iteration number = " << count << endl;
        // record the states reached within last step
        sp = s;
        // compute image of next step
#ifdef DEBUG
        cout << "reachable states in product: " << endl;
        bdd_print_set(cout, _state_vars.get_dict(), sp);
        cout << endl;
#endif
        s = sp | pre_image(sp);
        ++ count;
    }
    clock_t c_end = clock();
    cout << "Finished backward exploration in " << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms\n";
    return s;
}

buddy_bdd
dfwa_min_bdd::initialize_partition()
{

	clock_t c_start = clock();
    buddy_bdd reach = forward_explore();
    reach = reach & backward_explore();
    clock_t c_end = clock();
    cout << "Finished buddy exploration in " << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms\n";

    // buddy is a bit fast than cudd. Not sure whether it is due to the encoding or
    // that buddy is just faster than cudd
#ifdef DEBUG
	c_start = clock();
	buddy_bdd b_reach = _aut.explore();
	b_reach = b_reach & _aut.back_explore();
	c_end = clock();
	cout << "Finished buddy exploration in " << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms\n";
#endif

    cout << "reachable non-final states: " << endl;
    buddy_bdd non_finals = (! _finals) & reach;
    //vector<buddy_bdd> curr_states(_state_vars[0]);
    //generate_all_bits(curr_states, 0, non_finals, _manager->bddOne());
    // get block variables
    buddy_bdd block_0 = new_block_number(0, 0);
    // careful about nonreachable states
    buddy_bdd block_1 = new_block_number(1, 0);
    buddy_bdd partition = (block_0 & _finals & reach) | (block_1 & non_finals);
    // debug info
    cout << "Initial partition: " << endl;
    //curr_states.insert(curr_states.end(), _block_vars[0].begin(), _block_vars[0].end());
    //generate_all_bits(curr_states, 0, partition, _manager->bddOne());
    // check whether we need to add sink transitions
    // s <-> s'
    //buddy_bdd nonreach = ! cudd_reach;
#ifdef __ORIGINAL__
    buddy_bdd next_reach = bdd_replace(reach, _curr_to_next_state_pair);
#else
    buddy_bdd next_reach = bdd_replace(reach, _aut._curr_to_next_pairs);
#endif
    _trans = _trans & reach & next_reach;
    return partition;
}


/*-------------------------------------------------------------------*/
// symbolic DFA minimization
/*-------------------------------------------------------------------*/ 
void 
dfwa_min_bdd::minimize()
{
	clock_t c_start = clock();
	cout << "Starting DFA minimization..." << endl;
    // sigref based minimization for algorithms
    // 1. P(X, K) is the initial partition of states
    buddy_bdd partition = initialize_partition();
    unsigned prev_block_number = 2;
    unsigned iteration_num = 0;
    
    while(true)
    {
    	clock_t iter_start = clock();
        ++ iteration_num;
        cout << "The number of blocks at iteration " << iteration_num << " is " << prev_block_number << endl;
        // compute the signature sigf(s, a, k) = exists X'. T(X, A, X') and P(X', K)
        // first P(X, K) => P(X', K)
        cout << "Permute P(X, K) => P(X', K)  " << endl;
#ifdef __ORIGINAL__
        buddy_bdd state_next_partition = bdd_replace(partition, _curr_to_next_state_pair);
#else
        buddy_bdd state_next_partition = bdd_replace(partition, _aut._curr_to_next_pairs);
#endif
        // sigf(s, a, k) = exists X'. T(X, A, X') and P(X', K)
        cout << "Compute sigf(s, a, k) = exists X'. T(X, A, X') and P(X', K)  " << endl;
        buddy_bdd sigf = bdd_relprod(state_next_partition, _trans, _next_cube);

        cout << "Permute P(X, K) => P(X, K')  " << endl;
        buddy_bdd block_next_partition = bdd_replace(partition, _curr_to_next_block_pair);
        cout << "sigf(s, a, k) = sigf(s, a, k) & P(s, k')" << endl;
        cout << "#sigf = " << bdd_nodecount(sigf) << " #p = " << bdd_nodecount(block_next_partition) << endl;
        sigf = sigf & block_next_partition;
        //sigf(s, a, k) [X, A, K]
        //vector<buddy_bdd> new_state_label_block_1(_state_vars[0]);
        //new_state_label_block_1.insert(new_state_label_block_1.end(), _label_vars.begin(), _label_vars.end());
        //new_state_label_block_1.insert(new_state_label_block_1.end(), _block_vars[0].begin(), _block_vars[0].end());
        //generate_all_bits(new_state_label_block_1, 0, sigf, _manager->bddOne());
        // refine the partition, traverse the BDD to compute block numbers
        cout << "Compute new partition from sigf(s, a, k)   " << endl;
        //print(sigf);
        pair<buddy_bdd, unsigned> result = refine_partition(sigf);
        unsigned curr_block_number = result.second;
        // if the number of blocks is not changed, then get out of the loop 
        clock_t iter_end = clock();
        cout << "Finished iteration "<< iteration_num << " in " << 1000.0 * (iter_end - iter_start) / CLOCKS_PER_SEC << " ms\n";
        if(prev_block_number == curr_block_number)
        {
            break;
        }
        prev_block_number = curr_block_number;
        // set refined partition
        partition = result.first;
        // debug info
        //cout << "refined partition: " << endl << partition << endl;
        //vector<buddy_bdd> new_state_label_block(_state_vars[0]);
        //new_state_label_block.insert(new_state_label_block.end(), _block_vars[0].begin(), _block_vars[0].end());
        //generate_all_bits(new_state_label_block, 0, partition, _manager->bddOne());
    }

    _num_min_states = prev_block_number;
    cout << "Finished DFA minimization after "<< iteration_num 
         << " iterations. The number of states in minimal DFA is " << _num_min_states << endl;
    // construct new dfa
    // determine whether we can remove some variables
    reduce_block_variables(partition);
    
    cout << "minimized dfa: " << endl;
    for(unsigned i = 0; i < _block_vars.size(); i ++)
    {
        for(unsigned j = 0; j < _block_vars[i].size(); j ++)
        {
            cout << "var_" << i << " " << _block_vars[i][j] << endl;
        }
    }
    
    //cout << "label cube = " << _label_cube << endl;
    cout << "init: " << endl;
    //cout << "before: " << _init << endl;
    cout << "_curr_cube =" << _curr_cube << endl;
    // I(K) = exists X. I(X) and P(X, K)
    _init = bdd_relprod(_init, partition, _curr_cube);
   // cout << "after: " << _init << endl;
    
    cout << "finals: " << endl;
    // F(K) = exists X. F(X) and P(X, K)
    //cout << "before: " << _finals << endl;
    //cout << "and: " << (_finals & partition) << endl;
    _finals = bdd_relprod(_finals, partition, _curr_cube);
    //cout << "after: " << _finals << endl;
    
    cout << "trans: " << endl;
    // P(X', K')
    //vector<buddy_bdd> curr_states(_state_vars[0]);
    //curr_states.insert(curr_states.end(), _block_vars[0].begin(), _block_vars[0].end());
    //vector<buddy_bdd> next_states(_state_vars[1]);
    //next_states.insert(next_states.end(), _block_vars[1].begin(), _block_vars[1].end());
    
    cout << "P(X, K) => P(X', K')" << endl;
#ifdef __ORIGINAL__
    bddPair* curr_to_next = bdd_mergepairs(_curr_to_next_state_pair, _curr_to_next_block_pair);
    buddy_bdd next_partition = bdd_replace(partition, curr_to_next);
#else
    bddPair* curr_to_next = bdd_mergepairs(_aut._curr_to_next_pairs, _curr_to_next_block_pair);
    buddy_bdd next_partition = bdd_replace(partition, curr_to_next);
#endif
    bdd_freepair(curr_to_next);
    cout << "Computing the transition relation for minimal DFA..." << endl;
    clock_t t_start = clock();

    // R(s, a, k') := (∃t : T (s, a, t) ∧ P(t, k'))
    cout << "R(s, a, k') := ∃t : T (s, a, t) ∧ P(t, k')" << endl;
    cout << "#nodes in partition: " << bdd_nodecount(partition) << endl;
    cout << "#nodes in trans: " << bdd_nodecount(_trans) << endl;

    buddy_bdd minimized_trans = bdd_relprod(_trans, next_partition, _next_cube);
    // T P (s, a, t') := (∃s : R(s, a, k') ∧ P(s, k))
    cout << "T(k, a, k') := ∃s :R(s, a, k') ∧ P(s, k)" << endl;
    minimized_trans = bdd_relprod(minimized_trans, partition, _curr_cube);
    //cout << minimized_trans << endl;
    /*
    buddy_bdd minimized_trans = make_quotient_trans(_trans, partition, next_partition);
     */
    clock_t t_end = clock();
    cout << "Finished transition relation computation in " << 1000.0 * (t_end - t_start) / CLOCKS_PER_SEC << " ms\n";
    cout << "Output Transitions" << endl;
    //vector<buddy_bdd> label_block(_label_vars);
    //label_block.insert(label_block.end(), _block_vars[0].begin(), _block_vars[0].end());
    //label_block.insert(label_block.end(), _block_vars[1].begin(), _block_vars[1].end());
    //generate_all_bits(label_block, 0, minimized_trans, _manager->bddOne());
    cout << "The number of nodes in minimal transitions" << bdd_nodecount(minimized_trans) << endl;
    _trans = minimized_trans;
    cout << endl;
    clock_t c_end = clock();

    cout << "Finished DFA minimization in " << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms\n";
}

/*-------------------------------------------------------------------*/
// reduce the number of block variables as much as possible
/*-------------------------------------------------------------------*/
void
dfwa_min_bdd::reduce_block_variables(buddy_bdd_ptr partition)
{
    // determine whether we can remove some variables
	cout << "Reduced block variables..." << endl;
	unsigned num = 0;
    while(true)
    {
        buddy_bdd var_dd = _block_vars[0].back();
        buddy_bdd result = (! var_dd) & partition;
        if(result != partition || _block_vars[0].size() < 2)
        {
        	break;
        }
        ++ num;
        _block_vars[0].pop_back();
        cout << "number of vars: " << _aut._state_vars.get_var_num(0) << endl
             << "block[0].size() = " << _block_vars[0].size() << endl
			 << "block[1].size() = " << _block_vars[1].size() << endl;
        _block_vars[1].pop_back();
        // remove redundant variables
        partition = bdd_exist(partition, var_dd);
    }
    cout << "Finished reducing " << num << " block variables..." << endl;
}

/*-------------------------------------------------------------------*/
// move the minimal DFA from CUDD to BuDDy
// the ordering of dfa is (A, X, X')
/*-------------------------------------------------------------------*/
dfwa_ptr
dfwa_min_bdd::move_dfwa()
{
	cout << "labels in previous dfa: " << _aut._label_cube << endl;
	bdd_print_set(cout, _aut.get_dict(), _aut._label_cube);
	cout << endl;
	dfwa* result = new dfwa(_aut.get_dict(), _aut._label_cube);
	cout << "labels in minimal dfa: " << endl;
	bdd_print_set(cout, result->get_dict(), result->_label_cube);
	cout << endl;
	//mapping block variables back to buddy variables
	unordered_map<int, int> cudd_to_buddy_map;
	// K -> S -> BUDDY
	//vector<int> block_vars;
    result->_state_vars._copies = 2;
    // now we add variables from op1 and op2
    cout << "Moving minimal DFA to BuDDy..." << endl;
    // check whether this part can be improved

	// compute mapping for state variables
	// needs to reuse the variables in _aut and the ordering of variables
    // should be increasing to get a better performance
    result->_state_vars.add_ordered_bdd_vars(_aut._state_vars);
	assert(result->_state_vars._dd_vars[0].size() >= _block_vars[0].size());
	cout << "Now remove useless variables" << endl;
	int pop_num = result->_state_vars.get_bdd_vars(0).size() - _block_vars[0].size();
	while(pop_num > 0)
	{
		result->_state_vars.pop_back_vars();
		-- pop_num;
	}
	bddPair* block_to_states_pair = bdd_newpair();
	//_map.clear();
	for(unsigned i = 0; i < _block_vars[0].size(); i ++)
	{
		// current version of state variables
		int index = bdd_var(_block_vars[0][i]);
		int result_index = bdd_var(result->_state_vars.get_var(0, i));
		bdd_setpair(block_to_states_pair, index, result_index);
		//_map[index] = result_index;
		//cout << "map: " << index << "  -> " << bdd_var(result->_state_vars.get_var(0, i)) << endl;
		//curr_state_vars.push_back(bdd_ithvar(buddy_index));
		//++ num;
		// next version of state variables
		index = bdd_var(_block_vars[1][i]);
		//buddy_index = _cudd_to_buddy_vars[num];
		result_index = bdd_var(result->_state_vars.get_var(1, i));
		//_map[index] = result_index;
		bdd_setpair(block_to_states_pair, index, result_index);
	}
	//vector<vector<buddy_bdd>> state_vars;
	//state_vars.push_back(curr_state_vars);
	//state_vars.push_back(next_state_vars);
	//result._state_vars.add_bdd_vars(curr_state_vars);
	//result._state_vars.add_bdd_vars(next_state_vars);
	// compute mapping for labels
	bddPair* result_pair = bdd_mergepairs(block_to_states_pair, _anony_to_label_pair);
	cout << "Migrating dfa representation from CUDD to BuDDy..." << endl;
	clock_t c_start = clock();
	//cout << "cudd init: " << _init << endl;
	// compute the initial states
	result->_init = bdd_replace(_init, block_to_states_pair);
	//cout << "buddy init: " << result->_init << endl;
	//cout << "cudd _finals: " << _finals << endl;
	result->_finals = bdd_replace(_finals, block_to_states_pair);
	//cout << "buddy finals: " << result->_finals << endl;
    cout << "The number of nodes in minimal transition is " << bdd_nodecount(_trans)  << endl;
    //_map.insert(_anony_to_label.begin(), _anony_to_label.end());

    result->_trans = bdd_replace(_trans, result_pair);
    //result->_trans = move_to(_trans);
	//cout << "buddy trans: " << result->_trans << endl;
	clock_t c_end = clock();
	cout << "Finished migrating dfa representation in "
		 << 1000.0 * (c_end - c_start)/CLOCKS_PER_SEC << "ms ..." << endl;

	result->_curr_cube = result->_state_vars.get_cube(0);
	result->_next_cube = result->_state_vars.get_cube(1);
    // make pairs

	result->_curr_to_next_pairs = result->_state_vars.make_pair(0, 1);
	result->_next_to_curr_pairs = result->_state_vars.make_pair(1, 0);
	result->_reach = bddtrue;//result.explore();
    cout << "Finished computing dfa representation in BuDDy..." << endl;

    bdd_freepair(block_to_states_pair);
    bdd_freepair(result_pair);
    /*
    cout << "is_even(8) = " << is_even(8) << endl;
    cout << "is_even(9) = " << is_even(9) << endl;
    cout << "is_even(0) = " << is_even(0) << endl;
    cout << "is_even(1) = " << is_even(1) << endl;
    */
    return *result;
}

// ------------------ printing functions ------------------------------

void
dfwa_min_bdd::print(buddy_bdd bdd)
{
    cout << "formula : " << bdd << endl;
}

/*-------------------------------------------------------------------*/
// generate all possible truth values of a vector of variables
/*-------------------------------------------------------------------*/
void
dfwa_min_bdd::generate_all_bits(vector<buddy_bdd>& vars, unsigned index, buddy_bdd_ptr dd, buddy_bdd temp)
{ 
#ifdef DEBUG
    if (index == vars.size()) 
    {
        buddy_bdd inter = temp & dd;
        if(! inter.IsZero())
        {
            cout << "SAT: ";
            for(unsigned i = 0; i < vars.size(); i ++)
            {
                inter = temp & vars[i]; 
                if(inter.IsZero())
                {
                    cout << " !";
                }else
                {
                    cout << " ";
                }
                cout << vars[i] ;
            }
            cout <<  endl;
        }
        return; 
    }
    // true
    generate_all_bits(vars, index + 1, dd, temp & vars[index]); 
    // false
    generate_all_bits(vars, index + 1, dd, temp & (!vars[index])); 
#endif
}
