#include "dfwamin3.hh"

/*-------------------------------------------------------------------*/
// prepare state, label and block variables for DFA minimization
/*-------------------------------------------------------------------*/
void
dfwa_min_sylvan::prepare_variables()
{
    // 0. create duplicate label variables for dfa minimization
    // A' duplicate variables for labels
	// temporary variables
    vector<int> label_vars;
    get_list_var_indices(label_vars, _aut._label_cube);
    for(uint32_t i = 0; i < label_vars.size(); i ++)
    {
    	int buddy_index = label_vars[i];
        _sylvan_to_buddy_vars.push_back(buddy_index);
        // from buddy to cudd for migration
    	//_buddy_to_cudd_vars[buddy_index] = i;
    	        // new label variables
    	sylvan_bdd var_dd = sylvan_bdd::bddVar(i);
    }

    // 1. create state variables for dfa minimization
    // (X, X')
    //cout << "buddy state cube: " << endl;
    // current and next state variable vectors
    vector<uint32_t> curr_states;
    _state_vars.push_back(curr_states);
    vector<uint32_t> next_states;
    _state_vars.push_back(next_states);

    // cube bdd
    //_curr_cube = sylvan_bdd::bddOne();
    //_next_cube = sylvan_bdd::bddOne();
    
    _num_states_vars = _aut._state_vars.get_var_num(0);
    // for permuting variables (starting from label_vars.size())
    // ordered state variables will speed up migration
    uint32_t cudd_index = label_vars.size();
    for(unsigned i = 0; i < _num_states_vars; i ++)
    {
    	// current state
    	int buddy_index = _aut._state_vars.get_var_id(0, i);
        _sylvan_to_buddy_vars.push_back(buddy_index);
        _buddy_to_sylvan_vars[buddy_index] = cudd_index;
        //cout << "buddy: " << buddy_index << " <-> " << _buddy_to_cudd_vars[buddy_index] << endl;
        // new state variables
        sylvan_bdd var_dd = sylvan_bdd::bddVar(cudd_index);
        _state_vars[0].push_back(cudd_index);
        _curr_cube.add(cudd_index);

        ++ cudd_index;

        // next state
        buddy_index = _aut._state_vars.get_var_id(1, i);
        _sylvan_to_buddy_vars.push_back(buddy_index);
        _buddy_to_sylvan_vars[buddy_index] = cudd_index;
        //cout << "buddy: " << buddy_index << " <-> " << _buddy_to_cudd_vars[buddy_index] << endl;
                // new state variables
        var_dd = sylvan_bdd::bddVar(cudd_index);
        _state_vars[1].push_back(cudd_index);
        _next_cube.add(cudd_index);

        ++ cudd_index;
    }
    // end of state variables (number of state variables)
    _num_states_vars = 2 * _num_states_vars + label_vars.size();
    // 2. create label variables
    // be careful about not changing _aut labels
    // A: the real label variables in minimization

    //_label_cube = sylvan_bdd::bddOne();
    for(unsigned i = 0; i < label_vars.size(); i ++)
    {
    	int buddy_index = label_vars[i];
    	cudd_index = i + _num_states_vars;
        _sylvan_to_buddy_vars.push_back(buddy_index);
        // we do not use this copy of label variables for migration
        // from buddy to cudd
        _buddy_to_sylvan_vars[buddy_index] = cudd_index;
        // new label variables
        sylvan_bdd var_dd = sylvan_bdd::bddVar(cudd_index);
        _label_cube.add(cudd_index);
        _label_vars.push_back(cudd_index);
    }
    // end of label variables
    _last_pos_labels = label_vars.size() + _num_states_vars;
    
    _curr_label_cube.add(_label_cube);
    _curr_label_cube.add(_curr_cube);

    _next_label_cube.add(_label_cube);
    _next_label_cube.add(_next_cube);

    // 3. create block vectors
    // K block variables
    vector<uint32_t> curr_block;
    _block_vars.push_back(curr_block);
    vector<uint32_t> next_block;
    _block_vars.push_back(next_block);
    
    unsigned num_block_vars = _num_states_vars - label_vars.size();
    cudd_index = _last_pos_labels;
    // prepare block variables
    for(unsigned i =0; i < num_block_vars; i ++)
    {
        sylvan_bdd var_dd = sylvan_bdd::bddVar(cudd_index);
        //cout << "block variable " << var_dd << endl;
        if(is_even(i))
        {
            _block_vars[0].push_back(cudd_index);
        }else
        {
            _block_vars[1].push_back(cudd_index);
        }
        ++ cudd_index;
    }

    _num_vars = 4 * _num_states_vars + 2 * label_vars.size();
}

/*-------------------------------------------------------------------*/
// prepare Sylvan DFA minimization
/*-------------------------------------------------------------------*/
void
dfwa_min_sylvan::prepare()
{
	int workers = 6;
	lace_init(workers, 1024*1024*16);
	lace_startup(0, NULL, NULL);
	LACE_ME;

	size_t max = 16LL<<30;
	//if (max > getMaxMemory()) max = getMaxMemory()/10*9;
	//sylvan_set_sizes(1LL<<22, 1LL<<26, 1LL<<22, 1LL<<26);
	sylvan_set_sizes(1LL<<26, 1LL<<28, 1LL<<25, 1LL<<26);
	sylvan_init_package();

	// Initialize the BDD module with granularity 1 (cache every operation)
	// A higher granularity (e.g. 6) often results in better performance in practice
	sylvan_init_bdd();

    // create variables of dfwa_min_sylvan::
	//cout << "minimized before prepare " << _aut._label_cube << endl;
    prepare_variables();
    
    // disable dynamic variable ordering
#ifdef DEBUG
    cout << "buddy aut _init: " << _aut._init << endl;
    bdd_print_set(cout, _aut.get_dict(), _aut._init);
    cout << endl;
#endif
    cout << "Migrating BuDDy to Sylvan ..." << endl;
    // copy the initial state
    clock_t c_start = clock();
    _init = move_to_sylvan(_aut._init);

    //cout << "minimized init: " << _init.NodeCount() << endl;
    // copy the transition relation
#ifdef DEBUG
    cout << "buddy aut _trans: " << endl;
    bdd_print_set(cout, _aut.get_dict(), _aut._trans);
    cout << endl;

    cout << "The number of nodes in transition is " << bdd_nodecount(_aut._trans) << endl;
    buddy_bdd reach = _aut.explore();
    buddy_bdd all = bdd_replace(reach, _aut._curr_to_next_pairs);
    all = all & reach;
    _aut._trans = _aut._trans & all;
#endif
    //int count = bdd_nodecount(_aut._trans);
    //cout << "The number of nodes in transition after reach is " << bdd_nodecount(_aut._trans) << endl;
    _trans = move_to_sylvan(_aut._trans);
    // copy the final states
    //cout << "minimized trans: " << _trans << endl;
    _finals = move_to_sylvan(_aut._finals);
    //cout << "minimized finals " << _finals << endl;
    //cout << "minimized labels " << _aut._label_cube << endl;
    clock_t c_end = clock();
    cout << "Finished migrating BuDDy to Sylvan in "
        	 << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << "ms ..."
			 << "node for trans " << _trans.NodeCount() << endl;
    // change label variables to block variables
    vector<uint32_t> prev_label_vec;
    for(uint32_t i = 0; i < _label_vars.size(); i ++)
    {
    	prev_label_vec.push_back(i);
    }
    vector<uint32_t> curr_label_vec = _label_vars;
    assert (prev_label_vec.size() == curr_label_vec.size());

   // _trans = _trans.Permute(prev_label_vec, curr_label_vec);
    c_end = clock();
    cout << "Finished renaming Sylvan in "
        	 << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << "ms ..." << endl;
    //cout << "Now restore the ordering for minimization" << endl;
    //_manager->ShuffleHeap(order);
    c_end = clock();
    cout << "Finished migrating BuDDy to Sylvan in "
    	 << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << "ms ..." << endl;


}

dfwa_min_sylvan::~dfwa_min_sylvan()
{
	sylvan_quit();
}

void
dfwa_min_sylvan::output(ostream& os)
{
    os << "dfa: " << endl;
    os << "init: " << endl;
    os << _init.NodeCount() << endl;
    
    os << "trans: " << endl;
    os << _trans.NodeCount() << endl;
    
    os << "finals: " << endl;
    os << _finals.NodeCount() << endl;
}

// -------------- Symbolic DFA Minimization Algorithms ----------------

/*-------------------------------------------------------------------*/
// create num block variables
/*-------------------------------------------------------------------*/

sylvan_bdd
dfwa_min_sylvan::new_block_number(int block_number, int copy)
{
    // now check whether block_number can be represented with current
    // number of block variables
    //cout << "block_number = " << block_number << endl;
    sylvan_bdd dd = sylvan_bdd::bddOne();
    int bit = 1;
    for (uint32_t& bit_var : _block_vars[copy])
    {
    	sylvan_bdd bit_var_dd =  sylvan_bdd::bddVar(bit_var);
        dd = dd & ((block_number & bit) != 0 ? bit_var_dd : (! bit_var_dd));
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
sylvan_bdd
dfwa_min_sylvan::refine_partition_rec(sylvan_bdd_ptr dd, unsigned & block_number
, unordered_map<sylvan_BDD, sylvan_bdd> &computed_table)
{
	sylvan_BDD node = dd.GetBDD();
    //cout << "Compute new partition from sigf(s, a, k) in refine_partition_rec  " << endl;
    // check whether the node has been visited before
    unordered_map<sylvan_BDD, sylvan_bdd>::const_iterator it = computed_table.find(node);
    if(it != computed_table.end())
    {
        return it->second;
    }
    // false means cannot reach a block
    if(dd.isZero())
    {
    	return dd;
    }

    sylvan_bdd result;
    const uint32_t top_index = dd.TopVar();
    //cout << "topvar = " << top_index << _num_states_vars << endl;
    // check whether it is a state variable
    if(! is_state_variable(top_index) || dd.isOne())
    {
    	// This is a node which representing a block
        //cout << "block  " << block_number << endl;
        result = new_block_number(block_number, 0);
        ++ block_number;

        // debug info
        vector<uint32_t> label_block(_label_vars);
        label_block.insert(label_block.end(), _block_vars[0].begin(), _block_vars[0].end());
        //generate_all_bits(label_block, 0, dd, _manager->bddOne());
    }else
    {
        // Traverse the BDD further.
        //cout << "State variables identified  " << endl;
        sylvan_bdd dd_var = sylvan_bdd::bddVar(top_index);
        sylvan_bdd low_cofactor = dd.Else();
        sylvan_bdd low  = refine_partition_rec(low_cofactor, block_number, computed_table);
        sylvan_bdd high_cofactor = dd.Then();
        sylvan_bdd high = refine_partition_rec(high_cofactor, block_number, computed_table);
        result = dd_var.Ite(high, low);
    }
    computed_table[node] = result;
    return result;
}

pair<sylvan_bdd, unsigned>
dfwa_min_sylvan::refine_partition(sylvan_bdd_ptr sig)
{
	unsigned int block_number = 0;
    unordered_map<sylvan_BDD, sylvan_bdd> computed_table;
    sylvan_bdd partition = refine_partition_rec(sig, block_number, computed_table);
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
sylvan_bdd
dfwa_min_sylvan::next_image(sylvan_bdd_ptr curr)
{
	sylvan_bdd next = curr.AndAbstract(_trans, _curr_label_cube);
	// next to current
    next = next.Permute(_state_vars[1], _state_vars[0]);
    return next;
}

/*-------------------------------------------------------------------*/
// compute previous image of curr
/*-------------------------------------------------------------------*/
sylvan_bdd
dfwa_min_sylvan::pre_image(sylvan_bdd_ptr curr)
{
	sylvan_bdd next = curr.Permute( _state_vars[0], _state_vars[1]);
    // and exist
	sylvan_bdd pre = next.AndAbstract(_trans, _next_label_cube);
    return pre;
}
/*-------------------------------------------------------------------*/
// compute reachable state space
/*-------------------------------------------------------------------*/
sylvan_bdd
dfwa_min_sylvan::forward_explore()
{
	clock_t c_start = clock();
    sylvan_bdd s = _init;
    sylvan_bdd sp = sylvan_bdd::bddZero();
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
    cout << "Finished forward exploration in " << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " " << s.NodeCount() << " ms\n";
    return s;
}

sylvan_bdd
dfwa_min_sylvan::backward_explore()
{
	clock_t c_start = clock();
	sylvan_bdd s = _finals;
	sylvan_bdd sp = sylvan_bdd::bddZero();
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
    cout << "Finished backward exploration in " << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " " << s.NodeCount() << " ms\n";
    return s;
}

sylvan_bdd
dfwa_min_sylvan::initialize_partition()
{

	clock_t c_start = clock();
    sylvan_bdd reach = forward_explore();
    reach = reach & backward_explore();
    clock_t c_end = clock();
    cout << "Finished sylvan exploration in " << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC  << " ms " << reach.NodeCount()<< "\n";

    // buddy is a bit fast than cudd. Not sure whether it is due to the encoding or
    // that buddy is just faster than cudd
#ifdef DEBUG
	c_start = clock();
	buddy_bdd b_reach = _aut.explore();
	b_reach = b_reach & _aut.back_explore();
	//sylvan_bdd copy = move_to_sylvan(b_reach);
	c_end = clock();
	cout << "Finished buddy exploration in " << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms " << bdd_nodecount(b_reach)<< "\n";
#endif

    //cout << "reachable non-final states: " << endl;
    sylvan_bdd non_finals = (! _finals) & reach;
    vector<uint32_t> curr_states(_state_vars[0]);
    //generate_all_bits(curr_states, 0, non_finals, _manager->bddOne());
    // get block variables
    sylvan_bdd block_0 = new_block_number(0, 0);
    // careful about nonreachable states
    sylvan_bdd block_1 = new_block_number(1, 0);
    sylvan_bdd partition = (block_0 & _finals & reach) | (block_1 & non_finals);
    // debug info
    //cout << "Initial partition: " << endl;
    curr_states.insert(curr_states.end(), _block_vars[0].begin(), _block_vars[0].end());
    //generate_all_bits(curr_states, 0, partition, _manager->bddOne());
    // check whether we need to add sink transitions
    // s <-> s'
    //sylvan_bdd nonreach = ! cudd_reach;
    sylvan_bdd next_reach = reach.Permute( _state_vars[0], _state_vars[1]);
    _trans = _trans & reach & next_reach;
    return partition;
}


/*-------------------------------------------------------------------*/
// symbolic DFA minimization
/*-------------------------------------------------------------------*/ 
void 
dfwa_min_sylvan::minimize()
{
	clock_t c_start = clock();
	cout << "Starting DFA minimization..." << endl;
    // sigref based minimization for algorithms
    // 1. P(X, K) is the initial partition of states
    sylvan_bdd partition = initialize_partition();
    //cout << "Nodes in partition = " << partition.NodeCount() << endl;
    unsigned prev_block_number = 2;
    unsigned iteration_num = 0;
    
    //sylvan_bdd result_sigf = sylvan_bdd::bddOne();
    while(true)
    {
    	clock_t iter_start = clock();
        ++ iteration_num;
        cout << "Number of blocks at iteration " << iteration_num << " is " << prev_block_number << endl;
        // compute the signature sigf(s, a, k) = exists X'. T(X, A, X') and P(X', K)
        // first P(X, K) => P(X', K)
        //cout << "Permute P(X, K) => P(X', K)  " << endl;
        sylvan_bdd state_next_partition = partition.Permute( _state_vars[0], _state_vars[1]);
        //cout << "Nodes in partition = " << state_next_partition.NodeCount() << endl;
        // sigf(s, a, k) = exists X'. T(X, A, X') and P(X', K)
        //cout << "Compute sigf(s, a, k) = exists X'. T(X, A, X') and P(X', K)  " << endl;
        sylvan_bdd sigf = state_next_partition.AndAbstract(_trans, _next_cube);
        //cout << "Nodes in partition = " << sigf.NodeCount() << endl;
        //cout << "Permute P(X, K) => P(X, K')  " << endl;
        sylvan_bdd block_next_partition = partition.Permute( _block_vars[0], _block_vars[1]);
        //cout << "Nodes in partition = " << block_next_partition.NodeCount() << endl;
        //cout << "sigf(s, a, k) = sigf(s, a, k) & P(s, k') = exists X'. P(s, k') and T(X, A, X') and P(X', K) " << endl;
        /*
        sylvan_bdd temp = sylvan_bdd::bddZero();
        vector<sylvan_bdd> partition_vec;
        for(unsigned i = 0; i < prev_block_number; i ++)
        {
        	sylvan_bdd partition_index = new_block_number(i, 1);
        	partition_index = sigf & partition_index & block_next_partition;
        	temp |= partition_index;
        }*/
        //sigf = temp;
        sigf = sigf & block_next_partition;
        //result_sigf = sigf;
        //cout << "Nodes in partition = " << sigf.NodeCount() << endl;
        //sigf(s, a, k) [X, A, K]
        vector<uint32_t> new_state_label_block_1(_state_vars[0]);
        new_state_label_block_1.insert(new_state_label_block_1.end(), _label_vars.begin(), _label_vars.end());
        new_state_label_block_1.insert(new_state_label_block_1.end(), _block_vars[0].begin(), _block_vars[0].end());
        //generate_all_bits(new_state_label_block_1, 0, sigf, _manager->bddOne());
        // refine the partition, traverse the BDD to compute block numbers
        cout << "Computing new partition from signature..." << endl;
        //print(sigf);
        pair<sylvan_bdd, unsigned> result = refine_partition(sigf);
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
        //cout << "Nodes in partition = " << partition.NodeCount() << endl;
        // debug info
        //cout << "refined partition: " << endl << partition << endl;
        vector<uint32_t> new_state_label_block(_state_vars[0]);
        new_state_label_block.insert(new_state_label_block.end(), _block_vars[0].begin(), _block_vars[0].end());
        //generate_all_bits(new_state_label_block, 0, partition, _manager->bddOne());
    }

    //cout << "T_min(k', a, k) = exists X, X'. P(s, k') and T(X, A, X') and P(X', K) " << endl;
    //result_sigf = result_sigf.ExistAbstract(_curr_cube);
    //cout << "#T_min(k', a, k) = " << result_sigf.NodeCount() << endl;
    _num_min_states = prev_block_number;
    cout << "Finished DFA minimization after "<< iteration_num 
         << " iterations. \nThe number of states in minimal DFA is " << _num_min_states << endl;
    // construct new dfa
    // determine whether we can remove some variables
    reduce_block_variables(partition);
    /*
    cout << "minimized dfa: " << endl;
    for(unsigned i = 0; i < _block_vars.size(); i ++)
    {
        for(unsigned j = 0; j < _block_vars[i].size(); j ++)
        {
            cout << "var_" << i << " " << _block_vars[i][j] << endl;
        }
    }
    */
    //cout << "label cube = " << _label_cube << endl;
    //cout << "init: " << endl;
    //cout << "before: " << _init << endl;
    //cout << "_curr_cube =" << endl;
    // I(K) = exists X. I(X) and P(X, K)
    _init = _init.AndAbstract(partition, _curr_cube);
   // cout << "after: " << _init << endl;
    
    //cout << "finals: " << endl;
    // F(K) = exists X. F(X) and P(X, K)
    //cout << "before: " << _finals << endl;
    //cout << "and: " << (_finals & partition) << endl;
    _finals = _finals.AndAbstract(partition, _curr_cube);
    //cout << "after: " << _finals << endl;
    
    //cout << "trans: " << endl;
    // P(X', K')
    vector<uint32_t> curr_states(_state_vars[0]);
    curr_states.insert(curr_states.end(), _block_vars[0].begin(), _block_vars[0].end());
    vector<uint32_t> next_states(_state_vars[1]);
    next_states.insert(next_states.end(), _block_vars[1].begin(), _block_vars[1].end());
    
    //cout << "P(X, K) => P(X', K')" << endl;
    sylvan_bdd next_partition = partition.Permute(curr_states, next_states);
    cout << "Computing the transition relation for minimal DFA..." << endl;
    clock_t t_start = clock();
    //cout << "#T(s, a, t) = " << _trans.NodeCount() << " #P(s, k) = " << partition.NodeCount()
    //		<< " #P(t, k') = " << next_partition.NodeCount() << endl;

    // R(s, a, k') := (∃t : T (s, a, t) ∧ P(t, k'))
    //cout << "R(s, a, k') := ∃t : T (s, a, t) ∧ P(t, k')" << endl;
    //cout << "#nodes in partition: " << partition.NodeCount() << endl;
    //cout << "#nodes in trans: " << _trans.NodeCount() << endl;

    sylvan_bdd minimized_trans = _trans.AndAbstract(partition, _curr_cube);
    //cout << "#nodes in mini_trans: " << minimized_trans.NodeCount() << endl;
    // T P (s, a, t') := (∃s : R(s, a, k') ∧ P(s, k))
    //cout << "T(k, a, k') := ∃s :R(s, a, k') ∧ P(s, k)" << endl;
    minimized_trans = minimized_trans.AndAbstract(next_partition, _next_cube);
    //cout << minimized_trans << endl;

    //cout << "R(s, a, k') := ∃s ∃t : P(s, k) ∧ T (s, a, t) ∧ P(t, k')" << endl;
    //sylvan_bdd minimized_trans = make_quotient_trans(_trans, partition, next_partition);
    /*
    cout << "R(s, a, k') := P(s, k) ∧ P(t, k')" << endl;
    sylvan_bdd minimized_trans = partition & next_partition;
    cout << "#P(t, k') ∧ P(t, k') = " << minimized_trans.NodeCount() << endl;
    sylvan_bdd_set curr_next_cube;
    curr_next_cube.add(_curr_cube);
    curr_next_cube.add(_next_cube);
    minimized_trans = minimized_trans.AndAbstract(_trans, curr_next_cube);
    */
    clock_t t_end = clock();
    cout << "Finished transition relation computation in " << 1000.0 * (t_end - t_start) / CLOCKS_PER_SEC << " ms\n";
    //cout << "Output Transitions" << endl;
    vector<uint32_t> label_block(_label_vars);
    label_block.insert(label_block.end(), _block_vars[0].begin(), _block_vars[0].end());
    label_block.insert(label_block.end(), _block_vars[1].begin(), _block_vars[1].end());
    //generate_all_bits(label_block, 0, minimized_trans, _manager->bddOne());
    _trans = minimized_trans;
    //cout << endl;
    clock_t c_end = clock();
    cout << "Finished DFA minimization in " << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms\n";
}

// compute the transition relation for minimal DFA

/*-------------------------------------------------------------------*/
// making transition relation for quotient DFA
// probably work better for multi-core bdd
/*-------------------------------------------------------------------*/

sylvan_bdd
dfwa_min_sylvan::make_quotient_trans(sylvan_bdd_ptr trans, sylvan_bdd_ptr curr_partition, sylvan_bdd_ptr next_partition)
{
	// T(X, X', A), P(X, K), P(X', K')
	map_t_sylvan computed_table;
	set<uint32_t> curr_set;
	set<uint32_t> next_set;
	vector<set<uint32_t>> var_sets;
	var_sets.push_back(curr_set);
	var_sets.push_back(next_set);
	for(unsigned i = 0; i < var_sets.size(); i ++)
	{
		for(unsigned j = 0; j < _state_vars[i].size(); j ++)
		{
			var_sets[i].insert(_state_vars[i][j]);
		}
	}

	return make_quotient_trans(trans, curr_partition, next_partition, var_sets, computed_table);
}
sylvan_bdd
dfwa_min_sylvan::make_quotient_trans(sylvan_bdd trans, sylvan_bdd curr_partition, sylvan_bdd next_partition
		, vector<set<uint32_t>>& var_sets, map_t_sylvan& computed_table)
{
	//cout << "Entering recursive construction for transition relation..." << endl;
	sylvan_BDD first = trans.GetBDD();
	sylvan_BDD second = curr_partition.GetBDD();
	sylvan_BDD third = next_partition.GetBDD();
	//cout << "f0 = " << first << " f1 = " << second << " f2 = " << third << endl;
	map_t_sylvan::const_iterator it = computed_table.find(make_tuple(first, second, third));
	if(it != computed_table.end())
	{
		return it->second;
	}
	if(trans.isZero() || curr_partition.isZero() || next_partition.isZero())
	{
		return sylvan_bdd::bddZero();
	}
	//cout << "Needs to do construction for transition relation..." << endl;
	//cout << "f0 = " << trans << " f1 = " << curr_partition << " f2 = " << next_partition << endl;
	// computed variables
	uint32_t top_var = min(trans.TopVar(), curr_partition.TopVar());
	top_var = min(top_var, next_partition.TopVar());
	sylvan_bdd top_var_dd = sylvan_bdd::bddVar(top_var);
	//cout << "top_var = " << top_var << endl;
	sylvan_bdd result;
	BddMap then_map(top_var, sylvan_bdd::bddOne());
	BddMap else_map(top_var, sylvan_bdd::bddZero());
	if(var_sets[0].find(top_var) != var_sets[0].end()
	|| var_sets[1].find(top_var) != var_sets[1].end())
	{
		//cout << "current variable: " << endl;

		sylvan_bdd low = make_quotient_trans(
				  trans.Compose(else_map)
				, curr_partition.Compose(else_map)
				, next_partition
				, var_sets, computed_table);
		sylvan_bdd high = make_quotient_trans(
						  trans.Compose(then_map)
						, curr_partition.Compose(then_map)
						, next_partition, var_sets, computed_table);
		result = high | low;
	}else
	if(var_sets[1].find(top_var) != var_sets[1].end())
	{
		//cout << "next variable: " << endl;
		sylvan_bdd low = make_quotient_trans(
				  trans.Compose(else_map)
				, curr_partition
				, next_partition.Compose(else_map)
				, var_sets, computed_table);
		sylvan_bdd high = make_quotient_trans(
						  trans.Compose(then_map)
						, curr_partition
						, next_partition.Compose(then_map)
						, var_sets, computed_table);
		result = high | low;
	}else
	{
		// needs to pick one state then  A and (K, K')
		result = trans & curr_partition & next_partition;
	}
	//cout << "return from recursive function for transitin relation " << endl;
	computed_table[make_tuple(first, second, third)] = result;
	return result;
}

/*-------------------------------------------------------------------*/
// reduce the number of block variables as much as possible
/*-------------------------------------------------------------------*/
void
dfwa_min_sylvan::reduce_block_variables(sylvan_bdd_ptr partition)
{
    // determine whether we can remove some variables
	cout << "Reducing block variables..." << endl;
	unsigned num = 0;
    while(true)
    {
        sylvan_bdd var_dd = _block_vars[0].back();
        sylvan_bdd result = (! var_dd) & partition;
        if(result != partition || _block_vars[0].size() < 2)
        {
        	break;
        }
        ++ num;
        _block_vars[0].pop_back();
        //cout << "number of vars: " << _num_states_vars << endl
        //     << "block[0].size() = " << _block_vars[0].size() << endl
		//	 << "block[1].size() = " << _block_vars[1].size() << endl;
        _block_vars[1].pop_back();
        // remove redundant variables
        partition = partition.ExistAbstract(var_dd);
    }
    cout << "Finished reducing " << num << " block variables..." << endl;
}

/*-------------------------------------------------------------------*/
// move the minimal DFA from Sylvan to BuDDy
// the ordering of dfa is (A, X, X')
/*-------------------------------------------------------------------*/
dfwa_ptr
dfwa_min_sylvan::move_dfwa()
{
	//cout << "labels in previous dfa: " << _aut._label_cube << endl;
	//bdd_print_set(cout, _aut.get_dict(), _aut._label_cube);
	//cout << endl;
	dfwa* result = new dfwa(_aut.get_dict(), _aut._label_cube);
	//cout << "labels in minimal dfa: " << endl;
	//bdd_print_set(cout, result->get_dict(), result->_label_cube);
	//cout << endl;
	//mapping block variables back to buddy variables
	unordered_map<int, int> cudd_to_buddy_map;
	// K -> S -> BUDDY
	//vector<int> block_vars;
    result->_state_vars._copies = 2;
    // now we add variables from op1 and op2
    //cout << "Moving minimal DFA to BuDDy..." << endl;
    // check whether this part can be improved

	// compute mapping for state variables
	// needs to reuse the variables in _aut and the ordering of variables
    // should be increasing to get a better performance
    result->_state_vars.add_ordered_bdd_vars(_aut._state_vars);
	assert(result->_state_vars._dd_vars[0].size() >= _block_vars[0].size());
	for(unsigned i = 0; i < _block_vars[0].size(); i ++)
	{
		// current version of state variables
		unsigned int index = _block_vars[0][i];
		//int buddy_index = _cudd_to_buddy_vars[num];
		cudd_to_buddy_map[index] = bdd_var(result->_state_vars.get_var(0, i));
		//cout << "map: " << index << "  -> " << bdd_var(result->_state_vars.get_var(0, i)) << endl;
		//curr_state_vars.push_back(bdd_ithvar(buddy_index));
		//++ num;
		// next version of state variables
		index = _block_vars[1][i];
		//buddy_index = _cudd_to_buddy_vars[num];
		cudd_to_buddy_map[index] = bdd_var(result->_state_vars.get_var(1, i));
		//cout << "map: " << index << "  -> " << bdd_var(result->_state_vars.get_var(1, i)) << endl;
		//next_state_vars.push_back(bdd_ithvar(buddy_index));
		//++ num;
	}

	//cout << "Now remove useless variables" << endl;
	int pop_num = result->_state_vars.get_bdd_vars(0).size() - _block_vars[0].size();
	while(pop_num > 0)
	{
		result->_state_vars.pop_back_vars();
		-- pop_num;
	}
	//vector<vector<buddy_bdd>> state_vars;
	//state_vars.push_back(curr_state_vars);
	//state_vars.push_back(next_state_vars);
	//result._state_vars.add_bdd_vars(curr_state_vars);
	//result._state_vars.add_bdd_vars(next_state_vars);
	// compute mapping for labels
	for(unsigned i = _num_states_vars; i < _last_pos_labels; i ++)
	{
		cudd_to_buddy_map[i] = _sylvan_to_buddy_vars[i];
		//cout << "map: " << i << "  -> " << _cudd_to_buddy_vars[i] << endl;
	}
	cout << "Migrating DFA representation from Sylvan to BuDDy..." << endl;
	clock_t c_start = clock();
	//cout << "cudd init: " << _init << endl;
	// compute the initial states
	result->_init = move_to_buddy(_init, cudd_to_buddy_map);
	//cout << "buddy init: " << result->_init << endl;
	//cout << "cudd _finals: " << _finals << endl;
	result->_finals = move_to_buddy(_finals, cudd_to_buddy_map);
	//cout << "buddy finals: " << result->_finals << endl;
    cout << "Number of nodes in minimal transition is " << _trans.NodeCount() << endl;
	result->_trans = move_to_buddy(_trans, cudd_to_buddy_map);
	//cout << "buddy trans: " << result->_trans << endl;
	clock_t c_end = clock();
	cout << "Finished migrating DFA representation from Sylvan to BuDDy in "
		 << 1000.0 * (c_end - c_start)/CLOCKS_PER_SEC << "ms ..." << endl;

	result->_curr_cube = result->_state_vars.get_cube(0);
	result->_next_cube = result->_state_vars.get_cube(1);
    // make pairs

	result->_curr_to_next_pairs = result->_state_vars.make_pair(0, 1);
	result->_next_to_curr_pairs = result->_state_vars.make_pair(1, 0);
	result->_reach = bddtrue;//result.explore();
    cout << "Finished computing DFA representation in BuDDy..." << endl;
    /*
    cout << "is_even(8) = " << is_even(8) << endl;
    cout << "is_even(9) = " << is_even(9) << endl;
    cout << "is_even(0) = " << is_even(0) << endl;
    cout << "is_even(1) = " << is_even(1) << endl;
    */
    return *result;
}

/*-------------------------------------------------------------------*/
// copy a BDD from BuDDy to Sylvan
/*-------------------------------------------------------------------*/  
sylvan_bdd
dfwa_min_sylvan::move_to_sylvan(buddy_bdd dd)
{
    if(dd == bddtrue)
    {
        return sylvan_bdd::bddOne();
    }else
    if(dd == bddfalse)
    {
        return sylvan_bdd::bddZero();
    }else
    {
    	//nodes_num = 0;
        // traverse a BDD data structure
    	unordered_map<int, sylvan_bdd> computed_table;
    	sylvan_bdd result = move_to_cudd(dd, computed_table);
        return result;
    }
}

sylvan_bdd
dfwa_min_sylvan::move_to_cudd(buddy_bdd dd, unordered_map<int, sylvan_bdd>& computed_table)
{
	int dd_id = dd.id();
	unordered_map<int, sylvan_bdd>::const_iterator it = computed_table.find(dd_id);
	if(it != computed_table.end())
	{
		return it->second;
	}
	if(dd == bddtrue)
	{
		//cout << "true: node id = " << dd_id << endl;
		return sylvan_bdd::bddOne();
	}else
	if(dd == bddfalse)
	{
		//cout << "false: node id = " << dd_id << endl;
		return sylvan_bdd::bddZero();
	}else
	{
		//nodes_num = nodes_num + 1;
		//cout << "the number of nodes we have visited: " << nodes_num << " node id = " << dd_id << endl;
		// traverse a BDD data structure
		sylvan_bdd high = move_to_cudd(bdd_high(dd), computed_table);
		sylvan_bdd low  = move_to_cudd(bdd_low(dd), computed_table);
		int buddy_index = bdd_var(dd);
		int cudd_var_index = _buddy_to_sylvan_vars[buddy_index];
		sylvan_bdd var_dd = sylvan_bdd::bddVar(cudd_var_index);
		sylvan_bdd result = var_dd.Ite(high, low);
		computed_table[dd_id] = result;
		return result;
	}
}

/*-------------------------------------------------------------------*/
// copy a BDD from Sylvan to BuDDy
/*-------------------------------------------------------------------*/
buddy_bdd
dfwa_min_sylvan::move_to_buddy(sylvan_bdd dd, unordered_map<int, int>& cudd_to_buddy_map)
{
	if (dd.isOne())
	{
		return bddtrue;
	} else
	if (dd.isZero())
	{
		return bddfalse;
	} else
	{
		// traverse a BDD data structure
		unordered_map<sylvan_BDD, buddy_bdd> computed_table;
		buddy_bdd result = move_to_buddy(dd, cudd_to_buddy_map, computed_table);
		return result;
	}
}

buddy_bdd
dfwa_min_sylvan::move_to_buddy(sylvan_bdd_ptr dd, unordered_map<int, int>& cudd_to_buddy_map
		, unordered_map<sylvan_BDD, buddy_bdd>& computed_table)
{
	sylvan_BDD node = dd.GetBDD();
	unordered_map<sylvan_BDD, buddy_bdd>::const_iterator it = computed_table.find(node);
	if(it != computed_table.end())
	{
		return it->second;
	}
	if (dd.isOne())
	{
		return bddtrue;
	} else
	if (dd.isZero())
	{
		return bddfalse;
	} else
	{
		// traverse a BDD data structure
		const unsigned int top_index = dd.TopVar();
		sylvan_bdd cudd_var_dd = sylvan_bdd::bddVar(top_index);
		sylvan_bdd high_cofactor = dd.Then();
		buddy_bdd high = move_to_buddy(high_cofactor, cudd_to_buddy_map, computed_table);
		sylvan_bdd low_cofactor = dd.Else();
		buddy_bdd low = move_to_buddy(low_cofactor, cudd_to_buddy_map, computed_table);

		int buddy_var_index = cudd_to_buddy_map[top_index];
		//cout << "buddy var index: " << buddy_var_index << endl;
		buddy_bdd var_dd = bdd_ithvar(buddy_var_index);
		buddy_bdd result = bdd_ite(var_dd, high, low);
		computed_table[node] = result;
		return result;
	}
}

// ------------------ printing functions ------------------------------

void
dfwa_min_sylvan::print(sylvan_bdd bdd)
{
    //cout << "formula : " << bdd << endl;
    vector<sylvan_bdd> vec;
}

/*-------------------------------------------------------------------*/
// generate all possible truth values of a vector of variables
/*-------------------------------------------------------------------*/
void
dfwa_min_sylvan::generate_all_bits(vector<sylvan_bdd>& vars, unsigned index, sylvan_bdd_ptr dd, sylvan_bdd temp)
{ 
#ifdef DEBUG
    if (index == vars.size()) 
    {
        sylvan_bdd inter = temp & dd;
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
