#include "dfwamin.hh"

/*-------------------------------------------------------------------*/
// prepare state, label and block variables for DFA minimization
/*-------------------------------------------------------------------*/
void
dfwa_min::prepare_variables()
{
    // 0. create duplicate label variables for dfa minimization
    // A' duplicate variables for labels
	// temporary variables
    vector<int> label_vars;
    get_list_var_indices(label_vars, _aut._label_cube);
    for(unsigned i = 0; i < label_vars.size(); i ++)
    {
    	int buddy_index = label_vars[i];
        _cudd_to_buddy_vars.push_back(buddy_index);
        // from buddy to cudd for migration
    	_buddy_to_cudd_vars[buddy_index] = i;
    	        // new label variables
    	cudd_bdd var_dd = _manager->bddVar(i);
    }

    // 1. create state variables for dfa minimization
    // (X, X')
    cout << "buddy state cube: " << endl;
    // current and next state variable vectors
    vector<cudd_bdd> curr_states;
    _state_vars.push_back(curr_states);
    vector<cudd_bdd> next_states;
    _state_vars.push_back(next_states);

    // cube bdd
    _curr_cube = _manager->bddOne();
    _next_cube = _manager->bddOne();
    
    _num_states_vars = _aut._state_vars.get_var_num(0);
    // for permuting variables (starting from label_vars.size())
    // ordered state variables will speed up migration
    int cudd_index = label_vars.size();
    for(unsigned i = 0; i < _num_states_vars; i ++)
    {
    	// current state
    	int buddy_index = _aut._state_vars.get_var_id(0, i);
        _cudd_to_buddy_vars.push_back(buddy_index);
        _buddy_to_cudd_vars[buddy_index] = cudd_index;
        //cout << "buddy: " << buddy_index << " <-> " << _buddy_to_cudd_vars[buddy_index] << endl;
        // new state variables
        cudd_bdd var_dd = _manager->bddVar(cudd_index);
        _state_vars[0].push_back(var_dd);
        _curr_cube = _curr_cube & var_dd;

        ++ cudd_index;

        // next state
        buddy_index = _aut._state_vars.get_var_id(1, i);
        _cudd_to_buddy_vars.push_back(buddy_index);
        _buddy_to_cudd_vars[buddy_index] = cudd_index;
        //cout << "buddy: " << buddy_index << " <-> " << _buddy_to_cudd_vars[buddy_index] << endl;
                // new state variables
        var_dd = _manager->bddVar(cudd_index);
        _state_vars[1].push_back(var_dd);
        _next_cube = _next_cube & var_dd;

        ++ cudd_index;
    }
    // end of state variables (number of state variables)
    _num_states_vars = 2 * _num_states_vars + label_vars.size();
    // 2. create label variables
    // be careful about not changing _aut labels
    // A: the real label variables in minimization

    _label_cube = _manager->bddOne();
    for(unsigned i = 0; i < label_vars.size(); i ++)
    {
    	int buddy_index = label_vars[i];
    	int cudd_index = i + _num_states_vars;
        _cudd_to_buddy_vars.push_back(buddy_index);
        // we do not use this copy of label variables for migration
        // from buddy to cudd
        //_buddy_to_cudd_vars[buddy_index] = cudd_index;
        // new label variables
        cudd_bdd var_dd = _manager->bddVar(cudd_index);
        _label_cube = _label_cube & var_dd;
        _label_vars.push_back(var_dd);
    }
    // end of label variables
    _last_pos_labels = label_vars.size() + _num_states_vars;
    
    // 3. create block vectors
    // K block variables
    vector<cudd_bdd> curr_block;
    _block_vars.push_back(curr_block);
    vector<cudd_bdd> next_block;
    _block_vars.push_back(next_block);
    
    unsigned num_block_vars = _num_states_vars - label_vars.size();
    // prepare block variables
    for(unsigned i =0; i < num_block_vars; i ++)
    {
        cudd_bdd var_dd = _manager->bddVar();
        //cout << "block variable " << var_dd << endl;
        if(is_even(i))
        {
            _block_vars[0].push_back(var_dd);
        }else
        {
            _block_vars[1].push_back(var_dd);
        }
    }
}

/*-------------------------------------------------------------------*/
// prepare CUDD DFA minimization
/*-------------------------------------------------------------------*/
void
dfwa_min::prepare()
{
    // create variables of dfwa_min::
	//cout << "minimized before prepare " << _aut._label_cube << endl;
    prepare_variables();
    
    // disable dynamic variable ordering
    _manager->AutodynDisable();
    // now reorder things
    const unsigned var_num = _manager->ReadSize();
    int* order = new int[var_num];
    // first is labels
    for(unsigned i = 0; i < var_num; i ++)
    {
        order[i] = i;
    }
    _manager->ShuffleHeap(order);
#ifdef DEBUG
    cout << "buddy aut _init: " << _aut._init << endl;
    bdd_print_set(cout, _aut.get_dict(), _aut._init);
    cout << endl;
#endif
    cout << "Migrating BuDDy to CUDD ..." << endl;
    // copy the initial state
    clock_t c_start = clock();
    _init = move_to_cudd(_aut._init);

    //cout << "minimized init: " << _init << endl;
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
    //cout << "Number of nodes in transition after reach is " << bdd_nodecount(_aut._trans) << endl;
    _trans = move_to_cudd(_aut._trans);
    cout << "Number of nodes in transition after renaming is: " << bdd_nodecount(_trans) << endl;
    // copy the final states
    //cout << "minimized trans: " << _trans << endl;
    _finals = move_to_cudd(_aut._finals);
    //cout << "minimized finals " << _finals << endl;
    //cout << "minimized labels " << _aut._label_cube << endl;
    clock_t c_end = clock();
        cout << "Finished migrating BuDDy to CUDD in "
        	 << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << "ms ..."
			 << endl;
    // change label variables to block variables
    vector<cudd_bdd> prev_label_vec;
    for(unsigned i = 0; i < _label_vars.size(); i ++)
    {
    	prev_label_vec.push_back(_manager->bddVar(i));
    }
    vector<cudd_bdd> curr_label_vec;
    for(unsigned i = _num_states_vars; i < _last_pos_labels; i ++)
    {
    	curr_label_vec.push_back(_manager->bddVar(i));
    }
    assert (prev_label_vec.size() == curr_label_vec.size());

    _trans = cudd_permute(_trans, prev_label_vec, curr_label_vec);
    c_end = clock();
    cout << "Finished renaming CUDD in "
        	 << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << "ms ... " << _trans.nodeCount() << endl;
    //cout << "Now restore the ordering for minimization" << endl;
    //_manager->ShuffleHeap(order);
    c_end = clock();
    cout << "Finished migrating BuDDy to CUDD in "
    	 << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << "ms ..." << endl;


}

dfwa_min::~dfwa_min()
{
}

void
dfwa_min::output(ostream& os)
{
    os << "dfa: " << endl;
    os << "init: " << endl;
    os << _init << endl;
    vector<cudd_bdd> vec;
    vec.push_back(_init);
    _manager->DumpDot(vec);
    vec.pop_back();
    os << endl;
    
    os << "trans: " << endl;
    os << _trans << endl;
    vec.push_back(_trans);
    _manager->DumpDot(vec);
    vec.pop_back();
    os << endl;
    
    os << "finals: " << endl;
    os << _finals << endl;
    vec.push_back(_finals);
    _manager->DumpDot(vec);
    vec.pop_back();
    os << endl;
}
/*-------------------------------------------------------------------*/
// permute variables for from and to
/*-------------------------------------------------------------------*/ 
cudd_bdd
dfwa_min::cudd_permute(cudd_bdd_ptr dd, vector<cudd_bdd>& from, vector<cudd_bdd>& to)
{
    // get the number of variables in manager
    const unsigned var_num = _manager->ReadSize();
    int* permute = new int[var_num];
    for(unsigned i = 0; i < var_num; i ++)
    {
        permute[i] = i;
    }
    // swap the variables' indices
    for(unsigned i = 0; i < from.size(); i ++)
    {
        const unsigned from_index = from[i].NodeReadIndex();
        const unsigned to_index = to[i].NodeReadIndex();
        permute[from_index] = to_index;
        permute[to_index] = from_index;
    }
    cudd_node_ptr result = Cudd_bddPermute(_manager->getManager(), dd.getNode(), permute);
    delete[] permute;
    return cudd_bdd(*_manager, result);
}

// -------------- Symbolic DFA Minimization Algorithms ----------------

/*-------------------------------------------------------------------*/
// create num block variables
/*-------------------------------------------------------------------*/

cudd_bdd
dfwa_min::new_block_number(int block_number, int copy)
{
    // now check whether block_number can be represented with current
    // number of block variables
    //cout << "block_number = " << block_number << endl;
    cudd_bdd dd = _manager->bddOne();
    int bit = 1;
    for (cudd_bdd& bit_var : _block_vars[copy]) 
    {
        cudd_bdd bit_var_not = !bit_var;
        dd = dd & ((block_number & bit) != 0 ? bit_var : bit_var_not);
        bit <<= 1;
    }
    //cout << "returned block dd = " << dd << endl;
    
    return dd;
}

/*-------------------------------------------------------------------*/
// traverse BDD structure to assign block number
// NOTE: do not use regular node
// will need the number of states amount of times 
// NOT USED ANY MORE
/*-------------------------------------------------------------------*/ 
void
dfwa_min::compute_block_number_rec(cudd_bdd_ptr dd, unordered_map<cudd_node_ptr, int> &table, int & block_number)
{
    cudd_node_ptr node = dd.getNode();
    //cout << "Compute new partition from sigf(s, a, k) in compute_block_number_rec  " << endl;
    // check whether node has been visited before
    unordered_map<cudd_node_ptr, int>::const_iterator it = table.find(node);
    if(it != table.end())
    {
        //cout << "Found in table  " << endl;
        return ;
    }
    vector<cudd_bdd> label_block(_label_vars);
    label_block.insert(label_block.end(), _block_vars[0].begin(), _block_vars[0].end());
    if (dd.IsOne()) 
    {
      //cout << "Constant node computed  " << dd << endl;
      // never reach a constant node except the leaf-0.
      table[node] = block_number;
      //cout << "block number = " << block_number << endl;
      generate_all_bits(label_block, 0, dd, _manager->bddOne());
      ++ block_number;
    }else 
    // must not be constant zero
    if(! dd.IsZero())
    {
        //cout << "subformula: " << dd << endl;
        const unsigned int top_index = dd.NodeReadIndex();
        // check whether it is a state variable
        if (! is_state_variable(top_index)) 
        {
          //cout << "Non-state variables identified  " << endl;
          //Reached a block (signature). Compute its block number,
          table[node] = block_number;
          //cout << "block number = " << block_number << endl;
          generate_all_bits(label_block, 0, dd, _manager->bddOne());
          ++ block_number;
        } else 
        {
          // Traverse the BDD further.
          //cout << "State variables identified  " << endl;
          cudd_bdd dd_var = _manager->bddVar(top_index);
          cudd_bdd low_cofactor = dd.Cofactor(!dd_var);
          // low branch in CUDD
          compute_block_number_rec(low_cofactor, table, block_number);
          cudd_bdd high_cofactor = dd.Cofactor(dd_var);
          // high branch in CUDD
          compute_block_number_rec(high_cofactor, table, block_number);
        }
    }
}
/*-------------------------------------------------------------------*/
// traverse BDD structure to compute refined partition
// 1. increase the number of block variables on-the-fly: will need exponential
//    more times of BDD traversing
// 2. set to the number of state variables and then remove unessential variables:
//    much less times of BDD traversing
/*-------------------------------------------------------------------*/ 
cudd_bdd
dfwa_min::refine_partition_rec(cudd_bdd_ptr dd, unsigned & block_number
, unordered_map<cudd_node_ptr, cudd_bdd> &computed_table)
{
    cudd_node_ptr node = dd.getNode();
    //cout << "Compute new partition from sigf(s, a, k) in refine_partition_rec  " << endl;
    // check whether the node has been visited before
    unordered_map<cudd_node_ptr, cudd_bdd>::const_iterator it = computed_table.find(node);
    if(it != computed_table.end())
    {
        return it->second;
    }
    // false means cannot reach a block
    if(dd.IsZero())
    {
    	return dd;
    }

    cudd_bdd result;
    const unsigned int top_index = dd.NodeReadIndex();
    // check whether it is a state variable
    if(! is_state_variable(top_index) || dd.IsOne())
    {
    	// This is a node which representing a block
        //cout << "block  " << block_number << endl;
        result = new_block_number(block_number, 0);
        ++ block_number;

        // debug info
        vector<cudd_bdd> label_block(_label_vars);
        label_block.insert(label_block.end(), _block_vars[0].begin(), _block_vars[0].end());
        generate_all_bits(label_block, 0, dd, _manager->bddOne());
    }else
    {
        // Traverse the BDD further.
        //cout << "State variables identified  " << endl;
        cudd_bdd dd_var = _manager->bddVar(top_index);
        cudd_bdd low_cofactor = dd.Cofactor(!dd_var);
        cudd_bdd low  = refine_partition_rec(low_cofactor, block_number, computed_table);
        cudd_bdd high_cofactor = dd.Cofactor(dd_var);
        cudd_bdd high = refine_partition_rec(high_cofactor, block_number, computed_table);
        result = dd_var.Ite(high, low);
    }
    computed_table[node] = result;
    return result;
}

pair<cudd_bdd, unsigned>
dfwa_min::refine_partition(cudd_bdd_ptr sig)
{
	unsigned int block_number = 0;
    unordered_map<cudd_node_ptr, cudd_bdd> computed_table;
    cudd_bdd partition = refine_partition_rec(sig, block_number, computed_table);
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
cudd_bdd
dfwa_min::next_image(cudd_bdd_ptr curr)
{
	cudd_bdd next = curr.AndAbstract(_trans, _curr_cube & _label_cube);
	// next to current
    next = cudd_permute(next, _state_vars[1], _state_vars[0]);
    return next;
}

/*-------------------------------------------------------------------*/
// compute previous image of curr
/*-------------------------------------------------------------------*/
cudd_bdd
dfwa_min::pre_image(cudd_bdd_ptr curr)
{
	cudd_bdd next = cudd_permute(curr, _state_vars[0], _state_vars[1]);
    // and exist
	cudd_bdd pre = next.AndAbstract(_trans, _next_cube  & _label_cube);
    return pre;
}
/*-------------------------------------------------------------------*/
// compute reachable state space
/*-------------------------------------------------------------------*/
cudd_bdd
dfwa_min::forward_explore()
{
	clock_t c_start = clock();
    cudd_bdd s = _init;
    cudd_bdd sp = _manager->bddZero();
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
    cout << "Finished backward exploration in " << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " " << s.nodeCount() << " ms\n";
    return s;
}

cudd_bdd
dfwa_min::backward_explore()
{
	clock_t c_start = clock();
	cudd_bdd s = _finals;
	cudd_bdd sp = _manager->bddZero();
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
    cout << "Finished backward exploration in " << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " " << s.nodeCount() << " ms\n";
    return s;
}

cudd_bdd
dfwa_min::initialize_partition()
{

	clock_t c_start = clock();
    cudd_bdd reach = forward_explore();
    reach = reach & backward_explore();
    clock_t c_end = clock();
    //cout << "Finished cudd exploration in " << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms " << reach.nodeCount()<< "\n";
    // buddy is a bit faster than cudd. Not sure whether it is due to the encoding or
    // that buddy is just faster than cudd
#ifdef DEBUG
	c_start = clock();
	buddy_bdd b_reach = _aut.explore();
	b_reach = b_reach & _aut.back_explore();
	cudd_bdd copy = move_to_cudd(b_reach);
	c_end = clock();
	cout << "Finished buddy exploration in " << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms " << bdd_nodecount(b_reach)<< "\n";
#endif

    //cout << "reachable non-final states: " << endl;
    cudd_bdd non_finals = (! _finals) & reach;
    vector<cudd_bdd> curr_states(_state_vars[0]);
    generate_all_bits(curr_states, 0, non_finals, _manager->bddOne());
    // get block variables
    cudd_bdd block_0 = new_block_number(0, 0);
    // careful about nonreachable states
    cudd_bdd block_1 = new_block_number(1, 0);
    cudd_bdd partition = (block_0 & _finals & reach) | (block_1 & non_finals);
    // debug info
    //cout << "Initial partition: " << endl;
    curr_states.insert(curr_states.end(), _block_vars[0].begin(), _block_vars[0].end());
    generate_all_bits(curr_states, 0, partition, _manager->bddOne());
    // check whether we need to add sink transitions
    // s <-> s'
    //cudd_bdd nonreach = ! cudd_reach;
    cudd_bdd next_reach = cudd_permute(reach, _state_vars[0], _state_vars[1]);
    _trans = _trans & reach & next_reach;
    return partition;
}


/*-------------------------------------------------------------------*/
// symbolic DFA minimization
/*-------------------------------------------------------------------*/ 
void 
dfwa_min::minimize()
{
	clock_t c_start = clock();
	cout << "Starting DFA minimization..." << endl;
    // sigref based minimization for algorithms
    // 1. P(X, K) is the initial partition of states
    cudd_bdd partition = initialize_partition();
    //cout << "Nodes in partition = " << partition.nodeCount() << endl;
    unsigned prev_block_number = 2;
    unsigned iteration_num = 0;
    
    //cudd_bdd result_sigf = _manager->bddOne();
    while(true)
    {
    	clock_t iter_start = clock();
        ++ iteration_num;
        cout << "Number of blocks at iteration " << iteration_num << " is " << prev_block_number << endl;
        // compute the signature sigf(s, a, k) = exists X'. T(X, A, X') and P(X', K)
        // first P(X, K) => P(X', K)
        //cout << "Permute P(X, K) => P(X', K)  " << endl;
        cudd_bdd state_next_partition = cudd_permute(partition, _state_vars[0], _state_vars[1]);
        //cout << "Nodes in partition = " << state_next_partition.nodeCount() << endl;
        // sigf(s, a, k) = exists X'. T(X, A, X') and P(X', K)
        //cout << "Compute sigf(s, a, k) = exists X'. T(X, A, X') and P(X', K)  " << endl;
        cudd_bdd sigf = state_next_partition.AndAbstract(_trans, _next_cube);
        //cout << "Nodes in partition = " << sigf.nodeCount() << endl;
        //cout << "Permute P(X, K) => P(X, K')  " << endl;
        cudd_bdd block_next_partition = cudd_permute(partition, _block_vars[0], _block_vars[1]);
        //cout << "Nodes in partition = " << block_next_partition.nodeCount() << endl;
        //cout << "sigf(s, a, k) = sigf(s, a, k) & P(s, k')" << endl;
        sigf = sigf & block_next_partition;
        //result_sigf = sigf;
        //cout << "Nodes in partition = " << sigf.nodeCount() << endl;
        //sigf(s, a, k) [X, A, K]
        vector<cudd_bdd> new_state_label_block_1(_state_vars[0]);
        new_state_label_block_1.insert(new_state_label_block_1.end(), _label_vars.begin(), _label_vars.end());
        new_state_label_block_1.insert(new_state_label_block_1.end(), _block_vars[0].begin(), _block_vars[0].end());
        generate_all_bits(new_state_label_block_1, 0, sigf, _manager->bddOne());
        // refine the partition, traverse the BDD to compute block numbers
        cout << "Compute new partition from signature..." << endl;
        //print(sigf);
        pair<cudd_bdd, unsigned> result = refine_partition(sigf);
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
        //cout << "Nodes in partition = " << partition.nodeCount() << endl;
        // debug info
        //cout << "refined partition: " << endl << partition << endl;
        vector<cudd_bdd> new_state_label_block(_state_vars[0]);
        new_state_label_block.insert(new_state_label_block.end(), _block_vars[0].begin(), _block_vars[0].end());
        generate_all_bits(new_state_label_block, 0, partition, _manager->bddOne());
    }

    //cout << "T_min(k', a, k) = exists X, X'. P(s, k') and T(X, A, X') and P(X', K) " << endl;
    //result_sigf = result_sigf.ExistAbstract(_curr_cube);
    //cout << "#T_min(k', a, k) = " << result_sigf.nodeCount() << endl;

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
    }*/
    
    //cout << "label cube = " << _label_cube << endl;
    //cout << "init: " << endl;
    //cout << "before: " << _init << endl;
    //cout << "_curr_cube =" << _curr_cube << endl;
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
    vector<cudd_bdd> curr_states(_state_vars[0]);
    curr_states.insert(curr_states.end(), _block_vars[0].begin(), _block_vars[0].end());
    vector<cudd_bdd> next_states(_state_vars[1]);
    next_states.insert(next_states.end(), _block_vars[1].begin(), _block_vars[1].end());
    
    //cout << "P(X, K) => P(X', K')" << endl;
    cudd_bdd next_partition = cudd_permute(partition, curr_states, next_states);
    cout << "Computing the transition relation for minimal DFA..." << endl;
    clock_t t_start = clock();

    // R(s, a, k') := (∃t : T (s, a, t) ∧ P(t, k'))
    //cout << "R(s, a, k') := ∃t : T (s, a, t) ∧ P(t, k')" << endl;
    //cout << "#nodes in partition: " << partition.nodeCount() << endl;
    //cout << "#nodes in trans: " << _trans.nodeCount() << endl;

    cudd_bdd minimized_trans = _trans.AndAbstract(next_partition, _next_cube);
    // T P (s, a, t') := (∃s : R(s, a, k') ∧ P(s, k))
    //cout << "T(k, a, k') := ∃s :R(s, a, k') ∧ P(s, k)" << endl;
    minimized_trans = minimized_trans.AndAbstract(partition, _curr_cube);
    //cout << minimized_trans << endl;
    /*
    cudd_bdd minimized_trans = make_quotient_trans(_trans, partition, next_partition);
     */
    clock_t t_end = clock();
    cout << "Finished transition relation computation in " << 1000.0 * (t_end - t_start) / CLOCKS_PER_SEC << " ms\n";
    //cout << "Output Transitions" << endl;
    vector<cudd_bdd> label_block(_label_vars);
    label_block.insert(label_block.end(), _block_vars[0].begin(), _block_vars[0].end());
    label_block.insert(label_block.end(), _block_vars[1].begin(), _block_vars[1].end());
    generate_all_bits(label_block, 0, minimized_trans, _manager->bddOne());
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
cudd_bdd
dfwa_min::make_quotient_trans(cudd_bdd_ptr trans, cudd_bdd_ptr curr_partition, cudd_bdd_ptr next_partition)
{
	// T(X, X', A), P(X, K), P(X', K')
	map_t computed_table;
	set<unsigned> curr_set;
	set<unsigned> next_set;
	vector<set<unsigned>> var_sets;
	var_sets.push_back(curr_set);
	var_sets.push_back(next_set);
	for(unsigned i = 0; i < var_sets.size(); i ++)
	{
		for(unsigned j = 0; j < _state_vars[i].size(); j ++)
		{
			var_sets[i].insert(_state_vars[i][j].NodeReadIndex());
		}
	}

	return make_quotient_trans(trans, curr_partition, next_partition, var_sets, computed_table);
}
cudd_bdd
dfwa_min::make_quotient_trans(cudd_bdd trans, cudd_bdd curr_partition, cudd_bdd next_partition
		, vector<set<unsigned>>& var_sets, map_t& computed_table)
{
	//cout << "Entering recursive construction for transition relation..." << endl;
	cudd_node_ptr first = trans.getNode();
	cudd_node_ptr second = curr_partition.getNode();
	cudd_node_ptr third = next_partition.getNode();
	//cout << "f0 = " << first << " f1 = " << second << " f2 = " << third << endl;
	map_t::const_iterator it = computed_table.find(make_tuple(first, second, third));
	if(it != computed_table.end())
	{
		return it->second;
	}
	if(trans.IsZero() || curr_partition.IsZero() || next_partition.IsZero())
	{
		return _manager->bddZero();
	}
	//cout << "Needs to do construction for transition relation..." << endl;
	//cout << "f0 = " << trans << " f1 = " << curr_partition << " f2 = " << next_partition << endl;
	// computed variables
	unsigned top_var = min(trans.NodeReadIndex(), curr_partition.NodeReadIndex());
	top_var = min(top_var, next_partition.NodeReadIndex());
	cudd_bdd top_var_dd = _manager->bddVar(top_var);
	//cout << "top_var = " << top_var << endl;
	cudd_bdd result;
	if(var_sets[0].find(top_var) != var_sets[0].end())
	{
		//cout << "current variable: " << endl;
		cudd_bdd low = make_quotient_trans(
				  trans.Cofactor(! top_var_dd)
				, curr_partition.Cofactor(! top_var_dd)
				, next_partition
				, var_sets, computed_table);
		cudd_bdd high = make_quotient_trans(
						  trans.Cofactor(top_var_dd)
						, curr_partition.Cofactor(top_var_dd)
						, next_partition, var_sets, computed_table);
		result = low | high;
	}else
	if(var_sets[1].find(top_var) != var_sets[1].end())
	{
		//cout << "next variable: " << endl;
		cudd_bdd low = make_quotient_trans(
				  trans.Cofactor(! top_var_dd)
				, curr_partition
				, next_partition.Cofactor(! top_var_dd)
				, var_sets, computed_table);
		cudd_bdd high = make_quotient_trans(
						  trans.Cofactor(top_var_dd)
						, curr_partition
						, next_partition.Cofactor(top_var_dd)
						, var_sets, computed_table);
		result = low | high;
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
dfwa_min::reduce_block_variables(cudd_bdd_ptr partition)
{
    // determine whether we can remove some variables
	cout << "Reducing block variables..." << endl;
	unsigned num = 0;
    while(true)
    {
        cudd_bdd var_dd = _block_vars[0].back();
        cudd_bdd result = (! var_dd) & partition;
        if(result != partition || _block_vars[0].size() < 2)
        {
        	break;
        }
        ++ num;
        _block_vars[0].pop_back();
        /*
        cout << "number of vars: " << _num_states_vars << endl
             << "block[0].size() = " << _block_vars[0].size() << endl
			 << "block[1].size() = " << _block_vars[1].size() << endl;
        */
        _block_vars[1].pop_back();
        // remove redundant variables
        partition = partition.ExistAbstract(var_dd);
    }
    cout << "Finished reducing " << num << " block variables..." << endl;
}

/*-------------------------------------------------------------------*/
// move the minimal DFA from CUDD to BuDDy
// the ordering of dfa is (A, X, X')
/*-------------------------------------------------------------------*/
dfwa_ptr
dfwa_min::move_dfwa()
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
		unsigned int index = _block_vars[0][i].NodeReadIndex();
		//int buddy_index = _cudd_to_buddy_vars[num];
		cudd_to_buddy_map[index] = bdd_var(result->_state_vars.get_var(0, i));
		//cout << "map: " << index << "  -> " << bdd_var(result->_state_vars.get_var(0, i)) << endl;
		//curr_state_vars.push_back(bdd_ithvar(buddy_index));
		//++ num;
		// next version of state variables
		index = _block_vars[1][i].NodeReadIndex();
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
		cudd_to_buddy_map[i] = _cudd_to_buddy_vars[i];
		//cout << "map: " << i << "  -> " << _cudd_to_buddy_vars[i] << endl;
	}
	cout << "Migrating DFA representation from CUDD to BuDDy..." << endl;
	clock_t c_start = clock();
	//cout << "cudd init: " << _init << endl;
	// compute the initial states
	result->_init = move_to_buddy(_init, cudd_to_buddy_map);
	//cout << "buddy init: " << result->_init << endl;
	//cout << "cudd _finals: " << _finals << endl;
	result->_finals = move_to_buddy(_finals, cudd_to_buddy_map);
	//cout << "buddy finals: " << result->_finals << endl;
    cout << "Number of nodes in minimal transition is " << _trans.nodeCount() << endl;
	result->_trans = move_to_buddy(_trans, cudd_to_buddy_map);
	//cout << "buddy trans: " << result->_trans << endl;
	clock_t c_end = clock();
	cout << "Finished migrating DFA representation from CUDD to BuDDy in "
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
// copy a BDD from BuDDy to CUDD
/*-------------------------------------------------------------------*/  
cudd_bdd
dfwa_min::move_to_cudd(buddy_bdd dd)
{
    if(dd == bddtrue)
    {
        return _manager->bddOne();
    }else
    if(dd == bddfalse)
    {
        return _manager->bddZero();
    }else
    {
    	nodes_num = 0;
        // traverse a BDD data structure
    	unordered_map<int, cudd_bdd> computed_table;
    	cudd_bdd result = move_to_cudd(dd, computed_table);
        return result;
    }
}

cudd_bdd
dfwa_min::move_to_cudd(buddy_bdd dd, unordered_map<int, cudd_bdd>& computed_table)
{
	int dd_id = dd.id();
	unordered_map<int, cudd_bdd>::const_iterator it = computed_table.find(dd_id);
	if(it != computed_table.end())
	{
		return it->second;
	}
	if(dd == bddtrue)
	{
		//cout << "true: node id = " << dd_id << endl;
		return _manager->bddOne();
	}else
	if(dd == bddfalse)
	{
		//cout << "false: node id = " << dd_id << endl;
		return _manager->bddZero();
	}else
	{
		nodes_num = nodes_num + 1;
		//cout << "the number of nodes we have visited: " << nodes_num << " node id = " << dd_id << endl;
		// traverse a BDD data structure
		cudd_bdd high = move_to_cudd(bdd_high(dd), computed_table);
		cudd_bdd low  = move_to_cudd(bdd_low(dd), computed_table);
		int buddy_index = bdd_var(dd);
		int cudd_var_index = _buddy_to_cudd_vars[buddy_index];
		cudd_bdd var_dd = _manager->bddVar(cudd_var_index);
		cudd_bdd result = var_dd.Ite(high, low);
		computed_table[dd_id] = result;
		return result;
	}
}

/*-------------------------------------------------------------------*/
// copy a BDD from CUDD to BuDDy
/*-------------------------------------------------------------------*/
buddy_bdd
dfwa_min::move_to_buddy(cudd_bdd dd, unordered_map<int, int>& cudd_to_buddy_map)
{
	if (dd.IsOne())
	{
		return bddtrue;
	} else
	if (dd.IsZero())
	{
		return bddfalse;
	} else
	{
		// traverse a BDD data structure
		unordered_map<cudd_node_ptr, buddy_bdd> computed_table;
		buddy_bdd result = move_to_buddy(dd, cudd_to_buddy_map, computed_table);
		return result;
	}
}

buddy_bdd
dfwa_min::move_to_buddy(cudd_bdd_ptr dd, unordered_map<int, int>& cudd_to_buddy_map
		, unordered_map<cudd_node_ptr, buddy_bdd>& computed_table)
{
	cudd_node_ptr node = dd.getNode();
	unordered_map<cudd_node_ptr, buddy_bdd>::const_iterator it = computed_table.find(node);
	if(it != computed_table.end())
	{
		return it->second;
	}
	if (dd.IsOne())
	{
		return bddtrue;
	} else
	if (dd.IsZero())
	{
		return bddfalse;
	} else
	{
		// traverse a BDD data structure
		const unsigned int top_index = dd.NodeReadIndex();
		cudd_bdd cudd_var_dd = _manager->bddVar(top_index);
		cudd_bdd high_cofactor = dd.Cofactor(cudd_var_dd);
		cudd_bdd low_cofactor = dd.Cofactor(!cudd_var_dd);

		buddy_bdd high = move_to_buddy(high_cofactor, cudd_to_buddy_map, computed_table);
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
dfwa_min::print(cudd_bdd bdd)
{
    cout << "formula : " << bdd << endl;
    vector<cudd_bdd> vec;
    vec.push_back(bdd);
    _manager->DumpDot(vec);
    vec.pop_back();
    cout << endl;
}

/*-------------------------------------------------------------------*/
// generate all possible truth values of a vector of variables
/*-------------------------------------------------------------------*/
void
dfwa_min::generate_all_bits(vector<cudd_bdd>& vars, unsigned index, cudd_bdd_ptr dd, cudd_bdd temp)
{ 
#ifdef DEBUG
    if (index == vars.size()) 
    {
        cudd_bdd inter = temp & dd;
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
