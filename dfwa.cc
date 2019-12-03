#include "dfwa.hh"

/*-------------------------------------------------------------------*/
// initialize dd representation of a product DFA
/*-------------------------------------------------------------------*/
dfwa::dfwa(bdd_dict_ptr dict, bdd label_cube)
:_state_vars(dict), _label_cube(label_cube)
{
    _aut = nullptr;
}

struct GreaterThanByBddSize
{
  bool operator()(bdd& f1, bdd& f2) const
  {
    int size_1 = bdd_nodecount(f1);
    int size_2 = bdd_nodecount(f2);
    if(size_1 < size_2)
    {
        return false;
    }
    return size_1 >= size_2;
  }
};

/*-------------------------------------------------------------------*/
// binary operation construction for computing the disjunction of BDDs
/*-------------------------------------------------------------------*/
bdd
binary_disjunct_rec(vector<bdd> & trans, vector<bdd> & curr_st_bdd, unsigned left, unsigned right)
{
    if(left > right)
    {
        return bddfalse;
    }else
    if(left == right)
    // return (s, a, t)
    {
        return curr_st_bdd[left] & trans[left];
    }else
    // left < right
    {
        unsigned mid = (left + right) / 2;
        bdd op1 = binary_disjunct_rec(trans, curr_st_bdd, left, mid);
        bdd op2 = binary_disjunct_rec(trans, curr_st_bdd, mid + 1, right);
        return op1 | op2;
    }
}
/*-------------------------------------------------------------------*/
// initialize dd representation of a DFA
// make complete dfa
/*-------------------------------------------------------------------*/
dfwa::dfwa(twa_graph_ptr aut, bdd label_cube, set<unsigned>& finals, const char* name)
:_state_vars(aut->get_dict(), aut, 2, name, 0, aut->num_states() - 1), _label_cube(label_cube)
{
    _aut = aut;
    // label cube
    //cout << "construct state bdd #S = " << aut->num_states() << endl;
    vector<bdd> curr_st_bdd;
    vector<bdd> next_st_bdd;
    for(unsigned i = 0; i < aut->num_states(); i ++)
    {
        bdd dd = _state_vars.new_value(0, i);
        #ifdef DEBUG
        cout << "state i = " << i << endl;
        bdd_print_sat(cout, aut->get_dict(), dd);
        cout << endl;
        #endif
        curr_st_bdd.push_back(dd);
        dd = _state_vars.new_value(1, i);
        #ifdef DEBUG
        cout << "state i' = " << i << endl;
        bdd_print_sat(cout, aut->get_dict(), dd);
        cout << endl;
        #endif
        next_st_bdd.push_back(dd);
    }
    //cout << " curr_size = " << curr_st_bdd.size() << endl;
    _init = curr_st_bdd[aut->get_init_state_number()];
    #ifdef DEBUG
    //cout << "init  = " << aut->get_init_state_number() << endl;
    bdd_print_sat(cout, aut->get_dict(), _init);
    cout << endl;
    #endif
    _trans = bddfalse;
    _finals = bddfalse;
    // later if not needed, remove _reach
    _reach = bddfalse;
    cout << "Computing the transition relation..." << endl;
    // needs to be improved by huffman ?
    //priority_queue<bdd, std::vector<bdd>, GreaterThanByBddSize> pq;
    //vector<bdd> map2tr;
    for(unsigned s = 0; s < aut->num_states(); s ++)
    {
        _reach = _reach | curr_st_bdd[s];
        if(finals.find(s) != finals.end()) 
        {
			_finals = _finals | curr_st_bdd[s];
        }
        bdd sdd = curr_st_bdd[s];
        bdd tdd = bddfalse;
        //bdd outs = bddfalse;
        for(auto& tr : aut->out(s)) 
        {
            tdd = tdd | (tr.cond & next_st_bdd[tr.dst]); 
            //outs = outs | tr.cond;
        }
        sdd = sdd & tdd;
        /*
        cout << "tr state i = " << s << endl;
        bdd_print_sat(cout, aut->get_dict(), sdd);
        cout << endl;
        cout << endl;
        */
		//if(outs != bddtrue)
		//{
			// add a sink state
			//sdd = sdd | (curr_st_bdd[s] & (! outs) & next_st_bdd[aut->num_states()]);
		//}
        _trans = _trans | sdd;

        //map2tr.push_back(sdd);
        //cout << "encoding state = " << s << endl;
    }
    // add a sink state
    // _trans = _trans | (curr_st_bdd[aut->num_states()] & next_st_bdd[aut->num_states()]);
    //cout << "Binary construction for transition relation" << endl;
    // binary divide and conquer disjunction of bdds
    //_trans = binary_disjunct_rec(map2tr, curr_st_bdd, 0, map2tr.size() - 1);
    // compute transitions
    /*
    while(pq.size() > 1)
    {
        bdd f1 = pq.top();
        pq.pop();
        bdd f2 = pq.top();
        pq.pop();
        bdd f = f1 | f2;
        pq.push(f);
    }
    _trans = pq.top();
    pq.pop();
    */
    cout << "Finished computing the transition relation..." << endl;
    _curr_cube = _state_vars.get_cube(0);
    _next_cube = _state_vars.get_cube(1);
    // make pairs
    _curr_to_next_pairs = _state_vars.make_pair(0, 1);
    _next_to_curr_pairs = _state_vars.make_pair(1, 0);
    
    //bdd all = bdd_replace(_reach, _curr_to_next_pairs);
    //all = all & _reach;
    //_trans = _trans & all;
    #ifdef DEBUG
    cout << "trans = " << endl;
    bdd_print_sat(cout, aut->get_dict(), _trans);
    cout << endl;
    
    cout << "finals = " << endl;
    bdd_print_sat(cout, aut->get_dict(), _finals);
    cout << endl;
    #endif
}

dfwa::~dfwa()
{
    //cout << "HELLO start" << _curr_to_next_pairs << endl;
    if(_curr_to_next_pairs != nullptr)
    {
        bdd_freepair(_curr_to_next_pairs);
        _curr_to_next_pairs = nullptr;
    }
    //cout << "HELLO second" << _next_to_curr_pairs << endl;
    if(_next_to_curr_pairs != nullptr)
    {
        bdd_freepair(_next_to_curr_pairs);
        _next_to_curr_pairs = nullptr;
    }
    
    //cout << "HELLO end" << endl;
}

/*-------------------------------------------------------------------*/
// compute next step image of curr 
// FIXED (note the returned image contains propositions)
/*-------------------------------------------------------------------*/
bdd
dfwa::next_image(bdd curr)
{
    bdd next = bdd_relprod(_trans, curr, _curr_cube & _label_cube);
    next = bdd_replace(next, _next_to_curr_pairs);
    return next;
}

/*-------------------------------------------------------------------*/
// compute previous image of curr 
// FIXED (note the returned image contains propositions)
/*-------------------------------------------------------------------*/
bdd 
dfwa::pre_image(bdd curr)
{
    bdd next = bdd_replace(curr, _curr_to_next_pairs);
	bdd pre = bdd_relprod(_trans, next, _next_cube  & _label_cube);
    return pre;
}
/*-------------------------------------------------------------------*/
// compute reachable state space
/*-------------------------------------------------------------------*/
bdd
dfwa::explore()
{
    bdd s = _init;
    bdd sp = bddfalse;
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
    return s;
}

bdd
dfwa::back_explore()
{
    bdd s = _finals;
    bdd sp = bddfalse;
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
        cout << "The number of node in reverse R(" << count << ") is " << bdd_nodecount(s) << endl;
        ++ count;
    }
    return s;
}

bool
dfwa::is_empty()
{
	_reach = explore();
	return (_reach & _finals) == bddfalse;
}
/*-------------------------------------------------------------------*/
// get data from dfwa
/*-------------------------------------------------------------------*/
bdd
dfwa::get_init()
{
    return _init;
}

bdd
dfwa::get_trans()
{
    return _trans;
}

bdd
dfwa::get_finals()
{
    return _finals;
}

void
dfwa::output(ostream& os)
{
    os << "dfwa: " << endl;
    os << "init: " << endl;
    bdd_print_set(os, _state_vars.get_dict(), _init);
    os << endl;
    
    os << "trans: " << endl;
    bdd_print_set(os, _state_vars.get_dict(), _trans);
    os << endl;
    
    os << "finals: " << endl;
    bdd_print_set(os, _state_vars.get_dict(), _finals);
    os << endl;
    
}

void
dfwa::output_dfwa(ostream& os)
{
	os << "LISA DFA: " << endl;
	//os << "//ID VAR HIGH LOW " << endl;

	string alive_ap(ALIVE_AP);
	bdd_dict_ptr dict = get_dict();
	int index_alive = dict->varnum(formula::ap(alive_ap));
	bdd dd_alive = bdd_ithvar(index_alive);
	bdd label_cube = bdd_exist(_label_cube, dd_alive);

	vector<int> vars;
	get_list_var_indices(vars, label_cube);
	os << "LABEL VARS:";
	for(unsigned i = 0; i < vars.size(); i ++)
	{
		os << " " << vars[i];
	}
	os << endl;

	vars.clear();
	get_list_var_indices(vars, _curr_cube);
	os << "CURR STATE VARS:";
	for(unsigned i = 0; i < vars.size(); i ++)
	{
		os << " " << vars[i];
	}
	os << endl;

	vars.clear();
	get_list_var_indices(vars, _next_cube);
	os << "NEXT STATE VARS:";
	for(unsigned i = 0; i < vars.size(); i ++)
	{
		os << " " << vars[i] ;
	}
	os << endl;

	os << "INIT: " << _init.id() << endl;
	output_bdd(os, _init);

	os << "FINAL: " << _finals.id() << endl;
	output_bdd(os, _finals);

	bdd tr = bdd_exist(_trans, dd_alive);

	os << "TRANS: " << tr.id() << endl;
	output_bdd(os, tr);

}

void
dfwa::make_complete()
{
	bdd trans = bddtrue;
	for(unsigned i = 0; i < _state_vars._dd_vars[0].size(); i ++)
	{
		bdd var_0 = _state_vars._dd_vars[0][i];
		bdd var_1 = _state_vars._dd_vars[1][i];

		trans = trans & bdd_biimp(var_0, var_1);
	}

}
/*-------------------------------------------------------------------*/
// intersection dfwa:
// result is an empty dfwa
/*-------------------------------------------------------------------*/
void 
intersect_dfwa(dfwa_ptr result, dfwa_ptr op1, dfwa_ptr op2)
{
    // set the copies of state variables
    result._state_vars._copies = 2;
    // now we add variables from op1 and op2
    cout << "Computing the intersection product..." << endl;
    // check whether this part can be improved
    result._state_vars.add_bdd_vars(op1._state_vars);
    result._state_vars.add_bdd_vars(op2._state_vars);
    cout << "#AP1 = " << op1._state_vars._dd_vars[0].size() << " #AP2 = " <<  op2._state_vars._dd_vars[0].size() << endl;
    cout << "#PRO = " << result._state_vars._dd_vars[0].size() << endl;
    // 
    result._init = op1.get_init() & op2.get_init();
    // especially for the transition relation
    cout << "Computing transition relation in the intersection product..." << endl;
    result._trans = op1.get_trans() & op2.get_trans();
    cout << "Finished computing transition relation in the intersection product..." << endl;
    result._finals = op1.get_finals() & op2.get_finals();
    //result._reach = bddfalse;
    result._curr_cube = result._state_vars.get_cube(0);
    result._next_cube = result._state_vars.get_cube(1);
    // make pairs
    result._curr_to_next_pairs = result._state_vars.make_pair(0, 1);
    result._next_to_curr_pairs = result._state_vars.make_pair(1, 0);
    // compute reachable state space
    cout << "Computing reachable state space in the product..." << endl;
    result._reach = bddtrue;//result.explore();
    cout << "Finished computing reachable state space in the product..." << endl;
    //cout << "reachable states in product: " << endl;
    //bdd_print_set(cout, result._state_vars.get_dict(), result._reach);
    //cout << endl;
    // needs to whether this is useful
    /*
    bdd all = bdd_replace(result._reach, result._curr_to_next_pairs);
    all = all & result._reach;
    result._trans = result._trans & all;
    */
    cout << "Finished computing the intersection product..." << endl;
}

/*-------------------------------------------------------------------*/
// product for dfwa: func is passed for computing final states
// ASSUMPTION:
//  1. labels of op1 and op2 are the same
//  2. the state variables of op1 and op2 are different
/*-------------------------------------------------------------------*/
void
check_assumption(dfwa_ptr op1, dfwa_ptr op2)
{
	/*
	if(op1._label_cube != op2._label_cube)
	{
		cerr << "product: The propositions are not the same" << endl;
		cerr << "op1: ";
		bdd_print_set(cerr, op1._state_vars.get_dict(), op1._label_cube);
		cerr << endl << "op2: ";
		bdd_print_set(cerr, op1._state_vars.get_dict(), op2._label_cube);
		cerr << endl;
		exit(-1);
	}
	*/

	vector<bdd>& vars_1 = op1._state_vars.get_bdd_vars(0);
	vector<bdd>& vars_2 = op2._state_vars.get_bdd_vars(0);

	// check whether there are common variables
	set<int> vars;
	for(unsigned i = 0; i < vars_1.size(); i ++)
	{
		vars.insert(bdd_var(vars_1[i]));
	}
	for(unsigned i = 0; i < vars_2.size(); i ++)
	{
		if(vars.find(bdd_var(vars_2[i])) != vars.end())
		{
			cerr << "product: The state variables are not different -> ";
			bdd_print_set(cerr, op1._state_vars.get_dict(), vars_2[i]);
			cerr << endl;
			for(unsigned j = 0; j < vars_1.size(); j ++)
			{
				bdd_print_set(cerr, op1._state_vars.get_dict(), vars_1[j]);
			    cerr << endl;
			}
			exit(-1);
		}
	}
}
// keep the variables in increasing order
dfwa_ptr
product_dfwa(dfwa_ptr op1, dfwa_ptr op2, function<bdd(bdd&, bdd&)> func)
{
	check_assumption(op1, op2);

	bdd label_cube = op1._label_cube & op2._label_cube;
	cout << "labels in op1: " << endl;
	bdd_print_set(cout, op1.get_dict(), op1._label_cube);
	cout << endl;
	cout << "labels in op2: " << endl;
	bdd_print_set(cout, op1.get_dict(), op2._label_cube);
	cout << endl;
	dfwa* result = new dfwa(op1.get_dict(), label_cube);
	cout << "labels in product: " << endl;
	bdd_print_set(cout, op1.get_dict(), result->_label_cube);
	cout << endl;
    // set the copies of state variables
    result->_state_vars._copies = 2;
    // now we add variables from op1 and op2
    cout << "Computing the product..." << endl;
    // check whether this part can be improved
    result->_state_vars.add_bdd_vars(op1._state_vars);
    result->_state_vars.add_bdd_vars(op2._state_vars);
    cout << "state vars in op1: " << endl;
    bdd_print_set(cout, op1.get_dict(), op1._curr_cube);
    cout << endl;
    cout << "state vars in op2: " << endl;
    bdd_print_set(cout, op1.get_dict(), op2._curr_cube);
    cout << endl;
    cout << "#AP1 = " << op1._state_vars._dd_vars[0].size() << " #AP2 = " <<  op2._state_vars._dd_vars[0].size() << endl;
    cout << "#PRO = " << result->_state_vars._dd_vars[0].size() << endl;

    cout << "_init in op1: " << endl;
    //bdd_print_set(cout, op1.get_dict(), op1.get_init());
    cout << endl;
    cout << "_init in op2: " << endl;
    //bdd_print_set(cout, op1.get_dict(), op2.get_init());
    cout << endl;
    result->_init = op1.get_init() & op2.get_init();
    cout << "initial in product: " << endl;
    //bdd_print_set(cout, op1.get_dict(), result->_init);
    cout << endl;
    // especially for the transition relation
    cout << "#trans1 = " << bdd_nodecount(op1.get_trans()) << " #trans2 = " << bdd_nodecount(op2.get_trans()) << endl;
    cout << "Computing transition relation in the intersection product..." << endl;
    result->_trans = op1.get_trans() & op2.get_trans();
    cout << "trans in product: " << endl;
    //bdd_print_set(cout, op1.get_dict(), result->_trans);
    cout << endl;
    cout << "Finished computing transition relation in the intersection product..." << endl;
    bdd finals_1 = op1.get_finals();
    bdd finals_2 = op2.get_finals();
    cout << "_finals in op1: " << endl;
    //bdd_print_set(cout, op1.get_dict(), finals_1);
    cout << endl;
    cout << "_finals in op2: " << endl;
    //bdd_print_set(cout, op1.get_dict(), finals_2);
    cout << endl;
    result->_finals = func(finals_1, finals_2);
    cout << "finals in product: " << endl;
    //bdd_print_set(cout, op1.get_dict(), result->_finals);
    cout << endl;
    result->_reach = bddtrue;
    //result->_label_cube = op1._label_cube & op2._label_cube;
    result->_curr_cube = result->_state_vars.get_cube(0);
    cout << "state vars_0 in product: " << endl;
    bdd_print_set(cout, op1.get_dict(), result->_curr_cube);
    cout << endl;
    result->_next_cube = result->_state_vars.get_cube(1);
    cout << "state vars_1 in product: " << endl;
    bdd_print_set(cout, op1.get_dict(), result->_next_cube);
    cout << endl;
    // make pairs
    result->_curr_to_next_pairs = result->_state_vars.make_pair(0, 1);
    result->_next_to_curr_pairs = result->_state_vars.make_pair(1, 0);
    // compute reachable state space
    //cout << "Computing reachable state space in the product..." << endl;
    //result->_reach = result->explore();//bddtrue;//
    //cout << "reachable : " << result->_reach << endl;
    //cout << "Finished computing reachable state space in the product..." << endl;
    //cout << "reachable finals: " << (result->_reach & result->_finals) << endl;
    //cout << "reachable states in product: " << endl;
    //bdd_print_set(cout, result._state_vars.get_dict(), result._reach);
    //cout << endl;
    // needs to whether this is useful
    // only for test

    //bdd all = bdd_replace(result->_reach, result->_curr_to_next_pairs);
    //all = all & result->_reach;
    //cout << "reachable states two copies: " << all << endl;
    //result->_trans = result->_trans & all;

    cout << "Finished computing the product..." << endl;
    return *result;
}

dfwa_ptr
product_dfwa_and(dfwa_ptr op1, dfwa_ptr op2)
{
	return product_dfwa(op1, op2, local_bdd_and);
}

dfwa_ptr
product_dfwa_or(dfwa_ptr op1, dfwa_ptr op2)
{
	return product_dfwa(op1, op2, local_bdd_or);
}

dfwa_ptr
product_dfwa_minus(dfwa_ptr op1, dfwa_ptr op2)
{
	// first make sure op2 is complete
	return product_dfwa(op1, op2, local_bdd_not_and);
}

