#include "dfwanew.hh"
#include "debug.hh"

/*-------------------------------------------------------------------*/
// initialize dd representation of a product DFA
/*-------------------------------------------------------------------*/
dfwa_new::dfwa_new(bdd_dict_ptr dict, bdd& label_cube)
:_state_vars(dict), _label_cube(label_cube)
{
    _aut = nullptr;
    _curr_to_bdd_pairs = nullptr;
}
/*-------------------------------------------------------------------*/
// initialize dd representation of a DFA
/*-------------------------------------------------------------------*/
dfwa_new::dfwa_new(twa_graph_ptr aut, bdd& label_cube, set<unsigned>& finals, const char* name)
:_aut(aut), _state_vars(aut->get_dict(), aut, 1, name, 0, aut->num_states()), _label_cube(label_cube)
{
    _dict = aut->get_dict();
    // label cube
    //cout << "construct state bdd #S = " << aut->num_states() << endl;
    vector<bdd> curr_st_bdd;
    _finals = bddfalse;
    vector<bdd> pre_st_bdd;
    for(unsigned s = 0; s < aut->num_states(); s ++)
    {
        pre_st_bdd.push_back(bddfalse);
    }
    
    for(unsigned s = 0; s < aut->num_states(); s ++)
    {
    	int state = s + 1;
    	// state index from 1 to num_states
        bdd dd = _state_vars.new_value(0, state);
        //#ifdef DEBUG
        cout << "state i = " << s << endl;
        bdd_print_sat(cout, aut->get_dict(), dd);
        cout << endl;
        //#endif
        curr_st_bdd.push_back(dd);
        // final states
        if(finals.find(s) != finals.end()) 
        {
			_finals = _finals | curr_st_bdd[s];
			cout << "final: " << s << endl;
        }
        // need to store bits
        bdd trd = bddfalse;
        for(auto& tr : aut->out(s)) 
        {
            cout << "src = " << s << " dst = " << tr.dst << endl;
            bdd tr_dd = curr_st_bdd[s] & tr.cond;
            pre_st_bdd[tr.dst] = pre_st_bdd[tr.dst] | tr_dd;
        }
    }
    cout << " curr_size = " << curr_st_bdd.size() << endl;
    _init = curr_st_bdd[aut->get_init_state_number()];
    #ifdef DEBUG
    //cout << "init  = " << aut->get_init_state_number() << endl;
    bdd_print_sat(cout, aut->get_dict(), _init);
    cout << endl;
    #endif
    
    
    unsigned num_bits = _state_vars.get_var_num(0);
    cout << "Computing the transition relation for each bit... " << num_bits << endl;
    for(unsigned i = 0; i < num_bits; i ++)
    {
        _bitseq.push_back(bddfalse);
    }
    // needs to be improved by huffman ?
    //bdd all = bddfalse;
    for(unsigned s = 0; s < aut->num_states(); s ++)
    {
    	int state = s + 1;
        vector<bdd> nxt_bits;
        _state_vars.new_value(nxt_bits, 0, state);
        // list of bits encoding
        for(unsigned i = 0; i < nxt_bits.size(); i ++)
        {
            if(nxt_bits[i] == bddtrue)
            {
                // bits is true
                _bitseq[i] = _bitseq[i] | pre_st_bdd[s];
            }
            #ifdef DEBUG
            cout << " bit " << i << " truth = " << (nxt_bits[i] == bddtrue)<< endl;
            bdd_print_set(cout, _state_vars.get_dict(), _bitseq[i]);
            cout << endl;
            #endif
        }
        //all = all | pre_st_bdd[s];
    }
    cout << "Finished computing the transition relation for each bit..." << endl;
    // needs to add another bit for the rest of transitions
    /*
    bdd rest = !all;
    if(rest != bddfalse)
    {
        // needs to add another bit to store the rest of transitions
        unsigned var_num = _dict->register_anonymous_variables(1, _aut);
        // now add it to vars
        bdd dd = bdd_ithvar(var_num);
        _state_vars._dd_vars[0].push_back(dd);
        
        _init = _init & !dd;
        _finals = _finals & !dd;
        for(unsigned i = 0; i < _bitseq.size(); i ++)
        {
            _bitseq[i] = _bitseq[i] & !dd;
        }
        
        _bitseq.push_back(rest | dd);
    }
    */
    _curr_cube = _state_vars.get_cube(0);
    // make pairs
    _curr_to_bdd_pairs = bdd_newpair();
    const vector<bdd>& list_vars = _state_vars.get_bdd_vars(0);
    for(unsigned i = 0; i < list_vars.size(); i ++)
    {
        #ifdef DEBUG
        cout << "var = " << endl;
        bdd_print_set(cout, _state_vars.get_dict(), list_vars[i]);
        cout << endl;
        cout << "rep = " << endl;
        bdd_print_set(cout, _state_vars.get_dict(), _bitseq[i]);
        cout << endl;
        #endif
        bdd_setbddpair(_curr_to_bdd_pairs, bdd_var(list_vars[i]), _bitseq[i]);
        cout << "The number of nodes in B_{" << i << "} is " << bdd_nodecount(_bitseq[i]) << endl;
    }
    cout << "#list_vars = " << list_vars.size() << endl;
}

dfwa_new::~dfwa_new()
{
    delete _curr_to_bdd_pairs;
}

/*-------------------------------------------------------------------*/
// compute previous image of curr 
// FIXED (note the returned image contains propositions)
/*-------------------------------------------------------------------*/
bdd
dfwa_new::pre_image(bdd& curr)
{
    bdd pre = bdd_veccompose(curr, _curr_to_bdd_pairs);
    pre = bdd_exist(pre, _label_cube);
    return pre;
}

bdd
dfwa_new::back_explore()
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
        cout << "The number of node in R(" << count << ") is " << bdd_nodecount(s) << endl;
        ++ count;
    }
    return s;
}

void
dfwa_new::output(ostream& os)
{
    os << "dfwan: " << endl;
    os << "init: " << endl;
    bdd_print_set(os, _state_vars.get_dict(), _init);
    os << endl;
    
    os << "trans: " << endl;
    for(unsigned i = 0; i < _bitseq.size(); i ++)
    {
        os << "bit " << i << " :" << endl;
        bdd_print_set(os, _state_vars.get_dict(), _bitseq[i]);
        os << endl;
        
        os << "not bit " << i << " :" << endl;
        bdd_print_set(os, _state_vars.get_dict(), !_bitseq[i]);
        os << endl;
    }
    
    os << "finals: " << endl;
    bdd_print_set(os, _state_vars.get_dict(), _finals);
    os << endl;
    
    cout << "init predecessors: " << endl;
    bdd temp = pre_image(_init);
    bdd_print_set(os, _state_vars.get_dict(), temp);
    os << endl;
    
    cout << "finals predecessors: " << endl;
    temp = pre_image(_finals);
    bdd_print_set(os, _state_vars.get_dict(), temp);
    os << endl;
    
}
/*-------------------------------------------------------------------*/
// intersection dfwan:
// result is an empty dfwan
/*-------------------------------------------------------------------*/
void 
intersect_dfwan(dfwa_new_ref result, dfwa_new_ref op1, dfwa_new_ref op2)
{
    // now we add variables from op1 and op2
    cout << "Computing the intersection product..." << endl;
    // check whether this part can be improved
    result._state_vars.add_bdd_vars(op1._state_vars);
    result._state_vars.add_bdd_vars(op2._state_vars);
    // 
    result._init = op1.get_init() & op2.get_init();
    // especially for the transition relation
    cout << "Computing transition relation in the intersection product..." << endl;
    for(unsigned i = 0; i < op1._bitseq.size(); i ++)
    {
        result._bitseq.push_back(op1._bitseq[i]);
    }
    for(unsigned i = 0; i < op2._bitseq.size(); i ++)
    {
        result._bitseq.push_back(op2._bitseq[i]);
    }
    cout << "Finished computing transition relation in the intersection product..." << endl;
    result._finals = op1.get_finals() & op2.get_finals();
    //result._reach = bddfalse;
    result._curr_cube = result._state_vars.get_cube(0);
    // make pairs
    // compute reachable state space
    //result._curr_to_bdd_pairs = bdd_newpair();
    vector<bdd> list_vars = result._state_vars.get_bdd_vars(0);
    for(unsigned i = 0; i < list_vars.size(); i ++)
    {
        #ifdef DEBUG
        cout << "var = " << endl;
        bdd_print_set(cout, result._state_vars.get_dict(), list_vars[i]);
        cout << endl;
        cout << "rep = " << endl;
        bdd_print_set(cout, result._state_vars.get_dict(), result._bitseq[i]);
        cout << endl;
        #endif
       // bdd_setbddpair(result._curr_to_bdd_pairs, bdd_var(list_vars[i]), result._bitseq[i]);
    }
    cout << "Finished computing the intersection product..." << endl;
}

dfwa_new_ref
product_dfwa_new(dfwa_new_ref op1, dfwa_new_ref op2, function<bdd(bdd&, bdd&)> func)
{
    // now we add variables from op1 and op2
    cout << "Computing the intersection product..." << endl;
    bdd label_cube = op1._label_cube & op2._label_cube;
    dfwa_new_ptr result = new dfwa_new(op1.get_dict(), label_cube);
    // check whether this part can be improved
    result->_state_vars._copies = 1;
    result->_state_vars.add_bdd_vars(op1._state_vars);
    result->_state_vars.add_bdd_vars(op2._state_vars);
    //
    result->_init = op1.get_init() & op2.get_init();
    // especially for the transition relation
    cout << "Computing transition relation in the intersection product..." << endl;
    // probably this is not needed
    for(unsigned i = 0; i < op1._bitseq.size(); i ++)
    {
        result->_bitseq.push_back(op1._bitseq[i]);
    }
    for(unsigned i = 0; i < op2._bitseq.size(); i ++)
    {
        result->_bitseq.push_back(op2._bitseq[i]);
    }
    cout << "Finished computing transition relation in the intersection product..." << endl;
    bdd finals_1 = op1.get_finals();
    bdd finals_2 = op2.get_finals();
    result->_finals = func(finals_1, finals_2);

    result->_curr_cube = result->_state_vars.get_cube(0);
    // make pairs
    // compute reachable state space
    result->_curr_to_bdd_pairs = bdd_mergepairs(op1._curr_to_bdd_pairs, op2._curr_to_bdd_pairs);
    cout << "Finished computing the intersection product..." << endl;

    return *result;
}

dfwa_new_ref
product_dfwa_new_and(dfwa_new_ref op1, dfwa_new_ref op2)
{
	return product_dfwa_new(op1, op2, local_bdd_and);
}
dfwa_new_ref
product_dfwa_new_or(dfwa_new_ref op1, dfwa_new_ref op2)
{
	return product_dfwa_new(op1, op2, local_bdd_or);
}
dfwa_new_ref
product_dfwa_new_minus(dfwa_new_ref op1, dfwa_new_ref op2)
{
	return product_dfwa_new(op1, op2, local_bdd_not_and);
}
