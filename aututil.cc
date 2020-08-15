/*-------------------------------------------------------------------*/
// Yong Li (liyong@ios.ac.cn)
/*-------------------------------------------------------------------*/

#include "aututil.hh"

// quantify out all variables in @cube on transitions
// bdd_vars should include alive_ap
twa_graph_ptr
project(spot::bdd_dict_ptr dict, twa_graph_ptr aut, bdd cube, bool min)
{
    // now construct the dfa
    spot::twa_graph_ptr ret = make_twa_graph(aut->get_dict());

    ret->set_buchi();
    ret->prop_state_acc();
    // add one extra state for accepting state

    unsigned num_states = aut->num_states();
    ret->new_states(num_states);
    ret->set_init_state(aut->get_init_state_number());
    // for identifying accepting states
    string alive_ap(ALIVE_AP);
    int var_index = dict->varnum(formula::ap(alive_ap));
    bdd p = bdd_ithvar(var_index);
    // register propositions
    //cout << "var map: " << endl;
    map<formula, int> &var_map = ret->get_dict()->var_map;
    map<formula, int>::const_iterator iter = var_map.begin();
    while (iter != var_map.end())
    {
        formula key = iter->first;
        int value = iter->second;
        //cout << key << " -> " << value << endl;
        ret->register_ap(key);
        iter++;
    }
    // project on transitions
    for (unsigned s = 0; s < num_states; ++s)
    {
        // check if it is predecessor of accepting
        //out << "State " << s << ":\n";
        bool acc = false;
        // check whether current state s is accepting in DFA
        for (auto &t : aut->out(s))
        {
            if (t.acc.has(0) || (t.cond & !p) != bddfalse)
            {
                acc = true;
                break;
            }
        }
        if (acc)
        {
            // make accepting state s a sink accepting state
            ret->new_edge(s, s, bddtrue, {0});
            //cout << "state: " << s << endl;
        }
        else
        {
            for (auto &t : aut->out(s))
            {
                // quantify out output variables; transition t: t.cond is the bdd
                bdd new_label = bdd_exist(t.cond, cube);
                //cout << "projected label : " << new_label << " old label: " << t.cond << " cube: " << cube << endl;
                ret->new_edge(t.src, t.dst, new_label);
            }
        }
        // there may be a transition to the sink accepting state
    }
    // for sure original sink accepting state is no longer reachable
    // hopefully reduce some states
    clock_t post_start = clock();
    spot::postprocessor post;
    post.set_type(spot::postprocessor::BA);
    //post.set_pref(spot::postprocessor::Deterministic);
    post.set_level(spot::postprocessor::Low);
    ret = post.run(ret);
    clock_t post_end = clock();
    cout << "Total CPU time used for reducing projected DBA: "
         << 1000.0 * (post_end - post_start) / CLOCKS_PER_SEC << " ms\n";
    if (min)
    {
        clock_t c_start = clock();
        twa_graph_ptr tmp = spot::minimize_wdba(ret);
        clock_t c_end = clock();
        cout << "Total CPU time used for minimizing projected DBA: "
             << 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC << " ms\n";
        //ofstream outfile1("output3.hoa");
        //print_hoa(outfile1, C);
        return tmp;
    }
    else
    {
        return ret;
    }
}

// reverse a DFA： make accepting states initial states and reverse transitions
// bdd_vars should include alive_ap
twa_graph_ptr
reverse(spot::bdd_dict_ptr dict, twa_graph_ptr aut)
{
    clock_t rev_start = clock();
    // now construct the NFA
    spot::twa_graph_ptr ret = make_twa_graph(aut->get_dict());

    ret->set_buchi();
    ret->prop_state_acc();
    // add one extra state for accepting state

    // set the maximal state the unique initial state
    unsigned num_states = aut->num_states() + 1;
    
    ret->new_states(num_states);
    int init_state = aut->num_states();
    ret->set_init_state(init_state);
    // for identifying accepting states
    string alive_ap(ALIVE_AP);
    int var_index = dict->varnum(formula::ap(alive_ap));
    bdd p = bdd_ithvar(var_index);
    // register propositions
    //cout << "var map: " << endl;
    map<formula, int> &var_map = ret->get_dict()->var_map;
    map<formula, int>::const_iterator iter = var_map.begin();
    while (iter != var_map.end())
    {
        formula key = iter->first;
        int value = iter->second;
        //cout << key << " -> " << value << endl;
        ret->register_ap(key);
        iter++;
    }
    // reverse the transitions of the DFA
    set<unsigned> acc_states;
    unsigned sink = -1;
    for (unsigned s = 0; s < num_states - 1; ++s)
    {
        // check if it is predecessor of accepting
        //out << "State " << s << ":\n";
        bool acc = false;
        // check whether current state s is the sink accepting state in wDBA
        bool found_sink = false;
        if (!found_sink)
        {
            for (auto &t: aut->out(s))
            {
                // loop transition and accepting transition
                if (t.acc.has(0) && s == t.dst)
                {
                    sink = s;
                    found_sink = true;
                    break;
                }
            }
        }
        // if not sink state in wDBA, check whether s is accepting in DFA
        if(s != sink)
        {
            for (auto &t : aut->out(s))
            {
                if (t.acc.has(0) || (t.cond & !p) != bddfalse)
                {
                    acc = true;
                    break;
                }
            }
        }
        
        if (acc)
        {
            // this state may be the sink or accepting states in the DFA
            acc_states.insert(s);
            //cout << "state: " << s << endl;
        }
        else
        {
            // not accepting states, then just reverse transitions
            // no transitions to the sink state
            for (auto &t : aut->out(s))
            {
                // transition t: t.cond is the bdd
                ret->new_edge(t.dst, t.src, t.cond);
            }
        }
    }
    assert( sink != -1);
    // initial state to sink, as accepting state in NFA
    ret->new_edge(aut->get_init_state_number(), sink, !p);
    ret->new_acc_edge(sink, sink, !p, true);
    // reverse transitions for accepting states in reverse DFA
    for(unsigned s : acc_states)
    {
        for (auto &t : aut->out(s))
        {
            if(t.dst != sink)
            {
                // (s, a, t) -> (t, a, s)
                ret->new_edge(t.dst, t.src, t.cond); 
            }
        }
    }
    // now we need to add transitions from the accepting states for the initial state
    // copy existing transitions
    for(unsigned s : acc_states)
    {
        // note that here automaton is ret
        for (auto &t : ret->out(s))
        {
            if(t.dst != sink)
            {
                // (s, a, t) -> (t, a, s)
                ret->new_edge(init_state, t.dst, t.cond); 
            }
        }
    }
    // merge transitions with same source and desition states
    ret->merge_edges();
    cout << "Resulting automaton is " << (is_deterministic(ret)? "deterministic" : "nondeterministic") << endl;
    clock_t rev_end = clock();
    cout << "Total CPU time used for reversing DFA: "
         << 1000.0 * (rev_end - rev_start) / CLOCKS_PER_SEC << " ms\n";

    return ret;
}

