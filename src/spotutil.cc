/*-------------------------------------------------------------------*/
// Yong Li (liyong@ios.ac.cn)
/*-------------------------------------------------------------------*/

#include "spotutil.hh"

#include "debug.hh"
#include "mona.hh"


bool
is_even(unsigned value)
{
	return (value & 1) == 0;
}

/*-------------------------------------------------------------------*/
// get the list of var indices from a cube
/*-------------------------------------------------------------------*/
void
get_list_var_indices(vector<int>& vars, bdd cube)
{
    while(cube != bddtrue)
    {
        int index = bdd_var(cube);
        vars.push_back(index);
        cube = bdd_high(cube);
    }
}

/*-------------------------------------------------------------------*/
// check the language inclusion of two input NBAs
/*-------------------------------------------------------------------*/
string
is_twa_included(twa_graph_ptr aut_a, twa_graph_ptr aut_b) 
{
	auto complement_aut_b = dualize(aut_b);
	spot::twa_word_ptr word = aut_a->intersecting_word(complement_aut_b);
	std::stringstream ss;
	if(word != nullptr) {
		//spot::twa_graph_ptr result = word->as_automaton();
		word->simplify();
		ss << (*word);
		 // Print the resulting automaton.
		//print_hoa(std::cout, result);
    }else {
		ss << "";
	}
	return ss.str();
}

/*-------------------------------------------------------------------*/
// check the equivalence of two input NBAs
/*-------------------------------------------------------------------*/
string
is_twa_equivalent(twa_graph_ptr aut_a, twa_graph_ptr aut_b)
{
    string word = is_twa_included(aut_a, aut_b);
    if(word.size() != 0)
    {
        return word;
    }
    word = is_twa_included(aut_b, aut_a);
    return word;
}

/*-------------------------------------------------------------------*/
// compute the connectives and operators of a formula
/*-------------------------------------------------------------------*/
unsigned
get_size_formula(formula& f)
{
    // use id as the key to record the length and the number of props in a formula
    //std::cout << "operator" << f.id() << std::endl;
    unsigned count = 0;
    if(f.kind() == op::Not || f.kind() == op::X || f.kind() == op::strong_X || f.kind() == op::F || f.kind() == op::G )
    {
        formula f1 = f[0];
        count += 1 + get_size_formula(f1);
    }else
    if(f.kind() == op::U || f.kind() == op::R || f.kind() == op::M || f.kind() == op::W  || f.kind() == op::Implies || f.kind() == op::Equiv )
    {
        formula f1 = f[0], f2 = f[1];
        count += 1 + get_size_formula(f1) + get_size_formula(f2);
    }else
    if(f.kind() == op::And || f.kind() == op::Or)
    {
        count = 1;
        for(formula child: f)
        {
            count += get_size_formula(child);
        }
    }
    return count;
}

/*-------------------------------------------------------------------*/
// compute the number of propositions of a formula
/*-------------------------------------------------------------------*/
void
get_formula_aps(formula& f, set<formula>& aps)
{
    if(f.kind() == op::Not || f.kind() == op::X || f.kind() == op::strong_X || f.kind() == op::F || f.kind() == op::G )
    {
        formula f1 = f[0];
        return get_formula_aps(f1, aps);
    }else
    if(f.kind() == op::U || f.kind() == op::R || f.kind() == op::M || f.kind() == op::W  || f.kind() == op::Implies || f.kind() == op::Equiv )
    {
        formula f1 = f[0], f2 = f[1];
        get_formula_aps(f1, aps);
        get_formula_aps(f2, aps);
    }else
    if(f.kind() == op::And || f.kind() == op::Or)
    {
        for(formula child: f)
        {
            get_formula_aps(child, aps);
        }
    }else 
    if(f.kind() == op::ap)
    {
        aps.insert(f);
    }else if(f.kind() != op::tt && f.kind() != op::ff)
    {
        cerr << "Error formula in get_formula_aps(): " << f << endl;
        exit(-1);
    }
}


unsigned
get_size_formula_ap(formula& f)
{
    set<formula> aps;
    get_formula_aps(f, aps);
    #ifdef DEBUG
    for(auto& ap : aps)
    {
        cout << "AP = " << ap << endl;
    }
    #endif
    return aps.size();
}

/*-------------------------------------------------------------------*/
// check the equivalence of two input NBAs
/*-------------------------------------------------------------------*/
twa_graph_ptr
trans_formula(formula f, bdd_dict_ptr dict, unsigned num_ap_for_mona)
{
    twa_graph_ptr aut;
    unsigned num_aps = get_size_formula_ap(f);

    if(num_aps > num_ap_for_mona)
    {
        DEBUG_STDOUT( "mona translating ");
        aut = translate_ltlf_mona(f, dict);
        DEBUG_STDOUT( "done mona translating ");
    }else
    {
        DEBUG_STDOUT( "spot translating ");
        
        //cout << "spot translating" << endl;
        option_map m;
        translator trans(dict, &m);
        formula ltl = from_ltlf(f, ALIVE_AP);
        
        //spot::tl_simplifier simpl;
        //ltl = simpl.simplify(ltl);
        aut = trans.run(ltl);
        //spot::postprocessor post;
        //post.set_type(spot::postprocessor::BA);
        //post.set_pref(spot::postprocessor::Deterministic);
        //post.set_pref(spot::postprocessor::Deterministic); // or ::Deterministi
        //print_hoa(std::cout, aut) << '\n';
        DEBUG_STDOUT( "done spot translating ");
    }

    return aut;
}
/*-------------------------------------------------------------------*/
// compute the final states of a DFA
/*-------------------------------------------------------------------*/
void
get_final_states(twa_graph_ptr aut, set<unsigned>& finals)
{
  stack<unsigned> todo;
  unsigned init = aut->get_init_state_number();
  todo.push(init);
  set<unsigned> seen;
  seen.insert(init);
  while (!todo.empty())
    {
      unsigned s = todo.top();
      todo.pop();
      // transitions of s 
      for (auto& e: aut->out(s))
      {
        if (seen.insert(e.dst).second)
          todo.push(e.dst);
        if(e.acc.has(0)) 
        {
            //cout << "final state: " << e.dst << endl;
            finals.insert(e.dst);
        }
      }
    }
}

void
get_nonfinal_states(twa_graph_ptr A, set<unsigned>& nonfinals)
{
    for(unsigned s = 0; s < A->num_states(); s ++)
	{
		for (auto& e: A->out(s))
		{
			if(! e.acc.has(0))
			{
				nonfinals.insert(s);
			}
		}
	}
}

/*-------------------------------------------------------------------*/
// obtain the set of final states of a DFA 
// (has transition to accepting state of weak DBA on !alive )
/*-------------------------------------------------------------------*/
void
compute_final_states(twa_graph_ptr A, set<unsigned>& finals)
{
	bdd p2 = bdd_ithvar(A->register_ap(ALIVE_AP));
	for(unsigned s = 0; s < A->num_states(); s ++)
	{
		for (auto& e: A->out(s))
		{
			if((e.cond & !p2) != bddfalse && s != e.dst)
			{
				finals.insert(s);
			}
		}
	}
}

/*-------------------------------------------------------------------*/
// obtain the set of accepting states of a weak DBA
/*-------------------------------------------------------------------*/
void
compute_accepting_states(twa_graph_ptr A, set<unsigned>& acc)
{
    for(unsigned s = 0; s < A->num_states(); s ++)
	{
		for (auto& e: A->out(s))
		{
			if(e.acc.has(0))
			{
				acc.insert(s);
			}
		}
	}
}
/*-------------------------------------------------------------------*/
// convert a DFA to a weak DBA by adding alive propositions
/*-------------------------------------------------------------------*/
twa_graph_ptr
dfa_to_wdba(twa_graph_ptr aut, bdd_dict_ptr dict, set<unsigned> finals)
{
    return nullptr;
}

void
output_bdd(ostream& os, bdd dd, set<int>& computed_table)
{
	int dd_id = dd.id();
	if(computed_table.find(dd_id) != computed_table.end())
	{
		return ;
	}
	if(dd == bddtrue)
	{
		os << dd_id << " -1 t 0" << endl;
	}else
	if(dd == bddfalse)
	{
		os << dd_id << " -1 f 0" << endl;
	}else
	{
		bdd high = bdd_high(dd);
		output_bdd(os, high, computed_table);
		bdd low = bdd_low(dd);
		output_bdd(os, low, computed_table);
		os << dd_id << " " << bdd_var(dd) << " " << high.id() << " " << low.id() << endl;
	}
	computed_table.insert(dd_id);
}

void
output_bdd(ostream& os, bdd dd)
{
	set<int> computed_table;
	output_bdd(os, dd, computed_table);
}

bdd
local_bdd_and(bdd& op1, bdd& op2)
{
	//cout << "Inside local_bdd_and..." << endl;
	return op1 & op2;
}

bdd
local_bdd_or(bdd& op1, bdd& op2)
{
	//cout << "Inside local_bdd_or..." << endl;
	return op1 | op2;
}

bdd
local_bdd_not_and(bdd& op1, bdd& op2)
{
	//cout << "Inside local_bdd_not_and..." << endl;
	return op1 & !op2;
}
