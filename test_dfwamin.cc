#include <iostream>
#include <fstream>

#include <spot/twaalgos/hoa.hh>
#include <spot/twaalgos/minimize.hh>
#include <spot/twaalgos/powerset.hh>

#include <spot/twa/bddprint.hh>
#include <spot/twa/twagraph.hh>


#include <spot/tl/formula.hh>
#include <spot/tl/print.hh>
#include <spot/tl/parse.hh>

#include "dfwa.hh"
#include "dfwamin.hh"
#include "dfwavar.hh"
#include "mona.hh"
#include "spotutil.hh"
#include "minimize.hh"


using namespace std;
using namespace spot;

void
get_formulas(vector<formula>& lst, formula f)
{
    //cout << "split formula " << endl;
    if(f.kind() == op::And)
    {
        // needs to limit the number of conjunctions if no minimization is used before producting two FAs
        for(formula child: f)
          {
              lst.push_back(child);
              //cout << "subformula: " << child << endl;
          }
    }else
    {
        lst.push_back(f);
    }
    cout << "split formula " << endl;
}
void
compute_final_states(twa_graph_ptr A, set<unsigned>& finals)
{
	bdd p2 = bdd_ithvar(A->register_ap("alive"));
	for(unsigned s = 0; s < A->num_states(); s ++)
	{
		for (auto& e: A->out(s))
		{
			if((e.cond & !p2) != bddfalse)
			{
				finals.insert(e.dst);
			}
		}
	}
}

void
test_equivalent(twa_graph_ptr A, twa_graph_ptr small)
{
	//cout << "subset construction" << endl;
	//power_map pm;
	//A = tba_determinize(A); //tgba_powerset(A, pm);
	cout << "A automaton: " << endl;
	print_hoa(std::cout, A);
	cout << endl;
	cout << "B automaton: " << endl;
	print_hoa(std::cout, small);
	cout << endl;
	set<formula> set_of_props;
	for(auto &f : A->ap())
	{
	   set_of_props.insert(f);
	}
	for(auto &f : small->ap())
	{
	   set_of_props.insert(f);
	}
	bdd_dict_ptr dict = A->get_dict();

	// there is alive states
	bdd label_cube = bddtrue;
	for(auto& f : set_of_props)
	{
	   cout << "ap: " << f << endl;
	   label_cube = label_cube & bdd_ithvar(dict->varnum(f));
	}
	set<unsigned> finals_A;
	compute_final_states(A, finals_A);
    dfwa*  a_A = new dfwa(A, label_cube, finals_A);
    dfwa_ptr dfa_A = *a_A;
    dfa_A.output(cout);
    set<unsigned> finals_B;
    compute_final_states(small, finals_B);
    dfwa* a_B = new dfwa(small, label_cube, finals_B);
    dfwa_ptr dfa_B = *a_B;
    dfa_B.output(cout);
    {
    	cout << "A <=> B" << endl;
        dfwa_ptr pro = product_dfwa_minus(dfa_A, dfa_B);
        	  //bdd reach = pro.explore();
        if(pro.is_empty())
        {
          cout << "A is subset of B" << endl;
        }else
        {
        	cerr << "ERROR" << endl;
        	exit(-1);
        }
        //pro.free_variables();
        delete &pro;
        pro = product_dfwa_minus(dfa_B, dfa_A);
        if(pro.is_empty())
        {
           cout << "B is subset of A" << endl;
        }else
        {
        	cerr << "ERROR" << endl;
        	exit(-1);
        }
        //pro.free_variables();
        delete &pro;
    }
    {

		// minimize A
		cudd manager; // = new cudd();
		cout << "cudd minimization" << endl;
		cudd_ptr mgr = &manager;
		dfwa_min dfa_m(mgr, dfa_A);
		dfa_m.output(cout);
			  //delete manager;
		dfa_m.minimize();

		dfwa_ptr min = dfa_m.move_dfwa();
		//min.output(cout);
		//delete manager;

		dfwa_ptr pro = product_dfwa_minus(min, dfa_B);
			  //bdd reach = pro.explore();
		if(pro.is_empty())
		{
		  cout << "A is subset of B" << endl;
		}else
		{
			cerr << "ERROR" << endl;
			exit(-1);
		}
		delete &pro;
		pro = product_dfwa_minus(dfa_B, min);
		if(pro.is_empty())
		{
		   cout << "B is subset of A" << endl;
		}else
		{
			cerr << "ERROR" << endl;
			exit(-1);
		}
		delete &pro;

		delete &min;
    }
    delete a_A;
    delete a_B;

}
void
test(bdd_dict_ptr dict);

/**
 * the use of anonymous variables
 * */
void
test_anon(bdd_dict_ptr dict)
{
	int* first = new int[5];
	unsigned var_index = dict->register_anonymous_variables(5, first);
	bdd first_cube = bddtrue;
	for(unsigned i = 0; i < 5; i ++)
	{
		first_cube &= bdd_ithvar(var_index + i);
		first[i] = var_index + i;
		cout << "first: " << var_index + i << endl;
	}
	int* second = new int[3];
	var_index = dict->register_anonymous_variables(3, second);
	bdd second_cube = bddtrue;
	for(unsigned i = 0; i < 3; i ++)
	{
		second_cube &= bdd_ithvar(var_index + i);
		second[i] = var_index + i;
		cout << "second: " << var_index + i << endl;
	}

	cout << "first: " << endl;
	bdd_print_set(cout, dict, first_cube);
	cout << endl;

	cout << "second: " << endl;
	bdd_print_set(cout, dict, second_cube);
	cout << endl;
	dict->unregister_all_my_variables(first);
	dict->unregister_all_my_variables(second);
	delete first;
	delete second;
}

void test_2();
void
test3(bdd_dict_ptr dict);

int main(int argc, char** argv)
{
    if(argc < 1 || argc > 3)
		std::cout << "please input formula file" << std::endl;

	ifstream ltlfile(argv[1]);
	string line;
	bdd_dict_ptr dict = make_bdd_dict();

	// ltlfile.open(argv[1]);

	cout << "Opened file" << endl;

	if (ltlfile.is_open())
	{
	  cout << "We are inside the file" << endl;
	  int i = 0;
	  while (getline(ltlfile, line))
	  {
		cout << "formula: " << line << endl;

		auto pf1 = spot::parse_infix_psl(line.c_str());
		if (pf1.format_errors(std::cerr))
		{
			std::cerr << "error: " << line << std::endl;
			return -1;
		}
		cout << "spot translating..." << endl;
		twa_graph_ptr A = trans_formula(pf1.f, dict);
		option_map opt;
		// needs option map when spot version is less than 2.7.5
		spot::postprocessor post(&opt);
		post.set_type(spot::postprocessor::BA);
		post.set_pref(spot::postprocessor::Deterministic);
		//post.set_pref(spot::postprocessor::Deterministic); // or ::Deterministi
		A = post.run(A);
		cout << "mona translating..." << endl;
		twa_graph_ptr B = translate_ltlf_mona(pf1.f, dict);
		test_equivalent(A, B);
	  }
	}

  test(dict);
  test_anon(dict);
  test_2();
  test3(dict);
  return 0;
}

void test_2()
{
	 // The bdd_dict is used to maintain the correspondence between the
	  // atomic propositions and the BDD variables that label the edges of
	  // the automaton.
	  spot::bdd_dict_ptr dict = spot::make_bdd_dict();
	  // This creates an empty automaton that we have yet to fill.
	  spot::twa_graph_ptr aut = make_twa_graph(dict);

	  // Since a BDD is associated to every atomic proposition, the
	  // register_ap() function returns a BDD variable number
	  // that can be converted into a BDD using bdd_ithvar().
	  bdd p1 = bdd_ithvar(aut->register_ap("p1"));
	  bdd p2 = bdd_ithvar(aut->register_ap("alive"));

	  set<formula> set_of_props;
	  for(auto &f : aut->ap())
	  {
	    set_of_props.insert(f);
	  }
	  // there is alive states
	  bdd label_cube = bddtrue;
	  for(auto& f : set_of_props)
	  {
	    label_cube = label_cube & bdd_ithvar(dict->varnum(f));
	  }

	  // Set the acceptance condition of the automaton to Inf(0)&Inf(1)
	  aut->set_buchi();
	  //aut->prop_state_acc(true);
	  // States are numbered from 0.
	  aut->new_states(6);
	  // The default initial state is 0, but it is always better to
	  // specify it explicitely.
	  aut->set_init_state(0U);

	  // new_edge() takes 3 mandatory parameters: source state,
	  // destination state, and label.  A last optional parameter can be
	  // used to specify membership to acceptance sets.
	  //a 0, b 1, c 2, d 3, e 4, f 5,


	  aut->new_edge(0, 1, !p1);
	  aut->new_acc_edge(0, 3, p1, true);

	  aut->new_edge(1, 0, !p1);
	  aut->new_acc_edge(1, 2, p1, true);

	  aut->new_edge(2, 5, p1);
	  aut->new_acc_edge(2, 4, !p1, true);

	  aut->new_edge(3, 5, p1);
	  aut->new_acc_edge(3, 4, !p1, true);

	  aut->new_edge(4, 5, p1);
	  aut->new_acc_edge(4, 4, !p1, true);

	  aut->new_edge(5, 5, bddtrue);

	   /*
	  aut->new_edge(0, 1, p2 & !p1);
	  aut->new_edge(0, 3, p2 & p1);

	  aut->new_edge(1, 0, p2& !p1);
	  aut->new_edge(1, 2, p2&p1);

	  aut->new_edge(2, 5, p2& p1);
	  aut->new_edge(2, 4, p2 & !p1);

	  aut->new_edge(3, 5, p2&p1);
	  aut->new_edge(3, 4, p2&!p1);

	  aut->new_edge(4, 5, p2&p1);
	  aut->new_edge(4, 4, p2&!p1);

	  aut->new_edge(5, 5, p2);

	  aut->new_edge(2, 6, !p2);
	  aut->new_edge(3, 6, !p2);
	  aut->new_edge(4, 6, !p2);
	  aut->new_edge(6, 6, !p2, {0});
	  */
	  // Print the resulting automaton.
	  print_hoa(std::cout, aut);

	  //
	  cout << endl << "minimized "<< endl;
	  //aut = minimize_dfa(aut);
	  twa_graph_ptr autB = minimize_wdba(aut);
	  print_hoa(std::cout, autB);

	  //dfwa_var vars(dict, aut, 2, "s", 0, 6);

	  //cout << vars.get_lower() << endl;
	  //cout << vars.get_upper() << endl;

	  set<unsigned> finals;
	  finals.insert(0);
	  //finals.insert(3);
	  //finals.insert(4);
	  dfwa dfab(autB, label_cube, finals);
	  set<unsigned> finalsa;
	  finalsa.insert(2);
	  finalsa.insert(3);
	  finalsa.insert(4);
	  dfwa dfaa(aut, label_cube, finalsa);

	  dfaa.output(cout);

	  dfab.output(cout);

	  dfwa_ptr result = product_dfwa_and(dfab, dfaa);

	  result.output(cout);

	  // minimization
	  cudd manager; // = new cudd();
	  cout << "cudd" << endl;
	  cudd_ptr mgr = &manager;
	  dfwa_min df(mgr, result);
	  df.output(cout);
	  	  //delete manager;
	  df.minimize();

	  dfwa_ptr min =  df.move_dfwa();
	  min.output(cout);

	  delete &result;
	  delete &min;

}

void
test(bdd_dict_ptr dict)
{
	  // The bdd_dict is used to maintain the correspondence between the
	  // atomic propositions and the BDD variables that label the edges of
	  // the automaton.
	  //spot::bdd_dict_ptr dict = spot::make_bdd_dict();
	  // This creates an empty automaton that we have yet to fill.
	  spot::twa_graph_ptr aut = make_twa_graph(dict);

	  // Since a BDD is associated to every atomic proposition, the
	  // register_ap() function returns a BDD variable number
	  // that can be converted into a BDD using bdd_ithvar().
	  bdd p1 = bdd_ithvar(aut->register_ap("p1"));
	  //bdd p2 = bdd_ithvar(aut->register_ap("alive"));

	  set<formula> set_of_props;
	  for(auto &f : aut->ap())
	  {
	    set_of_props.insert(f);
	  }
	  // there is alive states
	  bdd label_cube = bddtrue;
	  for(auto& f : set_of_props)
	  {
	    label_cube = label_cube & bdd_ithvar(dict->varnum(f));
	  }

	  // Set the acceptance condition of the automaton to Inf(0)&Inf(1)
	  aut->set_buchi();
	  //aut->prop_state_acc(true);
	  // States are numbered from 0.
	  aut->new_states(6);
	  // The default initial state is 0, but it is always better to
	  // specify it explicitely.
	  aut->set_init_state(0U);

	  // new_edge() takes 3 mandatory parameters: source state,
	  // destination state, and label.  A last optional parameter can be
	  // used to specify membership to acceptance sets.
	  //a 0, b 1, c 2, d 3, e 4, f 5,


	  aut->new_edge(0, 1, !p1);
	  aut->new_acc_edge(0, 3, p1, true);

	  aut->new_edge(1, 0, !p1);
	  aut->new_acc_edge(1, 2, p1, true);

	  //aut->new_edge(2, 5, p1);
	  aut->new_acc_edge(2, 4, !p1, true);

	  //aut->new_edge(3, 5, p1);
	  aut->new_acc_edge(3, 4, !p1, true);

	  //aut->new_edge(4, 5, p1);
	  aut->new_acc_edge(4, 4, !p1, true);

	  //aut->new_edge(5, 5, bddtrue);

	   /*
	  aut->new_edge(0, 1, p2 & !p1);
	  aut->new_edge(0, 3, p2 & p1);

	  aut->new_edge(1, 0, p2& !p1);
	  aut->new_edge(1, 2, p2&p1);

	  aut->new_edge(2, 5, p2& p1);
	  aut->new_edge(2, 4, p2 & !p1);

	  aut->new_edge(3, 5, p2&p1);
	  aut->new_edge(3, 4, p2&!p1);

	  aut->new_edge(4, 5, p2&p1);
	  aut->new_edge(4, 4, p2&!p1);

	  aut->new_edge(5, 5, p2);

	  aut->new_edge(2, 6, !p2);
	  aut->new_edge(3, 6, !p2);
	  aut->new_edge(4, 6, !p2);
	  aut->new_edge(6, 6, !p2, {0});
	  */
	  // Print the resulting automaton.
	  print_hoa(std::cout, aut);

	  //
	  cout << endl << "minimized "<< endl;
	  //aut = minimize_dfa(aut);
	  twa_graph_ptr autB = minimize_wdba(aut);
	  print_hoa(std::cout, autB);

	  //dfwa_var vars(dict, aut, 2, "s", 0, 6);

	  //cout << vars.get_lower() << endl;
	  //cout << vars.get_upper() << endl;

	  set<unsigned> finals;
	  finals.insert(0);
	  //finals.insert(3);
	  //finals.insert(4);
	  dfwa dfab(autB, label_cube, finals);
	  set<unsigned> finalsa;
	  finalsa.insert(2);
	  finalsa.insert(3);
	  finalsa.insert(4);
	  dfwa dfaa(aut, label_cube, finalsa);

	  dfaa.output(cout);

	  dfab.output(cout);

	  dfwa result(dict, label_cube);
	  intersect_dfwa(result, dfab, dfaa);

	  result.output(cout);

	  cudd manager; // = new cudd();
	  cout << "cudd" << endl;
	  cudd_ptr mgr = &manager;
	  dfwa_min df(mgr, dfaa);
	  df.output(cout);
	  //delete manager;
	  df.minimize();

	  dfwa_ptr min =  df.move_dfwa();
	  min.output(cout);
	  //delete manager;

	  dfwa_ptr pro = product_dfwa_minus(min, dfab);
	  //bdd reach = pro.explore();
	  if(pro.is_empty())
	  {
		  cout << "A is subset of B" << endl;
	  }
	  delete &pro;
	  pro = product_dfwa_minus(dfab, min);
	  if(pro.is_empty())
	   {
	  	  cout << "B is subset of A" << endl;
	   }
	  delete &pro;
	  delete &min;
}

void
test3(bdd_dict_ptr dict)
{
	  // The bdd_dict is used to maintain the correspondence between the
	  // atomic propositions and the BDD variables that label the edges of
	  // the automaton.
	  //spot::bdd_dict_ptr dict = spot::make_bdd_dict();
	  // This creates an empty automaton that we have yet to fill.
	  spot::twa_graph_ptr aut = make_twa_graph(dict);

	  // Since a BDD is associated to every atomic proposition, the
	  // register_ap() function returns a BDD variable number
	  // that can be converted into a BDD using bdd_ithvar().
	  bdd p1 = bdd_ithvar(aut->register_ap("p1"));

	  set<formula> set_of_props;
	  for(auto &f : aut->ap())
	  {
	    set_of_props.insert(f);
	  }
	  // there is alive states
	  bdd label_cube = bddtrue;
	  for(auto& f : set_of_props)
	  {
	    label_cube = label_cube & bdd_ithvar(dict->varnum(f));
	  }

	  // Set the acceptance condition of the automaton to Inf(0)&Inf(1)
	  aut->set_buchi();
	  //aut->prop_state_acc(true);
	  // States are numbered from 0.
	  aut->new_states(3);
	  // The default initial state is 0, but it is always better to
	  // specify it explicitely.
	  aut->set_init_state(0U);

	  // new_edge() takes 3 mandatory parameters: source state,
	  // destination state, and label.  A last optional parameter can be
	  // used to specify membership to acceptance sets.
	  //a 0, b 1, c 2, d 3, e 4, f 5,


	  aut->new_edge(0, 1, p1);

	  aut->new_edge(1, 1, p1);

	  aut->new_acc_edge(1, 2, !p1, true);


	  //aut->new_edge(5, 5, bddtrue);

	   /*
	  aut->new_edge(0, 1, p2 & !p1);
	  aut->new_edge(0, 3, p2 & p1);

	  aut->new_edge(1, 0, p2& !p1);
	  aut->new_edge(1, 2, p2&p1);

	  aut->new_edge(2, 5, p2& p1);
	  aut->new_edge(2, 4, p2 & !p1);

	  aut->new_edge(3, 5, p2&p1);
	  aut->new_edge(3, 4, p2&!p1);

	  aut->new_edge(4, 5, p2&p1);
	  aut->new_edge(4, 4, p2&!p1);

	  aut->new_edge(5, 5, p2);

	  aut->new_edge(2, 6, !p2);
	  aut->new_edge(3, 6, !p2);
	  aut->new_edge(4, 6, !p2);
	  aut->new_edge(6, 6, !p2, {0});
	  */
	  // Print the resulting automaton.
	  print_hoa(std::cout, aut);

	  //
	  cout << endl << "minimized "<< endl;
	  //aut = minimize_dfa(aut);
	  twa_graph_ptr autB = minimize_dfa(aut);
	  print_hoa(std::cout, autB);

}

