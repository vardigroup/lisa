#include <iostream>

#include <spot/twaalgos/hoa.hh>
#include <spot/twaalgos/minimize.hh>

#include <spot/twa/bddprint.hh>
#include <spot/twa/twagraph.hh>


#include "dfwavar.hh"
#include "dfwa.hh"

using namespace std;
using namespace spot;

int main(void)
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

  result.output_dfwa(cout);

  // minimization

  delete &result;
  
  return 0;
}