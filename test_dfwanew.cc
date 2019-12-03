#include <iostream>
#include <chrono> 

#include <spot/twaalgos/hoa.hh>
#include <spot/twaalgos/minimize.hh>

#include <spot/twa/bddprint.hh>
#include <spot/twa/twagraph.hh>


#include "dfwavar.hh"
#include "dfwa.hh"
#include "dfwanew.hh"
#include "mona.hh"
#include "spotutil.hh"


using namespace std;
using namespace spot;

void 
test(bdd_dict_ptr dict, twa_graph_ptr aut, bool t)
{
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
  
  //dfwa_var vars(dict, aut, 2, "s", 0, 6);
  
  //cout << vars.get_lower() << endl;
  //cout << vars.get_upper() << endl;

  set<unsigned> finals;
  get_final_states(aut, finals);
  
  clock_t c_start = clock();
  if(! t)
  {
      dfwa dfa_1(aut, label_cube, finals);
      clock_t b_start = clock();
      dfa_1.back_explore();
      clock_t b_end = clock();
      cout << "reachability analysis : " << 1000.0 * (b_end-b_start) / CLOCKS_PER_SEC << " ms\n";
  }else
  {
      dfwa_new dfa_2(aut, label_cube, finals);
      //dfa_2.output(cout);
      clock_t b_start = clock();
      dfa_2.back_explore();
      clock_t b_end = clock();
      cout << "reachability analysis : " << 1000.0 * (b_end-b_start) / CLOCKS_PER_SEC << " ms\n";
  }
  
  clock_t c_end = clock();
  cout << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms\n";
  
}

void
test1();

int main(int argc, char**argv)
{
  // readfile
  bdd_dict_ptr dict = spot::make_bdd_dict();
  if(argc > 3)
  {
      cout << "input formula file" << endl;
      exit(-1);
  }
  
   ifstream ltlfile(argv[1]);
    //cout << "ltlfile = " << argv[1] << "partfile = " << argv[2] << endl;
    formula input_f;
    if (ltlfile.is_open()) 
    {
        string line;
        getline (ltlfile, line);
        cout << "formula: " << line << endl;
        auto pf1 = parse_infix_psl(line.c_str());
        if (pf1.format_errors(std::cerr))
        {
            std::cerr << "error: " << line << std::endl;
            return -1;
        }
        input_f = pf1.f;
    }
    bool t = false;
    if(argc > 2)
    {
        t = true;
    }
    // construct dfa with MONA
    twa_graph_ptr aut = translate_ltlf_mona(input_f, dict);
    test(dict, aut, t);
  //test1();
  return 0;
}

void 
test1()
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
  bdd p2 = bdd_ithvar(aut->register_ap("p2"));
  
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
  dfwa_new dfab(autB, label_cube, finals);
  dfab.output(cout);
  
  set<unsigned> finalsa;
  finalsa.insert(2);
  finalsa.insert(3);
  finalsa.insert(4);
  dfwa_new dfaa(aut, label_cube, finalsa);
  
  dfaa.output(cout);
}
