#include <iostream>
#include <utility>
#include <iostream>
#include <string>
#include <chrono> 

#include <spot/twaalgos/hoa.hh>
#include <spot/twaalgos/minimize.hh>

#include <spot/twa/bddprint.hh>
#include <spot/twa/twagraph.hh>

#include <spot/tl/formula.hh>
#include <spot/tl/print.hh>
#include <spot/tl/parse.hh>


#include "dfwavar.hh"
#include "dfwa.hh"
#include "synt.hh"
#include "mona.hh"
#include "spotutil.hh"
#include "debug.hh"

using namespace std;
using namespace spot;

// read mona file


/*-------------------------------------------------------------------*/
// construct the DFA from mona output file
/*-------------------------------------------------------------------*/
/*
twa_graph_ptr
read_dfwa_from_mona_file(dfwa& dfa, const char * file_name, bdd_dict_ptr dict)
{
    // dfa file
    ifstream dfa_file(file_name);
    // ordered sequence of propositions
    vector<string> atom_props;
    // final states in the DFA
    vector<bool> final_states;
    // stores the node number of the state
    vector<int> behaviour;
    // node tuples for the transitions
    vector<tuple<int, int, int>> node_data;
    // number of vars
    int num_vars;
    // number of states
    int num_states;
    // init_state
    int init_state;
    // number of bdd nodes
    int num_bdd_node;
    
    if (dfa_file.is_open()) 
    {
        bool flag = false;
        string line;
        while(getline(dfa_file, line))
        {
            vector<string> temp;
            if(flag)
            {
                if(str_contain(line, "end") )
                {
                    break;
                }
                // not end then parsing structure of bdd
                str_split(line, temp, ' ');
                #ifdef DEBUG
                for(int i = 0; i < temp.size(); i ++)
                {
                    cout << temp[i] << " ";
                }
                cout <<  endl;
                cout << "#bdd = " << temp.size() << endl;
                #endif
                tuple<int, int, int> tp = make_tuple(stoi(temp[0]), stoi(temp[1]), stoi(temp[2]));
                
                node_data.push_back(tp);
            }
            if(str_contain(line, "MONA DFA")) 
            {
                DEBUG_STDOUT( "parsing starts" );
            }else
            // now we parse the output
            //number of variables: 5
            if(str_contain(line, "number of variables"))
            {
                string delimiter = ":";
                string number = line.substr(line.find(delimiter) + 1);
                num_vars = stoi(number);
                DEBUG_STDOUT("#AP=" + to_string( num_vars));
            }else
            // variables: P149 P170 P172 P53 P93
            if(str_contain(line, "variables") && !str_contain(line, "number"))
            {
                // split by white space and delete white space
                string delimiter = ":";
                line = line.substr(line.find(delimiter) + 1);
                str_split(line, atom_props, ' ');
                #ifdef DEBUG
                for(int i = 0; i < atom_props.size(); i ++)
                {
                    cout << "AP: " << atom_props[i] << endl;
                }
                #endif
            }else
            //states: 19
            if(str_contain(line, "states"))
            {
                string delimiter = ":";
                string number = line.substr(line.find(delimiter) + 1);
                num_states = stoi(number);
                DEBUG_STDOUT("#S=" + to_string( num_states));
            }else
            //initial: 0
            if(str_contain(line, "initial"))
            {
                string delimiter = ":";
                string number = line.substr(line.find(delimiter) + 1);
                init_state = stoi(number);
                DEBUG_STDOUT("I = " + to_string( init_state));
            }else
            // bdd nodes: 76
            if(str_contain(line, "bdd nodes"))
            {
                string delimiter = ":";
                string number = line.substr(line.find(delimiter) + 1);
                num_bdd_node = stoi(number);
                DEBUG_STDOUT("#bdd = " + to_string( num_bdd_node));
            }else 
            // final: -1 -1 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            if(str_contain(line, "final"))
            {
                string delimiter = ":";
                line = line.substr(line.find(delimiter) + 1);
                str_split(line, temp, ' ');
                if(temp.size() != num_states)
                {
                    DEBUG_STDERR("error final states");
                }else
                {
                    for(int i = 0;i < temp.size(); i ++)
                    {
                        final_states.push_back(temp[i] == "1");
                        #ifdef DEBUG
                        if(temp[i] == "1")
                        {
                            cout << "final : " << i << endl;
                        }
                        #endif
                    }
                }
            }else 
            // behaviour: 0 1 1 10 17 23 2 29 36 42 48 54 60 62 65 68 70 72 74
            // map state to the bdd node
            if(str_contain(line, "behaviour"))
            {
                string delimiter = ":";
                line = line.substr(line.find(delimiter) + 1);
                str_split(line, temp, ' ');
                DEBUG_STDOUT("behaviour: ");
                for(int i = 0; i < temp.size(); i ++)
                {
                    DEBUG_STDOUT(" " + temp[i]);
                    behaviour.push_back(stoi(temp[i]));
                }
                cout << endl;
                
            }else
            // bdd: start of bdd 
            if(str_contain(line, "bdd:"))
            {
                flag = true;
            }
            //end
        }
        dfa_file.close();
    }
    // now construct the dfa
    spot::twa_graph_ptr aut = make_twa_graph(dict);
    vector<bdd> bdd_props;
    // get bdd repr for propositions
    for(int i = 0; i < atom_props.size(); i ++)
    {
        bdd p = bdd_ithvar(aut->register_ap(atom_props[i]));
        bdd_props.push_back(p);
    }
    // add another state
    //bdd alive = bdd_ithvar(aut->register_ap(ALIVE_AP));
    //bdd_props.push_back(alive);
    
    //aut->set_buchi();
    //aut->prop_state_acc();
    // add one extra state for accepting state
    //aut->new_states(num_states + 1);
    //aut->set_init_state(init_state);
    
    // now construct transition system of aut
    DEBUG_STDOUT("behaviour size = " + to_string( behaviour.size()));
    // need a map to store computed results
    unordered_map<int, vector<tuple<bdd, int>>> node2bdd;
    dfwa_var vars(dict, aut, 2, "s", 0, behaviour.size() - 1);
    dfa._state_vars = vars;
    dfa._trans = bddfalse;
    dfa._finals = bddfalse;
    // new state variables
    for(int i = 0; i < behaviour.size(); i ++)
    {
        // construct the transition of a state
        //, int state , bdd label , int node_id , vector<bdd> props
        //, vector<tuple<int, int, int>> node_data
        DEBUG_STDOUT("state " + to_string(i) + " behaviour " + to_string(behaviour[i]));
        vector<tuple<bdd, int>> succs;
        if(node2bdd.find(behaviour[i]) != node2bdd.end())
        {
            succs = node2bdd[behaviour[i]];
        }else
        {
            construct_twa_trans(succs, behaviour[i], bdd_props, node_data, node2bdd);
        }
        bdd tr = bddfalse;
        bdd src = dfa._state_vars.new_value(0, i);
        for(tuple<bdd, int> &t : succs)
        {
            // from state i to get<1>(t) via get<0>(t)
            //aut->new_edge(i, get<1>(t), get<0>(t) & alive);
            bdd tgt = dfa._state_vars.new_value(1, i);
            tr = tr | (src & tgt & get<0>(t));
        }
        dfa._trans = dfa._trans | tr;
        // accepting to extra state
        if(final_states[i])
        {
            dfa._finals = dfa._finals | src;
        }
    }
    // traverse intial state
    tuple<int, int, int> node = node_data[behaviour[0]];
    dfa._init = dfa._state_vars.new_value(0, get<1>(node));
    
    // compute cubes
    dfa._curr_cube = dfa._state_vars.get_cube(0);
    dfa._next_cube = dfa._state_vars.get_cube(1);
    dfa._curr_to_next_pairs = dfa._state_vars.make_pair(0, 1);
    dfa._next_to_curr_pairs = dfa._state_vars.make_pair(1, 0);
        // compute reachable state space
    cout << "Computing reachable state space in the product..." << endl;
    dfa._reach = dfa.explore();
    cout << "Finished computing reachable state space in the product..." << endl;
    //cout << "reachable states in product: " << endl;
    //bdd_print_set(cout, result._state_vars.get_dict(), result._reach);
    //cout << endl;
    bdd all = bdd_replace(dfa._reach, dfa._curr_to_next_pairs);
    all = all & dfa._reach;
    dfa._trans = dfa._trans & all;
    cout << "Finished computing the intersection product..." << endl;
    #ifdef DEBUG
        print_hoa(std::cout, aut);
        cout << endl;
    #endif
    
    return aut;
}
*/
/*-------------------------------------------------------------------*/
// execute mona to construct DFA for the input ltlf formula
/*-------------------------------------------------------------------*/
/*
twa_graph_ptr
translate_ltlf_mona_dfwa(dfwa& dfa, formula f, bdd_dict_ptr dict)
{
    string fol_file_name = "./fol.ltlf";
    ofstream ofs (fol_file_name, ofstream::out);
    // output in a parsable formula
    ofs << str_psl(f, true);
    ofs.close();
    //cout << "formula: " << str_psl(f, true) << endl;
    // now call ltlf2fol
    string mona_file_name = "./ltlf.mona";
    string command = "./ltlf2fol NNF " + fol_file_name + " > " + mona_file_name;
    system(command.c_str());
    string dfa_file_name = "./mona.dfa";
    command = "mona -u -xw " + mona_file_name+ " >" + dfa_file_name;
    system(command.c_str());
    // if this turns to be a bottleneck, we need pthread to read from pipe
    return read_dfwa_from_mona_file(dfa, dfa_file_name.c_str(), dict);
}
*/
// read part file

void read_from_part_file(const char * file_name, vector<string>& input, vector<string>& output)
{
    // const char * file_name
    ifstream part_file(file_name);
    if (part_file.is_open()) 
    {
        bool flag = false;
        string line;
        while(getline(part_file, line))
        {
            if(str_contain(line, "inputs"))
            {
                string delimiter = ":";
                line = line.substr(line.find(delimiter) + 1);
                str_split(line, input, ' ');
            }else
            if(str_contain(line, "outputs"))
            {
                string delimiter = ":";
                line = line.substr(line.find(delimiter) + 1);
                str_split(line, output, ' ');
            }else
            {
                cout << "read partfile error!" <<endl;
                cout << file_name <<endl;
                cout << line <<endl;
            }
        }
    }
}


int main(int argc, char** argv)
{
  
    if(argc <= 1 || argc > 4)
    {
          cout << "usage: <program> <ltlf_file> <part_file>" << endl;
          exit(-1);
    }
      
    bdd_dict_ptr dict = make_bdd_dict();
    ifstream ltlfile(argv[1]);
    //cout << "ltlfile = " << argv[1] << "partfile = " << argv[2] << endl;
    formula input_f;
    if (ltlfile.is_open()) 
    {
        string line;
        getline (ltlfile, line);
        //cout << "formula: " << line << endl;
        auto pf1 = parse_infix_psl(line.c_str());
        if (pf1.format_errors(std::cerr))
        {
            std::cerr << "error: " << line << std::endl;
            return -1;
        }
        input_f = pf1.f;
    }
    twa_graph_ptr aut = translate_ltlf_mona(input_f, dict);
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

    // partition aps
    vector<string> input;
    vector<string> output;
    
    read_from_part_file(argv[2], input, output);
    bdd input_cube = bddtrue;
    bdd output_cube = bddtrue;
    for(string & in : input)
    {
        bdd p = bdd_ithvar(aut->register_ap(in.c_str()));
        input_cube = input_cube & p;
    }
    for(string & out : output)
    {
        bdd p = bdd_ithvar(aut->register_ap(out.c_str()));
        output_cube = output_cube & p;
    }
    bdd p = bdd_ithvar(aut->register_ap(ALIVE_AP));
    output_cube = output_cube & p;
    //print_hoa(std::cout, aut);
  
      //dfwa_var vars(dict, aut, 2, "s", 0, 6);
      
      //cout << vars.get_lower() << endl;
      //cout << vars.get_upper() << endl;
      
      set<unsigned> finals;
      get_final_states(aut, finals);
      dfwa dfa(aut, label_cube, finals);
      
      if(argc > 3)
      {
          dfa.output(cout);
          cout << "input cube " << endl;
          bdd_print_set(cout, dfa._state_vars.get_dict(), input_cube);
          cout << endl;
          cout << "output cube " << endl;
          bdd_print_set(cout, dfa._state_vars.get_dict(), input_cube);
          cout << endl;
      }
      clock_t c_start = clock();
      auto t_start = chrono::high_resolution_clock::now();
      synt syn(dfa, input_cube, output_cube);
      
      syn.is_realizable();
       clock_t c_end = clock();
        auto t_end = chrono::high_resolution_clock::now();
     cout << "Total CPU time used: "
              << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms\n"
              << "Total wall clock time passed: "
              << std::chrono::duration<double, std::milli>(t_end-t_start).count()
              << " ms\n";
      return 0;
}

void
test()
{
    spot::bdd_dict_ptr dict = spot::make_bdd_dict();
  // This creates an empty automaton that we have yet to fill.
  spot::twa_graph_ptr aut = make_twa_graph(dict);

  // Since a BDD is associated to every atomic proposition, the
  // register_ap() function returns a BDD variable number
  // that can be converted into a BDD using bdd_ithvar().
  bdd p1 = bdd_ithvar(aut->register_ap("i"));
  bdd p2 = bdd_ithvar(aut->register_ap("o"));
  
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
  
  bdd input_cube = p1;
  bdd output_cube = p2;
  
  

  // Set the acceptance condition of the automaton to Inf(0)&Inf(1)
  aut->set_buchi();
  //aut->prop_state_acc(true);
  // States are numbered from 0.
  aut->new_states(4);
  // The default initial state is 0, but it is always better to
  // specify it explicitely.
  aut->set_init_state(0U);

  // new_edge() takes 3 mandatory parameters: source state,
  // destination state, and label.  A last optional parameter can be
  // used to specify membership to acceptance sets.
  //a 0, b 1, c 2, d 3, e 4, f 5, 
  
  
  aut->new_edge(0, 1, (!p1 | p2) & (!p2 | p1));
  aut->new_edge(0, 2, ! ((!p1 | p2) & (!p2 | p1)));

  aut->new_edge(1, 1, p1 & p2);
  aut->new_edge(1, 2, !p1 & p2);
  aut->new_edge(1, 3, !p2);
  
  aut->new_edge(2, 2, !p1 | p2);
  aut->new_edge(2, 3, p1 & (!p2));

  aut->new_edge(3, 3, bddtrue);
    
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
  
  //dfwa_var vars(dict, aut, 2, "s", 0, 6);
  
  //cout << vars.get_lower() << endl;
  //cout << vars.get_upper() << endl;

  set<unsigned> finals;
  finals.insert(3);
  //finals.insert(3);
  //finals.insert(4);
  dfwa dfa(aut, label_cube, finals);
  
  dfa.output(cout);
  cout << "input cube " << endl;
  bdd_print_set(cout, dfa._state_vars.get_dict(), input_cube);
  cout << endl;
  cout << "output cube " << endl;
  bdd_print_set(cout, dfa._state_vars.get_dict(), input_cube);
  cout << endl;
  
  synt syn(dfa, input_cube, output_cube);
  
  syn.is_realizable();
  

  
}
