
/*****************************************************************************
Yong Li
* 
*****************************************************************************/
#include "mona.hh"

#include "debug.hh"
#include "spotutil.hh"

/*-------------------------------------------------------------------*/
// check whether str contains match
/*-------------------------------------------------------------------*/
bool
str_contain(string str, const char* match)
{
    return str.find(match) != string::npos;
}
/*-------------------------------------------------------------------*/
// split str by whitespace and return the list of nonspace strings
/*-------------------------------------------------------------------*/
void
str_split(string str, vector<string>& result, char delim)
{
    std::stringstream ss(str);
    std::string token;
    while (getline(ss, token, delim)) 
    {
        if(! token.empty())
        {
            result.push_back(token);
        }
    }
}
/*-------------------------------------------------------------------*/
// construct the BDD from mona node file
// used a map to store the nodes which have been visited before
//TODO: construct MTBDD for the successors of the state, this should be more efficient
/*-------------------------------------------------------------------*/
void
construct_twa_trans(
  vector<tuple<bdd, int>>& succs
, int node_id
, vector<bdd>& props
, vector<tuple<int, int, int>>& node_data
, unordered_map<int, vector<tuple<bdd, int>>>& node2bdd)
{
    if(node2bdd.find(node_id) != node2bdd.end())
    {
        // first check whether already computed
        succs = node2bdd[node_id];
        return ;
    }
    // tuple data (x, l, r) where x is the proposition index
    // l is the left node
    // r is the right node
    tuple<int, int, int> node = node_data[node_id];
    // reached the leaf node
    int prop = get<0>(node);
    if(prop == -1)
    {
        // next state is get<1>(node), label is true
        succs.push_back(make_tuple(bddtrue, get<1>(node)));
    }else
    {
        // get the variable dd representation
        bdd prop_dd = props[prop];
        // low branch
        vector<tuple<bdd, int>> lsuccs;
        construct_twa_trans(lsuccs, get<1>(node), props, node_data, node2bdd);
        for (vector<tuple<bdd, int>>::iterator it = lsuccs.begin() ; it != lsuccs.end(); ++it)
        {
            tuple<bdd, int> tp = *it;
            succs.push_back(make_tuple(!prop_dd & get<0>(tp), get<1>(tp)));
        }
        // high branch
        vector<tuple<bdd, int>> rsuccs;
        construct_twa_trans(rsuccs, get<2>(node), props, node_data, node2bdd);
        for (vector<tuple<bdd, int>>::iterator it = rsuccs.begin() ; it != rsuccs.end(); ++it)
        {
            tuple<bdd, int> tp = *it;
            succs.push_back(make_tuple(prop_dd & get<0>(tp), get<1>(tp)));
        }
    }
    node2bdd[node_id] = succs;
}


/*-------------------------------------------------------------------*/
// construct the DFA from mona output file
/*-------------------------------------------------------------------*/
twa_graph_ptr
read_from_mona_file(const char * file_name, bdd_dict_ptr dict)
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
    bdd alive = bdd_ithvar(aut->register_ap(ALIVE_AP));
    //bdd_props.push_back(alive);
    
    aut->set_buchi();
    aut->prop_state_acc();
    // add one extra state for accepting state
    aut->new_states(num_states + 1);
    //aut->set_init_state(init_state);
    
    // now construct transition system of aut
    DEBUG_STDOUT("behaviour size = " + to_string( behaviour.size()));
    // need a map to store computed results
    unordered_map<int, vector<tuple<bdd, int>>> node2bdd;
    
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
        for(tuple<bdd, int> &t : succs)
        {
            aut->new_edge(i, get<1>(t), get<0>(t) & alive);
        }
        
        // accepting to extra state
        if(final_states[i])
        {
            aut->new_edge(i, num_states, !alive);
        }
    }
    
    // now set accepting states
    aut->new_edge(num_states, num_states, !alive, {0});
    aut->merge_edges();
    // traverse intial state
    tuple<int, int, int> node = node_data[behaviour[0]];
    aut->set_init_state(get<1>(node));
    
    #ifdef DEBUG
        print_hoa(std::cout, aut);
        cout << endl;
    #endif
    
    
    return aut;
}

/*-------------------------------------------------------------------*/
// execute mona to construct DFA for the input ltlf formula
// NOT depend on ltlf2fol of Syft anymore
/*-------------------------------------------------------------------*/
twa_graph_ptr
translate_ltlf_mona(formula f, bdd_dict_ptr dict)
{
    // code depending on ltlf2fol of Syft
    /*
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
    */
    string mona_file_name = "./ltlf.mona";
    ofstream ofs(mona_file_name, ofstream::out);
    ofs << "#LTLf formula" << endl;
    ofs << "#" << str_psl(f, true) << endl;
    formula bnf = get_bnf(f);
    ofs << "# Backus normal form" << endl;
    ofs << "#" << str_psl(bnf, true) << endl;
    // the BNF form, and then convert it to fol formula
    trans_ltlf2fol(ofs, bnf);
    ofs.close();
    string dfa_file_name = "./mona.dfa";
    string command = "mona -u -xw " + mona_file_name+ " >" + dfa_file_name;
    int r = system(command.c_str());
    // if this turns to be a bottleneck, we need pthread to read from pipe
    return read_from_mona_file(dfa_file_name.c_str(), dict);
}
