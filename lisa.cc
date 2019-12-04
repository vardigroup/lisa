
// standard 
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <unordered_map>
#include <string>
#include <vector>
#include <queue>
#include <list>
#include <set>
#include <chrono> 

// spot 
#include <spot/twaalgos/hoa.hh>
#include <spot/twaalgos/translate.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/product.hh>
#include <spot/twaalgos/word.hh>
#include <spot/twaalgos/degen.hh>
#include <spot/twaalgos/remprop.hh>
#include <spot/twaalgos/minimize.hh>

#include <spot/twa/twagraph.hh>

#include <spot/tl/formula.hh>
#include <spot/tl/print.hh>
#include <spot/tl/parse.hh>
#include <spot/tl/relabel.hh>
#include <spot/tl/ltlf.hh>
#include <spot/tl/simplify.hh>

#include <spot/misc/optionmap.hh>
#include <spot/misc/timer.hh>

#include "mona.hh"
#include "spotutil.hh"
#include "dfwavar.hh"
#include "dfwa.hh"
#include "debug.hh"
#include "synt.hh"
#include "dfwamin.hh"
#include "dfwamin2.hh"
#include "dfwamin3.hh"
#include "minimize.hh"
#include "synt.hh"

using namespace spot;
using namespace std;

#define info(o) cout << "[INFO] " << ( o ) << endl
#define erro(o) cerr << "[ERRO] " << ( o ) << endl

// options

static struct opt_t
{
	const char* _ltlfile_name = nullptr;
	const char* _parfile_name = nullptr;

	bool _symbolic = true;
	bool _minimization = false;

	unsigned _num_ap_for_mona = 7;
	unsigned _num_product = 6;
	unsigned _num_st_for_single = 800;
	unsigned _num_st_for_product = 2500;
	int _num_last_automata = -1;

	bool _synthesis = false;
	bool _out_start = false;
	bool _env_first = false;

	uint8_t _bdd = 0;

}* opt;


void 
get_formulas(vector<formula>& lst, formula f)
{
    cout << "Breaking formula into small pieces..." << endl;
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
    //cout << "split formula " << endl;
}

class dfwa_pair
{
public:
	unsigned _num_states;
	bool _is_explicit;
	twa_graph_ptr _twa;
	dfwa* _dfa = nullptr;
	unsigned _num_propduct = 0;

	formula _formula;

	dfwa_pair(twa_graph_ptr aut, unsigned num_states, bool is_explicit, formula& f)
	: _num_states(num_states), _is_explicit(is_explicit), _formula(f)
	{
		_twa = aut;
	}
	dfwa_pair(dfwa* aut, unsigned num_states, bool is_explicit, formula& f)
		: _num_states(num_states), _is_explicit(is_explicit), _formula(f)
	{
		_dfa = aut;
		//_twa = new shared_ptr<twa_graph>(new twa_graph);
	}
};

struct GreaterThanByDfwaSize
{
  bool operator()(dfwa_pair& p1, dfwa_pair& p2) const
  {
    if(p1._num_states < p2._num_states)
    {
        return false;
    }else
    if(p1._num_states == p2._num_states)
    {
    	return !p1._is_explicit;
    }
    return p1._num_states >= p2._num_states;
  }
};

bool 
compare_aut_size(twa_graph_ptr p1, twa_graph_ptr p2)
{
    if(p1->num_states() == p2->num_states())
    {
        return false;
    }
    return p1->num_states() > p2->num_states();
}

twa_graph_ptr minimize_explicit(twa_graph_ptr A)
{
	twa_graph_ptr C = spot::minimize_wdba(A);
	//A = spot::minimize_obligation(A);
	// check equivalence of two automata
#ifdef DEBUG
	string word = is_twa_equivalent(A, C);
	if(word.size() == 0)
	{
		cout << "A: equivalent two automata" << endl;
	}
#endif
	return C;
}

tuple<dfwa*, unsigned, bool>
minimize_symbolic(dfwa_ptr P)
{
	unsigned num_states;
	tuple<dfwa*, unsigned, bool> result;
	if(opt->_bdd == 1)
	{
		cudd manager;
		cudd_ptr mgr = &manager;
		// call dfa minimization
		dfwa_min min(mgr, P);
		min.minimize();
		dfwa_ptr res = min.move_dfwa();
		// check wehter the emptiness is changed
		num_states = min.get_num_min_states();
		result = make_tuple<>(&res, num_states, true);
	}else
	if(opt->_bdd == 2)
	{
		dfwa_min_bdd min(P);
		min.minimize();
		dfwa_ptr res = min.move_dfwa();
		// check wehter the emptiness is changed
		num_states = min.get_num_min_states();
		result = make_tuple<>(&res, num_states, true);
	}else
	{
		dfwa_min_sylvan min(P);
		min.minimize();
		dfwa_ptr res = min.move_dfwa();
		// check wehter the emptiness is changed
		num_states = min.get_num_min_states();
		result = make_tuple<>(&res, num_states, true);
	}
	return result;
}

tuple<dfwa*, unsigned, bool>
make_product(bdd_dict_ptr dict, dfwa* A, dfwa* B, unsigned num_prod)
{
	unsigned num_states;

	dfwa_ptr P = product_dfwa_and(*A, *B);
	//cout << "labels in product and: " << P._label_cube << endl;
	//cout << "state_0 in product and: " << P._curr_cube << endl;
	//cout << "state_1 in product and: " << P._next_cube << endl;
	//cout << "product: " << endl;
	//P.output(cout);

	//unsigned var_num = 1;
	//cout << "condition: " << (P._state_vars.get_var_num(0) > var_num) << endl;

	if(num_prod > (opt->_num_product))
	{
		tuple<dfwa*, unsigned, bool> result = minimize_symbolic(P);
		delete &P;
		cout << "return from minimal product..." << endl;
		return result;
	}else
	{
		num_states = bdd_nodecount(P._trans);
		return make_tuple<>(&P, num_states, false);
	}

}



void print_usage()
{
	cout << "Usage: lisa [OPTION...] [FILENAME[/COL]...]" << endl;
	cout << "Read a formula file and output the number of states of the constructed DFA" << endl << endl;
	cout << " Input options:" << endl;
	cout << " -h  " << "                  show this help page" << endl;
	cout << " -exp" << "                  use only explicit method (default false)" << endl;
	cout << " -min" << "                  minimize the last symbolic DFA (default false)" << endl;
	cout << " -syn" << "                  synthesize after DFA construction (default false)" << endl;
	cout << " -bdd" << "                  use buddy for DFA minimization" << endl;
	cout << " -syl" << "                  use sylvan for DFA minimization (default)" << endl;
	cout << " -cdd" << "                  use cudd for DFA minimization" << endl;
	cout << " -nap" << "  <int>           number of atomic propositions for calling mona (default 7)" << endl;
	cout << " -npr" << "  <int>           number of products for calling minimization (default 6)" << endl;
	cout << " -nia" << "  <int>           number of states of individual DFA for calling symbolic approach (default 800)" << endl;
	cout << " -npa" << "  <int>           number of states of product DFA for calling symbolic approach (default 2500)" << endl;
	cout << " -lst" << "  <int>           number of last automata for calling symbolic approach (default -1)" << endl;
	cout << " -out" << "                  print out the wining strategy if realizable" << endl;
	cout << " -part" << " <file>          the file specifying the input and output propositions" << endl;
	cout << " -ltlf" << " <file>          the file specifying the input LTLf formula" << endl;
	cout << " -env" << "                  environment plays first" << endl;
}

void parse_opt(int argc, char** argv)
{
	// first one is lisa, separated by space
	if(argc == 1)
	{
		print_usage();
	}
	for(int i = 1; i < argc; i ++)
	{
		string s(argv[i]);
		//cout << argv[i] << endl;
		if(s.size() == 0)
		{
			continue;
		}
		if(s == "-exp")
		{
			opt->_symbolic = false;
			//cout << "hello" << s << endl;
			continue;
		}
		if(s == "-min")
		{
			opt->_minimization = true;
			continue;
		}
		if(s == "-syn")
		{
			opt->_synthesis = true;
			continue;
		}
		if(s == "-out")
		{
			opt->_out_start = true;
			continue;
		}
		if(s == "-env")
		{
			opt->_env_first = true;
			continue;
		}
		if(s == "-nap" && i + 1 < argc)
		{
			opt->_num_ap_for_mona = stoi(argv[i + 1]);
			i ++;
			continue;
		}
		if(s == "-npr" && i + 1 < argc)
		{
			opt->_num_product = stoi(argv[i + 1]);
			i ++;
			continue;
		}
		if(s == "-nia" && i + 1 < argc)
		{
			opt->_num_st_for_single = stoi(argv[i + 1]);
			i ++;
			continue;
		}
		if(s == "-lst" && i + 1 < argc)
		{
			opt->_num_last_automata = stoi(argv[i + 1]);
			i ++;
			continue;
		}
		if(s == "-npa" && i + 1 < argc)
		{
			opt->_num_st_for_product = stoi(argv[i + 1]);
			i ++;
			continue;
		}
		if(s == "-ltlf" && i + 1 < argc)
		{
			opt->_ltlfile_name = argv[i + 1];
			//cout << "hello" << argv[i+1] << endl;
			i ++;
			continue;
		}
		if(s == "-part" && i + 1 < argc)
		{
			opt->_parfile_name = argv[i + 1];
			i ++;
			continue;
		}
		if(s == "-cdd")
		{
			opt->_bdd = 1;
			continue;
		}
		if(s == "-bdd")
		{
			opt->_bdd = 2;
			continue;
		}
		if(s == "-syl")
		{
			opt->_bdd = 0;
			continue;
		}
		if(s == "-h")
		{
			print_usage();
			exit(0);
		}else
		{
			erro("wrong input options: " + s);
			print_usage();
			exit(-1);
		}
	}
	// validity checking
	if(opt->_ltlfile_name == nullptr )
	{
		erro( "missing LTLf file name");
		exit(-1);
	}
	if(opt->_synthesis && ( opt->_parfile_name == nullptr))
	{
		erro("missing proposition partition file name");
		exit(-1);
	}

}

dfwa*
symbolize_twa(bdd_dict_ptr dict, twa_graph_ptr aut)
{
	// there is alive states
	bdd label_cube = bddtrue;
	for (auto f : aut->ap())
	{
		bdd f_var = bdd_ithvar(aut->register_ap(f));
		label_cube = label_cube & f_var;
		//cout << "formula : " << f << " index: " << aut->register_ap(f) << endl;
	}

	set<unsigned> finals_aut;
	compute_final_states(aut, finals_aut);
	/*
	cout << "final states: " << endl;
	for(unsigned k : finals_aut)
	{
		cout << "final: " << k << endl;
	}*/
	/*
	if(aut->num_states() < 20)
	{
		ofstream outfile("output" + to_string(number) + ".hoa");
		print_hoa(outfile, aut);
		for(unsigned k : finals_aut)
		{
			cout << "output: " << k << endl;
		}
		number ++;
	}*/
	// now compute dfwa
	dfwa* A = new dfwa(aut, label_cube, finals_aut);
	return A;
}
void
read_from_part_file(const char *file_name, vector<string>& input, vector<string>& output)
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
                exit(-1);
            }
        }
    }
}


int 
main(int argc, char** argv)
{
    opt_t o;
    opt = &o;
    parse_opt(argc, argv);

    ifstream ltlfile(opt->_ltlfile_name);
    string line;
    priority_queue<dfwa_pair, std::vector<dfwa_pair>, GreaterThanByDfwaSize> autlist;
    clock_t c_start = clock();
    spot::bdd_dict_ptr dict = spot::make_bdd_dict();
    formula input_f;
	cout << "Starting the decomposition phase" << endl;
    if (ltlfile.is_open()) 
    {
        getline (ltlfile, line);
        //cout << "formula: " << line << endl;
        auto pf1 = spot::parse_infix_psl(line.c_str());
        if (pf1.format_errors(std::cerr))
        {
            std::cerr << "Error: " << line << std::endl;
            return -1;
        }
        // formula 
        input_f = pf1.f;
        vector<formula> lst;
        //cout << "parsed: " << pf1.f << endl;
        get_formulas(lst, pf1.f);
        /*
        set<formula> ap_set;
        get_formula_aps(input_f, ap_set);
        for(formula f : ap_set)
        {
        	int index = dict->register_acceptance_variable(f, opt);
        	cout << f << "->" <<  index << endl;
        }*/
        //reorganize_formulas(lst);
        while(lst.size() > 0)
        {
            // translating automata
            formula f = lst.back();
            lst.pop_back();
            // cout << str_psl(f, true) << endl;
            twa_graph_ptr aut = trans_formula(f, dict, opt->_num_ap_for_mona);
            // cout << aut->num_states() << endl;
            dfwa_pair pair(aut, aut->num_states(), true, f);
            pair._num_propduct = 0;
            // cout << "st = " << aut->num_states() << endl;
            autlist.push(pair);
        }
        ltlfile.close();
    }

    //cout << "splited formulas" << endl;
    // do products 
    bdd_autoreorder(BDD_REORDER_WIN2ITE);
	cout << "Starting the composition phase" << endl;

    set<twa_graph_ptr> optimized;
    while(autlist.size() > 1) 
    {
        cout << "Number of DFAs in the set: " << autlist.size() << endl;
        dfwa_pair first = autlist.top();
        autlist.pop();
        dfwa_pair second = autlist.top();
        autlist.pop();
        cout << "Number of states or nodes in M1 and M2: " << first._num_states
        		<< ",  " << second._num_states << endl;
        formula result_formula = formula::And({first._formula, second._formula});
        //cout << result_formula << endl;
        bool must_symbolic = opt->_num_last_automata > 0 && autlist.size() + 2 <= opt->_num_last_automata;
        if(first._is_explicit && second._is_explicit)
        {
        	twa_graph_ptr A = first._twa;
        	twa_graph_ptr B = second._twa;
        	if(optimized.find(A) == optimized.end())
        	{
        		A = minimize_explicit(A);
        		optimized.insert(A);
        	}
        	if(optimized.find(B) == optimized.end())
        	{
				B = minimize_explicit(B);
				optimized.insert(B);
			}

        	if( !opt->_symbolic || (! must_symbolic && (A->num_states() < opt->_num_st_for_single && B->num_states() < opt->_num_st_for_single
        	&& (A->num_states() * B->num_states() < opt->_num_st_for_product))))
        	{
        		// explict representation used
        		twa_graph_ptr P = spot::product(A, B);
        		//cout << "explicit minimization starts..." << endl;
        		P = spot::minimize_wdba(P);
        		optimized.insert(P);
        		dfwa_pair pair(P, P->num_states(), true, result_formula);
        		pair._num_propduct = 1;
        		cout << "Number of states in explicit product is: " << P->num_states() << endl;
        		autlist.push(pair);
        	}else
        	{
        		dfwa* fst = symbolize_twa(dict, A);
        		dfwa* snd = symbolize_twa(dict, B);
        		tuple<dfwa*, unsigned, bool> result = make_product(dict, fst, snd, 2);
        		dfwa_pair pair(get<0>(result), get<1>(result), false, result_formula);
        		if(get<2>(result))
        		{
        			pair._num_propduct = 1;
        		}else
        		{
        			pair._num_propduct = 2;
        		}
        		cout << "Number of nodes in symbolic product is: " << get<1>(result) << endl;
        		autlist.push(pair);
        		delete fst;
        		delete snd;
        	}
        }else if(first._is_explicit || second._is_explicit)
        {
        	// needs symbolic automata
        	dfwa* A = nullptr;
        	dfwa* B = nullptr;

        	if(first._is_explicit)
        	{
        		twa_graph_ptr aut = first._twa;
        		B = second._dfa;
        		// make sure it is weak DBA
				if (optimized.find(aut) == optimized.end())
				{
					aut = minimize_explicit(aut);
				}
				// now compute dfwa
				A = symbolize_twa(dict, aut);
        	}else
        	{
        		twa_graph_ptr aut = second._twa;
        		B = first._dfa;
        		if (optimized.find(aut) == optimized.end())
        		{
        			aut = minimize_explicit(aut);
        		}
        		// now compute dfwa
        		A = symbolize_twa(dict, aut);
        	}
        	unsigned num = first._num_propduct + second._num_propduct + 1;
			tuple<dfwa*, unsigned, bool> result = make_product(dict, A, B, num);
			dfwa_pair pair(get<0>(result), get<1>(result), false, result_formula);
			if(get<2>(result))
			{
				pair._num_propduct = 1;
			}else
			{
				pair._num_propduct = num;
			}
			cout << "Number of nodes in symbolic product is: " <<  get<1>(result) << endl;
			autlist.push(pair);
			delete B;
        }else
        {
        	// two symbolic automata
        	dfwa* A = first._dfa;
        	dfwa* B = second._dfa;
        	unsigned num = first._num_propduct + second._num_propduct;
        	tuple<dfwa*, unsigned, bool> result = make_product(dict, A, B, num);
        	dfwa_pair pair(get<0>(result), get<1>(result), false, result_formula);
        	cout << "Number of nodes in symbolic product is: " << get<1>(result) << endl;
        	autlist.push(pair);
        	if(get<2>(result))
        	{
        		pair._num_propduct = 1;
        	}else
        	{
        		pair._num_propduct = num;
        	}

        	delete A;
        	delete B;
        }
    }
    clock_t c_end = clock();
    cout << "Finished constructing minimal dfa in "
    		 << 1000.0 * (c_end - c_start)/CLOCKS_PER_SEC << "ms ..." << endl;
    dfwa_pair pair = autlist.top();
	cout << "Number of states (or nodes) is: " << pair._num_states << endl;
    if(pair._is_explicit && optimized.find(pair._twa) == optimized.end())
    {
    	// in case we only have one DFA and it is not minimized
    	pair._twa = minimize_explicit(pair._twa);
    }
    cout << "Final result (or number of nodes): " << pair._num_states << endl;
    if(! pair._is_explicit && ! opt->_synthesis)
    {
    	if(opt->_minimization)
    	{
    		minimize_symbolic(*pair._dfa);
    	}
        delete pair._dfa;
        exit(0);
    }
    if(pair._is_explicit && ! opt->_synthesis)
    {
    	if(opt->_out_start)
    	{
    		// output
    		ofstream outfile("output.hoa");
    		print_hoa(outfile, pair._twa);
    	}
    	exit(0);
    }
    dfwa* aut = nullptr;
    if(! pair._is_explicit)
	{
    	aut = pair._dfa;
	}else
	{
		aut = symbolize_twa(dict, pair._twa);
	}
    // synthesis
	vector<string> input;
	vector<string> output;
	//cout << "read part file " << endl;
	if(opt->_parfile_name != nullptr)
	{
		read_from_part_file(opt->_parfile_name, input, output);
	}else
	{
		cerr << "Please input the file name for inputs and outputs" << endl;
		exit(-1);
	}
	//cout << "The number of nodes in transition is " << bdd_nodecount(aut->_trans) << endl;
	//cout << "finished reading part file " << endl;
	// NOTE that some propositions may not be used in DFA
	bdd input_cube = bddtrue;
	bdd output_cube = bddtrue;
    //set<formula> set_aps;
    //get_formula_aps(input_f, set_aps);
	//set<formula>::iterator it;
	map<formula, int>& var_map = dict->var_map;
#ifdef DEBUG
	map<formula, int>::const_iterator iter = var_map.begin();
	while (iter != var_map.end())
	{
		formula key = iter->first;
		int value = iter->second;
		cout << key << " -> " << value << endl;
		iter ++;
	}
#endif
	//cout << dict->var_map << endl;
	//cout << "partition of propositions" << endl;
	for(string& in : input)
	{
		formula f = formula::ap(in);
		// not in var map
		if(var_map.count(f) == 0)
		{
			continue;
		}
		//cout << "in: " << in << " ";
		int var_index = dict->varnum(f);
		//cout << var_index << endl;
		bdd p = bdd_ithvar(var_index);
		input_cube = input_cube & p;
	}
	//cout << "done with input propositions" << endl;
	for(string& out : output)
	{
		formula f = formula::ap(out);
		// not in var map
		if(var_map.count(f) == 0)
		{
			continue;
		}
		//cout << "out: " << out  << " ";
		int var_index = dict->varnum(f);
		//cout << var_index << endl;
		bdd p = bdd_ithvar(var_index);
		output_cube = output_cube & p;
	}
	string alive_ap(ALIVE_AP);
	int var_index = dict->varnum(formula::ap(alive_ap));
	bdd p = bdd_ithvar(var_index);
	output_cube = output_cube & p;
	//cout << "out: " << alive_ap << " " << var_index << endl;

	//aut->output(cout);
	{
		clock_t c_start = clock();
		auto t_start = chrono::high_resolution_clock::now();
		//aut->output_dfwa(cout);
		synt syn(*aut, input_cube, output_cube);
		cout << "Starting to synthesize " << endl;
		if(opt->_env_first)
		{
			syn.env_play_first();
			cout << "Environment will play first" << endl;
		}else{
			cout << "System will play first" << endl;
		}
		syn.is_realizable();
		if(opt->_out_start)
		{
			syn.synthesize();
		}
		cout << "Finished synthesizing" << endl;

		clock_t c_end = clock();
		cout << 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC << " ms\n";
		auto t_end = chrono::high_resolution_clock::now();
		cout << "Total CPU time used: "
			  << 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC << " ms\n"
			  << "Total wall clock time passed: "
			  << std::chrono::duration<double, std::milli>(t_end-t_start).count()
			  << " ms\n";
	}
	if(aut != nullptr)
	{
		delete aut;
	}
	opt = nullptr;
	//dict->unregister_all_my_variables(opt);
	//dict->unregister_acceptance_variable(opt);
}
