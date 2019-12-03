
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
#include <spot/twaalgos/powerset.hh>

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

	uint8_t _bdd = 0;

}* opt;


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
	cout << " -min" << "                  use minimization explicit method (default false)" << endl;
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
		cout << "formula : " << f << " index: " << aut->register_ap(f) << endl;
	}

	set<unsigned> finals_aut;
	compute_final_states(aut, finals_aut);
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
    if (ltlfile.is_open()) 
    {
        getline (ltlfile, line);
        cout << "formula: " << line << endl;
        auto pf1 = spot::parse_infix_psl(line.c_str());
        if (pf1.format_errors(std::cerr))
        {
            std::cerr << "error: " << line << std::endl;
            return -1;
        }
        // formula 
        input_f = pf1.f;
        vector<formula> lst;
        cout << "parsed: " << pf1.f << endl;
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
            cout << "subformula: " << str_psl(f, true) << endl;
            twa_graph_ptr aut = trans_formula(f, dict, opt->_num_ap_for_mona);
            cout << "NFA #S = " << aut->num_states() << endl;
            
            // determinize
            bool input_is_det = is_deterministic(aut);
            twa_graph_ptr det_aut;
            if(input_is_det) 
            {
                det_aut = std::const_pointer_cast<twa_graph>(aut);
            }else
            {
                det_aut = tgba_powerset(aut);
            }
            cout << "DFA #S = " << det_aut->num_states() << endl;
            // minimize
        	twa_graph_ptr P = spot::minimize_wdba(det_aut);
            cout << "mDFA #S = " << P->num_states() << endl;
        }
        ltlfile.close();
    }
    return 0;
}
