
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
#include <spot/tl/print.hh>
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
#include "minimize.hh"

using namespace spot;
using namespace std;

static void
print_raw_formula(const char *f)
{
  std::cout << "formula: " << f << std::endl;
}

int num_aps = 0;

unsigned
traverse_formula(formula f)
{
    // use id as the key to record the length and the number of props in a formula
    //std::cout << "operator" << f.id() << std::endl;
    if(f.kind() == op::Not || f.kind() == op::X || f.kind() == op::F || f.kind() == op::G )
    {
        return 1 + traverse_formula(f[0]);
    }else
    if(f.kind() == op::U || f.kind() == op::R || f.kind() == op::M || f.kind() == op::W  || f.kind() == op::Implies || f.kind() == op::Equiv )
    {
        return 1 + traverse_formula(f[0]) + traverse_formula(f[1]);
    }else
    if(f.kind() == op::And || f.kind() == op::Or)
    {
        unsigned count = 1;
        for(formula child: f)
        {
            count += traverse_formula(child);
        }
        return count;
    }else 
    {
        return 1;
    }
}

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
bool 
compare_ltl_size(formula& f1, formula& f2)
{
    unsigned size_1 = list_formula_props(f1).size() + traverse_formula(f1);
    unsigned size_2 = list_formula_props(f2).size() + traverse_formula(f1);
    if(size_1 == size_2)
    {
        return false;
    }
    return size_1 > size_2;
}

struct GreaterThanBySize
{
  bool operator()(formula& f1, formula& f2) const
  {
    unsigned size_1 = list_formula_props(f1).size() + traverse_formula(f1);
    unsigned size_2 = list_formula_props(f2).size() + traverse_formula(f2);
    if(size_1 < size_2)
    {
        return false;
    }
    return size_1 >= size_2;
  }
};
// change it for 
void
reorganize_formulas(vector<formula> & lst)
{
	const int num = 100;
    if(lst.size() < num)
    {
        return ;
    }
    priority_queue<formula, std::vector<formula>, GreaterThanBySize> pq;
    // priority queue for number of propositions
    // until we have 100 subfomulas
    while(lst.size() > 0)
    {
        pq.push(lst.back());
        lst.pop_back();
    }
    /*
    while(pq.size() > 0)
    {
        formula f = pq.top();
        cout << list_formula_props(f).size() + traverse_formula(f) << " ";
        pq.pop();
    }*/
    // reorganize formulas
    int cc = 0;
    while(pq.size() > num)
    {
        formula f1 = pq.top();
        pq.pop();
        formula f2 = pq.top();
        pq.pop();
        // connect them together
        cout << "f1 = " << traverse_formula(f1) << endl;
        cout << "f2 = " << traverse_formula(f2) << endl;
        formula f = formula::And({f1, f2});
        pq.push(f);
        cc ++;
    }
    cout << "number = " << cc << endl;
    //exit(-1);
    
    while(pq.size() > 0)
    {
        lst.push_back(pq.top());
        pq.pop();
    }
}

// input one formula
twa_graph_ptr
parse_formula(const char* f, bdd_dict_ptr dict)
{
    print_raw_formula(f);

  {
    // parse the input formula
    auto pf1 = spot::parse_infix_psl(f);
    if (pf1.format_errors(std::cerr))
    {
      std::cerr << "error: " << f << std::endl;
      return nullptr;
    }
    // formula 
    auto f1 = pf1.f;
    std::cout << f1 << std::endl;
    // make a dictionary
    //spot::bdd_dict_ptr dict = spot::make_bdd_dict();
    //traverse_formula(f1, dict);
    //std::cout<< num_aps << std::endl;
    
    // translate to ltlf
    //auto ltlf = spot::from_ltlf(f1, "alive");
    spot::translator trans(dict);
    //trans.set_type(spot::postprocessor::BA);
    trans.set_pref(spot::postprocessor::Small);
    spot::twa_graph_ptr aut = trans.run(spot::from_ltlf(pf1.f));

    // removing ap alive is not correct
    //spot::remove_ap rem;
    //rem.add_ap("alive");
    //aut = rem.strip(aut);

    spot::postprocessor post;
    //post.set_type(spot::postprocessor::BA);
    post.set_pref(spot::postprocessor::Small); // or ::Deterministic
    aut = post.run(aut);

    //print_hoa(std::cout, aut) << '\n';
    
    return aut;
  }
}

class dfwa_pair
{
public:
	unsigned _num_states;
	bool _is_explicit;
	twa_graph_ptr _twa;
	dfwa* _dfa = nullptr;

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

pair<dfwa*, unsigned>
make_product(bdd_dict_ptr dict, dfwa* A, dfwa* B)
{
	unsigned num_states;

	dfwa_ptr P = product_dfwa_and(*A, *B);
	cout << "labels in product and: " << P._label_cube << endl;
	cout << "state_0 in product and: " << P._curr_cube << endl;
	cout << "state_1 in product and: " << P._next_cube << endl;
	cout << "product: " << endl;
	//P.output(cout);

	unsigned var_num = 1;
	cout << "condition: " << (P._state_vars.get_var_num(0) > var_num) << endl;
	if(P._state_vars.get_var_num(0) > var_num)
	{
		cudd manager;
		cudd_ptr mgr = &manager;
		// call dfa minimization
		dfwa_min min(mgr, P);
		min.minimize();
		// make sure it is in heap
		dfwa_ptr res = min.move_dfwa();
		// check wehter the emptiness is changed

		num_states = min.get_num_min_states();
		delete &P;
		cout << "return from minimal product..." << endl;
		return make_pair<>(&res, num_states);
	}else
	{
		num_states = bdd_nodecount(P._trans);
		return make_pair<>(&P, num_states);
	}

}

static unsigned number = 0;

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

	if(aut->num_states() < 20)
	{
		ofstream outfile("output" + to_string(number) + ".hoa");
		print_hoa(outfile, aut);
		for(unsigned k : finals_aut)
		{
			cout << "output: " << k << endl;
		}
		number ++;
	}
	// now compute dfwa
	dfwa* A = new dfwa(aut, label_cube, finals_aut);
	return A;
}

int 
main(int argc, char** argv)
{
    if(argc < 1 || argc > 3)
        std::cout << "please input formula file" << std::endl;
    ifstream ltlfile(argv[1]);
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
        reorganize_formulas(lst);
        /*
        cout << "formulas splited: " << lst.size() << endl;
        ofstream subfile("subformula.txt");
        while(lst.size() > 0)
        {
            formula f = lst.back();
            lst.pop_back();
            subfile << str_psl(f, true) << endl;
        }
        return 1;
        */
        while(lst.size() > 0)
        {
            // translating automata
            formula f = lst.back();
            lst.pop_back();
            cout << str_psl(f, true) << endl;
            twa_graph_ptr aut = trans_formula(f, dict);
            cout << aut->num_states() << endl;
            dfwa_pair pair(aut, aut->num_states(), true, f);
            cout << "st = " << aut->num_states() << endl;
            autlist.push(pair);
        }
        ltlfile.close();
    }
    //return 1;
    cout << "splited formulas" << endl;
    // do products 
    set<twa_graph_ptr> optimized;
    const unsigned int SIG_NUM = 800;
    const unsigned int PRO_NUM = 2500;
    while(autlist.size() > 1) 
    {
        cout << "loop starts: sorting " << autlist.size() << endl;
        cout << "sorted ..." << endl;
        dfwa_pair first = autlist.top();
        autlist.pop();
        dfwa_pair second = autlist.top();
        autlist.pop();
        cout << "poped two elements #fst = " << first._num_states
        		<< " #snd = " << second._num_states << endl;
        formula result_formula = formula::And({first._formula, second._formula});
        cout << result_formula << endl;
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

        	if(A->num_states() < SIG_NUM && B->num_states() < SIG_NUM
        	&& (A->num_states() * B->num_states() < PRO_NUM))
        	{
        		// explict representation used
        		twa_graph_ptr P = ::product(A, B);
        		P = spot::minimize_wdba(P);
        		optimized.insert(P);
        		dfwa_pair pair(P, P->num_states(), true, result_formula);
        		cout << "explicit product finished, result has " << P->num_states() << " states" << endl;
        		autlist.push(pair);
        	}else
        	{
        		dfwa* fst = symbolize_twa(dict, A);
        		dfwa* snd = symbolize_twa(dict, B);
        		cout << "first: " << endl;
        		//fst->output(cout);
        		cout << "second: " << endl;
        		//snd->output(cout);
        		pair<dfwa*, unsigned> result = make_product(dict, fst, snd);
        		/*if(get<0>(result)->is_empty())
        		        	{

        		        		cout << "empty: " << str_psl(result_formula, true) << endl;
        		        		exit(-1);
        		        	}*/
        		dfwa_pair pair(get<0>(result), get<1>(result), false, result_formula);
        		cout << "symbolic product finished, result has " << get<1>(result) << " states" << endl;
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
			pair<dfwa*, unsigned> result = make_product(dict, A, B);
			/*if(get<0>(result)->is_empty())
			        	{
			        		cout << "empty: " << str_psl(result_formula, true) << endl;
			        		exit(-1);
			        	}*/
			dfwa_pair pair(get<0>(result), get<1>(result), false, result_formula);
			cout << "symbolic product finished, result has " << get<1>(result) << " states" << endl;
			autlist.push(pair);
			delete B;
        }else
        {
        	// two symbolic automata
        	dfwa* A = first._dfa;
        	dfwa* B = second._dfa;
        	pair<dfwa*, unsigned> result = make_product(dict, A, B);
        	/*if(get<0>(result)->is_empty())
        	{
        		cout << "empty: " << str_psl(result_formula, true) << endl;
        		exit(-1);
        	}*/
        	dfwa_pair pair(get<0>(result), get<1>(result), false, result_formula);
        	cout << "symbolic product finished, result has " << get<1>(result) << " states" << endl;
        	autlist.push(pair);
        	delete A;
        	delete B;
        }
        if(false)
        {
            dfwa_pair pair = autlist.top();
            twa_graph_ptr mona_aut = translate_ltlf_mona(pair._formula, dict);
            if(pair._is_explicit)
                {
                	string word = is_twa_equivalent(pair._twa, mona_aut);
            		if (word.size() == 0)
            		{
            			cout << "Equivalent" << endl;
            		}
                }else
                {

                	dfwa* mona_dfwa = symbolize_twa(dict, mona_aut);
            		{
            			cout << "L(A) <= L(B)" << endl;
            			dfwa_ptr pro = product_dfwa_minus(*pair._dfa, *mona_dfwa);
            			//bdd reach = pro.explore();
            			if (pro.is_empty()) {
            				cout << "A is subset of B" << endl;
            			} else {
            				cerr << "ERROR, not equivalent" << endl;
            				exit(-1);
            			}
            			//pro.free_variables();
            			delete &pro;
            			cout << "L(B) <= L(A)" << endl;
            			pro = product_dfwa_minus(*mona_dfwa, *pair._dfa);
            			if (pro.is_empty()) {
            				cout << "B is subset of A" << endl;
            			} else {
            				cerr << "ERROR, not equivalent" << endl;
            				//exit(-1);
            			}
            			//pro.free_variables();
            			delete &pro;
            		}
               }
        }

    }
    clock_t c_end = clock();
    cout << "Finished constructing minimal dfa in "
    		 << 1000.0 * (c_end - c_start)/CLOCKS_PER_SEC << "ms ..." << endl;
    dfwa_pair pair = autlist.top();
    cout << "Final result: " << pair._num_states << endl;
    if(! pair._is_explicit)
	{
		delete pair._dfa;
	}
    exit(0);
    // test equivalence
    twa_graph_ptr mona_aut = translate_ltlf_mona(input_f, dict);
    if(pair._is_explicit)
    {
    	string word = is_twa_equivalent(pair._twa, mona_aut);
		if (word.size() == 0)
		{
			cout << "Equivalent" << endl;
		}
    }else
    {
    	/*
    	dfwa* mona_dfwa = symbolize_twa(dict, mona_aut);
		{
			cout << "L(A) <= L(B)" << endl;
			dfwa_ptr pro = product_dfwa_minus(*pair._dfa, *mona_dfwa);
			//bdd reach = pro.explore();
			if (pro.is_empty()) {
				cout << "A is subset of B" << endl;
			} else {
				cerr << "ERROR" << endl;
				exit(-1);
			}
			//pro.free_variables();
			delete &pro;
			cout << "L(B) <= L(A)" << endl;
			pro = product_dfwa_minus(*mona_dfwa, *pair._dfa);
			if (pro.is_empty()) {
				cout << "B is subset of A" << endl;
			} else {
				cerr << "ERROR" << endl;
				exit(-1);
			}
			//pro.free_variables();
			delete &pro;
		}
		*/
    }

    set<unsigned> finals;
    mona_aut = minimize_dfa(mona_aut);
    cout << "minimized result : " << mona_aut->num_states() << endl;


}
