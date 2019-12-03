
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

#include "debug.hh"
#include "dfwa.hh"
#include "dfwavar.hh"
#include "mona.hh"
#include "spotutil.hh"
#include "synt.hh"

using namespace spot;
using namespace std;

static struct opt_t
{
  spot::bdd_dict_ptr dict = spot::make_bdd_dict();
}* opt;

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
    if(lst.size() < 1000)
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
    while(pq.size() > 1000)
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

bool 
compare_aut_size(twa_graph_ptr p1, twa_graph_ptr p2)
{
    if(p1->num_states() == p2->num_states())
    {
        return false;
    }
    return p1->num_states() > p2->num_states();
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

/*-------------------------------------------------------------------*/
// synthesize with an explicit DFA 
/*-------------------------------------------------------------------*/      
void
explict_synthesize(twa_graph_ptr aut, char* partfile)
{
    set<formula> set_of_props;
    for(auto &f : aut->ap())
    {
        set_of_props.insert(f);
    }
    bdd_dict_ptr dict = aut->get_dict();
    // there is alive states
    bdd label_cube = bddtrue;
    for(auto& f : set_of_props)
    {
        label_cube = label_cube & bdd_ithvar(dict->varnum(f));
    }

    // partition aps
    vector<string> input;
    vector<string> output;
    
    read_from_part_file(partfile, input, output);
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
}

/*-------------------------------------------------------------------*/
// synthesizing 
/*-------------------------------------------------------------------*/ 

void
explict_synthesize(bdd_dict_ptr dict
, twa_graph_ptr A
, twa_graph_ptr B
, vector<string>& input
, vector<string>& output)
{
    // the set of propositions appears in two DFAs
    set<formula> set_of_props;
    for(auto &f : A->ap())
    {
        set_of_props.insert(f);
    }
    for(auto &f : B->ap())
    {
        set_of_props.insert(f);
    }
    // there is alive states
    bdd label_cube = bddtrue;
    for(auto& f : set_of_props)
    {
        label_cube = label_cube & bdd_ithvar(dict->varnum(f));
    }
    set<unsigned> finals_A;
    get_final_states(A, finals_A);
    set<unsigned> finals_B;
    get_final_states(B, finals_B);
    // now construct automata
    dfwa s_A(A, label_cube, finals_A);
    dfwa s_B(B, label_cube, finals_B);
  
    //dfaa.output(cout);
  
    //dfab.output(cout);
  
    dfwa result(dict, label_cube);
    intersect_dfwa(result, s_A, s_B);
    //result.output(cout);
    cout << "Finished product" << endl;
    cout << " #path = " << bdd_pathcount(result._reach) << endl;
    
    // synthesis
    bdd input_cube = bddtrue;
    bdd output_cube = bddtrue;
    for(string & in : input)
    {
        bdd p = bdd_ithvar(A->register_ap(in));
        input_cube = input_cube & p;
    }
    for(string & out : output)
    {
        bdd p = bdd_ithvar(A->register_ap(out));
        output_cube = output_cube & p;
    }
    string alive_ap(ALIVE_AP);
    bdd p = bdd_ithvar(A->register_ap(alive_ap));
    output_cube = output_cube & p;
    
    //print_hoa(std::cout, aut);
     cout << "A: " << endl;
    for(spot::formula ap: A->ap())
    {
        cout << ap.ap_name() << endl;
    }
    cout << "B: " << endl;
    for(spot::formula ap: B->ap())
    {
        cout << ap.ap_name() << endl;
    }
    
    cout << "now construct synt " << endl;
      
    clock_t c_start = clock();
    auto t_start = chrono::high_resolution_clock::now();
    synt syn(result, input_cube, output_cube);
    cout << "Starting to synthesize " << endl;
    syn.is_realizable();
    cout << "Finished synthesizing" << endl;

    clock_t c_end = clock();
    cout << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms\n";
    auto t_end = chrono::high_resolution_clock::now();
    cout << "Total CPU time used: "
              << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms\n"
              << "Total wall clock time passed: "
              << std::chrono::duration<double, std::milli>(t_end-t_start).count()
              << " ms\n";
}

int 
main(int argc, char** argv)
{
    if(argc < 1 || argc > 3)
    {
        cout << "please formula file and part file" << endl;
        cout << "<usage>: " << endl;
        cout << "         <program> <ltlf_file> <part_file>" << endl;
        exit(-1);
    }
    
    // ltlf to a list of automata
    ifstream ltlfile(argv[1]);
    string line;
    vector<twa_graph_ptr> autlist;
    opt_t o;
    spot::bdd_dict_ptr dict = spot::make_bdd_dict();
    formula input_f;
    if (ltlfile.is_open()) 
    {
        getline (ltlfile, line);
        cout << "formula: " << line << endl;
        auto pf1 = spot::parse_infix_psl(line.c_str());
        if (pf1.format_errors(std::cerr))
        {
            cerr << "error: " << line << endl;
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
            twa_graph_ptr aut = trans_formula(f, dict, 7);
            cout << aut->num_states() << endl;
            autlist.push_back(aut);
        }
        ltlfile.close();
    }
    //return 1;
    cout << "Done splitting formulas" << endl;
    // do products 
    bool opt = true;
    
    // record the pointers of all minimized DFA
    set<twa_graph_ptr> optimized;
    while(autlist.size() > 1) 
    {
        cout << "loop starts: sorting " << autlist.size() << endl;
        for(int i = 0; i < autlist.size(); i ++)
        {
            cout << i << " st = " << autlist[i]->num_states() << endl;
        }
        sort (autlist.begin(), autlist.end(), compare_aut_size);
        cout << "sorted ..." << endl;
        twa_graph_ptr A = autlist.back();
        autlist.pop_back();
        twa_graph_ptr B = autlist.back();
        autlist.pop_back();
        cout << "poped two elements" << endl;
        /*
        cout << "A : ";
        for(spot::formula ap: A->ap())
        {
           cout << ap.ap_name() << endl;
           B->register_ap(ap.ap_name());
        }
        cout << "B : ";
        for(spot::formula ap: B->ap())
        {
           cout << ap.ap_name() << endl;
           A->register_ap(ap.ap_name());
        }
        */
        cout << "A: " << endl;
        for(spot::formula ap: A->ap())
        {
           cout << ap.ap_name() << endl;
        }
        cout << "B: " << endl;
        for(spot::formula ap: B->ap())
        {
           cout << ap.ap_name() << endl;
        }
        cout << "A#states=" << A->num_states() << " isdet= " << is_deterministic(A)
             << " B#states=" << B->num_states() << " isdet= " << is_deterministic(B) << endl;
        
        if(opt)
        {
            //spot::postprocessor post;
            //post.set_type(spot::postprocessor::BA);
            //post.set_pref(spot::postprocessor::Small); // or ::Deterministic
            //post.set_pref(spot::postprocessor::Low);
            //post.set_pref(spot::postprocessor::Deterministic);
            //A = post.run(A);  
            //B = post.run(B);
            //A = spot::minimize_wdba(A);
            twa_graph_ptr C;
            if(optimized.find(A) == optimized.end())
            {
               C = spot::minimize_wdba(A);
               optimized.insert(C);
            }else 
            {
                C = A;
                DEBUG_STDOUT( "No need to minimize A");
            }
            //A = spot::minimize_obligation(A); 
            // check equivalence of two automata
            #ifdef DEBUG
            string word = is_twa_equivalent(A, C);
            if(word.size() == 0)
            {
                cout << "A: equivalent two automata" << endl;
            }
            #endif
            A = C;
            if(optimized.find(B) == optimized.end())
            {
               C = spot::minimize_wdba(B);
               optimized.insert(C);
            }else 
            {
                C = B;
                DEBUG_STDOUT("No need to minimize B");
            }
            #ifdef DEBUG
            word = is_twa_equivalent(B, C);
            if(word.size() == 0)
            {
                cout << "B: equivalent two automata" << endl;
            }
            #endif
            B = C;
        }
        const unsigned int NUM = 3000;
        if(autlist.size() == 0)
        {
            // do symbolic computation for the last product
            // collect the propositions in formula
            set<formula> set_of_props;
            for(auto &f : A->ap())
            {
                set_of_props.insert(f);
            }
            for(auto &f : B->ap())
            {
                set_of_props.insert(f);
            }
            // there is alive states
            bdd label_cube = bddtrue;
            for(auto& f : set_of_props)
            {
                label_cube = label_cube & bdd_ithvar(dict->varnum(f));
            }
            set<unsigned> finals_A;
            //compute_final_states(A, finals_A);
            get_final_states(A, finals_A);
            set<unsigned> finals_B;
            //compute_final_states(B, finals_B);
            get_final_states(B, finals_B);
            // now construct automata
            dfwa s_A(A, label_cube, finals_A);
            dfwa s_B(B, label_cube, finals_B);
          
            //dfaa.output(cout);
          
            //dfab.output(cout);
          
            dfwa result(dict, label_cube);
            intersect_dfwa(result, s_A, s_B);
            //result.output(cout);
            cout << "Finished product" << endl;
            cout << "#states=" << bdd_satcount(result._reach) << " #path = " << bdd_pathcount(result._reach) << endl;
            
            // synthesis
            vector<string> input;
            vector<string> output;
            cout << "read part file " << endl;
            read_from_part_file(argv[2], input, output);
             cout << "finished reading part file " << endl;
            bdd input_cube = bddtrue;
            bdd output_cube = bddtrue;
            for(string & in : input)
            {
                bdd p = bdd_ithvar(A->register_ap(in));
                input_cube = input_cube & p;
            }
            for(string & out : output)
            {
                bdd p = bdd_ithvar(A->register_ap(out));
                output_cube = output_cube & p;
            }
            string alive_ap(ALIVE_AP);
            bdd p = bdd_ithvar(A->register_ap(alive_ap));
            output_cube = output_cube & p;
            cout << "now construct synt " << endl;
            //print_hoa(std::cout, aut);
             cout << "A: " << endl;
            for(spot::formula ap: A->ap())
            {
                cout << ap.ap_name() << endl;
            }
            cout << "B: " << endl;
            for(spot::formula ap: B->ap())
            {
                cout << ap.ap_name() << endl;
            }
            //return 0;
              //dfwa_var vars(dict, aut, 2, "s", 0, 6);
              
              //cout << vars.get_lower() << endl;
              //cout << vars.get_upper() << endl;
              
              clock_t c_start = clock();
              auto t_start = chrono::high_resolution_clock::now();
              synt syn(result, input_cube, output_cube);
              cout << "Starting to synthesize " << endl;
              syn.is_realizable();
              syn.synthesize();
              cout << "Finished synthesizing" << endl;

              clock_t c_end = clock();
              cout << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms\n";
                auto t_end = chrono::high_resolution_clock::now();
             cout << "Total CPU time used: "
                      << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms\n"
                      << "Total wall clock time passed: "
                      << std::chrono::duration<double, std::milli>(t_end-t_start).count()
                      << " ms\n";
            //synt syn();
            return 0;
          
        }
        A = ::product(std::move(A), B);
        if(opt )
        {
            //spot::sat_minimize(A, "")
            // minimize a weak dba 
            #ifdef DEBUG 
            cout << "optimize weak DBA" << endl;
            if(A->prop_weak() && is_deterministic(A))
            {
                cout << "weak DBA" << endl;
            }else
            {
                cout << "not weak DBA" << endl;
            }
            #endif
            twa_graph_ptr C = spot::minimize_wdba(A);
            optimized.insert(C);
            #ifdef DEBUG 
            //A = spot::minimize_obligation(A); 
            // check equivalence of two automata
            string word = is_twa_equivalent(A, C);
            if(word.size() == 0)
            {
                cout << "equivalent two automata" << endl;
            }
            #endif 
            A = C;
        }else
        {
            cout << "minimizing nbas ..." << endl;
            spot::postprocessor post;
            //post.set_type(spot::postprocessor::BA);
            post.set_pref(spot::postprocessor::Small); // or ::Deterministic
            //post.set_pref(spot::postprocessor::Low);
            //post.set_pref(spot::postprocessor::Deterministic);
            A = post.run(A);  
            cout << "done minimizing nbas ..." << endl;
            //B = post.run(B);
        }
        autlist.push_back(A);
        cout << " P#states=" << A->num_states() << endl;
    }
    
    twa_graph_ptr result = autlist.back();
    cout << "#states: "<< result->num_states() << " det= " << is_deterministic(result) << endl;
    explict_synthesize(result, argv[2]);
    //ofstream outfile("output.hoa");
    //print_hoa(outfile, result);
    
}
