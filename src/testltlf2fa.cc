/*-------------------------------------------------------------------*/
// Yong Li (liyong@ios.ac.cn)
/*-------------------------------------------------------------------*/

// standard 
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <unordered_map>
#include <string>
#include <vector>

// spot 
#include <spot/twaalgos/hoa.hh>
#include <spot/twaalgos/translate.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/product.hh>
#include <spot/twaalgos/word.hh>
#include <spot/twaalgos/minimize.hh>
#include <spot/twaalgos/contains.hh>

#include <spot/twa/twagraph.hh>

#include <spot/tl/formula.hh>
#include <spot/tl/print.hh>
#include <spot/tl/parse.hh>
#include <spot/tl/print.hh>
#include <spot/tl/ltlf.hh>

#include <spot/misc/optionmap.hh>
#include <spot/misc/timer.hh>

#include "spotutil.hh"
#include "mona.hh"
#include "ltlf2fol.hh"
#include "aututil.hh"


using namespace spot;
using namespace std;

int
main(int argc, char** argv)
{
    if(argc < 1 || argc > 3)
        std::cout << "please input formula file" << std::endl;
    
    ifstream ltlfile(argv[1]);
    string line;
    bdd_dict_ptr dict = make_bdd_dict();
    // ltlfile.open(argv[1]);

    //cout << "Opened file" << endl;
    
    if (ltlfile.is_open())
    {
      //cout << "We are inside the file" << endl;
      while (getline(ltlfile, line))
      {
        //cout << "formula: " << line << endl;
        auto pf1 = spot::parse_infix_psl(line.c_str());
        if (pf1.format_errors(std::cerr))
        {
            std::cerr << "error: " << line << std::endl;
            return -1;
        }
        //cout << "parsed formula: " << endl;
        formula f = pf1.f;
        cout << str_psl(pf1.f, true) << endl;
        formula f1 = from_ltlf(f, "alive");
        // formula 1
        formula r = get_bnf(f);
        formula r1 = from_ltlf(r, "alive");
        //spot::formula f = spot::parse_formula("(a U b) U a");
        //spot::formula g = spot::parse_formula("b U a");
        std::cout << " BNF and Original LTLf formulas:"<< (spot::are_equivalent(f1, r1) ?
                "Equivalent\n" : "Not equivalent\n");
                
        formula rr = get_nnf(f);
        formula rr1 = from_ltlf(r, "alive");
        std::cout << " NNF and Original LTLf formulas: " << (spot::are_equivalent(f1, rr1) ?
                "Equivalent\n" : "Not equivalent\n");
        // automata construction
        option_map m;
        translator trans(dict, &m);
        formula ltl = from_ltlf(f, "alive");
        //spot::tl_simplifier simpl;
        //ltl = simpl.simplify(ltl);
        twa_graph_ptr aut1 = trans.run(ltl);
        // construct from MONA
        string mona_file_name = "./ltlf.mona";
        ofstream ofs (mona_file_name, ofstream::out);
        ofs << "#LTLf formula" << endl;
        ofs << "#" << str_psl(f, true) << endl;
        ofs << "# Backus normal form" << endl;
        // output in a parsable formula
        // the BNF form, and then convert it to fol formula
        ltlf_to_fol(ofs, r);
        ofs.close();
        //cout << "formula: " << str_psl(f, true) << endl;
        // now call ltlf2fol
        //string mona_file_name = "./ltlf.mona";
        //string command = "./ltlf2fol NNF " + fol_file_name + " > " + mona_file_name;
        //system(command.c_str());
        string dfa_file_name = "./mona.dfa";
        string command = "mona -u -xw " + mona_file_name+ " >" + dfa_file_name;
        system(command.c_str());
    // if this turns to be a bottleneck, we need pthread to read from pipe
        twa_graph_ptr aut2 = read_from_mona_file(dfa_file_name.c_str(), dict);
        cout << "#Spot DFA: " << aut1->num_states() << endl;
        cout << "#Mona DFA: " << aut2->num_states() << endl;
        // check equivalence
        string word = is_twa_equivalent(aut1, aut2);
        cout << "Spot and MONA DFA:";
        if(word.size() == 0)
        {
            cout << "Equivalent" << endl;
        }else
        {
            cout << "CE: " << word << endl;
            cout << "Inequivalent" << endl;
            exit(-1);
        }
        // another automaton
        //rans_prefixltlf2fol(ostream &os, formula &f)
        //string mona_file_name1 = "./ltlf1.mona";
        //ofstream ofs1(mona_file_name1, ofstream::out);
        //ofs1 << "#LTLf formula" << endl;
        //ofs1 << "#" << str_psl(f, true) << endl;
        //ofs1 << "# Backus normal form" << endl;
        // output in a parsable formula
        // the BNF form, and then convert it to fol formula
        //trans_prefixltlf2fol(ofs1, r);
        //ofs1.close();
        //cout << "formula: " << str_psl(f, true) << endl;
        // now call ltlf2fol
        //string mona_file_name = "./ltlf.mona";
        //string command = "./ltlf2fol NNF " + fol_file_name + " > " + mona_file_name;
        //system(command.c_str());
        //string dfa_file_name = "./mona.dfa";
        //command = "mona -u -xw " + mona_file_name1 + " >" + dfa_file_name;
        //system(command.c_str());
    // if this turns to be a bottleneck, we need pthread to read from pipe
        //twa_graph_ptr aut3 = read_from_mona_file(dfa_file_name.c_str(), dict);
        //cout << "Aut3: " << aut3->num_states() << endl;
        //ofstream outfile("outputaut.hoa");
		//print_hoa(outfile, aut3);
        //bool isIn = contains(aut3, aut2);
        /*
        if(isIn)
        {
            cout << "Contains" << endl;
        }else
        {
            cout << "Not contained" << endl;
        }
        isIn = contains(aut2, aut3);
        if(isIn)
        {
            cout << "Contains" << endl;
        }else
        {
            cout << "Not contained" << endl;
        }
        */
        string mona_file_name2 = "./ltlf2.mona";
        ofstream ofs2(mona_file_name2, ofstream::out);
        ofs2 << "#LTLf formula" << endl;
        ofs2 << "#" << str_psl(f, true) << endl;
        ofs2 << "# Backus normal form" << endl;
        // output in a parsable formula
        // the BNF form, and then convert it to fol formula
        ltlf_to_pfol(ofs2, r);
        ofs2.close();
        //cout << "formula: " << str_psl(f, true) << endl;
        // now call ltlf2fol
        //string mona_file_name = "./ltlf.mona";
        //string command = "./ltlf2fol NNF " + fol_file_name + " > " + mona_file_name;
        //system(command.c_str());
        //string dfa_file_name = "./mona.dfa";
        command = "mona -u -xw " + mona_file_name2 + " >" + dfa_file_name;
        system(command.c_str());
    // if this turns to be a bottleneck, we need pthread to read from pipe
        twa_graph_ptr aut4 = read_from_mona_file(dfa_file_name.c_str(), dict);
        cout << "#Past DFA: " << aut4->num_states() << endl;
        //aut4->merge_edges();
        //ofstream outfile1("output4.hoa");
		//print_hoa(outfile1, aut4);
        //outfile1.close();
        // reverse
        //cout << "Output aut5" << endl;
        twa_graph_ptr aut5 = reverse(dict, aut4);
        cout << "#Reverse NFA: " << aut5->num_states() << endl;
        //aut5->merge_edges();
        //ofstream outfile2("output5.hoa");
		//print_hoa(outfile2, aut5);
        //outfile2.close();
        word = is_twa_equivalent(aut1, aut5);
        cout << "Spot DFA and Reverse NFA: ";
        if(word.size() == 0)
        {
            cout << "Equivalent" << endl;
        }else
        {
            cout << "CE: " << word << endl;
            cout << "Inequivalent" << endl;
            exit(-1);
        }
      }
    }
    ltlfile.close();
    return 0;
}
    

