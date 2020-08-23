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
#include "dfwa.hh"

using namespace spot;
using namespace std;

int main(int argc, char **argv)
{
    if (argc < 1 || argc > 3)
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
            cout << "formula: " << f << endl;
            string mona_file_name2 = "./ltlf2.mona";
            ofstream ofs2(mona_file_name2, ofstream::out);
            ofs2 << "#LTLf formula" << endl;
            ofs2 << "#" << str_psl(f, true) << endl;
            ofs2 << "# Backus normal form" << endl;
            formula r = get_bnf(f);
            // output in a parsable formula
            // the BNF form, and then convert it to fol formula
            ltlf_to_pfol(ofs2, r);
            ofs2.close();
            //cout << "formula: " << str_psl(f, true) << endl;
            // now call ltlf2fol
            //string mona_file_name = "./ltlf.mona";
            //string command = "./ltlf2fol NNF " + fol_file_name + " > " + mona_file_name;
            //system(command.c_str());
            string dfa_file_name = "./mona.dfa";
            string command = "mona -u -xw " + mona_file_name2 + " >" + dfa_file_name;
            system(command.c_str());
            // if this turns to be a bottleneck, we need pthread to read from pipe
            twa_graph_ptr dfa = read_from_mona_file(dfa_file_name.c_str(), dict);
            cout << "#Past DFA: " << dfa->num_states() << endl;
            //aut4->merge_edges();
            //ofstream outfile1("output4.hoa");
            //print_hoa(outfile1, aut4);
            //outfile1.close();
            // reverse
            //cout << "Output aut5" << endl;
            twa_graph_ptr nfa = reverse(dict, dfa);
            cout << "#Reverse NFA: " << nfa->num_states() << endl;
            //aut5->merge_edges();
            //ofstream outfile2("output5.hoa");
            //print_hoa(outfile2, aut5);
            //outfile2.close();
            // symbolic encoding
            nfa = reduce(nfa, 1);
            cout << "#Reduced FA: " << nfa->num_states() << endl;
            set<unsigned> nonfinals;
            get_nonfinal_states(nfa, nonfinals);
            map<formula, int> &var_map = nfa->get_dict()->var_map;
            map<formula, int>::const_iterator iter = var_map.begin();
            bdd labelcube = bddtrue;
            while (iter != var_map.end())
            {
                formula key = iter->first;
                int value = iter->second;
                //cout << key << " -> " << value << endl;
                labelcube = labelcube & bdd_ithvar(value);
                iter++;
            }
            clock_t start = clock();
            if(is_deterministic(nfa))
            {
                cout << "Deterministic" << endl;
                dfwa *ret = new dfwa(nfa, labelcube, nonfinals);
                delete ret;
            }else
            {
                cout << "Nondeterministic" << endl;
                dfwa *ret = new dfwa(nfa, labelcube, nonfinals, true);
                delete ret; 
            }
            
            clock_t end = clock();
            cout << "Total CPU time used for encoding NFA: " << 1000.0 * (end - start) / CLOCKS_PER_SEC << " ms\n";
            
        }
    }
    ltlfile.close();
    return 0;
}
