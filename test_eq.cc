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

#include <spot/twa/twagraph.hh>

#include <spot/tl/formula.hh>
#include <spot/tl/print.hh>
#include <spot/tl/parse.hh>
#include <spot/tl/print.hh>

#include <spot/misc/optionmap.hh>
#include <spot/misc/timer.hh>

#include "mona.hh"
#include "spotutil.hh"

using namespace spot;
using namespace std;

twa_graph_ptr
translate_ltlf_spot(formula f, bdd_dict_ptr dict)
{
    twa_graph_ptr A = trans_formula(f, dict);
    A = minimize_wdba(A);
    return A;
}

int
main(int argc, char** argv)
{
    if(argc < 1 || argc > 3)
        std::cout << "please input formula file" << std::endl;
    
    ifstream ltlfile(argv[1]);
    string line;
    bdd_dict_ptr dict = make_bdd_dict();

    // ltlfile.open(argv[1]);

    cout << "Opened file" << endl;
    
    if (ltlfile.is_open())
    {
      cout << "We are inside the file" << endl;
      int i = 0;
      while (getline(ltlfile, line))
      {
        cout << "formula: " << line << endl;

        auto pf1 = spot::parse_infix_psl(line.c_str());
        if (pf1.format_errors(std::cerr))
        {
            std::cerr << "error: " << line << std::endl;
            return -1;
        }
        
        twa_graph_ptr A = translate_ltlf_spot(pf1.f, dict);
        
        twa_graph_ptr B = translate_ltlf_mona(pf1.f, dict);
        
        cout << "#A = " << A->num_states() << " isdet= "<< is_deterministic(A) << " #B = " << B->num_states() << " isdet= " << is_deterministic(B) << endl;
        string word = is_twa_equivalent(A, B);
        if(word.size() == 0)
        {
            cout << "Equivalent" << endl;
            unsigned init = B->get_init_state_number();
            bdd tr = bddfalse;
            cout << "curr=" << init << endl;
            for (auto& e: B->out(init))
            {
                // cout 
                cout << "succ=" << e.dst << endl;
                tr |= e.cond;
            }
            string ap("alive");
            bdd ap_dd = bdd_ithvar(B->register_ap(ap));
            if(tr == ap_dd)
            {
                cout << "universal tr" << endl;
            }
            ofstream outfile("output.hoa");
            print_hoa(outfile, B);
        }else
        {
            cerr << "Inquivalent formula " << pf1.f << endl;
            string name = "A" + to_string(i) + ".hoa";
            ofstream afile(name.c_str());
            print_hoa(afile, A);
            afile.close();
            name = "B" + to_string(i) + ".hoa";
            ofstream bfile(name.c_str());
            print_hoa(bfile, B);
            bfile.close();
            return -1;
        }
        ++ i;
      }
    }
    ltlfile.close();
    return 0;
}
    

