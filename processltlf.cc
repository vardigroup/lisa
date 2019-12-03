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

#include "spotutil.hh"
#include "mona.hh"

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
        spot::formula f = spot::parse_formula(line.c_str());
        std::cout << "parsed: " << f << '\n';
        auto pf1 = spot::parse_infix_psl(line.c_str());
        if (pf1.format_errors(std::cerr))
        {
            std::cerr << "error: " << line << std::endl;
            return -1;
        }
        //cout << "parsed formula: " << endl;
        cout << str_psl(pf1.f, true) << endl;
        // formula 1
        
      }
    }
    ltlfile.close();
    return 0;
}
    

