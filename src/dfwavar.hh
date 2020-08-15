/*-------------------------------------------------------------------*/
// Yong Li (liyong@ios.ac.cn)
/*-------------------------------------------------------------------*/

#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

#include <spot/twaalgos/hoa.hh>

#include <spot/twa/twagraph.hh>

#include <spot/tl/formula.hh>

#include <spot/misc/bddlt.hh>

using namespace std;
using namespace spot;

static unsigned index_for_vars = 0;

unsigned
get_num_bits(int value);

/**
 * Use dict to register anonymous variables as state variables.
 * remember to register the variables whenever use the state variables
 * and release the variables in destructor.
 * */
class dfwa_var
{
    private:
        void prepare_vars();
        const char* UNDERSCORE = "_";
        void order_vars();

    public:
        dfwa_var(bdd_dict_ptr dict, twa_graph_ptr aut, unsigned copies, string name, int lower, int upper)
        : _dict(dict), _aut(aut), _copies(copies), _name(name), _lower(lower), _upper(upper)
        {
        	assert(lower <= upper);
            prepare_vars();
        }

        unsigned _copies = 0;
        vector<vector<bdd>> _dd_vars;
        bdd_dict_ptr _dict;
        vector<vector<string>> _dd_names;
        //vector<unsigned> _dd_nums;
        string _name;
        int _lower;
        int _upper;
        twa_graph_ptr _aut;
        

        // for product
        dfwa_var(bdd_dict_ptr dict);
        
        ~dfwa_var();
        int get_lower()
        {
            return _lower;
        }
        int get_upper()
        {
            return _upper;
        }
        
        string get_name()
        {
            return _name;
        }
        
        unsigned get_copies()
        {
            return _copies;
        }
        
        unsigned get_var_num(int copy)
        {
        	return _dd_vars[copy].size();
        }

        bdd_dict_ptr get_dict()
        {
            return _dict;
        }
        
        bdd& get_var(int copy, int index)
        {
        	assert(index >= 0
        			&& copy < _copies
					&& _dd_vars[copy].size() > index);
        	return _dd_vars[copy][index];
        }

        int get_var_id(int copy, int index)
        {
            return bdd_var(get_var(copy, index));
        }

        void pop_back_vars();

        void add_bdd_vars(dfwa_var& vars);
        
        // when the ordering of the variables does not matter
        void add_ordered_bdd_vars(dfwa_var& vars);

        //void add_bdd_vars(vector<bdd>& vars_vec);

        vector<bdd>& get_bdd_vars(unsigned copy)
        {
            return _dd_vars[copy];
        }
        // only contains variables
        int to_integer(bdd dd);
        
        bdd new_value(int copy, int value);
        
        void new_value(vector<bdd>& bits, int copy, int value);
        
        bdd get_cube(int copy);
        
        bddPair* make_pair(unsigned fst_copy, unsigned snd_copy);
};

typedef dfwa_var& dfwa_var_ptr;
