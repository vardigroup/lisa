#ifdef DEBUG

#include "dfwavar.hh"

/*-------------------------------------------------------------------*/
// get the number of bits needed to represent an integer value
/*-------------------------------------------------------------------*/
unsigned
get_num_bits(int value)
{
    if(value <= 0)
    {
        return 0;
    }
    if(value == 1)
    {
        return 1;
    }
    -- value;
    unsigned count = 0;
    while (value) 
    {
        cout << " val = " << value << endl;
        value = value >> 1;
        ++ count;
    }
    return count;
}
/*-------------------------------------------------------------------*/
// new empty bdd variables
/*-------------------------------------------------------------------*/
dfwa_var::dfwa_var(bdd_dict_ptr dict)
: _dict(dict)
{
    _aut = nullptr;
    _copies = 0;
}

dfwa_var::~dfwa_var()
{
	_dict->unregister_all_my_variables(this);
}
/*-------------------------------------------------------------------*/
// new bdd variables
/*-------------------------------------------------------------------*/

void
dfwa_var::prepare_vars()
{
    // numbers
    int num_values = get_upper() - get_lower() + 1;
    unsigned num_bits = get_num_bits(num_values);
    //cout << "needs " << num_bits << " bits" << endl;
    for (unsigned copy = 0; copy < _copies; copy++) 
    {
        vector<bdd> vec;
        _dd_vars.push_back(vec);
        vector<string> vec_str;
        _dd_names.push_back(vec_str);
    }

    //int varNum = bdd.getNumVars();
    // n consecutive BDD variables which will be used only by for_me.
    unsigned var_num = _dict->register_anonymous_variables(_copies * num_bits, this);
    //unsigned var_num =_dict->register_anonymous_variables(_copies * num_bits, this);
    for (unsigned bit = 0; bit < num_bits; bit++)
    {
        //_dd_nums.push_back(var_num);
        for (unsigned copy = 0; copy < _copies; copy++) 
        {
            bdd dd = bdd_ithvar(var_num);
            string dd_name = _name + UNDERSCORE + to_string(bit) + UNDERSCORE + to_string(copy);
            _dd_vars[copy].push_back(dd);
            _dd_names[copy].push_back(dd_name);
            //cout << "dd_name = " << dd_name << " var_num = " << var_num << endl;
            ++ var_num;
        }
    }
}
/*-------------------------------------------------------------------*/
// a bdd to an integer
/*-------------------------------------------------------------------*/
int 
dfwa_var::to_integer(bdd dd)
{
    int value = 0;
    int bit = 1;
    for (bdd bit_var : get_bdd_vars(0)) 
    {
        bdd result = bit_var & dd;
        if (result != bddfalse) 
        {
            value |= bit;
        }
        bit <<= 1;
    }
    value += _lower;
    return value;
}

/*-------------------------------------------------------------------*/
// an integer to a bdd
/*-------------------------------------------------------------------*/
bdd
dfwa_var::new_value(int copy, int value)
{
    value -= get_lower();
    bdd dd = bddtrue;
    int bit = 1;
    for (bdd bit_var : get_bdd_vars(copy)) 
    {
        bdd bit_var_not = !bit_var;
        dd = dd & ((value & bit) != 0 ? bit_var : bit_var_not);
        bit <<= 1;
    }
    return dd;
}

void
dfwa_var::new_value(vector<bdd>& bits, int copy, int value)
{
    value -= get_lower();
    int bit = 1;
    for (bdd bit_var : get_bdd_vars(copy)) 
    {
        if((value & bit) != 0)
        {
            bits.push_back(bddtrue);
        }else
        {
            bits.push_back(bddfalse);
        }
        bit <<= 1;
    }
}

/*-------------------------------------------------------------------*/
// get the list of variables of copy
/*-------------------------------------------------------------------*/
bdd
dfwa_var::get_cube(int copy)
{
    bdd cube = bddtrue;
    for (bdd dd : get_bdd_vars(copy)) 
    {
        bdd temp = cube & dd;
        cube = temp;
    }
    return cube;
}

/*-------------------------------------------------------------------*/
// add more variables
/*-------------------------------------------------------------------*/
void
dfwa_var::add_bdd_vars(dfwa_var& vars)
{
    unsigned copies = max(_copies, vars.get_copies());
    cout << "copies: " << copies << endl;
    for(unsigned i = 0; i < copies; i ++)
    {
        if(i < _dd_vars.size())
        {
            // add variables to the copy i vector
            vector<bdd>& temp = _dd_vars[i];
            vector<bdd>& add = vars.get_bdd_vars(i);
            for(unsigned j = 0; j < add.size(); j ++)
            {
                temp.push_back(add[j]);
            }
        }else
        {
            //create a vector and add variables
        	cout << "add state variables" << endl;
            vector<bdd> temp;
            vector<bdd>& add = vars.get_bdd_vars(i);
            for(unsigned j = 0; j < add.size(); j ++)
            {
                temp.push_back(add[j]);
            }
            _dd_vars.push_back(temp);
        }
    }
    // make sure the order is increasing

}

/*-------------------------------------------------------------------*/
// make bdd pair by variable indices
/*-------------------------------------------------------------------*/
bddPair*
dfwa_var::make_pair(unsigned fst_copy, unsigned snd_copy)
{
    bddPair* pair = bdd_newpair();
    /*
    for(unsigned i = 0; i < _dd_nums.size(); i ++)
    {
        bdd_setpair(pair, (int)(_dd_nums[i] + fst_copy), (int)(_dd_nums[i] + snd_copy));
        cout << "fst = " << (_dd_nums[i] + fst_copy) << " snd = " <<  (_dd_nums[i] + snd_copy) << endl;
    }
    */
    vector<bdd> curr = get_bdd_vars(fst_copy);
    vector<bdd> next = get_bdd_vars(snd_copy);
    for(unsigned i = 0; i < curr.size(); i ++)
    {
        bdd_setpair(pair, bdd_var(curr[i]), bdd_var(next[i]));
        //cout << "fst = " << bdd_var(curr[i]) << " snd = " <<  bdd_var(next[i]) << endl;
    }
    return pair;
}

#endif
