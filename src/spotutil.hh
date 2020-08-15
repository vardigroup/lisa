/*-------------------------------------------------------------------*/
// Yong Li (liyong@ios.ac.cn)
/*-------------------------------------------------------------------*/

#pragma once

// Enable the strong Next operator for LTLf formula
#define SPOT_WANT_STRONG_X 1
#define SPOT_HAS_STRONG_X 1

#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <stack>
#include <set>

#include <spot/twaalgos/dualize.hh>
#include <spot/twaalgos/hoa.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/product.hh>
#include <spot/twaalgos/word.hh>
#include <spot/twaalgos/degen.hh>
#include <spot/twaalgos/translate.hh>
#include <spot/twa/twagraph.hh>
#include <spot/parseaut/public.hh>

#include <spot/tl/formula.hh>
#include <spot/tl/ltlf.hh>
#include <spot/tl/simplify.hh>

#include <spot/misc/optionmap.hh>
#include <spot/misc/hash.hh>
#include <chrono>

#include <bddx.h>


using namespace std;
using namespace spot;

const char* const ALIVE_AP = "alive";

string
is_twa_included(twa_graph_ptr aut_a, twa_graph_ptr aut_b);

string
is_twa_equivalent(twa_graph_ptr aut_a, twa_graph_ptr aut_b);

twa_graph_ptr
trans_formula(formula f, bdd_dict_ptr dict, unsigned num_ap_for_mona);

twa_graph_ptr
dfa_to_wdba(twa_graph_ptr aut, bdd_dict_ptr dict, set<unsigned> finals);

void
get_final_states(twa_graph_ptr aut, set<unsigned>& finals);

void
get_nonfinal_states(twa_graph_ptr A, set<unsigned>& nonfinals);

void
compute_final_states(twa_graph_ptr A, set<unsigned>& finals);

void
compute_accepting_states(twa_graph_ptr A, set<unsigned>& acc);

unsigned
get_size_formula_ap(formula& f);

void
get_formula_aps(formula& f, set<formula>& aps);

unsigned
get_size_formula(formula& f);

bool
is_even(unsigned value);

void
get_list_var_indices(vector<int>& vars, bdd cube);

void
output_bdd(ostream& os, bdd dd);

bdd
local_bdd_and(bdd& op1, bdd& op2);

bdd
local_bdd_or(bdd& op1, bdd& op2);

bdd
local_bdd_not_and(bdd& op1, bdd& op2);
