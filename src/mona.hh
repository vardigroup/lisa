/*-------------------------------------------------------------------*/
// Yong Li (liyong@ios.ac.cn)
/*-------------------------------------------------------------------*/

#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <tuple>
#include <sstream>
#include <unordered_map>

#include <spot/twaalgos/hoa.hh>
#include <spot/twa/twagraph.hh>
#include <spot/tl/formula.hh>
#include <spot/tl/parse.hh>
#include <spot/tl/print.hh>

#include "ltlf2fol.hh"

using namespace spot;
using namespace std;

bool
str_contain(string str, const char* match);

void
str_split(string str, vector<string>& result, char delim = ' ');
// translate ltlf formula to DFA via mona
twa_graph_ptr
translate_ltlf_mona(formula ltl, bdd_dict_ptr dict);

twa_graph_ptr
read_from_mona_file(const char * file_name, bdd_dict_ptr dict);
