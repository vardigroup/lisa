/*-------------------------------------------------------------------*/
// Yong Li (liyong@ios.ac.cn)
/*-------------------------------------------------------------------*/

#pragma once


#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/minimize.hh>
#include <spot/twa/twagraph.hh>
#include <spot/twaalgos/postproc.hh>

#include "spotutil.hh"

/*-------------------------------------------------------------------*/
// Project out all the variables specified in @cube the automaton @aut
// if @min is true, determinize and minimize the resulting automaton 
/*-------------------------------------------------------------------*/

twa_graph_ptr
project(spot::bdd_dict_ptr dict, twa_graph_ptr aut, bdd cube, bool min);


/*-------------------------------------------------------------------*/
// Reverse a finite automaton @aut, i.e.,
// 1. (s, a, t) -> (t, a, s)
// 2. s is accepting then make it nonaccepting
/*-------------------------------------------------------------------*/
twa_graph_ptr
reverse(spot::bdd_dict_ptr dict, twa_graph_ptr aut);


/*-------------------------------------------------------------------*/
// Reduce a finite automaton @aut, i.e.,
// @level is the amount of effort for reduction
// 0. low level
// 1. medium level
// 2. high level
/*-------------------------------------------------------------------*/
twa_graph_ptr
reduce(twa_graph_ptr aut, int level);
