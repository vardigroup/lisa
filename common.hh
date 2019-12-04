#pragma once

#include <bddx.h>
#include <cudd.h>
#include <cuddObj.hh>

// CUDD and BuDDy both use BDD, we need new names for CUDD::BDD
typedef CUDD::BDD cudd_bdd;

typedef CUDD::BDD& cudd_bdd_ptr;

typedef CUDD::Cudd* cudd_ptr;

typedef CUDD::Cudd cudd;

typedef DdNode* cudd_node_ptr;

typedef bdd buddy_bdd;

typedef bdd& buddy_bdd_ptr;
