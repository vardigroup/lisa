#pragma once
//#define TRACE
#ifdef TRACE
#  define trace std::cerr
#else
#  define trace while (0) std::cerr
#endif

#include <queue>
#include <deque>
#include <set>
#include <list>
#include <vector>
#include <sstream>
#include <spot/misc/hash.hh>
#include <spot/misc/bddlt.hh>
#include <spot/twaalgos/product.hh>
#include <spot/twaalgos/powerset.hh>
#include <spot/twaalgos/gtec/gtec.hh>
#include <spot/twaalgos/strength.hh>
#include <spot/twaalgos/sccfilter.hh>
#include <spot/twaalgos/sccinfo.hh>
#include <spot/twaalgos/ltl2tgba_fm.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/dualize.hh>
#include <spot/twaalgos/remfin.hh>
#include <spot/tl/hierarchy.hh>

using namespace spot;
using namespace std;

// minimize twa by Hopcroft's algorithm
twa_graph_ptr
minimize_dfa(twa_graph_ptr aut);
