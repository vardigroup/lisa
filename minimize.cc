# include "minimize.hh"

// This is called hash_set for historical reason, but we need the
// order inside hash_set to be deterministic.
typedef std::set<unsigned> hash_set;
static std::ostream&
dump_hash_set(const hash_set* hs,
              std::ostream& out)
{
  out << '{';
  const char* sep = "";
  for (auto i: *hs)
    {
      out << sep << i;
      sep = ", ";
    }
  out << '}';
  return out;
}

static std::string
format_hash_set(const hash_set* hs)
{
  std::ostringstream s;
  dump_hash_set(hs, s);
  return s.str();
}

// Find all states of an automaton.
static void
build_state_set(const const_twa_graph_ptr& a, hash_set* fin, hash_set* seen)
{
  std::stack<unsigned> todo;
  unsigned init = a->get_init_state_number();
  todo.push(init);
  seen->insert(init);
  while (!todo.empty())
    {
      unsigned s = todo.top();
      todo.pop();
      // transitions of s 
      for (auto& e: a->out(s))
      {
        if (seen->insert(e.dst).second)
          todo.push(e.dst);
        if(e.acc.has(0) && e.dst != s)
        {
            //cout << "final state: " << endl;
            fin->insert(s);
        }
      }
    }
}

// From the base automaton and the list of sets, build the minimal
// resulting automaton
static twa_graph_ptr
build_result(const const_twa_graph_ptr& a,
             std::list<hash_set*>& sets,
             hash_set* final)
{
    //cout << " quotient computing..." << endl;
    // get dict
    auto dict = a->get_dict();
    auto res = make_twa_graph(dict);
    // copy the ap of a
    res->copy_ap_of(a);
    //res->prop_state_acc(true);
    
    // For each set, create a state in the output automaton.  For an
    // input state s, state_num[s] is the corresponding the state in
    // the output automaton.
    std::vector<unsigned> state_num(a->num_states(), -1U);
    { 
        unsigned num = res->new_states(sets.size());
        //cout << "first state = " << num << endl;
        for (hash_set* h: sets)
        {
            //cout << "num = " << num << " " << format_hash_set(h) << endl;
            for (unsigned s: *h)
            {
                //cout << " state s= " << s << endl;
                state_num[s] = num;
            }
            num = num + 1;
        }
    }
    /*
    for(int i = 0; i < state_num.size(); i ++)
    {
        cout << i << " -> " << state_num[i] << endl;
    }
    
    cout << "final: ";
    for (set<unsigned>::iterator it = final->begin(); it!= final->end(); ++it)
    {
        cout << *it << " ";
    }
    cout << endl;
	*/
    if (!final->empty())
        res->set_buchi();

    // For each transition in the initial automaton, add the
    // corresponding transition in res.
  for (hash_set* h: sets)
    {
      // Pick one state.
      //cout << format_hash_set(h) << endl;
      unsigned src = *h->begin();
      unsigned src_num = state_num[src];
      bool accepting = (final->find(src) != final->end());

      // Connect it to all destinations.
      for (auto& e: a->out(src))
        {
          unsigned dn = state_num[e.dst];
          if ((int)dn < 0)  // Ignore useless destinations.
            continue;
          res->new_acc_edge(src_num, dn, e.cond, accepting);
        }
    }
    //cout << "done.. -> " << res->num_states() << endl;
  res->merge_edges();
  if (res->num_states() > 0)
    res->set_init_state(state_num[a->get_init_state_number()]);
  else
    res->set_init_state(res->new_state());
  return res;
}

// main function to compute the partition
// copied the function from spot
static twa_graph_ptr minimize_dfa_inner(const const_twa_graph_ptr& det_a,
                                  hash_set* final, hash_set* non_final)
{
  typedef std::list<hash_set*> partition_t;
  // current partition
  partition_t cur_run;
  // next partition
  partition_t next_run;

  // The list of equivalent states.
  partition_t done;

  // state to set map (initial capacity = num_states())
  std::vector<unsigned> state_set_map(det_a->num_states(), -1U);

  // Size of det_a
  unsigned size = final->size() + non_final->size();
  // Use bdd variables to number sets.  set_num is the first variable
  // available.
  // n consecutive BDD variables which will be used only by for_me.
  unsigned set_num =
    det_a->get_dict()->register_anonymous_variables(size, det_a);
  // variables which should be freed later
  std::set<int> free_var;
  for (unsigned i = set_num; i < set_num + size; ++i)
  {
    free_var.insert(i);
  }
  //
  std::map<int, int> used_var;

  hash_set* final_copy;
  // final is not empty
  if (!final->empty())
    {
      unsigned s = final->size();
      // first used variable
      used_var[set_num] = s;
      free_var.erase(set_num);
      if (s > 1)
      {
        // one partition
        cur_run.emplace_back(final);
      }
      else
      {
        // s <= 1
        done.emplace_back(final);
      }
      // state to set_num (index)
      for (auto i: *final)
      {
        state_set_map[i] = set_num;
      }
      // copy
      final_copy = new hash_set(*final);
    }
  else
    {
      // empty set
      final_copy = final;
    }

  if (!non_final->empty())
    {
      unsigned s = non_final->size();
      // nonfinal states -> set index num
      unsigned num = set_num + 1;
      used_var[num] = s;
      free_var.erase(num);
      if (s > 1)
      {
          // larger than 1
          cur_run.emplace_back(non_final);
      }
      else
      { 
          // size == 1 -> cannot split any more
          done.emplace_back(non_final);
      }
      // state -> set index
      for (auto i: *non_final)
      {
          state_set_map[i] = num;
      }
    }
  else
    {
      delete non_final;
    }

  // A bdd_states_map is a list of formulae (in a BDD form)
  // associated with a destination set of states.
  typedef std::map<bdd, hash_set*, bdd_less_than> bdd_states_map;

  bool did_split = true;

  while (did_split)
    {
        // record whether this run splits a set
      did_split = false;
      while (!cur_run.empty())
        {
          // Get a set to process.
          // cur is splitter ? is the set to split?
          hash_set* cur = cur_run.front();
          cur_run.pop_front();

          trace << "processing " << format_hash_set(cur)
                << std::endl;

          bdd_states_map bdd_map;
          // traverse the states in cur set
          // cur is the set of states to split
          for (unsigned src: *cur)
            {
              bdd f = bddfalse;
              // check the output of states in curr
              for (auto si: det_a->out(src))
                {
                  unsigned i = state_set_map[si.dst];
                  if ((int)i < 0)
                  {
                    // The destination state is not in our
                    // partition.  This can happen if the initial
                    // FINAL and NON_FINAL supplied to the algorithm
                    // do not cover the whole automaton (because we
                    // want to ignore some useless states).  Simply
                    // ignore these states here.
                    continue;
                  }
                  f |= (bdd_ithvar(i) & si.cond);
                }
               // construct cond & set of succs bdd representation

              // Have we already seen this formula ?
              bdd_states_map::iterator bsi = bdd_map.find(f);
              if (bsi == bdd_map.end())
                {
                  // No, create a new set.
                  hash_set* new_set = new hash_set;
                  new_set->insert(src);
                  // cond & succ -> src
                  bdd_map[f] = new_set;
                }
              else
                {
                  // Yes, add the current state to the set.
                  // add source state
                  bsi->second->insert(src);
                }
            }
          // cond & succ -> source map
          // for cur set (should be split)
          auto bsi = bdd_map.begin();
          if (bdd_map.size() == 1)
            {
              // The set was not split.
              // should not be split
              trace << "set " << format_hash_set(bsi->second)
                    << " was not split" << std::endl;
              next_run.emplace_back(bsi->second);
            }
          else
            {
                // needs to split cur
              did_split = true;
              // traverse bdd_map and split src?
              for (; bsi != bdd_map.end(); ++bsi)
                {
                  // states in the same partition
                  hash_set* set = bsi->second;
                  // Free the number associated to these states.
                  unsigned num = state_set_map[*set->begin()];
                  assert(used_var.find(num) != used_var.end());
                  unsigned left = (used_var[num] -= set->size());
                  // Make sure LEFT does not become negative (hence bigger
                  // than SIZE when read as unsigned)
                  assert(left < size);
                  if (left == 0)
                    {
                      used_var.erase(num);
                      free_var.insert(num);
                    }
                  // Pick a free number
                  assert(!free_var.empty());
                  num = *free_var.begin();
                  free_var.erase(free_var.begin());
                  used_var[num] = set->size();
                  // all states should be mapped to the same number
                  for (unsigned s: *set)
                  {
                      state_set_map[s] = num;
                  }
                  // Trivial sets can't be split any further.
                  if (set->size() == 1)
                    {
                      trace << "set " << format_hash_set(set)
                            << " is minimal" << std::endl;
                      done.emplace_back(set);
                    }
                  else
                    {
                      trace << "set " << format_hash_set(set)
                            << " should be processed further" << std::endl;
                      next_run.emplace_back(set);
                    }
                }
            }
          delete cur;
        }
      if (did_split)
        trace << "splitting did occur during this pass." << std::endl;
      else
        trace << "splitting did not occur during this pass." << std::endl;
      std::swap(cur_run, next_run);
    }
    // list of sets, insert sets from cur_run to done
  done.splice(done.end(), cur_run);

#ifdef TRACE
  cout << "Final partition: ";
  for (hash_set* hs: done)
  {
    cout << format_hash_set(hs) << ' ';
  }
  trace << std::endl;
#endif

  // Build the result.
  auto res = build_result(det_a, done, final_copy);

  // Free all the allocated memory.
  delete final_copy;

  for (hash_set* hs: done)
    delete hs;

  return res;
}

twa_graph_ptr
minimize_dfa(twa_graph_ptr aut)
{
    hash_set* final = new hash_set;
    hash_set* seen = new hash_set;
    hash_set* non_final = new hash_set;
    
    // compute final and non_final states
    build_state_set(aut, final, seen);
    //cout << "nonfinal: ";
    for (set<unsigned>::iterator it = seen->begin(); it!= seen->end(); ++it)
    {
        unsigned s = *it;
        if(final->find(s) == final->end())
        {
            //cout << s << " ";
            non_final->insert(s);
        }
    }
    //cout << endl;
    delete seen;
    return minimize_dfa_inner(aut, final, non_final);
}
