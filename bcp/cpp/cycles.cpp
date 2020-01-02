
class Cycle {
  Vector< Vector< int > > as; // sparse patterns, index j; s
  Vector< Vector< int > > ps; // sparse items, index i; t
   // as[j][s] is the s-th item in pattern j
   // ps[i][t] is the t-th pattern containing item i
   // ps[i] is sorted REVERSELY for each i
  Vector< Vector< int > > a; // dense patterns, index j; i
/// Cycle:
  Vector<int> s, p; // as[p[k]][s[k]] is the item index
    // currently chosen in pattern p[k],
    // i.e. patterns are indexed dense and items sparse.
    // In each level k, before searching for p[k],
    // it is set to -1 (i.e. -1 is pushed), the same for s[].

  void Preprocess(); // Eliminate patterns which are a subset of another????
};

bool Cycle::Find() {
  s.reserve(m); p.reserve(m); // no more than 2m nodes
  s.resize(1);  p.resize(1);
  for (p[0]=0;p[0]<n;++p[0]) { // for all patterns after p (incl.)
    s[0] = -1; // new i-search
    while (FindNextCyclableItem(p,s)) {
        // except the last such item in p and except those marked
      Labyrinth(p,s);
      if (p.size() > 1) { // a cycle found
        assert(p.size() == s.size()); // ???
	return 1;
      }
      else assert(p.size()==1 and s.size()==1);
    }
  }
  return 0;
}

/// p[0] and i[0] are the start:
void Cycle::Labyrinth(Vector<int> &p, Vector<int> &s) {
  p.push_back(-1); // init pattern search
NextPat:
  if (not FindNextPat(p,s)) {
    // so that p contains i[k-1] and some items not in p[k-1]
    if (1 == p.size()) return; // search completed
    // else go find Next Item at the prev. level
  }
  else
    if (NewPatternCanClose(p,s)) // pushes this item into s
      return;
NextItem:
  FindNextItem(p,s) // item not cont. in p[k-1] and not marked
  goto NextPat; // also marking outgoing items
}

void Cycle::FindNextItem(Vector<int> &p, Vector<int> &s) {
  assert(p.size() == s.size()); // invariant here
NextItem:
  if (++s.back() == as[p.back()].size()) { // no more items in p[k]
    s.pop_back(); // remove the place in stack
    return;
  }
  if (Mark[p.back()][s.back()])
    goto NextItem;
  if (ps[as[p.back()][s.back()]].size() < 2)
    // no 2 patterns with item as[p.back()][s.back()]
    goto NextItem;
  if (a[p[p.size()-2]][as[p.back()][s.back()]])
    goto NextItem;
  Mark[p.back()][s.back()] = 1;
 /// Marking items from the ONE-BEFORE-LAST pattern (cumulative for all p.):
  MarkItemsOfTheOneBeforeLastPattern(); // care of position
 // - this marking should be canceled when level up (see FindNextPat())
  p.push_back(-1); // init new pattern search
}

bool Cycle::FindNextPat(Vector<int> &p, Vector<int> &s) {
  assert(p.size() == s.size()+1); // invariant here
  if (not ...) {
    p.pop_back();
    UnmarkItemsOfTheOneBeforeLastPattern(); // care of position
    return;
  }
  s.push_back(-1);
  return 1;
}

/// AFTER CHOOSING A PATTERN, SHOULD WE FIRST LOOK IF IT CONTAINS ANY ITEMS FROM EARLIER PATTERNS??

