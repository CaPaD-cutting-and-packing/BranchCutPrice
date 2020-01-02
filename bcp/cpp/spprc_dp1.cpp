// spprc_dp1.cpp
// G. Belov, TU Dresden, 2005

#include "stdafx.h"
#include "spprc_dp1.h"
#include "lasthdr.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

namespace spprc_dp1
{

/**
 Dynamic programming for SPPRC,
 Shortest Path Problem with Resource Constraints
 with upper bounds on $x_{ii}$ ($x_{ij}$ are binary)
 and a single resource.
 Only a few arcs make the problem different from knapsack;
 exactly $k$ arcs, given in d1.

 The table filling procedure follows the "Improved DP"
 from
    U. Pferschy. Dynamic Programming revisited: improving knapsack algorithms.
    Computing, 63:419--430, 1999.
 also discussed in the book
    H. Kellerer, U. Pferschy, D. Pisinger. Knapsack Problems. Springer-Verlag, 2004.

 The overall COMPLEXITY is O(kmL)

 IMPLEMENTATION
 No raster points, particularly no reduced RP.
 Thus, all capacities <L are calculated optimally too.

 USAGE:
 For given m & L, call InitBasicData(m,L,&l_first).
 Before each run, you must set/change the profits d and d1
  (see below.)
 Then call Run(&b_first), where b is the array of
 upper bounds. To get the results, use
 GetObjValue(L) and RestoreSolution(L,x).
 In the end you may call Clear()
 (e.g. in your main destructor)
*/

/// INTERFACE /////////////////////////

/// 1. MANUAL INPUT - THE PROFITS /////
Vector<double> d; // dimension m+1
  // i.e. indexing [1..m]
//Vector<PVCont> d1; // must be 1..m
Vector<map<int,double> > dI;
typedef map<int,double> PVMap;
  // d1[i] is the list of all modified arcs
  // INCOMING FROM j (.first) plus the value
  // ASSUMPTION: b[j] > 0, b[i] is handled correctly

/// 2. INTERFACE FUNCTIONS ////////////
/// l_first points to the length of the first item
/// in an array of int's:
void InitBasicData(int m_, int L_, int* l_first);
void Run(int* b_first);
double GetObjValue(const int L0);
/// x[1..m] will be calculated:
void RestoreSolution(const int L0, Vector<int>& x);
void Clear();

/// INTERNALS: ////////////////////////

/// 1. INTERNAL TYPEDEFS //////////////

// typedef int value_t; // need scaling to handle real-valued profits
typedef double value_t; // for a while
typedef short int index_t; // item index
class SortedList;

/// 2. INTERNAL VARIABLES /////////////

/// 2.1. CONSTANTS (INPUT): ///////////
int m; // number of product types
int L; // knapsack length
/// These product data arrays have
/// the size m+1. The reason is that the products
/// are indexed 1...m because 0 is the source node
Vector<int> l; // product lengths
Vector<int> b; // upper bounds, some may be zero

/// 2.2. WORKING VARIABLES: ///////////
Vector< Vector<value_t> > f1; // the dynamic tables for each $i$, x[i]>=1.
 // in fact, we can forget these scales for positions l <= l[i].
Vector<value_t> f0; // the table for the case x[i]=0; in the end, the final table
Vector< Vector<index_t> > jji;
  // the indices $j$ of incoming $x_{ji}$; \not= i for bounded i!
Vector<index_t> ii0; // the indices $i$ of the best in each position $d$ (for f0)
Vector< Vector<index_t> > kki;
  // the values of $x_{ii}$ for the bounded $i$ and all $d$ (for unbounded, =0)
Vector<value_t>& p = d; // the integralized profits
  // (attention: summing up till L) - NOT DONE YET. Then use only p
Vector<bool> fTarget; // if $i$ receives a modified arc
typedef list<pair<int,double> > PVCont;
Vector<PVCont> outg_j, outg_not_j; // outgoing arcs in the modified KP graph
int callCnt=0; // debug counter

/// 3. INTERNAL FUNCTIONS /////////////

void InitRun(int* b_first) // before each run
{
  ++callCnt;
  int i;
  b.resize(m+1);
  copy(b_first, b_first+m, b.begin()+1);
  assert(f1.size() > m && f0.size()>L && jji.size() > m && ii0.size()>L
   && l.size() > m && kki.size() > m
    && b.size() > m && d.size() > m && d.size()>m && dI.size()>m
    && fTarget.size()>m && outg_j.size()>m && outg_not_j.size()>m && m>0);
  fill_n(f0.begin(), L+1, 0.0); // FILL F0.
  fill_n(ii0.begin(), L+1, 0);
  outg_not_j[0].clear();
  for (i=1;i<=m;++i) { // <=m !!!!! For {ii} arcs
    outg_not_j[i].clear(); outg_j[i].clear();
  }
  fill_n(fTarget.begin(),m+1,0);

  for (i=1;i<=m;++i) { // FOR EACH PRODUCT i
    assert(l[i] > 0 && l[i] <= L && b[i] >= 0);
    assertm(f1[i].size() > L && jji[i].size() > L && kki[i].size() > L,
      "spprc_dp: Call InitBasicData(m,L,l_first) before usage");
    fill_n(f1[i].begin(), L+1, -1e100); // FILL F[i].
//    fill_n(jji[i].begin(), l[i]+1, 0); // needed ???
    fill(kki[i].begin(), kki[i].end(), 0); // needed for unbounded to restore x.
    /// MODIFIED ARCS: FILLING THE OUTGOING ARRAYS
//    dI[i].sort(); // increasing. Already.
    for_each_in(dI[i], it, PVMap::iterator) {
//      log_n_(2.5," mdf ["<<it->first<<','<<i<<"]="<<it->second
      assert(it->first <= i); // "==": for x{ii}
      fTarget[i] = true;
      if (fabs(it->second) < 1e-6)
        log_n_(1, " OOPS! SPPRC1: No diff? ");
      if (it->first) // not 0
        assert(0 < b[it->first]);
      if (dI[i].begin() == it) // only for the earliest outgoing
        outg_not_j[it->first].push_back(make_pair(i,
	  (0==it->first) ? p[i]+it->second : p[i]));
	    // the normal value if it->first>0
      if (it->first>0 and it->first != i) // not for {ii}
        outg_j[it->first].push_back(make_pair(i,p[i] + it->second)); // +p[i]
      PVMap::iterator it1=it; int jUpper;
      if (dI[i].end() != ++it1) {
        assert(it1->first > it->first); // after sorting
	jUpper = it1->first; // till the next outg.
      }
      else jUpper = i;
      for (int j=it->first+1; j< jUpper; ++j)
        outg_j[j].push_back(make_pair(i,p[i]));
    }
  }
  /// No profit transformation to integers yet
} // EVERYTHING IS SET?

// Add 1 copy of item i based on Not_i
void Shift_Not_i_To_i(int i) { // in fact, propagate01
  if (fTarget[i] or b[i] < 1) return;
  const value_t pi = p[i];
  Vector<value_t>& fi = f1[i];
  Vector<index_t>& jjji = jji[i];
  for (int d=0, dt=l[i]; dt<=L; ++d, ++dt)
    // STARTING FROM 0. // only if really better than prev.
     // values possibly filled by longer arcs:
    if (f0[d] + pi > fi[dt])
    { // or "register value_t& v1=
      fi[dt] = f0[d] + pi;
      jjji[dt] = ii0[d]; // pointing to that i
    }
}

// Propagate from f1[j] to f1[i]
void FillScale11(const int j, const int i, const value_t pji)
{
  Vector<value_t>& fi = f1[i], &fj = f1[j];
  Vector<index_t>& jjii = jji[i];
  for (int d=l[i], dt=l[i]*2; // to ensure >=1
      dt<=L; ++d, ++dt)
    if (fi[dt] < fj[d] + pji) { // or "register value_t& v1=
      fi[dt] = fj[d] + pji;
			jjii[dt] = (short int)j;
    }
}

// Propagate from f0 to f1[i]
void FillScale01(const int i, const value_t pj0i)
{
  Vector<value_t>& f2 = f1[i];
  Vector<index_t>& jjii = jji[i];
  for (int d=0, dt=l[i]; dt<=L; ++d, ++dt)
    if (f2[dt] < f0[d] + pj0i) { // or "register value_t& v1=
      f2[dt] = f0[d] + pj0i;
      jjii[dt] = ii0[d];
    }
}

/// Propagate labels form node j to some nodes i>j.
/// If b[i] == 0 then no propagation to i.
void Propagate01(const int j)
{
  for_each_in(outg_not_j[j], oit, PVCont::iterator) {
    int iT = oit->first;
    assert(dI[iT].size());
    assert(dI[iT].begin()->first == j); // only the first outg to that iT
    if (b[iT] > 0)   // with the default p[iT] only for j>0
      FillScale01(iT, (0==j) ? oit->second : p[iT]);
  }
}
void Propagate11(const int j)
{
  if (j>0 and j<m)
  if (b[j]>0)
  for_each_in(outg_j[j], oit, PVCont::iterator) {
    int iT = oit->first;
    assert(j < iT);
    if (b[iT] > 0)
      FillScale11(j, iT, oit->second); // with the special p[j,i]
  }
}

// Combine 2 scales, a[i]==0 and a[i]>=1, into a[i+1]==0
void Combine_i_and_Not_i(int i) {
  // FillScale1(f1[i],f0,0,0) possible but not effective
  for (int d=0;d<=L;++d)
    if (f0[d] < f1[i][d]) {
      f0[d] = f1[i][d];
      ii0[d] = (short int)i;
    }
}

inline bool BoundedItem(int i) // + EDUARDO
{
  return (l[i] * (b[i]+1) <= L);
}

/// FillTable is a kind of propagation but
/// in the arc {ii}.
void FillTableU(const int i) // for unbounded products
{
  const value_t pi = dI[i].size() ? (
    ((--dI[i].end())->first == i) // {ii} arc
    ? (--dI[i].end())->second + p[i] : p[i] ) : p[i];
  Vector<value_t>& fi = f1[i];
  Vector<index_t>& jjii = jji[i];

  for (int d=l[i], dt=l[i]*2; dt<=L; ++d, ++dt)
    // STARTING FROM li; <li, empty!
    if (fi[dt] < fi[d] + pi) { // or "register value_t& v1=
      fi[dt] = fi[d] + pi;
			jjii[dt] = (short int)i;
    }
}

/// class SortedList //////////////////
/// A list of capacities d with equal d mod l[i]
/// and increasing d and delta (see the references)
class SortedList {
  struct Item {
    value_t delta, fi_d;
    int d;
    Item(int d_, value_t v1, value_t v2) {
      d=d_; delta=v1; fi_d=v2;
    }
  };
  list<Item> lst;
public:
  void insert(int d, value_t v1, value_t v2);
  void getmin(int& d,value_t& fi_d);
  void clear() {
    lst.clear();
  }
};

void SortedList::insert(int d, value_t v1, value_t v2)
{
  while (!lst.empty() && lst.back().delta >= v1)
    lst.pop_back();
  lst.push_back(Item(d,v1,v2));
}

inline void SortedList::getmin(int& d,value_t& fi_d)
{
  d = lst.front().d;
  fi_d = lst.front().fi_d;
  lst.pop_front();
}
/// end class SortedList //////////////

void FillTableB(const int i)
{ // For bounded items: runs only if b[i]>1. What if ==2?
  if (b[i] < 2 /*|| l[i]*2 > L*//*-follows from boundedness*/)
    return; // the latter case done in PropagateArcs()

  const int bi_1 = b[i] - 1; // using bi-1 after propagation
  const int li = l[i];
  const value_t pi = dI[i].size() ? (
    ((--dI[i].end())->first == i) // {ii} arc
    ? (--dI[i].end())->second + p[i] : p[i] ) : p[i];
  Vector<value_t>& fi = f1[i];
  Vector<index_t>& kkii = kki[i];

  int d0, d, k;
  value_t fi_d0;
  SortedList q;

  for (int r=0;r<li;++r) {
    d0 = r+li; // the starting capacity for a sequence
      // of possible consecutive updates.
      // As opposed to pure knapsack, we start not from r
    d = d0+li; // the capacity for an update operation
    k = 1; // k copies of item i are added to capacity d0
    fi_d0 = fi[d0];
    q.clear();

    while (d <= L) {
      register const value_t tv = fi_d0 + k*pi;
      if (fi[d] < tv) { // better value
        q.insert(d, pi*d - fi[d]*li, fi[d]);
	fi[d] = tv;
	kkii[d] = (short int)k; // add item i
//	if (d==1000 and i==94 and tv>1) {
//	  log__(tv);
//	}
	d += li;
	if (k < bi_1)
	  ++ k; // further update in this sequence possible
	else { // k == bi_1
	  q.getmin(d0,fi_d0);
	  k = (d-d0) / li;
	  assert(k <= bi_1);
	}
      }
      else { // no update
        d0= d;
	fi_d0 = fi[d0];
	k = 1;
	d += li;
	q.clear();
      }
    } // end while (d <= L)
  }
}

/// PUBLIC FUNCTIONS //////////////////
/// l_first points to the length of the first item:
void InitBasicData(int m_, int L_, int* l_first)
{
  Clear();
  m=m_; L=L_;
  assert(m>0);
  l.resize(m+1);
  copy(l_first, l_first+m, l.begin()+1);

  f1.resize(m+1);
  jji.resize(m+1);
  kki.resize(m+1);
  d.resize(m+1);
  dI.resize(m+1);
  f0.resize(L+1);
  ii0.resize(L+1);
  fTarget.resize(m+1);
  outg_j.resize(m+1);
  outg_not_j.resize(m+1);
  for (int i=1;i<=m;++i) {
    f1[i].resize(L+1);
    jji[i].resize(L+1);
    kki[i].resize(L+1);
  }
}

void DoneRun() {
  for (int i=1;i<=m;++i)
    dI[i].clear();
}

void Run(int* b_first)
{
  InitRun(b_first);
  for (int i=0;i<=m;++i)
  {
    if (i>0)
      if (0 == b[i]) continue; // ASSUMPTION ON d1 !!!
    Propagate01(i); // before FillTable() - needed for {ii} arcs
    if (i>0) {
      Shift_Not_i_To_i(i); // but only if i not a target
      if (BoundedItem(i))
        FillTableB(i);
      else
        FillTableU(i);
    }
    if (i<m) // after copmuting scale<a_i\ge1>[i] optimally
      Propagate11(i); // even for b[i]==0
    if (i>0) // after propagation from not_i:
      Combine_i_and_Not_i(i); // into Not_i[];
        // for i==m, this gives the final scale
  }
  DoneRun();
}

double GetObjValue(const int L0)
{
  assert(L0 >= 0 && L0 <= L);
  return f0[L0]; // CHANGE IF INTEGER PROFITS
}

/// x[1..m] will be calculated:
void RestoreSolution(const int L0, Vector<int>& x)
{
  assert(L0 >= 0 && L0 <= L);
  int ll = L0;
  value_t mval = -1, epsilon;
  int k, i0, mind = 0;
  x.resize(1);
  x.resize(m+1); // also fill_n(,,0) ?

  mval = f0[L0];
  mind = ii0[L0];
  epsilon = mval * 1e-6;  // CHANGE IF INTEGER PROFITS
  if (epsilon < 1e-8)
    epsilon = 1e-8; // no 0 multipliers
  while (mind > 0) { // while in a product node
    x[mind] += (k=kki[mind][ll]) + 1;
    assert(x[mind] <= b[mind]);
    ll -= l[mind] * k;
    i0 = jji[mind][ll];
    ll -= l[mind];
//    mval -= p[mind][mind] * k; // OF checking not implemented
//    mval -= p[i0][mind];
    assert(i0 <= mind);
    if (BoundedItem(mind))
      {assert(i0 < mind);}
    else {assert(0 == k);}
    mind = i0;
  }
  assert(mind == 0); // we are at the source
//  assert(fabs(mval) <= epsilon);
  assert(ll >= 0);
}

void Clear()
{
  m = L = 0;
  l.clear();
  b.clear();
  d.clear();
  dI.clear();
  p.clear();
  f1.clear();
  f0.clear();
  jji.clear();
  kki.clear();
  ii0.clear();
  fTarget.clear();
  outg_j.clear();
  outg_not_j.clear();
}

} // namespace spprc_dp;

}                                                // SS_END_NAMESPACE__
