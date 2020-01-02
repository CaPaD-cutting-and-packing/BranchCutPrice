// bkp_dp.cpp
// G. Belov, TU Dresden, 2004

#include "stdafx.h"
#include "bkp_dp.h"
#include "lasthdr.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

namespace bkp_dp
{

/**
 Dynamic programming for BKP, Bounded Knapsack Problem
 with positional values arising from cuts on variables of
 the Arc Flow Formulation of V. de Carvalho

 The table filling procedure follows the "Improved DP"
 from
    U. Pferschy. Dynamic Programming revisited: improving knapsack algorithms.
    Computing, 63:419--430, 1999.
 also discussed in the book
    H. Kellerer, U. Pferschy, D. Pisinger. Knapsack Problems. Springer-Verlag, 2004.

 The overall COMPLEXITY results as O(mL)

 IMPLEMENTATION
 No raster points, particularly no reduced RP.
 Thus, all capacities <L are calculated optimally too.
 No dominance (even w/o positional values)
 No sorting - external order accepted

 USAGE:
 For given m & L, call InitBasicData(m,L,&l_first) once.
 Items must be sorted so as you compute the positional values.
 Before each run, you must set/change general
 and positional item profits
  (see below.)
 Then call Run(&b_first), where b is the array of
 upper bounds. To get the results, use
 GetObjValue(L) and RestoreSolution(L,x).
 In the end you may call Clear()
 (e.g. in your main destructor)
*/

/// INTERFACE /////////////////////////

/// 1. MANUAL INPUT - THE PROFITS /////
Vector<double> d; // dimension m+1, indexes 1..m
// typedef list<pair<int,float> > PVICont;
 /// each list will be sorted internally
Vector<PVCont> pv; // for each product type 1..m,
 // pv[i] is the list of (p, increase v[i][p])

double deps;

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
Vector<value_t> f; // the dynamic table
//Vector<int> lw; // the 0-distance, i.e. lw[d] is how much waste has
  // the pattern in f[d]. Should be initialized? What with pvi?
Vector< Vector<index_t> > kk;
  // for unbounded i, kk[i][ll]==1 if present to pos. ll
Vector<value_t>& p = d; // the integralized profits
  // (attention: summing up till L) - NOT DONE YET

/// 3. INTERNAL FUNCTIONS /////////////

void InitRun(int* b_first) // before each run
{
  b.resize(m+1);
  copy(b_first, b_first+m, b.begin()+1);
  assert(f.size() > m && l.size() > m && kk.size() > m
    && b.size() > m && d.size() > m && m>0 && pv.size() > m);
  fill(f.begin(), f.end(), -1e+100); // filling "empties"
  f[0] = 0;
  for (int i=1;i<=m;++i) { // for each product i
    assert(l[i] > 0 && l[i] <= L
      && b[i] >= 0);
    assertm(kk[i].size() > L,
      "spprc_dp: Call InitBasicData(m,L,l_first) before usage");
    fill(kk[i].begin(), kk[i].end(), 0);
    pv[i].sort();
    if (not pv[i].empty())
      if (pv[i].back().first < L)
        pv[i].push_back(make_pair(L,0.0f));  // dummy end
    if (pv[i].empty())
      pv[i].push_back(make_pair(L,0.0f));  // dummy end
  }
  /// No profit transformation to integers yet
}

inline bool BoundedItem(int i)
{
  return (l[i] * (b[i]+1) <= L);
}

void FillTableU(const int i) // for unbounded products
{
  const value_t pi0 = p[i];
  value_t pi;
  Vector<index_t>& kki = kk[i];
  PVCont::iterator pvi = pv[i].begin();

  for (int d=0, dt=l[i]; dt<=L; ++d, ++dt)
  if (f[d] > -1e+50) {
    if (pvi->first > d) // only 1 comparison normally
      pi = pi0; // and 1 assignment
    else {
      while (pvi->first < d) ++ pvi;
      if (pvi->first == d)
        pi = pi0 + pvi->second; // no ++pvi because of lw[]
      else pi = pi0;
    }
    if (f[dt] < f[d] + pi) { // or "register value_t& v1=
      f[dt] = f[d] + pi;
      kki[dt] = 1;
    }
  }
  // for d
}

/// class SortedList //////////////////
/// A list of capacities d with equal d mod l[i]
/// and increasing d and delta (see the references)
class SortedList {
  struct Item {
    value_t delta, fi_d;
    int d;
    Item(int d_, value_t v1, value_t v2)
      : d(d_), delta(v1), fi_d(v2)
    { }
  };
  list<Item> lst;
public:
  void insert(int d, value_t v1, value_t v2);
  inline void getmin(int& d,value_t& fi_d);
  void clear() {
    lst.clear();
  }
  void add2all(const value_t );
};

void SortedList::insert(int d, value_t v1, value_t v2)
{
  while (!lst.empty() && lst.back().delta >= v1)
    lst.pop_back();
  lst.push_back(Item(d,v1,v2));
}

inline void SortedList::getmin(int& d,value_t& fi_d)
{
//  assert(!lst.empty());
  d = lst.front().d;
  fi_d = lst.front().fi_d;
  lst.pop_front();
}

void SortedList::add2all(const value_t v) {
  for_each_in(lst, il, list<Item>::iterator)
    il->fi_d += v;
}
/// end class SortedList //////////////

void FillTableB(const int i)
{ // For bounded items: runs only if b[i]>0. What if ==1?
  if (b[i] < 1) return;

  const int bi = b[i];
  const int li = l[i];
  const value_t pi0 = p[i];
  Vector<index_t>& kki = kk[i];

  int d0, d, k;
  value_t f_d0;
  SortedList q;

  for (int r=0;r<li;++r) {
    d0 = r; // the starting capacity for a sequence
      // of possible consecutive updates.
    d = d0+li; // the capacity for an update operation
    k = 1; // k copies of item i are added to capacity d0
    f_d0 = f[d0];
    q.clear();
    PVCont::iterator pvi = pv[i].begin();

    /// Here d is as dt in FillU, d0 as d there
    while (d <= L) {
      /// Checking the start position of the last item to be put:
      if (pvi->first <= d-li) {
        while (pvi->first < d-li) ++ pvi;
        if (pvi->first == d-li) { // updating level for this arc
          f_d0 += pvi->second;
          q.add2all(pvi->second); // no ++pvi because of lw
        }
      }
      register const value_t tv = f_d0 + k*pi0;
      if (f_d0 > -1e+50 and f[d] < tv) { // better value
        q.insert(d, pi0*d - f[d]*li, f[d]); // is this ok with pv? NUMERICAL ROUND-OFF HERE? TODO
	f[d] = tv;
	kki[d] = (short int)k; // add item i
	d += li;
	if (k < bi)
	  ++ k; // further update in this sequence possible
	else { // k == bi
	  q.getmin(d0,f_d0); // there are some - e.g. just inserted
	  k = (d-d0) / li;
	  assert(k <= bi);
	}
      }
      else { // no update
        d0= d;
	f_d0 = f[d0];
	k = 1;
	d += li;
	q.clear();
      }
    } // end while (d <= L)
  }
}


void FillTableB_SimpleAndLong(const int i)
{ // For bounded items: runs only if b[i]>0. What if ==1?
	if (b[i] < 1) return;
	assert(pv[i].size()<=1);    //  cout << '^' << flush;
	
	const int bi = b[i];
	const int li = l[i];
	const value_t pi0 = p[i];
	Vector<index_t>& kki = kk[i];
	
// 	int d0, d, k;
// 	value_t f_d0;
	
	for (int d=L-li; d>=0; --d) {
		if (f[d] > -1e+50) {
			for (int k=min(bi, (L-d)/li); k>=1; --k)
				if (f[d+k*li] < f[d] + pi0) { // or "register value_t& v1=
					f[d+k*li] = f[d] + pi0;
					kki[d+k*li] = (short int)k;
				}
		}
	} // for d
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

  f.resize(L+1);
  pv.resize(m+1); // init from caller
  kk.resize(m+1);
  d.resize(m+1);
  for (int i=1;i<=m;++i)
    kk[i].resize(L+1);
}

void Run(int* b_first)
{
  InitRun(b_first);
  for (int i=1;i<=m;++i)
  if (d[i] > deps or pv[i].size() > 1) // the simplest tst
  {
    if (BoundedItem(i)) {
			if (true || pv[i].size() > 1)
				FillTableB(i);
			else 
				FillTableB_SimpleAndLong(i);
		}
    else
      FillTableU(i);
  }
}

double GetObjValue(const int L0)
{
  assert(L0 >= 0 && L0 <= L);
  int LLm=L0; // monotonize:
  for (int LL = L0;LL>=0;--LL)
    if (f[LL] > f[LLm])
      LLm=LL;
//   cout << " BKP_DP: val-1 =" << (f[LLm]-1) << ' ' << flush;
  return f[LLm]; // CHANGE IF INTEGER PROFITS
}

/// x[1..m] will be calculated:
void RestoreSolution(const int L0, Vector<int>& x)
{
  assert(L0 >= 0 && L0 <= L);
  int LLm = L0;
  for (int LL = L0;LL>=0;--LL)
    if (f[LL] > f[LLm])
      LLm=LL;
  int ll = LLm;
//   value_t // mval = -1,
//     epsilon;
  int k, mind = m;
//  int lw0 = lw[L];
  x.resize(1);
  x.resize(m+1); // also fill_n(,,0) ?
  // + consider pos. values of forget mval

//   epsilon = f[L0] * 1e-6;  // CHANGE IF INTEGER PROFITS
//   if (epsilon < 1e-8)
//     epsilon = 1e-8; // no 0 multipliers
  for ( ; mind > 0; -- mind) {
  if (BoundedItem(mind)) { // while in a product node
    x[mind] += (k=kk[mind][ll]);
    assert(x[mind] <= b[mind]);
    ll -= l[mind] * k;
//    mval -= p[mind] * k;
  }
  else
    while(kk[mind][ll]) {
      ++ x[mind];
      ll -= l[mind];
//      mval -= p[mind]
    }
//    assert(lw[d])
  }
//  assert(lw[LLm] == ll);
  assert(mind == 0); // we are at the source
//  assert(fabs(mval) <= epsilon);
  assert(ll == 0); // equality when only at raster points
}

void Clear()
{
  m = L = 0;
  l.clear();
  b.clear();
  d.clear();
  p.clear();
  f.clear();
  kk.clear();
  pv.clear();
}

} // namespace bkp_dp;

}                                                // SS_END_NAMESPACE__
