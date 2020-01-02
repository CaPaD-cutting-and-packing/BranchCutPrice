// spprc_dp.cpp
// G. Belov, TU Dresden, 2004

#include "stdafx.h"
#include "spprc_dp.h"
#include "lasthdr.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

namespace spprc_dp
{

/**
 Dynamic programming for SPPRC,
 Shortest Path Problem with Resource Constraints
 with upper bounds on $x_{ii}$ ($x_{ij}$ are binary)
 and single resource

 The table filling procedure follows the "Improved DP"
 from
    U. Pferschy. Dynamic Programming revisited: improving knapsack algorithms.
    Computing, 63:419--430, 1999.
 also discussed in the book
    H. Kellerer, U. Pferschy, D. Pisinger. Knapsack Problems. Springer-Verlag, 2004.

 The overall COMPLEXITY results as O(mmL)

 IMPLEMENTATION
 No raster points, particularly no reduced RP.
 Thus, all capacities <L are calculated optimally too.

 USAGE:
 For given m & L, call InitBasicData(m,L,&l_first).
 Before each run, you must set/change the profits d
  (see below.)
 Then call Run(&b_first), where b is the array of
 upper bounds. To get the results, use
 GetObjValue(L) and RestoreSolution(L,x).
 In the end you may call Clear()
 (e.g. in your main destructor)
*/

/// INTERFACE /////////////////////////

/// 1. MANUAL INPUT - THE PROFITS /////
Vector< Vector<double> > d; // dimensions m+1,m+1
  // i.e. indexing [from 0..m][to 1..m]
  // e.g. d[0][i] is the profit of arc {0i}
  // i.e. the products are indexed 1...m

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
Vector< Vector<value_t> > f; // the dynamic tables for each $i$
Vector< Vector<index_t> > jj;
  // the indices $j$ of incoming $x_{ji}$
Vector< Vector<index_t> > kk;
  // the values of $x_{ii}$ for all $i$ and $L$
Vector< Vector<value_t> >& p = d; // the integralized profits
  // (attention: summing up till L) - NOT DONE YET

/// 3. INTERNAL FUNCTIONS /////////////

void InitRun(int* b_first) // before each run
{
  b.resize(m+1);
  copy(b_first, b_first+m, b.begin()+1);
  assert(f.size() > m && jj.size() > m && l.size() > m && kk.size() > m
    && b.size() > m && d.size() > m && m>0);
  for (int i=1;i<=m;++i) { // for each product i
    assert(d[i].size() > m && l[i] > 0 && l[i] <= L
      && b[i] >= 0);
    assertm(f[i].size() > L && jj[i].size() > L && kk[i].size() > L,
      "spprc_dp: Call InitBasicData(m,L,l_first) before usage");
    fill_n(f[i].begin(), l[i], 0); // filling zeros
    fill_n(jj[i].begin(), l[i], -1); // up to length l[i]
    fill(kk[i].begin(), kk[i].end(), 0);
  } // -1: "no arc coming to i with this resource position"
  /// No profit transformation to integers yet
}

/// Propagate labels form node i to nodes j>i.
/// If b[j] == 0 then no propagation to j.
void PropagateArcs(const int i)
{
  if (i == 0) // i.e. from the source node: simple
    for (int j=1;j<=m;++j) {
      const value_t vp = b[j] > 0 ? p[0][j] : 0;
      const index_t vi = b[j] > 0 ? 0 : -1;
      fill(f[j].begin() + l[j], f[j].end(), vp);
      fill(jj[j].begin() + l[j], jj[j].end(), vi);
    } // 0 : incoming from the source
  else // when i>0, use the existent table f[i]:
    for (int j=i+1;j<=m;++j)
    if (b[j] > 0) // !!!
    {
      const value_t pij = p[i][j];
      Vector<value_t>& fi = f[i];
      Vector<value_t>& fj = f[j];
      Vector<index_t>& jjj = jj[j];

      for (int d=l[i], dt=l[i]+l[j]; dt<=L; ++d, ++dt)
        if (fj[dt] < fi[d] + pij) { // or "register value_t& v1=
          fj[dt] = fi[d] + pij;
					jjj[dt] = (short int)i;
        }
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
  const value_t pi = p[i][i];
  Vector<value_t>& fi = f[i];
  Vector<index_t>& jji = jj[i];

  for (int d=l[i], dt=l[i]*2; dt<=L; ++d, ++dt)
    // STARTING FROM li; <li, empty!
    if (fi[dt] < fi[d] + pi) { // or "register value_t& v1=
      fi[dt] = fi[d] + pi;
			jji[dt] = (short int)i;
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
  const value_t pi = p[i][i];
  Vector<value_t>& fi = f[i];
  Vector<index_t>& kki = kk[i];

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
	kki[d] = (short int)k; // add item i
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

  f.resize(m+1);
  jj.resize(m+1);
  kk.resize(m+1);
  d.resize(m+1);
  d[0].resize(m+1);
  for (int i=1;i<=m;++i) {
    f[i].resize(L+1);
    jj[i].resize(L+1);
    kk[i].resize(L+1);
    d[i].resize(m+1);
  }
}

void Run(int* b_first)
{
  InitRun(b_first);
  for (int i=1;i<=m;++i)
  {
    PropagateArcs(i-1);
    if (BoundedItem(i))
      FillTableB(i);
    else
      FillTableU(i);
  }
}

double GetObjValue(const int L0)
{
  assert(L0 >= 0 && L0 <= L);
  value_t mval = 0;
  for (int i=1;i<=m;++i)
    if (f[i][L0] > mval)
      mval = f[i][L0];
  return mval; // CHANGE IF INTEGER PROFITS
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

  for (int i=1;i<=m;++i)
    if (f[i][L0] > mval) {
      mval = f[i][L0];
      mind = i;
    }
  epsilon = mval * 1e-6;  // CHANGE IF INTEGER PROFITS
  if (epsilon < 1e-8)
    epsilon = 1e-8; // no 0 multipliers
  while (mind > 0) { // while in a product node
    x[mind] += (k=kk[mind][ll]) + 1;
    assert(x[mind] <= b[mind]);
    ll -= l[mind] * k;
    i0 = jj[mind][ll];
    ll -= l[mind];
    mval -= p[mind][mind] * k;
    mval -= p[i0][mind];
    mind = i0;
  }
  assert(mind == 0); // we are at the source
  assert(fabs(mval) <= epsilon);
  assert(ll >= 0);
}

void Clear()
{
  m = L = 0;
  l.clear();
  b.clear();
  d.clear();
  p.clear();
  f.clear();
  jj.clear();
  kk.clear();
}

} // namespace spprc_dp;

}                                                // SS_END_NAMESPACE__
