#ifndef __SPPRC_DP1_H
#define __SPPRC_DP1_H

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
extern Vector<double> d; // dimension m+1
  // i.e. indexing [1..m]
//typedef list<pair<int,double> > PVCont;
//extern Vector<PVCont> d1; // must be 1..m
extern Vector<map<int,double> > dI;
  // dI[i] is the list of all modified arcs,
  // the double: the modification, increase relative to p[i]
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

} // namespace spprc_dp1

}                                                // SS_END_NAMESPACE__

#endif // __SPPRC_DP1_H
