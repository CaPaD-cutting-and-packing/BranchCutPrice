#ifndef __SPPRC_DP_H
#define __SPPRC_DP_H

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
 Then fill b[1..m] and call Run(). To get the results,
 use GetObjValue(L) and RestoreSolution(L,x).
 In the end you may call Clear()
 (e.g. in your main destructor)
*/

/// INTERFACE /////////////////////////

/// 1. MANUAL INPUT - THE PROFITS /////
extern Vector< Vector<double> > d; // dimensions m+1,m+1
  // i.e. indexing [from 0..m][to 1..m]
  // e.g. d[0][i] is the profit of arc {0i}
  // i.e. the products are indexed 1...m
  // d.empty() => basic data not initialized

/// 2. INTERFACE FUNCTIONS ////////////
/// l_first points to the length of the first item
/// in an array of int's:
void InitBasicData(int m_, int L_, int* l_first);
void Run(int* b_first);
double GetObjValue(const int L0);
/// x[1..m] will be calculated:
void RestoreSolution(const int L0, Vector<int>& x);
void Clear();

} // namespace spprc_dp

}                                                // SS_END_NAMESPACE__

#endif // __SPPRC_DP_H
