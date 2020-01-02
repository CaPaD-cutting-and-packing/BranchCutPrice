#ifndef __BKP_DP__
#define __BKP_DP__


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
extern Vector<double> d; // dimension m+1, indexes 1..m.
  // d.empty() => basic data not initialized
typedef list<pair<int,double> > PVCont;   // THIS SEEMS TO BE FOR AFF arc prices!
 /// each list will be sorted internally
extern Vector<PVCont> pv; // for each product type

extern double eps, deps;

/// 2. INTERFACE FUNCTIONS ////////////
/// l_first points to the length of the first item
/// in an array of int's:
void InitBasicData(int m_, int L_, int* l_first);
void Run(int* b_first);
double GetObjValue(const int L0);
/// x[1..m] will be calculated:
void RestoreSolution(const int L0, Vector<int>& x);
void Clear();

} // bkp_dp

}                                                // SS_END_NAMESPACE__

#endif // __BKP_DP__
