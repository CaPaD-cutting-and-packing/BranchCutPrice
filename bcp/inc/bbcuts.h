#ifndef __BBCUTS_H
#define __BBCUTS_H


// Use the solution of DP for bound as incumbent.
// How data storage ? 0-based arrays. Auto-dimensioning
// Cuts: case alpha==1 !!
// Slacks elimination ?
// 1. Level cut 2. Dual simpl 3. +Cuts. Constraint pool
// Stand over the toes.
// Not full paths may be better for the rnd. heuristic?
// Or change the latter.
// FOR 2D: NEED THE WHOLE d-Vector, not only for pieces,
// also for the width constr.


#include "problem.h"
#include "cutgen.h"


#define BBCuts__NotOnlyFullPatterns// MUST BE


namespace bgn {                     // SS_BEGIN_NAMESPACE__


// Algorithm B&B for col.gen. with superadditive cuts.
// Implemented static, so don't instantiate the class.
// The basic procedures for 1D, 1D MSL, 2D KP.
// The caller must set data and obey some rules.
/* E.g. nStepsMax should increase with time.
 cuts__ should contain the cuts' weights in u's.
 Before using the class, call Reallocate(m,mc).
 Before Init(), input vars must be set.
 Before Run(), Init() must have been called.
  */
class BBCuts {
public:
  struct Piece {
    size l;
    int
      b,
      i0;       // the original index before sorting
    double
      d,        // the simplex mult. or weight
      dAppr,
      q,
      dSApprA,
      qma,      // q max after (inclusive)
      dApprMA, // max dAppr after (inclusive)
      w;        // weight for sorting
    size
      lma,      // l min after (not inclusive)
//      lma1,   // l min after (inclusive)
      lsa;      // l sum after
    bool operator > (const Piece & pc2) const
    { return w > pc2.w; } // for sorting
  };//__________________________________________________
  typedef Vector<Piece> PieceContainer;
  //////////////////////////////////////////////////////
  /// INPUT:: (not to be modified inside)
  //____________________________________________________
/// Number of pieces types, addi constr, and cuts:
  static int m/*,mc,mu*/;
  static int mc; // = N addi constr.
  static size L; // Stock length
  static PieceContainer pieces__; // Each knows l,b,d
/// Cutting planes:
  static CutSet cuts__;    // with weights incaps.
  static d_vec d__;   // the weigths with those of addi
  static Pool<Column> * forbidden__;
/// Initial, eg incumbent, lower bound:
  static double zLowerInitial;
/// A priori lower bound, eg for primal simplex or
/// obtained from Lagrange relaxation.
/// E.g. CSP1: 1+eps
  static double zMin;
  static bool foundInit; // whether a good sol. already kn
/// The tolerance by which a better solution must be >=:
	static double eps; // =smth*bbpi_eps
// 	static double deps;
/// Early termination controls:
  static double nStepsMin; // =16384*8=nStepsMin0
// Use global prms nStepsMaxIncrRatio + interval
  static bool fRandomizeWeights;
// + How -- params
  static bool
    fConsiderNotOnlyFullPatterns,
    fTryAlsoNotOnlyFullPatterns;
  static double outputLevel;
  static double nStepsTooMuch;
  static bool fCheckBnd, fUpdateBnd;
  //////////////////////////////////////////////////////
  /// OUTPUT::
  //____________________________________________________
/// Result flag. =true e.g. if (z > zMin).
  static bool found;
/// The value of the best found solution
  static double z;
  static bool fETerm;
/// A corresponding solution in compact form
//  typedef ss::Pattern Pattern;
  static Column * colBest;
  static ColSet * colsetRes;
  //////////////////////////////////////////////////////
  /// STATISTICS::
  //////////////////////////////////////////////////////
private:
  /// VARIABLES::
  /// (now all static, no parallelization)
  /// (even the latter can use static with templates)
  /// or use separate executables
  //____________________________________________________
/// A lower bound, initally max(zMin,zLowerInitial)+eps,
/// thereafter (bestSolution + eps):
  static double lb;
/// An upper bound for the partial solution:
  static double ub1; // = dAppr*a - alfL;
  static double sad; // the accumulated o.f. value
  static PieceContainer pieces;
  static size dL;
  static int ms,
    k,
    pk;
  static Vector<IX> ix;
  static double nSteps;
  static CutSet cuts;
//  static set<Vector<IX> > forbidden;
  static Column colTmp, colTmpBetter;
  static double
    apprError,
    bndConstTerm,
    zfConstTerm;
  static Vector<double> dAppr; // appr multipliers
  static bool fPrevPieceType;
  static opt::OptSection opt;
  //////////////////////////////////////////////////////
  /// Methods::
  //____________________________________________________
public:
  static void Reallocate(int m);
// Precondition: All input vars set!
  static void Init();
// Precondition: Init() must have been called
  static void Run(); // need this ??
  // Randomization? Restarts? Time control ?
  // Perturbation: slightly change weights
protected:
  static void Optimize();
  static void InitEnumeration();
  static void RecalcSums();
  static double UB2(int ); //the bound for the free length
  static void AddMaxItemsOfThePiece();
  static void ConsiderCalcObj1();
  static void ConsiderCalcObj2();
  static void CalcObjective();
  static void SignalBetterSolution();
  static bool EnumerationOver();
  static void RemoveLastItem();
  static bool CheckForEarlyTermination();
  static void CopyInput();
  static void ReduceNCuts();
  static void CalcBounds();
  static void CalcConstTerms();
  static void SortPieces();
  static void InitServiceData();
  static void RandomizePieceWeightRatio(Piece *pc);
  static void SaveSolutionToTemp();
  static void SaveBetterSolutionFromTemp();
  static void CopyResult();
  static opt::OptContainer Options();
  static void Print();
  static void DoStatistics(int what, ostream &);
}; // BBCuts
  //____________________________________________________


inline double BBCuts::UB2(int k) {
  /// No dynamic programming yet!
  /// The simplest version of the lin. bound: cont. rel.
  register Piece * pc = &pieces[k];
  return // qma = max. ratio (dAppr[i]/l[i]), i>=k
    pc->qma >= 0 ?
    ( (dL < pc->lsa) ? // but add only non-neg in dSApprA
      (pc->qma * (double)dL)
      : pc->dSApprA )
      : pc->dApprMA;
  ;
} //____________________________________________________
inline void BBCuts::ConsiderCalcObj1() {
//#ifdef BBCuts__NotOnlyFullPatterns
//  if (fConsiderNotOnlyFullPatterns)
    CalcObjective(); ///////////////////////////////
//#endif
} //____________________________________________________
inline void BBCuts::ConsiderCalcObj2() {
//#ifdef BBCuts__NotOnlyFullPatterns
//  if (fConsiderNotOnlyFullPatterns)
    if (not fPrevPieceType)
    // if (pieces[k].dsl < 0) // if the pieces can dec. z
        CalcObjective(); ///////////////////////////////
//#endif
} //____________________________________________________
inline bool BBCuts::EnumerationOver() { // + some cleanup
  if (0 <= pk)
    return false;
  return true;
} //____________________________________________________


}                                                // SS_END_NAMESPACE__


#endif // __BBCUTS_H
