#ifndef __PROBLEM_H
#define __PROBLEM_H

/** Separation of algorithmic and problem-specific features.
 Class Problem defines problem- and model-dependent
 operations. Used by algorithms: BCP, SVC & deriv.
 BUT algorithm-specifics should go there to use it in full */

#include "lpcol.h"
#include "cuts.h"
//#include "bcp.h"
#include "pool.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

typedef long long size; // product length(, width, etc)
const size SIZE_MAX__ = LONG_LONG_MAX;

typedef int demand;
const demand DEMAND_MAX__ = INT_MAX;

/** A pair of integers. Should be pair<int ,int > */
struct IX {
  int i,x;
  void set(const int i_,const int x_) { i=i_; x=x_; }
  IX(const int i_=0,const int x_=0) :i(i_), x(x_) { }
  bool operator<(const IX& ix2) const
  { return i==ix2.i ? x<ix2.x : i<ix2.i; }
};

/** Internal representation of patterns.
 Contains the intensity $x_j$. */
struct Pattern {
  Vector<IX> ix; // The original constraints
//  Vector<IX> addi; // addi constraints
  typedef Vector<IX>::iterator iterator;
  typedef Vector<IX>::iterator ix_iterator;
  void PushIX(const int i,const int x)
  { ix.push_back(IX(i,x)); } // All to be PUSHed !

  double ofc; // the O.F. coeff
  mutable double x; // the intensity
  void SetObj(double o) {ofc=o;}
  double GetObj() { return ofc; }

  void * addi_info;
  void * & GetAddiInfo() { return addi_info; }

  explicit Pattern(const int m) : ofc(0)
  { ix.reserve(m); }
  Pattern() : ofc(0) { }
  void clear() { ix.clear(); }

  bool operator<(const Pattern& p) const
  { return ix < p.ix; } // sort before !!!!
}; // struct Pattern -------------------------

struct PatternList {
  list<Pattern> pl;
  Pattern &Add(const Pattern &p=Pattern())
  {pl.push_back(p); return pl.back();}
};

/// define the valid signs of dual multipliers:
typedef Vector<signed char> constrtypevec;

/// Abstract cutting&packing problem
class Problem : public Alg {
public:

  const char * TypeName() { return "[some problem]"; }
  /// (DIS-)ABILITY DEFINITIONS
  virtual bool SACutsPossible() { return false; }
  virtual bool BranchingPossible() { return false; }
  virtual bool CanBeInfeas() { return false; } // for most
  virtual bool CanBeInfeasOnRestrPool() { return true; } // for most
  virtual bool NoForbiddenColsWithHyperplanes() { return true; }
    // for most, because no b&b col.gen.
  virtual bool NoCutsWithHyperplanes() { return true; }

  /// Problem info:
  const char * infile;
  int inst;
  const char * prName;

  /// Some info which is to be kept consistent by the derived:
  constrtypevec validsign;
  Column clBest; // after col gen
  double redCostBest;
  d_vec b; // the b of Ax <>= b,
    // the (in)equality type being defined by validsign.
    // Only initial constraints.
  virtual Column * GetBestCol()
  { assert(redCostBest<1e100); return &clBest; }
  virtual double GetBestRedCost()
  { return redCostBest; }

  virtual void BrOnHP(const Vector<Column*>& colpool,
    const Vector<ColId>& cols, const Vector<double>& lpx,
    LPCut* &hpU, LPCut* &hpL)
  { // Branching on Hyperplanes
    throw runtime_error("BrOnHP Not implemented for this problem type.");
  }
  /// Provides temp. changes in ba, incl. existing cuts:
  virtual bool TempModifyLP(const Vector<Column*>& colpool,
    const Vector<ColId>& cols, const Vector<double>& lpx,
    const d_vec& lpd, list<LPCut*> & cuts1)
  { // Cutting on Hyperplanes: need to make a test
    // with a temp. modified LP? LP value known in Problem
    return false;
  }
  /// restores original ba:
  virtual void SeparateHP(const Vector<Column*>& colpool,
    const Vector<ColId>& cols, const Vector<double>& lpx,
    const d_vec& lpd, list<LPCut*> & cuts1)
  { // Cutting on Hyperplanes
    throw runtime_error("Problem specific cuts not implemented for this problem type.");
  }

/// This should be kept up2date by the using class
/// ini. filled by Problem, then by user (cuts' rhs):
  Vector<double> ba; // the rhs of all constr.,
    // class Problem needs this ??? Lagrange? but with cuts?
  CutList * cuts; // The cuts present in formulation
  d_vec b_cg; // the rhs for col. gen.,
    // may be reduced by branching (->proper columns)
  Pool<Column> * forbidden; // The dictionary of forbidden

/// Bounds info, should not be doubled this way!
  void SetLPBnd(double v) { lpb = v; }
  void SetLocalLPBnd(double v) { llpb = v; }
  void SetCurrentLPBnd(double v) { lpb1 = v; }
  void SetCurrentLPValue(double v) { lpv1 = v; }
  void SetHeurBnd(double v) { zi = v; } // when better MIP
  void UpdateHeurBnd(double v)
    { if ((v)<zi) zi = (v); } // when better MIP

  bool fLPOpt; // whether a true LP opt
    // is reached in last LP solution (no dual bnd worked)
  bool fLocalUB; // set by BCP

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
  virtual int Dim() =0; // the original dimension w/o cuts
  virtual int the_m() { return -1; }
    // The number of constr; not really used
  double GetLPBnd() {return lpb;}
  double GetLocalLPBnd() {return llpb;}
  double GetCurrentLPBnd() {return lpb1;}
  double GetCurrentLPValue() {return lpv1;}
  double GetHeurBnd() { return (zi); }

  double GetLocalLagrValue() {return lrv1;}// must be: llrv;
  double GetCurrentLagrValue() {return lrv1;} // for subgr
   // - current: not the best, use local for bound
  void UpdateLagrBnd(const d_vec &d);
  void CalcLagrRel(const d_vec &d);

  Problem(Env* penv,
    const char*in,const int i,const char*pn)
    : Alg(penv), infile(in), inst(i), prName(pn),
  cuts(0), lpb(-1e100), lpb1(1e100), lpv1(1e100),
  llrv(-1e100), lrv1(1e100), redCostBest(1e100) { }
  virtual ~Problem() { }

  /// Input problem using the starting value (2nd param)
  virtual int Read(istream&,long long) =0;
  virtual void Init() =0;
  // to init internal algorithms after data input
  virtual void Done() { }
  virtual void PrintLog(ostream &ofs) { }
  virtual void InitLP();
  // after adding cuts, before solving (modified) LP

// if exact (Ax=b), then no slacks needed in CSP1:
  virtual void HeuristicLPBasis(ColumnList &)=0;

  double lpb, llpb,
    lpb1, // last (current), in this iteration
    lpv1,
    llrv, // Lagrange bounds: local
    lrv1, // at the moment
    zi; // best IP known

// getting bound from a lp/lagr value:
  virtual double raster_ceil(const double v)
  { return ceil(v - 1e-8); }
  virtual double raster_below(const double v)
  { /*assert(floor(v) == v);*/ return v-1; }

  virtual void GenCol(ColSet&,const d_vec&) =0;

  // All tolerances better 1e-6 ... sometimes need even more!!!
  virtual double GetVEps() { return max(1e-6, fabs(GetLPBnd()) * 1e-6); }
   // "value-epsilon" -- relative to O.F. value
  virtual double GetDEps() {return 1e-6;}  // for dual vars
  virtual double GetRCEps() {return 1e-6;}  // for red. costs

  virtual double GetBBEps() {return 1e-6;} // for bb-colgen
  virtual double GetXEps() {return 1e-6;}
    // for the column intensities (used in rounding)

  virtual bool Optimum()
  { assert(lpb <= zi); return lpb == zi; }

  virtual LPCut * GetLevelCut()=0;
  // and sets the level to llb - CURRENTLY DISABLED
  virtual double GetLevelCutRHS() { return GetLocalLPBnd(); }
    // - when not negating // when is it adjusted?
  virtual double CalcRedCost
    (Column*, CutList *, const d_vec&); // for checking

////////////////////////////  ROUNDING + INTEGER SOLUTION
/// INPUT:
  bool fFastRounding; // when after dual simplex,
    // intergal LP solution is not opt!!!
  void SetFastRounding(bool f) { fFastRounding = f; }
//  Vector<Column*> lpcol; // but no slacks
  Vector<Pattern> pat; // receiving already patterns - not ok,
  // should be just columns
  Vector<double> lpx;
  d_vec lpd;  // simplex multipliers
// :END INPUT

/// Integer rounding:
/// feasibility (=integrality) must have been tested before
  virtual bool ConstructIntSol()=0; // returns1 if opt
/// For PMP, we need to unite equal patterns:
  virtual void ExtractBestSolution(Problem*) { } // calls this:
  virtual void ExtractSolution(Vector<Pattern>&,d_vec&) { }
  // returns 1 if opt
  virtual bool ExtractSolutionPart(Vector<Pattern>&,d_vec&)
  { return 0; } // used in SVC to combine with the rounded part

  /// OUTPUT:
  /// The best known solution:
  Vector<Pattern> patBest; Vector<double> xiBest;
  Vector<Pattern>& GetBestPatList() { return patBest; }
  Vector<double>& GetBestX() { return xiBest; }

  virtual void PrintColumn(ostream&,Column*);
  virtual void MakePattern(Pattern *p,  Column *c);
  virtual void MakeColumn(Column *c,Pattern *p);

  virtual void PrintBestSolution(ostream &) { }
  virtual void PrintProblem(ostream&) = 0;

  /// STATISTICS:
  virtual void InitStat() { }
  virtual void PrintStat(ostream&) { }
};

}                                                // SS_END_NAMESPACE__

#endif // __PROBLEM_H
