#ifndef __PROBL_CP22_H
#define __PROBL_CP22_H

// 2D 2-staged (=>guillotine) cutting (=>knapsack) problem
// item/problem rotation (1. stage horiz/vert)

// Level cut
// opt test: effective length => LPBnd better
// col gen
// solution printing

#include "problem.h"
//#include "lp.h"
#include "bb.h"
#include "bbcuts.h"
#include "cutgen.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

class CP22LevelCut : public SACutSE {
public:
  int Type() const { return -2; }
  int GetSlackCoef() { return -1; } // - gcd !!!
  double lhs;
  CP22LevelCut() { alfa = 0; }
  double LHS () const { return lhs; }
  double RHS () const { return 1e+100; }
  double GetRHS () const { return LHS(); }
  double CalcRHS__ (Column *c)
  { return lhs = c->GetObj(); }
//  double Calculate__(Column * c); // the real obj. coef.
  bool IntegerSlack() { return true; } // only for int gcd valid
};

class CP22 : public Problem {
public:

  const char * TypeName() { return "CP22"; }
  /// (DIS-)ABILITY DEFINITIONS
  bool SACutsPossible() { return true; }
  bool BranchingPossible() { return true; }
  bool CanBeInfeasOnRestrPool() { return false; }
protected:

struct Piece {
  size l,w;
  double c; // the weight
  int b;
  float wgt;
//  int i_; // Number in the input
  // - no, use its size
  int i0; // Number before sorting

  // Init w before sorting!
  bool operator > (const Piece & pc2) const
    { return wgt > pc2.wgt; } // for sorting
};

  typedef Vector<Piece> PieceContainer;

// Algorithms:
  auto_ptr<BB> bb;

// public:
  size L,L0, W,W0;
  int m,m0;
  PieceContainer pc,pc0;

  double cMax, cMin;
  double areaPriceMin;
  double gcdC; /// gcd of Ci's
  size wMin;
/// OPTTEST:
  int opttest;
  int lastPos; // raster point
  Vector<int> rpc;
  double rpM; // multiplier: (raster point) * rpM = real,
    // used to decrease discretization. Maaybe == gcdC.
    // if != gcdC, then somewhat greater

  //const v_float veps;

  virtual void FFSBasis
    (const PieceContainer &,ColumnList &); // ???
  virtual void SetObj(Column*);
  virtual void SetWConstr(Column*);
  virtual double GetW(Pattern *p);
  virtual void GenColPure(ColSet& cs,const d_vec& d);
  virtual void GenColWithCuts(ColSet& cs,const d_vec& d);

  double GetDEps() { return deps; /*1*eps*/ }
  double GetVEps() { return max(1e-6, -GetLPBnd()*deps); }  // ?????
  double GetBBEps() { return GetDEps(); }
  double GetXEps();

  virtual double raster_ceil(const double v);
  virtual double raster_below(const double v);
  virtual void InitOptTest();

public:
  CP22(Env *penv,
    const char *in,const int i,const char *n)
    :Problem(penv,in,i,n), bb(new BB()), cMax(1), m(0) { }
  virtual ~CP22() { }
  virtual void Init(); // Check consistency + initialized?
  virtual int Read(istream&,long long);
  virtual void InitProblem();
//  virtual void Reallocate();
  virtual void HeuristicLPBasis(ColumnList &bas)
  { FFSBasis(pc, bas); }

  virtual void FillSigns();

  virtual void GenCol(ColSet& cs,const d_vec& d);

//  virtual void CalcLagrRel(const d_vec &);
//  virtual void UpdateLagrBnd(const d_vec &);

/*  virtual void UpdateCurrentLPBnd(v_float );
  virtual void UpdateLPBnd();
  virtual void UpdateHeurBnd(v_float );

  virtual void SetCurrentLPBnd(v_float );
  virtual void SetLPBnd0(); */

  virtual int Dim() { return m+1; } // ???
  virtual int the_m() { return m; }

//  bool Optimum();

  virtual void PrintColumn(ostream&,Column*);
  void MakePattern(Pattern * p, Column * c);
  void MakeColumn(Column * c,Pattern * p);
//  CP22LevelCut levelCut;
  CSPLevelCut lc1;
  LPCut * GetLevelCut() { return &lc1; }

  void CopyForbCols(list<Column*> & fc);

////////////////////////////  ROUNDING + INTEGER SOLUTION
/// Input: see Problem
/// Output: see Problem
protected:
  Vector<double> xi;
  double zr;   // obj
  Vector<double> xf;  // frac parts
  Vector<double> br;  // reduced rhs
  Vector<int> fracused; // indices of rnd-up cols
  int nRPE;   // n res pr ext
  static int RPEMax;
  size WLeft;

  void GenColForHeur
    (ColSet *,const d_vec &d,
    Vector<double> &b, size WLeft);
  bool TestIntegrality();
  void SaveRoundedPart();
  bool Round();
  void GetRHS(Vector<double> &br);
  class SVC;
  friend class SVC;
  bool CompleteIntSol();
  bool VaryResProblem();
  void AddResidualSolution(SVC & svc);
  void ControlSolution();

public:
  bool ConstructIntSol();
  void PrintProblem(ostream &);

public:
  static opt::OptContainer Options();
  static opt::OptSection opt;
  static bool fFirstCut1stD;
  static bool fSortPieces, fMergePieces;
  static int nStartBasis;
  static double nStepsMin0;
  static double nStepsMinInc;
  static double deps, bb_eps;
  static double outputLevel;
};//___class_CSP1_______________________________________

}                                                // SS_END_NAMESPACE__

#endif // __PROBL_CP22_H
