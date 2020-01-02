#ifndef __PROBL_CSP2_H
#define __PROBL_CSP2_H

#include "probl_csp1.h"
//#include "lp.h"
#include "bb.h"
#include "bbcuts.h"

// 2D:
#include "ggk2.h"
//#include "wang.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

class CSP2 : public CSP1 {
public:

  const char * TypeName() { return "CSP2"; }
  /// (DIS-)ABILITY DEFINITIONS
  bool SACutsPossible() { return false; }

private:
struct Piece {
  size l,w;
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

// ALGORITHMS:
//  auto_ptr<BB> bb;
  GGK2 ggk2;
//  Wang wang; // ALL UNCONSTRAINED YET

// public:
  size L,L0,W,W0;
  int m,m0;
  PieceContainer pc,pc0;

// COORDINATE INFO OF THE COLUMNS:
  list<Vector<GGK2::IXY> > crd;

  //const v_float veps;

  virtual void FFSBasis
    (const PieceContainer &,ColumnList &); // ???
  virtual void GenColPure(ColSet& cs,const d_vec& d);
//  virtual void GenColWithCuts(ColSet& cs,const d_vec& d);

  double GetVEps() { return deps; /*1*eps*/ }
  double GetDEps() { return deps; }
//  v_float GetBBEps() { return bb_eps; }
//  double GetXEps();

//  virtual v_float raster_ceil(const v_float v)
//  { return ceilEps(v, deps); }

public:
  CSP2(Env *penv,
    const char *in,const int i,const char *n)
    :CSP1(penv,in,i,n),
    ggk2()//, wang(ggk2.val,ggk2.rpl,ggk2.rpw) 
  { }
  virtual ~CSP2() { }
  virtual void Init(); // Check consistency + initialized?
  virtual int Read(istream&,long);
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

  virtual int Dim() { return m; } // ???

//  bool Optimum();

  virtual void PrintColumn(ostream&,Column*);
//  void MakePattern(Pattern & p, Column & c);
  CSPLevelCut levelCut;
  LPCut * GetLevelCut();

////////////////////////////  ROUNDING + INTEGER SOLUTION
/// Input: see Problem
/// Output: see Problem
/// Working vars: see CSP1
private:
/*  Array<Pattern> pat; // store rounded part
  Array<double> xi;
  double zr;   // obj
  Vector<double> xf;  // frac parts
  Vector<double> br;  // reduced rhs
  Vector<int> fracused; // indices of rnd-up cols
  int nRPE;   // n res pr ext
  static int RPEMax;
*/
  void GenColForHeur
    (Column *col,const d_vec &d, Vector<double> &b);
//  void ConvertLPCols();
//  bool TestIntegrality();
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
  void PrintBestSolution(ostream &);
  void PrintProblem(ostream &);

private:
  static opt::OptContainer Options();
  static opt::OptSection opt;
  static bool fSortPieces, fMergePieces;
  static double WangGammaStep;
  static double WangNSol;
//  static double nStepsMinInc;
  static double deps;
  static int RPEMax;
  static double outputLevel;
};//___class_CSP2_______________________________________

}                                                // SS_END_NAMESPACE__

#endif // __PROBL_CSP2_H
