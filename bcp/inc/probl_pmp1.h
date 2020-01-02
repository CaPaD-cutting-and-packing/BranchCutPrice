#ifndef __PROBL_PMP1_H
#define __PROBL_PMP1_H

#include "probl_csp1.h"
#include "bb.h"
#include "bbcuts.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

// Cuts could accelerate finding integer?
// Every generated col: properly set obj
// attention: with artificial column much changes
// functionality definitions also here

class PMP1 : public CSP1 {
public:

// OPTIONS:
  double deltaK; // rndup(dK%*cspIP)
  static double deltaKpercent;
    // upper bound on addi stock rel to cspIP. def: 1%
  static double pF; // fixed (setup) cost
  static double pV; // variable cost (per stock sheet)

  double gcdPP; // = gcd (pF, pV)
  double q0max; // max multiplicity of a pattern.
    // Set in a node

  double
    lbD, // MIN number of different patterns
    ubD, // MAX
    cspIP,
    cspLP
  ;
  CSP1 csp; // used for rounding
  SSVC ssvc; // for sequencing a given sol.
  Random rndSeq; // frequency of appl

  const char * TypeName() { return "PMP1"; }
  /// (DIS-)ABILITY DEFINITIONS
  bool SACutsPossible() { return false; } // yet
  bool BranchingPossible() { return true; }
//  bool CanBeInfeasOnRestrPool() { return deltaK<1e50; }
    // SHOULD RET. 1 if deltaK < inf
    // but not with artificial column
  // BUT in CSP returns true => here also.
//protected:

  virtual void FFSBasis
    (const PieceContainer &,ColumnList &); // add (m,1)
//  virtual void GenCol(ColSet& cs,const d_vec& d);
  virtual void GenColPure(ColSet& cs,const d_vec& d);
//  virtual void GenColWithCuts(ColSet& cs,const d_vec& d);

  virtual double raster_ceil(const double v) // opt. test
  { return gcdPP * ceilEps(v/gcdPP, GetDEps()/gcdPP); }

public:
  PMP1(Env *_penv, const char *_in,int _i,const char *_n);
  virtual ~PMP1() { }
  virtual void Init(); // Check consistency + initialized?
//  virtual void InitProblem();
  virtual void PrintLog(ostream &);
  virtual void PrintStat(ostream &);

  virtual void FillSigns();

//  virtual void CalcLagrRel(const d_vec &);
//  virtual void UpdateLagrBnd(const d_vec &);

  virtual int Dim() { return m+1; }
    // !!!!!!!!!!!  ??? how in CSP? Set all to m?

//  virtual void PrintColumn(ostream&,Column*);
////////////////////////////  ROUNDING + INTEGER SOLUTION
/// Input: see Problem
/// Output: see Problem
//protected:
/*  class SVC : public CSP1::SVC {
  public:
    SVC(...);
  protected:
    void Execute();
  }

  friend class SVC; */
  virtual bool CompleteIntSol();

public:
  bool ConstructIntSol();
  virtual void ExtractBestSolution(Problem*);
  virtual void ExtractSolution(Vector<Pattern>&,d_vec&);
  // returns 1 if opt
  virtual bool ExtractSolutionPart(Vector<Pattern>&,d_vec&);

  virtual void PrintColumn(ostream&,Column *c);
  virtual void MakePattern(Pattern *p,  Column *c);
  virtual void MakeColumn(Column *c,Pattern *p);

// OPTIONS:
public:
  static opt::OptContainer Options();
  static opt::OptSection opt;

  static double CSPOutpLev;
  static double CSPTiLim;
  static double BPPTiLim;

  static bool fSequence;
  static double seqFreq;

//  static bool fModelEquality;
//  static bool fSortPieces, fMergePieces; // stay in CSP
  static int nStartBasis;
//  static double nStepsMin0;
//  static double nStepsMinInc;
//  static double deps, bb_eps;
  static double outputLevel;
};//___class_CSP1_______________________________________

}                                                // SS_END_NAMESPACE__

#endif // __PROBL_PMP1_H

