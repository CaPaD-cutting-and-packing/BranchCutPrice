#ifndef __PROBL_MCSP1_H
#define __PROBL_MCSP1_H

/**
  Header file for 1D-MCSP, i.e. with multiple stock lengths
*/

#include "probl_csp1.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

class MCSP1 : public CSP1 {
public:
  const char * TypeName() { return "MCSP1"; }
  /// (DIS-)ABILITY DEFINITIONS
  virtual bool CanBeInfeas() { return true; } // ..
//protected:

struct Stock {
  size L;
  int B; // upper bound
  double R; // the price
  float w; // weight for sorting. Init before!
//  int i_; // Number in the input
  int i0; // Number before sorting
  bool operator > (const Piece & pc2) const
    { return w > pc2.w; } // for sorting
}; // struct Piece;

  typedef Vector<Stock> StockContainer;

// Algorithms: //  auto_ptr<BB> bb;
//  BB bb; // no auto_ptr, because of object copying

// public:
  size Lmin, Lmax;
  double Rmin, Rmax; // material prices
  int M, M0;
  StockContainer sL,sL0;

	double GetVEps() { return max(1e-6, GetLPBnd() * deps); /*1*eps*/ } // ???
  double GetDEps() { return deps; }
  double GetBBEps() { return bb_eps; }
  double GetXEps();

  virtual double raster_ceil(const double v)
  { return ceilEps(v, deps); }

//public:
  CSP1(Env *penv,
    const char *in,const int i,const char *n);
  CSP1(int ,Env *penv, // this just for test
    const char *in,const int i,const char *n);
  virtual ~CSP1() { }
  virtual void Init(); // Check consistency + initialized?
//  virtual void Done();
  virtual void PrintLog(ostream &);
  virtual int Read(istream&,long);
  virtual void InitProblem();
  virtual void HeuristicLPBasis(ColumnList &bas)
  { FFSBasis(pc, bas); }
protected:
  virtual void FFSBasis
    (const PieceContainer &,ColumnList &); // ???
  virtual void GenColPure(ColSet& cs,const d_vec& d);
  virtual void GenColWithCuts(ColSet& cs,const d_vec& d);
  virtual void GenColWithHP(ColSet& cs,const d_vec& d);
  virtual void FillSigns(); // called from Init()

public:
  virtual void GenCol(ColSet& cs,const d_vec& d);

  /// Branching on Hyperplanes
  virtual void BrOnHP(const Vector<Column*>& colpool,
    const Vector<ColId>& cols, const Vector<double>& lpx,
    LPCut* &hpU, LPCut* &hpL);
  /// BRANCHING on VRP variables:
  Vector<Vector<float> > xVRP;

//  virtual void CalcLagrRel(const d_vec &);
//  virtual void UpdateLagrBnd(const d_vec &);

  virtual int Dim() { return m; } // ???
  virtual int the_m() { return m; }

  virtual void PrintColumn(ostream&,Column*);
  CSPLevelCut levelCut;
  LPCut * GetLevelCut();

////////////////////////////  ROUNDING + INTEGER SOLUTION
/// Input: see Problem
/// Output: see Problem
//protected:
  Vector<double> xi;
  double zr;   // obj
  Vector<double> xf;  // frac parts
  Vector<double> br;  // reduced rhs
  Vector<int> fracused; // indices of rnd-up cols
  int nRPE;   // n res pr ext
  static int RPEMax;

  virtual void GenColForHeur
    (Column *col,const d_vec &d, Vector<double> &b);
//  virtual void ConvertLPCols();
//  virtual bool TestIntegrality();
  virtual void SaveRoundedPart();
  virtual bool Round();
  virtual void GetRHS(Vector<double> &br);

  /// STATISTICS:
  virtual void InitStat();
  virtual void PrintStat(ostream&);

  static Vector<double> sta;

  //namespace {
class SVC: public Alg {
protected:
/////// INPUT via constructor:
  CSP1 * pr;
public:
  /*const*/ Vector<double> & b0;
  /*const*/ Vector<double> & d0;
  int m;
protected:
  const double lb, ub; // a lower/upper bnd
//////// VARS:
  Vector<double> bc;  // the remaining demands
  d_vec d;
  unsigned int kkk; // total iteration number
  int k; // iteration number
  double zz; // current value
  double result; // best value
  int xMin; // intensity of current pattern
  size waste;
  double SW;  // randomizer


  typedef Vector<Pattern> PatArray;
  PatArray pat;
  Vector<double> x;
public:
//////////// OUTPUT:
  PatArray patBest;
  Vector<double> xBest;


  SVC(CSP1* p_, Vector<double> & b_,
      Vector<double> & d_, const double lb_,
      const double ub_);
  virtual double Run();
protected:
  virtual void Execute();
  virtual void ConvertValues();
  virtual bool GenPat();
  virtual void CorrectValues();
  virtual void ControlSolution();
  virtual void SaveSolution();
  virtual void FillLastPattern();
  static opt::OptContainer Options();
  static opt::OptSection opt;
public:
  static int iterMax;
  static double outputLevel;
  static double patternUseRatio;
  static double pow_l;
}; ////////////////////// class CSP1::SVC

  friend class SVC;
  virtual bool CompleteIntSol();
  virtual bool VaryResProblem();
  virtual void AddResidualSolution(SVC & svc);
  virtual void ControlSolution();

public:
  virtual bool ConstructIntSol();
  virtual void PrintProblem(ostream &);

//private:
  static opt::OptContainer Options();
  static opt::OptSection opt;
  static bool fModelEquality;
  static bool fSortPieces, fMergePieces;
  static bool fEffectiveL;
  static int nStartBasis;
  static double nStepsMin0;
  static double nStepsMinInc;
  static bool fTestMSVC, fSP2Relaxation;
  static double deps, bb_eps;
  static double outputLevel;
  static bool fCompareBySPPRC_DP;

class MSVC: public SVC {
  // Always register max N open, max spread and min diff
  // for spread, also N different patterns intere.
  // but not VERY intere. because measure in intensities
  // BUT MINIMIZING N __FURTHER__ OPEN
  // Save Pareto-opt
  // ms = min spread, mos = min open stacks
  // + the way values are corrected? ie the weights?
  // a couple of initial runs w/o modified values
  // to estimate starting
protected:
  Random rnd;
  double rndCM; // then di in [1/rndCM*di,rndCM*di] uniformly
  d_vec dd; // the 'old' values
  int nPatLast;
  double h; // sum li b'i
  d_vec spread;
  int nOpen; // in and after a pattern
  int nStillOpen; // after a pattern
  int nOpenMax;
  double spreadMax;
  double nOpenMaxTarget, spreadMaxTarget;
  int nOpenMaxNext; // when restricting
  double mosR, msR;
  BBMos bb;
public:
  double nOpenMaxMin, spreadMaxMin, zzMin, nDiffPatMin;
  struct Solution {
    Solution() :no(0), ns(0), nz(0), nd(0), iter(0) { }
    double no, // N OPEN
      ns, // SPREAD
      nz; // MATERIAL INPUT
    // keeping all 3 indicators
    double nd; // N DIFFERENT PATTERNS
    int iter; // the iteration where it was found
    bool operator<(const Solution& s) const {
      return ((nz<s.nz)?1:
      ((nz==s.nz)?
        ((no<s.no)?1:
      ((no==s.no)?
        ((nd<s.nd)?1:0)
        :0)
        )
        :0)
        );
    }
  } sol0; //sOS, sS, sZ;
  typedef set<Solution> SolContainer;
  SolContainer sols;

  MSVC(CSP1* p_, Vector<double> & b_,
      Vector<double> & d_, const double lb_,
      const double ub_) : SVC(p_,b_,d_,lb_,ub_) {
    nOpenMaxMin = 1e100;
    nDiffPatMin = 1e100;
    spreadMaxMin = 1e100;
    zzMin = 1e100;
  }
//  virtual double Run();
protected:
  virtual void Execute();
//  virtual void ConvertValues();
  virtual bool GenPat();
  virtual void GenRestrCol(Column *col,const d_vec &d,Vector<double> &b);
  virtual void CorrectValues();
//  virtual void ControlSolution(); // to check&save every new
//  virtual void SaveSolution();
//  virtual void FillLastPattern();
  virtual void CheckSolution();
  virtual void CheckPattern();
  virtual void StartNewPlan();
  static opt::OptContainer Options();
  static opt::OptSection opt;
public:
  static int iterMax;
static bool CSP1::MSVC::fRestrictNOpen;
static int CSP1::MSVC::restrictDelay;
static int CSP1::MSVC::dNOpen;
static int CSP1::MSVC::nOpenMaxInitial;
  static double patternUseRatio;
  static double outputLevel;

  static double pow_l;
  static int weighScheme;

  static double nObjective; // 0=open stacks min,
     // 1=spread min, between?
  static double msReduc;
  static double msP;
  static double msR__;
  static double msRR;
  static double mosReduc;
  static double mosP;
  static double mosR__;
  static double mosRR;

// randomization for each plan?
  static double rndCMCM; // rndCM in [1,1+rnd*(rndCMCM-1)]
}; ////////////////////// class CSP1::MSVC

  friend class MSVC;


// Sequencing SVC
class SSVC: public MSVC {
    // a couple of initial runs w/o modified values
  // to estimate starting
public:
  //////////////// INPUT: ////////////////////
  //// The solution to be permutated.
  /// No different patterns inside
  Vector<Pattern> pat0;
  d_vec x0;
  //////////////// OUTPUT: //////////////////
  /// sols of MSVC
  SSVC(CSP1* p_, Vector<double> & b_,
      Vector<double> & d_, const double lb_,
      const double ub_) : MSVC(p_,b_,d_,lb_,ub_) { }
protected:
  int iNext; // th next pattern
  double zkpBest; // = d*pat0[iNext]
  i_vec fUsed;
  Random rndStart; // the starting pattern
//  virtual double Run();
protected:
  virtual void Execute();
//  virtual void ConvertValues();
  virtual bool GenPat();
//  virtual void GenRestrCol(Column *col,const d_vec &d,Vector<double> &b);
//  virtual void CorrectValues();
  virtual void ControlSolution() { } // fn empty -- permut
//  virtual void SaveSolution();
//  virtual void FillLastPattern();
  virtual void CheckSolution();
//  virtual void CheckPattern();
  virtual void StartNewPlan();
//  static opt::OptContainer Options();
//  static opt::OptSection opt;
public:
}; ////////////////////// class CSP1::SSVC

  friend class SSVC;


};//___class_CSP1_______________________________________


}                                                // SS_END_NAMESPACE__

#endif // __PROBL_MCSP1_H
