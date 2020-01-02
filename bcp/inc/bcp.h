#ifndef __BCP_H__32
#define __BCP_H__32

#include "pool.h"
#include "lp.h"
#include "lpcol.h"
#include "cutgen.h"
#include "branch.h"
#include "node.h"
#include "subgr.h"
#include "problem.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

// CHECK CUT COEFS
// #define CHECK_CUTS // Now option

#define N_RUNS_MAX_DEF 1 //20 -- do it later, rnd...
#define DEF_LP_ONLY false

//#undef OUTP_LEV__
//#define OUTP_LEV__ (0>=OutputInterval?1:   (0==cntNode%kNodeOutput) ? FMin(outputLevel,opt::GlobalOutputLevel()) : 0)

class BCP : public Alg {
public:
  auto_ptr<Problem> pr;
  auto_ptr<LP> lp, mip;
  auto_ptr<Subgr> subgr;

// BOUNDS:
  double
    gub, // global upper bound
    glb, llb, // global/local lower (LP) bound
    glb_others, // lpbound of all other active nodes
       // (glb = min(glb_others,llb))
    glb0, // initial glb
    clb, // current lower bound (not nec. true)
    gLb, lLb // global/local Lagrange bound
  ; // CORRESPONDING LP VALUES WHO IMPLY THESE:
  double glv,llv,glv0,clv,glv_others,gLv,lLv;
  double llbLast, llvLast; // in last CPA iteration
  int iterLLBLast;

// RESULTS && STATISTICS:
  enum StatusValues
  { none, opt, LPBnd, feas, LPErr, infeas, error };
  StatusValues status;
  static const char * statusName[];

  int lpCols, ipCols;

  Timer timer;
  MyTime timeLP0, timeIP,
    timeIPBest, // to find the best int sol
    timeLPBest;
  Timer tmCG1, tmCG2, tmRnd;
  MyTime timeCG1, timeCG2, timeRnd;
  bool fMIPSolBest; // whether MIP solver made it best
  int nTooLongColGen; // in InitFile() !!! or caller acc.
  double NnvAll, NnvTst; // NewValue() tests (red.cost bnd)
  int nErr2; // non-fatal errs

  int nextRNDNCols;
  double nCG, nCG_TOff;

// VARIABLES:
  bool firstinstance; // processed in a file
  bool fSubproblem;
  char outfile1[1024];
  char outfile2[1024];

  int nRuns;

/////////////////// LP: //////////////////////
// BETTER INCAPSULATION
  int Dim() { return pr->Dim() + (int)cuts.size(); }

  typedef Pool<Column> ColContainer;
  Pool<Column> allcols;
  typedef Vector<Column*> ColMainPool;
  Vector<Column*> colpool;
  typedef Vector<ColId> ColIDPool;
  Vector<ColId> cols; // the local cols in the LP
  Vector<LP::ColStatus> cstat, cstat__;
  Vector<double> lpx, lpx__;
  d_vec redCosts; // get when opt sol.
  d_vec lpd; // last simplex multipliers
  bool fLPOpt, fLPOpt0; // whether correct LP value was g.
    // lagrange bound not helped
  bool fLPIntOpt; // whether last LP integer and locally optimal
  bool fInfeasBasic; // whether an infeas slack is in
  bool fInfeasColsIn; // whether infeas cols are not set to 0
  double maxObjCoef;
	double minRC; // minimal reduced cost after pricing

  i_vec matind;
  d_vec matval; // for adding cols/cuts
  ColSet newCols;
  int Nnv; // number of new values;
	
	///////////////////////// MIP: ///////////////////////////
	Vector<Pattern> patMIP;
	int nMIPColsLast;                       // when solved

/////////////////// CUTS: ////////////////////
  typedef set<LPCut*,CmpPtrByVal<LPCut*> > CutRefSet;
    // used after del nodes ?
  CutRefSet cutpool; // all cuts known for the node
  CutList cuts; // all those added to the formulation,
    // bad: better incapsulate
  Vector<LPCut*> invCuts; // involved cuts (util)
  list<CGIV> cutViol;
  Vector<bool> fNTD; // if not to del cut

  int iter; // N of iter-s with ColGen
  int iterAll; // all iters, also w/o colgen
    // (simple pool reopt. with new cuts)
  int sum_iter, sum_iterAll;
	int nColsPrintLPVal;

  int delResetIterNext; // when to del all untight cuts
  double  delReset;
  Random rndCuts;
  bool fLevelCut; // turned off
  int maxCuts; // the max number of PURE cuts, not branching HP
  int MaxAllCuts() { return maxCuts + nHP; }

/////////////////// BRANCHING: ////////////////////
  int nHP; // number of BRANCHING hyperplanes (==depth!!!) added
  SolutionTree nodes;
//  list<VarBnd> ubnds; // all upper bounds
  SolutionTree::node_iterator theNode,
    lastNode,
    nextNodeDFS, // the left/right son of the lastNode,
       // if lastNode was split
    sonRight, sonLeft; // Right=Lower/Left=Upper bnd
  bool fLastNodeSplit;
  list<SolutionTree::node_iterator> backtrack;
  bool fReoptTree;
  int cntNode, // N selected nodes
    cntInitNode; // N nodes initialized for bounding proc.
  int depth;
  double nLeft, nRight; // summed skew values of all initialized nodes
  i_vec lb, ub; // the var. bounds set in the node.
    // However, not for new cols generated there

  Random rndBr;
  Random RndBrFrac; // random fractional part

// SERVICE
  static volatile bool fOutputTimer; // set by interrupt
  static void AlarmHandler(int =SIGVTALRM);
  Pool<Column> givsol; // a help solution
  double givsolobj;

// CALLS:
public:
  BCP(Problem *p_);
  virtual ~BCP() { }
  virtual void Run();
  virtual char * Version();

protected:
  virtual void Optimize();
  virtual bool Stop();
  virtual bool Bounding(); // true if the node to fathom
  virtual bool SolveLP();
  virtual void StartReoptTree();
//  virtual bool StopLP();
  virtual bool CPATailingOff();
  virtual bool StopCPA();
  virtual void Branch();
  virtual void SaveBasis();
  virtual void CopyBasTo(Node * );
  virtual void Create2Sons();
  virtual void BranchOnVar();
  virtual void BranchOnHyperplanes();

  // return index/ub value of a branching var (lb=ub+1)
  virtual pair<int,double> SelectBrVar_MostInfeas();
  virtual pair<int,double> SelectBrVar_PsCosts();
  virtual void InitPsCost(int iCol, bool upper, double x);
  virtual double PseudoCost(Column *c, double f);
  virtual void UpdatePseudoCost(); // After solving node

  virtual void Fathom(SolutionTree::node_iterator );
  virtual bool Select();
  virtual bool SelectNextNode();
  virtual bool ContradictionsInNewNode();
// + clearance? or in Branch(); or Select();?
  virtual bool InitNewNode();
  virtual void CompileColumns();
  virtual void CompileVarBounds();
  virtual char TightenVarBounds(); // returns 0 if contrad.
  virtual void CompileCuts();
  virtual void UpdateVarBounds();
  virtual void UpdateCuts();
  virtual void RestoreBasis();
  virtual void RecalcGLB();
  virtual void DeleteFathomed();

  virtual int LocalUpperBound(const Column *);
	virtual bool Rounding(int fFast=0);
	virtual bool RoundingMIP();
	virtual void ConvertLPCols(Vector<double> &lpx,
    Vector<Pattern> &pat,Vector<double> &patx);
	virtual void StoreMIPPatterns(Vector<double>& xMIP, Vector<Pattern> &pat, Vector<double> &patx);
  virtual bool FixAndSet();
  virtual bool NewValues();
  virtual void Fix(); // if glb or gub changed
  virtual void Separate(); // we do this cycle always
  virtual bool Optimum();
  virtual bool LocalOptimum();
  virtual void InitRun();
  virtual void InitOptimize();
  virtual void InitRootNode();
  virtual void InitBounding(); // init CPA
  virtual void DoneBounding();
  virtual void DoneOptimize();
  virtual void DoneRun();
  virtual bool NextOptimize();

  virtual void RecalcCutsRHS();
  virtual void ClearCoefs();
  virtual void  MarkTightCuts(); ////////////////////////
  virtual void DeleteSomeCuts(int mode);
  virtual void DelCut(int iCut, bool =false);
  virtual bool FindViolCutsInPool();
  virtual void ConstructNewCuts();      ////// + level cuts
  virtual bool ModifyAndResolveLP(); // may be needed - RCC
  virtual void RestoreLP();
  virtual void ConstructSACut
    (Vector<double> &bii, double xfi,
      Vector<LP::ColStatus> &cstat, d_vec &lpx);
  virtual void AddLevelCut(); // always ?
  virtual void AddSomeCuts(); // FROM THE POOL
  virtual void AddCut(LPCut *pcut);
  virtual void PrintCuts();
  virtual void PrintNCuts();
  virtual void CheckCuts();

  virtual double GetLPValue() { return lp->GetValue(); }
  virtual void SetCurrentLPV(double v);
  virtual void SetLocalLPV(double v);
  virtual void UpdateLocalLPV(); // from current
  virtual void UpdateGlobalLPV();
  virtual double GetLocalLagrBnd();
  int GetColIndex(int i);
  Column * GetCol(int i);
  Column * GetMainCol(int i);
  bool IsSlackCol(int i); // the cols of LP
  virtual void InitLP();
  virtual bool SolvePrimal();
  virtual bool SolveDual();
  virtual bool SolveDual_GetMultsOnly(); // for PsC init
  virtual void GetLPSolution();
  virtual bool StartPhaseI();
  virtual void EndPhaseI();
  virtual bool Price(); // Subgradient Optimization when no cuts
  virtual bool GenCol(d_vec& d, bool fSubgr=0);
  virtual bool AddCols();
  virtual void PricingOver();
  virtual void InitBasis();
  virtual void AddSlacks();
  virtual bool AddColToLP(LP * , Column * c, int j, bool fMIP = false);
  virtual void DelCol(int i);
  virtual bool AddCols(ColSet &cs);
  virtual bool AddCols(ColList &cl);

  virtual void ReadGivenSolution();
  virtual void CheckGivenSolution(); // + cuts?
  virtual void PrintIter();
  virtual bool TimeLimit();
  virtual void PrintLog();
  virtual void InitLog1();
  virtual void PrintContSol();
  virtual void PrintSolutions();
  virtual void InitLog2(); // None ?
  virtual void InitLog3();
//  void PrintOptions(ostream &ofs);


// OPTIONS:
public:
  virtual void PrintStatistics(ostream &os) {
    pr->PrintStat(os);
  }
  static opt::OptContainer Options();
  static opt::OptSection opt__;

  static double TimeLimit__;
  static bool fNoOpt;
  static int nRunsMax;
  static int maxDepth;

  static int nBFSDFS;
  static double nDive;
  static double BrStratRL;
  static bool fCGTailOff;

  static int nBranch;
  static float BrVarFrac;
  static bool fBranchPsCosts;
  static double dInfeasFracPart;
  static double Alfa1; // =2
  static double wTillUB;

  static bool fRedCostBnd;

  static bool fLocalReduce;
  static bool fLocalUB;
  static int nReducePool;
  static int kNodeNewCuts;
  static double MIPTiLim;
  static int RNDInterval;
  static int RNDNodeInterval;
  static int RNDColInterval;
  static int MIPInterval;
  static int MIPNodeInterval;
  static int kNodeOutput;
  static bool fPrintContSol0,
    fPrintContSolBest, fPrintIntSolBest,
    fPrintProblem, fPrintUnsolved;
  static bool fOriginalRedCost;
  static bool fUseLagrange;
  static double M__;
  // static double fInitBasis;
  static double outputLevel;
  static double outputInterval;

  static bool fLPOnly; // ----------------------------

  static bool fCheckGivenSol;

// CPA options
  static opt::OptContainer OptCuts();
  static opt::OptSection optCuts;

  static int firstCutAfterNode;
  static int maxCuts__;
  static int iterMax;
  static int iterRootMax;
  static int iterTailOff;
  static bool fSkipColGenAtLPBound;
  static int cutType;
  static int nPoolIterMax;
  static int nIterkMax;
  static int cutSelect;
  static int CGNormMax;
  static int CGFracParts;
  static double rndViolDev;
  static int cntDel0;
  static double cntDel0inc;
  static double delReset0;
  static double delResetInc;
  static bool fCheckCuts;  // ------------------------

// OBSOLETE:
  static int fSE; // ------------------------
  static int fPoolSearch;  // ------------------------
//  static int nMarkedMax;
  static bool fAddLevelCut;  // ------------------------

  static void PrintOptions(ostream&);
};

}                                                // SS_END_NAMESPACE__

#endif // __BCP_H__32
