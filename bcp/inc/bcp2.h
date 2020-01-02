#ifndef __BCP2_H__32
#define __BCP2_H__32

#include "bcp.h"
#include "probl_pmp1.h"


namespace bgn {                     // SS_BEGIN_NAMESPACE__

//#undef OUTP_LEV__
//#define OUTP_LEV__ (0>=OutputInterval?1:   (0==cntNode%kNodeOutput) ? FMin(outputLevel,opt::GlobalOutputLevel()) : 0)

class BCP2 : public BCP {
public:

// BOUNDS:

// RESULTS && STATISTICS:

// VARIABLES:
     // INIT THEM!!!
  double lpzfAdd; // the addition for lower-bounded cols
  bool fCSPOpt,
    fBPPOpt;
  PMP1 * pr;

/////////////////// LP: //////////////////////

/////////////////// CUTS: ////////////////////

/////////////////// BRANCHING: ////////////////////

// SERVICE

// CALLS:
public:
  BCP2(Problem *p_) : BCP(p_),
    pr(dynamic_cast<PMP1*>(p_)) // geil
  { }
  virtual ~BCP2() { }
//  virtual void Run();
  virtual char * Version();

protected:
  virtual void SolveCSP();
  virtual void SolveBPP();
  virtual bool CheckData();
  virtual pair<int,double> SelectBrVar_MostInfeas();
  virtual void UpdateVarBounds(); // here change coefs?
  virtual double GetLPValue() { return lp->GetValue() + lpzfAdd; }

//  virtual int LocalUpperBound(Column *); // changed?
  virtual bool Rounding(int fFast=0); // no MIP, feas != opt

  virtual void InitRun(); // invoking CSP& BPP
  virtual void InitOptimize();
  virtual void InitBounding(); // init CPA

  virtual void DoneRun();

//  virtual void PricingOver(); //playing with zf in root
//  virtual void InitBasis(); // which?
  virtual void PrintIter();

  virtual void PrintLog();
  virtual void InitLog1();
//  virtual void PrintContSol();
//  virtual void PrintSolutions();
  virtual void PrintStatistics(ostream &);
  static double NCSP, ND0, ND;
  static int nFeas;


// ALL ADDI OPTIONS IN PMP1
  static int nVarSel; // var selection

//  static void PrintOptions(ostream&); // virtual ?
    static opt::OptContainer Options();
  static opt::OptSection opt__;
  static double outputLevel;
};

}                                                // SS_END_NAMESPACE__

#endif // __BCP2_H__32

