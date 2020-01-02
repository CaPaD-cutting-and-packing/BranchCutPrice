#ifndef __LP_H
#define __LP_H

// An abstract LP solver interface

namespace bgn {                     // SS_BEGIN_NAMESPACE__

//  enum LPStatusValues
//  { none, ok, infeas, unbounded, LPCYCLE, error };
const char* const LPStatusName[] =
  { "none", "opt", "feas", "infeas", "unbnd", "error" };

/// Abstract LP Solver wrapper e.g. for Soplex, CPLEX usw.
class LP {
public:
  virtual void Open() =0;
  virtual void Close() =0;
  virtual void SolvePrimal() =0;
  virtual void SolveDual() =0;

  enum StatusValues
  { none, opt, feas, infeas, unbounded, error };
  StatusValues status;
  enum ColStatus { atLower, basic, atUpper, free_super };

///////////////////////////////////////////////////////
//////////// The constructor //////////////////
////// should set minimization sense to min
///// s=0: CPLEX
  static LP * CreateLP(int w=0);

  int NRows() { return Dim(); }
  virtual int Dim()=0;
  virtual int NCols()=0;

  virtual ~LP() { }

  virtual void AddCol(double obj,int nnz,
    int *cmatind,double *cmatval,double lb,double ub) =0;
  virtual void GetStatus(StatusValues &s) =0;
  virtual double GetValue() =0;
  virtual void GetMultipliers(Vector<double> &d) =0;
  virtual void GetRedCosts(d_vec &d) =0;
  virtual void InitConstr(Vector<double> &b, const char* sense=NULL) =0;
  virtual void GetSolution
    (Vector<ColStatus> *cs,Vector<double> *x) =0;
  // x.size() must be NCols() 'cause (x[i]!=0) => var non-int
  virtual void LoadBasis(Vector<ColStatus> &cs) =0;
  virtual double SolveMIP
  (Vector<double> &x, double tm, double UB=1e75) =0;
  virtual void GetBasisInverse
  (Vector<Vector<double> > &bi,Vector<double> &xx) =0;
  virtual void AddRow // equality
    (int nnz,int *rmatind,double *rmatval,double lub)=0;
  virtual void DelRow(int ir) =0;
  virtual void DelCol(int ic) =0;
  virtual void ChangeRHS(int ir,double lub) =0;
  virtual void ChangeVarBnds(int ic,double lb,double ub) =0;
  virtual void GetVarBnds(int ic,double &lb,double &ub) =0;
  virtual void GetObjCoef(int ic,double & c) =0;
  virtual void SetObjCoef(int ic,double c) =0;
  virtual void GetFullCol
    (int j,Vector<double> &col,double &obj) =0; // +rhs?
  virtual double GetLPCoef(int i, int j) =0;

  virtual void WriteModel(char*) =0;
};//___class_LPWrp______________________________________


}                                                // SS_END_NAMESPACE__

#endif // __LP_H
