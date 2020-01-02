#ifndef __WRPCPLEX_H
#define __WRPCPLEX_H

/* Bring in the CPLEX function declarations and the C library
   header file stdio.h with the following single include. */

#include "ilcplex/cplex.h"

/* Bring in the declarations for the string functions */

#include <stdlib.h>
#include <string.h>


//ILOSTLBEGIN

#include "lp.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

class WrpCpx : public LP {
public:

  int status; // for all operations
//    static 
   CPXENVptr     env;
   CPXLPptr      lp;
//   int iScale; // scaling mode, init to 0
   d_vec lb, ub, theX; // getting inverse
   i_vec cstat;
	 
	 // PARAMS
	 static double epsOpt, epsFeas;
	 static int nMIPThreads;
	 static opt::OptSection opt;
	 static opt::OptContainer Options();
	 
  WrpCpx();
  ~WrpCpx(); // IMPORTANT: Cplex must be closed manually
  virtual int Dim();
  virtual int NCols();

  virtual void Open();
  virtual void Close();
void SolvePrimal();
void SolveDual();
void AddCol(double obj,int nnz,
  int *cmatind,double *cmatval,double lb,double ub);
void GetStatus(StatusValues &s);
double GetValue();
void GetMultipliers(Vector<double> &d);
void GetRedCosts(d_vec &d);
void InitConstr(Vector<double> &b, const char* pSense=NULL);
void GetSolution
  (Vector<LP::ColStatus> *cs,Vector<double> *x);
void LoadBasis(Vector<LP::ColStatus> &cstat);
// x.size() must be NCols() 'cause x[i] => var non-int
double SolveMIP
  (Vector<double> &x, double tm, double UB=1e75);
void GetBasisInverse
  (Vector<Vector<double> > &bi,Vector<double> &xx);
void AddRow // equality
(int nnz,int *rmatind, double *rmatval,double lub);
void DelRow(int ir);
void DelCol(int ic);
void ChangeRHS(int ir,double lub);
void ChangeVarBnds(int ic,double lb,double ub);
void GetVarBnds(int ic,double &lb,double &ub);
  virtual void GetObjCoef(int ic,double & c);
  virtual void SetObjCoef(int ic,double c);
void GetFullCol
  (int j,Vector<double> &col,double &obj); // +rhs?
double GetLPCoef(int i, int j);

void WriteModel(char*);
};

}                                                // SS_END_NAMESPACE__

#endif // __WRPCPLEX_H
