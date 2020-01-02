#ifndef __SUBGR_H__32
#define __SUBGR_H__32


#include "lpcol.h"


namespace bgn {                     // SS_BEGIN_NAMESPACE__


// The subgradient algorithm
class Subgr {
  int m; // the dim
  int k; // iteration
  Vector<double> s_0,s_1,s_2,s_3; // prev. subgradients


  bool CalcSubgr(double lpBnd,Column *bestCol);


public:
  /////////////////////////////////////////////////////
  /////////////// INPUT ///////////////////
  //___________________________________________________
  Vector<double> d, // the modified multipliers
    d0; // the initial ?? why?
  int iterMax;
  double rho, rho_k;
  d_vec b; // the rhs. reduced by branching?
  Vector<signed char> vs; // valid sign of d[i]


  static double iterMaxRatio;
  static double weight1;
  static double rhoDef;
  static double normToThe;
  static double outputLevel;
  static opt::OptSection opt;
  static opt::OptContainer Options();


public:
  Subgr() :rho(rhoDef) { }
  void SetRHS(d_vec &b_,Vector<signed char> &vs_)
    { b=b_; vs = vs_; }
  void Init(const Vector<double> &,const int);
  bool UpdateMultipliers
    (double lagrval,double lpbCurr,Column *bestCol);
};


}                                                // SS_END_NAMESPACE__


#endif // __SUBGR_H__32
