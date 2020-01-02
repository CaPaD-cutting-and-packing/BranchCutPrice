// subgr.cpp  -- the subgradient method
// Author: Gleb <Belov@math.tu-dresden.de>

#include "stdafx.h"
#include "subgr.h"
#include "lasthdr.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

// + Farley!
/// with cuts: not nec. opt col gen!

bool Subgr::UpdateMultipliers
  (double lgvCurr,double lpbCurr,Column *bestCol) {
  int i;

  rho_k *= rho; // a param initially from rhoDef
  double forLagr = lgvCurr;
    // CP22: No lagr bnd yet
  double stepLen
    = (lpbCurr // Value ?
      - forLagr);
  if (forLagr < -1e+50) stepLen = 1; // = 0;
  stepLen *= rho_k;
  assert(d.size() && s_0.size());

  if (not CalcSubgr(lpbCurr,bestCol)) // s_0 is zero
    return false;
//  PrintVec(cout,s_0); log_ln(""); cin.get();
  double nrm =(normToThe<2) ?
    Norm(&s_0[0],&s_0[0]+m) : Norm2(&s_0[0],&s_0[0]+m);
//  PrintVec(cout,s_0); log_ln(""); cin.get();
  for (i=0;i<m;++i) s_0[i] /= nrm;

/*  if (OUTP_LEV__>=4) {
    log_ln("current subgr:");
    PrintVec(cout,s_0); log_ln(""); cin.get();
  }*/

  if (0==k)
  { s_3=s_2=s_1=s_0; }

  double w_1 = 1 - weight1;

  for (i=0; i<m; ++i)

    d[i] += stepLen
      * (weight1*s_0[i]
      + (w_1/2)*s_1[i]
      + (w_1/4)*s_2[i] + (w_1/4)*s_3[i]);
  for (i=0; i<m; ++i)
    if (d[i]*vs[i] < 0) d[i]=0;

  s_3 = s_2; s_2 = s_1; s_1 = s_0;
  ++ k; // !!!!!!!!!!!!
  return true;
} //___UpdateMultipliers________________________________

bool Subgr::CalcSubgr(double lpbCurr,Column *bestCol) {
  d_vec col(m); // Filled with 0's ?
  fill_n(&col[0],col.size(),0);
  for_each_in(bestCol->id,iti,Column::iterator)
    col[iti->i] = iti->d; // uncompacting the sparse
  int i;
  for (i=0;i<m; ++i)
    s_0[i] = b[i] - lpbCurr // Current !!?
      * col[i];
  for (i=0;i<m; ++i)
    if (fabs(s_0[i]) > 1e-6)
      return true;
  return false;
}

void Subgr::Init(const d_vec & d_,const int m_)
{
  assert(b.size() && vs.size());
  m = m_;
  rho_k = 1; // rho init from rhoDef
  iterMax = int(iterMaxRatio * m);
	int s = (int)d_.size();
  d.resize(s); d0.resize(s);
  s_1.resize(s); s_2.resize(s); s_3.resize(s);
  d0 = d = d_; s_0.resize(d_.size()); k=0;
}

opt::OptContainer Subgr::Options() {
  opt::OptContainer oc;
  oc
    << opt::MakeOpt(&iterMaxRatio, 0.15,
      "iterMaxRatio", "Relative to the dimension. Actually this doesn't help since I changed something..."
			" And seems bad for DP-col gen")
    << opt::MakeOpt(&weight1, 1,
      "weight1", "0..1: the weight of the last subgr.\n"
      "'The subgr of the 3 prev steps get the rest (2:1:1)")
    << opt::MakeOpt(&rhoDef, 0.99,
      "rhoDef",
      "Multiply step len by rho^k (k=step)")
    << opt::MakeOpt(&normToThe, 2,
      "normToThe",
      "Divide subgr by its norm^normToThe; =1,2")
    << opt::MakeOpt(&outputLevel, DEF_OUTP_LEVEL,
      "outputLevel", "");
  return oc;
} //____________________________________________________

opt::OptSection Subgr::opt
  ("Subgr", "The subgradient procedure",
  Subgr::Options(), opt::SolverCfg(), 5000);

double Subgr::iterMaxRatio = 0.15;
double Subgr::weight1 = 1;
double Subgr::rhoDef = 0.99;
double Subgr::normToThe=2;
double Subgr::outputLevel=DEF_OUTP_LEVEL;

}                                                // SS_END_NAMESPACE__
