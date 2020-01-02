#ifndef __CUTGEN_H
#define __CUTGEN_H

#include "pool.h"
#include "lpcol.h" // column
#include "cuts.h"
#include "branch.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

//class CutGen : public LPCut {
//};//____________________________________________________

class CSPLevelCut : public LPCut {
public:
  int Type() const { return -1; }
  int is_const() const { return 1; }
  bool IntegerSlack() const { return true; }
  double ofc; // const terms for col gen
  double lhs; // the level
  double Calc__(Column*c) {
    return 0;/*c->IsSlack() ? // Calc own slack:
      (c->GetCutSlackCut() == this ?
      c->GetCutSlackCoef():0): c->GetObj();*/
  }
  void CalcIntermSums(Column * c) { ofc = Calc__(c); }
  // Utilities for col gen:
  double CalcUsingIntermSums() { return ofc; }

  int GetSlackCoef() const { return -1; }
  double GetRHS () const { return lhs; }
// for the foloowing, lhs must be pre-set:
  double CalcRHS__(d_vec&) { return GetRHS(); }
  void CalcConstTerms(Column *c) { ofc = c->GetObj(); }
  void AssignConstTerms() { }

  double CalcApprCoef(int k) { return 0; }
  double CalcApprErrorL() { return ofc; } // for MCSP
  double CalcApprErrorU() { return ofc; }

  void Print(ostream&os)
  { os << "Level cut: rhs="<<lhs; }
  bool operator<(const LPCut& ) const
    { return true; } // level cut always first
  double CalcSlackValue
  (i_vec &iNZ, d_vec &xNZ, map<LPCut*,double> &slVal)
  { assert(false); return 0; }
};

// SE : Slack elimination
// SA : SuperAdditive

class SACutSE : public LPCut {
public:
  /// Dependence on original constraints:
  Vector<double>
    u, uSE;
  /// Dependence on other SuperAddi cuts:
  struct Dependence {
    LPCut * c;
    mutable double u, uSE;
    Dependence(double u1,double u2, LPCut*cg)
      : u(u1), uSE(u2), c(cg) { }
    Dependence(double u1, LPCut*cg)
      : u(u1), uSE(0), c(cg) { }
    Dependence(LPCut*cg)
      : u(0), uSE(0), c(cg) { }
    bool operator<(const Dependence& d) const
      { return c < d.c; }
  };
  typedef list<Dependence> DepContainer;
  DepContainer dep; // remove with zero weights
  Pool<Dependence> invCuts; // involved cuts

  /// The Alfa of the superadditive fct.
  double alfa; // =0 for CutSet
  int cutType; // 0: MI cut, 1: CG cut
  int fSE; // slack elim

//  int nUse; // How many use it

  double
    sua, suaSE,  // the intermediate sum for col gen
    sua0, suaSE0; // its ini. value
  double
    value, // some value for recursive calc.
    value2;// another one for alfa U/L
  double rhs;

  /*static*/ double F_Alfa(/*double alfa,*/double v);
  /*static*/ double F_Bar(/*double alfa,*/double v);
  static double
    F_Alfa(double alfa,double v,int cutType);
  static double
    F_Bar(double alfa,double v,int cutType);

//  bool fNodeVisited;

public:
  int Type() const { return 1; }
  bool IntegerSlack() { return (cutType!=0); } // alfa == 1
     // then for all cuts must be == 1 what means
     // Gomory fractional

  double GetRHS () const { return rhs; }
  double CalcRHS__ (d_vec&);

  SACutSE() :
    value(0), value2(0), alfa(0), cutType(0), fSE(0)
    { }
///////////////////// CALCULATION ////////////////////////
  void ClearNonRec(); // non-recursive
  double CalcCutSlackCoef__(LPCut *pcut);
  double GetCutSlackCoef__(LPCut *pcut);
  double Calc__(Column *);
  void CalcIntermSums(Column *);
  // Utilities for col gen:
  /// Obj func calculation:
  double CalcUsingIntermSums(); // Calc Pis
  /// Bound calculation:
  double CalcApprCoef(int k),       // Calc uAppr
         CalcApprErrorL(),
         CalcApprErrorU();

  void ClearSums();
  void AddToSums(int k, int x), Sub1FromSums(int k);

  void CalcConstTerms(Column *),
       AssignConstTerms();

  void resize(int m) { u.resize(m); uSE.resize(m); }

  void Print(ostream& =cout);
  void Number(int &nn);
  void ProduceListOfInvolved(Vector<LPCut*> &) const;

  void Sort() { dep.sort(); }
  bool operator<(const LPCut&) const;

///////////// STORING COEFS: ////////////////
/////////////////////////////////////////////
  Vector<double> rawCoef;
  int GetNCoefs();
  double Calc__(Column *c, int iCol);
  double GetCoef(int iCol);
  double CalcRawCoef(Column *c,int iCol);
// Dependence on variable bounds:
  typedef Pool<VarBnd> BndContainer;
  BndContainer bnds;
  double CalcSlackValue
    (i_vec &iNZ, d_vec &xNZ, map<LPCut*,double> &slVal);
};//____________________________________________________


class CutSet { // passed to col gen
public:
  typedef SACutSE::DepContainer DepContainer;
    // -- simpl. multipliers inst. of u's
  DepContainer dep;
  CutList invCuts;
  // involved cuts, detected by AddInvolved() when
  // called with a set of all possible

public:
  template <class vec>
    void CalcCoefs(const vec&); // ???
  // Utilities for col gen:
  /// Obj func calculation:
  double CalcCoefsUsingIntermSums();
  /// Bound calculation:
  double CalcApprCoefs(int k),
  // -- these return corresp. coefs of the cuts
  // multiplied by weights
         CalcApprError();

  void ClearValues();
  void PrepareRecursion();

  void AddToSums(int k, int x), Sub1FromSums(int k);

  void CalcConstTerms(Column*),
       AssignConstTerms();

  void ReduceWithZeroWeights();

  // Produce the exact list of cuts involved in col gen
  // as a subset of invCuts
  void AddInvolved(CutList *);

};//____________________________________________________

}                                                // SS_END_NAMESPACE__

#endif // __CUTGEN_H
