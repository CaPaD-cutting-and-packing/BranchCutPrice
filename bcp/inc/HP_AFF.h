#ifndef __HP_AFF_H
#define __HP_AFF_H

/** The branching hyperplanes
  corresponding to bounds on variables
  of the Arc Flow Formulation of 1D CSP
*/

//#include "pool.h"
#include "lpcol.h" // column
#include "cuts.h"
//#include "branch.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

/// Branching:
class BrAFF : public LPCut {
public: // private:
  /// 1-based product indexation:
  const int * l; // product lengths IN 0..m-1 INDEXATION
  const int i, p; // product i (INDEX. 1..m) at position p
  const int rhs; // the bound;
  const int slackSign; // defines "<=" or ">="

public:
  int Type() const { return -2; } // negative: for hyperplanes
  bool IntegerSlack() { return true; } // alfa == 1
  virtual bool CanBeDeleted() const { return false; }

  double GetRHS() const { return rhs; }
  double CalcRHS__ (d_vec&) {return rhs;}

  BrAFF(const int *L, int I=0, int P=0, int RHS=0, int SLACKSIGN=1) :
    l(L), i(I), p(P), rhs(RHS), slackSign(SLACKSIGN)
    { }
///////////////////// CALCULATION ////////////////////////
  int GetSlackCoef() const { return slackSign; }
  double Calc__(Column *);

/// Nothing for col.gen. by B&B:
  void CalcIntermSums(Column *) { }
  // Utilities for col gen:
  /// Obj func calculation:
  double CalcUsingIntermSums() { return 0; } // Calc Pis

  void Print(ostream& os=cout) {
    os << " xAFF["<< i << ',' << p << ']'
      << ( slackSign>0 ? "<=" : slackSign<0 ? ">=" : "=" )
      << rhs;
  }
  bool operator<(const LPCut& c) const {
    if (Type() != c.Type())
      return Type()<c.Type();
    const BrAFF& php = dynamic_cast<const BrAFF&>(c);
//    return i<php.i ? true : i==php.i ? j<php.j : false;
// this is not enough!! Because may be diff rhs
    if (i<php.i) return true;
    if (i>php.i) return false;
    if (p<php.p) return true;
    if (p>php.p) return false;
    if (rhs<php.rhs) return true;
    if (rhs>php.rhs) return false;
    if (slackSign<php.slackSign) return true;
    assert(slackSign != php.slackSign);
    return false;
    /* : i==php.i ? j<php.j ? true
      : j==php.j ? rhs<php.rhs ? true
      : rhs==php.rhs ? slack : false; */
  }

///////////// STORING COEFS: ////////////////
/////////////////////////////////////////////
  double Calc__(Column *c, int iCol)
  { return Calc__(c); }
  double GetCoef(int iCol) { assert(0); return 0; }
  double CalcSlackValue
    (i_vec &iNZ, d_vec &xNZ, map<LPCut*,double> &slVal)
  { assert(0); /* no pool restoration */ return 0; }
};// class BrAFF _______________________________


}                                                // SS_END_NAMESPACE__

#endif // __HP_AFF_H
