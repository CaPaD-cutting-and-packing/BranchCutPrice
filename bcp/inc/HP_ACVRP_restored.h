#ifndef __HP_ACVRP_H
#define __HP_ACVRP_H
/** The branching hyperplanes & cutting planes  corresponding to bounds on the variables  of the ACVRP-like formulation of 1D CSP*/
//#include "pool.h"
#include "lpcol.h" // column
#include "cuts.h"
//#include "branch.h"
namespace bgn {                     // SS_BEGIN_NAMESPACE__
/// Branching:
class BrVRP : public LPCut {public: // private:
/// 1-based product indexation:
  int i, j; // the arc, i:0..m, j:1..m
  int rhs; // the bound;
  int slackSign; // defines "<=" or ">="
  public:  int Type() const { return -1; } // negative: for hyperplanes
  bool IntegerSlack() { return true; } // alfa == 1
  virtual bool CanBeDeleted() const { return false; }
  double GetRHS() const { return rhs; }
  double CalcRHS__ (d_vec&) {return rhs;}
  BrVRP(int I=0, int J=0, int RHS=0, int SLACKSIGN=1) :    i(I), j(J), rhs(RHS), slackSign(SLACKSIGN)    { }
  ///////////////////// CALCULATION ////////////////////////
  int GetSlackCoef() const { return slackSign; }
  double Calc__(Column *);
  /// Nothing for col.gen. by B&B:
  void CalcIntermSums(Column *) { }  // Utilities for col gen:
    /// Obj func calculation:
  double CalcUsingIntermSums() { } // Calc Pis
  void Print(ostream& os=cout) {    os << "x["<< i << ',' << j << ']'      << ( slackSign>0 ? "<=" : slackSign<0 ? ">=" : "=" )      << rhs;  }
  bool operator<(const LPCut& c) const {    if (Type() != c.Type())      return Type()<c.Type();    const BrVRP& php = dynamic_cast<const BrVRP&>(c);//    return i<php.i ? true : i==php.i ? j<php.j : false;// this is not enough!! Because may be diff rhs
   if (i<php.i) return true;    if (i>php.i) return false;    if (j<php.j) return true;    if (j>php.j) return false;    if (rhs<php.rhs) return true;    if (rhs>php.rhs) return false;    if (slackSign<php.slackSign) return true;    assert(slackSign != php.slackSign);    return false;    /* : i==php.i ? j<php.j ? true      : j==php.j ? rhs<php.rhs ? true      : rhs==php.rhs ? slack : false; */
  }///////////// STORING COEFS: /////////////////////////////////////////////////////////////
  double Calc__(Column *c, int iCol)  { return Calc__(c); }
  double GetCoef(int iCol) { assert(0); }
  double CalcSlackValue    (i_vec &iNZ, d_vec &xNZ, map<LPCut*,double> &slVal)
    { assert(0); /* no pool restoration */}
};     // class BrVRP _______________________________/// Capacity cut

class Capacity : public LPCut {public: // private:
    /// 1-based product indexation:
    Vector<int> S; // S[i]=bool(i \in S), $ i \in \{1..m\}$
    // can change if we split item types
    double rhs; // the bound;
    int slackSign; // as parameter for branching
    public:  int Type() const { return 345; } // negative: for hyperplanes
    bool IntegerSlack() { return true; } // alfa == 1
    virtual bool CanBeDeleted() const { return false; }
    double GetRHS() const { return rhs; }
    double CalcRHS__ (d_vec&) {return rhs;}
    Capacity(const Vector<int> S_, double RHS=0, int SS=-1) :
      S(S_), rhs(RHS), slackSign(SS) {  }
      ///////////////////// CALCULATION ////////////////////////
      int GetSlackCoef() const { return slackSign; }
      double Calc__(Column *);
      /// Nothing for col.gen. by B&B:
      void CalcIntermSums(Column *) { }  // Utilities for col gen:
      /// Obj func calculation:
      double CalcUsingIntermSums() { } // Calc Pis
      void Print(ostream& os=cout) {    os << " capacity[...] >= " << rhs;  }
      bool operator<(const LPCut& c) const {    if (Type() != c.Type())      return Type() < c.Type();    const Capacity& php = dynamic_cast<const Capacity&>(c);//    return i<php.i ? true : i==php.i ? j<php.j : false;// this is not enough!! Because may be diff rhs
      if (rhs<php.rhs) return true;    if (rhs>php.rhs) return false;    if (slackSign<php.slackSign) return true;    if (slackSign>php.slackSign) return false;    if (S<php.S) return true;    assert(S != php.S);    return false;  }
      ///////////// STORING COEFS: /////////////////////////////////////////////////////////////
      double Calc__(Column *c, int iCol)  { return Calc__(c); }
      double GetCoef(int iCol) { assert(0); } // not to be deleted
      double CalcSlackValue    (i_vec &iNZ, d_vec &xNZ, map<LPCut*,double> &slVal)  { assert(0); /* no pool restoration */
    }};// class Capacity _______________________________
    }                                                // SS_END_NAMESPACE__
#endif // __brVRP_H
