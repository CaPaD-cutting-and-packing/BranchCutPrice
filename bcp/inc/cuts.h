/****************************************************
                          cuts.h  -  description
                             -------------------
    begin                : Mon Dec 9 2002
    copyright            : (C) 2002 by Gleb Belov
    email                : belov@math.tu-dresden.de
 ****************************************************/
#ifndef __CUTS_H__32
#define __CUTS_H__32

namespace bgn {                     // SS_BEGIN_NAMESPACE__

class Column;

class LPCut;
typedef Vector<LPCut*> CutList;

class LPCut {
public:
  int signature;
  mutable bool fNodeVisited; // for tree processing
  mutable bool fPresent;
  LPCut() :fNodeVisited(false), signature(1234567), fPresent(false) { }
  virtual ~LPCut() { signature = 654321; } // ?????
  void CheckID() const {
    assert(1234567 == signature);
  }

  mutable int no; // the ordered number
  mutable int __Mark;
  void Mark(int M) {CheckID();__Mark=M;}
  int Mark() {CheckID();return __Mark;}
  double cntDel, // -- how many iterations till deletion
  cntDel0; // -- initial for cntDel

  virtual int Type() const =0;
    // >0: cuts, <0: branching hyperplanes, z.B.
    // 1: SA Cut
    // -1: VRP branch
  virtual int is_const() const { return 0; }
  virtual bool IntegerSlack() const { return false; }
  virtual bool CanBeDeleted() const { return true; }

  virtual void ClearNonRec() { CheckID(); fNodeVisited = false; }
// Calc coef of cut slacks:
  virtual double CalcCutSlackCoef__(LPCut *pcut)
  { return (pcut == this) ? GetSlackCoef() : 0; }
  virtual double GetCutSlackCoef__(LPCut *pcut)
  { return (pcut == this) ? GetSlackCoef() : 0; }

// RHS calculation
  virtual int GetSlackCoef() const { return 1; } // <=
// Only rhs or lhs is valid for a cut and that is returned,
// must have been precalculated:
  virtual double GetRHS() const =0;
  virtual double CalcRHS__(d_vec&) =0; // each time ?

// Utilities for col gen:
// Coef calculation:
// (use ClearNonRec() on all involved before)
  virtual double Calc__(Column *) =0;
  virtual void CalcIntermSums(Column *)=0;
  virtual double CalcUsingIntermSums()=0;
  /// Bound calculation:
  virtual double CalcApprCoef(int k) { return 0; }
  virtual double CalcApprErrorL() { return 0; }
  virtual double CalcApprErrorU() { return 0; }

  virtual void AddToSums(int k, int x) { }
  virtual void Sub1FromSums(int k) { }
  virtual void CalcConstTerms(Column *) { }
  virtual void AssignConstTerms() { }

  virtual void Print(ostream& =cout)=0;
  virtual void Number(int &nn)
  { if (-1==no) no = nn++; }
  virtual void ProduceListOfInvolved(Vector<LPCut*> &cuts)
  const
    {
      if (fNodeVisited) return;
      fNodeVisited = true;
      cuts.push_back((LPCut*)this);
    }
  virtual void Sort() { }
  virtual bool operator<(const LPCut&) const =0;

///////////// STORING COEFS: ////////////////
/////////////////////////////////////////////
  virtual int GetNCoefs() { return 0; }
    // iCol <= GetNCoefs() !!!
  virtual double Calc__(Column *c, int iCol) { return 0; }
  // These 2 are needed only in SACut:
  virtual double GetCoef(int iCol) { return 0; }
  virtual double CalcRawCoef(Column *c,int iCol)
    { return 0; }

  // But this is needed for all:
  virtual double CalcSlackValue
  (i_vec &iNZ, d_vec &xNZ, map<LPCut*,double> &slVal)=0;
// deprecated:
/*  virtual
    double Calculate(Column*c)=0;
  double Calc(Column*c) { return Calculate(c); }
    virtual double CalcLRHS(Column *) =0;
  virtual
    void Clear() { fNodeVisited = false; }*/
};

struct CGIV { // violated cut descriptor
    const LPCut *c; // iterator in cutpool
    float wgt; // weight for sorting
    double v; // violation grade
    CGIV(const LPCut * pc,double vv)
      :c(pc), v(vv) { pc->CheckID(); }
    bool operator<(const CGIV& cgiv) const
      { return wgt > cgiv.wgt; } // reverse
};


}                                                // SS_END_NAMESPACE__

#endif // __CUTS_H__32
