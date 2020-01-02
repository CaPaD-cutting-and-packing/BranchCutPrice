#ifndef __LPCS_H__32
#define __LPCS_H__32


#include "lp.h"


namespace bgn {                     // SS_BEGIN_NAMESPACE__


template <class PtrD>
class CmpPtrByVal {
public: bool operator() (const PtrD v1, const PtrD v2) const
        {return (*v1) < (*v2);}
};


class ColTree
  //:public ID // -- GNU has problems
{
  unsigned isVar :1;
  unsigned infeas :1;
  unsigned isSlack :1;
  double obj;
  ColTree *parent;
public:
  // struct ID {
  int i;
  d_float d;
  operator ID () { return ID(i,d); }
  bool operator < (const ColTree& id) const
  { return (i<id.i) ? true : (i>id.i) ? false : d>id.d;}
  // };
  typedef set<ColTree*,CmpPtrByVal<ColTree*> > ColBranches;
  ColBranches sons;


//  bool operator==(const ID &id)
//  { return (i==id.i)&&(d==id.d); }


  ColTree(ID *id=0)
    :isVar(0), parent(0), infeas(0), isSlack(0), obj(0),
    i(0), d(0)
  { if (id) {i=id->i;d=id->d;} }
  bool IsVar() { return isVar; }
  void MarkAsVar(bool y=true) {isVar = y;}
  bool GetInfeas() {return infeas; }
  void SetInfeas(bool i) {infeas = i;}
  bool IsSlack() {return isSlack;}
  void SetSlack(bool y=true) {isSlack=y;}
  double GetObj() { return obj; }
  void SetObj(double o) {obj=o;}
  ColTree *GetParent() { return parent; }
  void SetParent(ColTree *p) {parent=p;}


  int Idx() {return i;}
  d_float Val() {return d;}
};


class LPColSet {
  list<ColTree> nodes; // list of vars?
  int nCols;
public:
  int idxMax; // n rows
  double objMax; // max objective coef
  ColTree colRoot;


  LPColSet() :idxMax(0), nCols(0), objMax(0) { }
  int NCols() {return nCols;}
protected:
  ColTree * AddNZTo(ColTree *ct, ID *id) {
    ColTree *ct1=add2(nodes,ColTree(id));
    ct->sons.insert(ct1);
    ct1->SetParent(ct);
    return ct1;
  }
public:
  ColTree * AddCol(Column &c) {
    sort(c.id.begin(),c.id.end()); // + check duplicates somewhere
    // + max dim
    ColTree * ct = &colRoot;
    ColTree * ct1;
    ColTree::ColBranches::iterator ct1it;
    for_each_in(c.id,iid,Column::iterator) {
      ColTree ctt(&*iid);
      ct1it = ct->sons.find
        (&ctt); // ?
      if (ct->sons.end() == ct1it)
        ct1 = AddNZTo(ct, &*iid);
      else ct1 = *ct1it;
      ct = ct1;
    }
// A really new column (may be non-maximal):
    if (not ct->IsVar()) {
      ct->MarkAsVar();
      ct->SetObj(c.GetObj()); // + UB,LB
      ++ nCols; // correct only here?


      if (ct->i > idxMax) idxMax = ct->i; // where ch?
      if (fabs(ct->GetObj()) > objMax)
        objMax = fabs(ct->GetObj());
    }
    ct->SetInfeas(false); // added explicitly
    ct->SetSlack(false);
    return ct;
  }
  ColTree * AddSlack(int i,int coef,bool infeas) {
    Column c(1);
    c.id.push_back(ID(i,coef));
    c.SetObj(0);
    ColTree * ct=FindCol(c);
    if (ct)
      return ct; // if added by user explicitly
    ct = AddCol(c);
    ct->SetSlack(true);
    ct->SetInfeas(infeas); // what if explicitely added?
     // a flag ?
    return ct;
  }
  // returns 0 if not found
  ColTree * FindCol(Column &c) {
    sort(c.id.begin(),c.id.end()); // + check duplicates somewhere
    // + max dim
    ColTree * ct = &colRoot;
    ColTree * ct1;
    ColTree::ColBranches::iterator ct1it;
    for_each_in(c.id,iid,Column::iterator) {
      ColTree ctt(&*iid);
      ct1it = ct->sons.find
        (&ctt); // ? ??!!
      if (ct->sons.end() == ct1it)
        return 0;
      else ct1 = *ct1it;
      ct = ct1;
    }
    if (ct->IsVar()) // !!!
      return ct;
    return 0;
  }
  // Produce column from tree representation,
  // reversely sorted
  static void MakeColumn(Column &c, ColTree *ct) {
    assert(ct->GetParent());
    assert(ct->IsVar());
    c.SetObj(ct->GetObj());
    c.id.clear();
    while (ct->GetParent()) { // w/o colRoot
      c.id.push_back((ID)(*ct));
      ct = ct->GetParent();
    }
  }
};//__LPColSet__________________________________________


}                                                // SS_END_NAMESPACE__


#endif // __LPCS_H__32


