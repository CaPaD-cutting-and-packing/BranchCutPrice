#ifndef __LPCOL_H__32
#define __LPCOL_H__32

// #include "cuts.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

class LPCut;

//typedef double Real;

typedef int colentry; // only integers in orig. constr
// level cuts are not those

struct ID {
  int i; 
  colentry d;
  void set(const int i_=0,const colentry d_=0)
  { i=i_; d=d_; }
  ID(const int i_=0,const colentry d_=0) :i(i_), d(d_) {}

  bool operator < (const ID& id) const
  { return (i==id.i) ? (d<id.d) : (i<id.i); }
};

class Column {
public:
  explicit Column(const int m=0) :ofc(0), /*c(0),*/ d(0),
    addi_info(0), nLB(0), nUB(0), fInfeas(0), nHidden(0) /*, lb(0), ub(1e100)*/
  { id.reserve(m); } // Now [] goes but values lost then
  Column (d_vec &v):ofc(0), /*c(0),*/ d(0),
      addi_info(0), nLB(0), nUB(0), fInfeas(0), nHidden(0)
  {
    id.reserve(v.size());
    for (unsigned i=0;i<v.size();++i)
      if (v[i]) PushID(i,(colentry)v[i]);
  }
  mutable Vector<ID> id;
  bool operator < (const Column & c) const {
    return (d==c.d) ? (id < c.id) : (d<c.d);
  }
  typedef Vector<ID>::iterator iterator;
  typedef Vector<ID>::const_iterator const_iterator;
  void PushID(const int i,const colentry d)
  { id.push_back(ID(i,d)); }

  double ofc;  // OF coeff
  void SetObj(const double c) { ofc=c; }
  double GetObj() const { return ofc; }

  // PSEUDO_COSTS:
  int nLB, nUB;
  double sumPsCLB, sumPsCUB;
//  double lb,ub;
/*  void SetLB(const double b) { lb = b; }
  double GetLB() const { return lb; }
  void SetUB(const double b) { ub = b; }
  double GetUB() const { return ub; }*/

// BEING CUT SLACK:
//  LPCut *c; // ==0 if not
  colentry d; // != 0 => column is SLACK!
  void MakeSlack(const int i, const colentry v)
  { clear(); PushID(i,v); d=v; ofc = 0; }
/*  void MakeSlack(LPCut * const p, const colentry v)
    { c = p; d = v; }*/
//  LPCut * GetCutSlackCut() const { return c; }
//  colentry GetCutSlackCoef() { return d; }
  bool IsSlack() const { return (d!=0); }
  colentry GetSlackCoef() const { return d; }

  bool fInfeas;
  bool fInfeasible() { return fInfeas; }

  mutable int nHidden;
  int Hidden() { return nHidden; }

  int index;
  int GetIndex() { return index; }
  
  void * addi_info;
  void * & GetAddiInfo() { return addi_info; }

  void clear()
  { id.clear(); d=0; ofc=0; }

template <class vec>
  double VecProd(vec & v) {
    double sum=0;
    for_each_in(id,iid,iterator) sum += v[iid->i] * iid->d;
    return sum;
  }

  void Sort() const { // before comparison
    sort(id.begin(),id.end());
    CheckRepetitions();
  } // Now:
  void SortWithMerging1st() const { // if one piece was pre-set
    sort(id.begin(),id.end());
    if (id.size()>1)
      if (id[0].i == id[1].i) {
        id[1].d += id[0].d;
        id.erase(id.begin());
      }
    CheckRepetitions();
  }
  void CheckRepetitions() const {
    unsigned i;
    for (i=1;i<id.size();++i)
      assert(id[i].i != id[i-1].i);
  }
  // Now:
/*  bool operator < (const Column & c) const
  { return id.size() ? (c.id.size() ?
    lexicographical_compare(id.begin(),id.end(),
    c.id.begin(),c.id.end()) : true) :
     (c.id.size() ? false
      : GetCutSlackCut() < c.GetCutSlackCut()); }*/
  
  ostream& print(ostream& os) {
    for_each_in(id,it,iterator)
      os<<(it->i)<<':'<<(it->d)<<" o";
    os<<GetObj();//<<" l"<<GetLB()<<" u"<<GetUB();
    return os;
  }
};
/* // MSVC has errors with this.
inline ostream& operator<< (ostream& os,Column& c) {
  for_each_in(c.id,it,Column::iterator)
    os<<(it->i)<<':'<<(it->d)<<' ';
  os<<c.GetObj();
  return os;
} */

  struct ColId {
    int j; // col index, == -1 if cut Slack, then slackCut != NULL
    LPCut *slackCut; // the cut whose slack this is
    ColId(const int j_=0,LPCut * const p_=NULL)
      :j(j_), slackCut(p_) { }
  };


struct ColumnSet { // orders sorted cols
  set<Column> cs;
  typedef set<Column>::iterator iterator;
  typedef set<Column>::const_iterator const_iterator;

  void Add(const Column& c) { c.Sort(); cs.insert(c); }
  void clear() { cs.clear(); }
};
typedef ColumnSet ColSet;

struct ColumnList {
  list<Column> cl;
  typedef list<Column>::iterator iterator;
  typedef list<Column>::const_iterator const_iterator;

  void Add(const Column& c) { cl.push_back(c); }
  Column& Add()
  { cl.push_back(Column()); return cl.back(); }
  void clear() { cl.clear(); }
};
typedef ColumnList ColList;

}                                                // SS_END_NAMESPACE__

#endif // __LPCOLS_H__32
