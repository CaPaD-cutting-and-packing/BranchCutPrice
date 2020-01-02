#ifndef __LPDEFS_H__32
#define __LPDEFS_H__32

#error Do not include me

namespace bgn {                     // SS_BEGIN_NAMESPACE__

typedef double Real;

typedef double d_float; // dual multipliers
typedef double v_float; // OF values;
typedef double   a_float; // restr matr coefs

typedef Vector<d_float> d_vec;
typedef Vector<a_float> a_Vector;
typedef Vector<signed char> constrtypevec; // !!!

  class LPCut;

typedef int colentry; // only integers in orig. constr
// level cuts are not those

struct ID {
  int i; 
  colentry d;
  void set(const int i_=0,const colentry d_=0)
  { i=i_; d=d_; }
  ID(const int i_=0,const colentry d_=0) :i(i_), d(d_) {}

  bool operator < (const ID& id) const
  { return (i<id.i) ? true : (i>id.i) ? false : d>id.d;}
};

struct Column {
  mutable Vector<ID> id;
  typedef Vector<ID>::iterator iterator;
  typedef Vector<ID>::const_iterator const_iterator;
  typedef iterator id_iterator;

// CUT SLACK:
  LPCut *c; // ==0 if not
  colentry d; // != 0 => column is SLACK!

  d_float ofc;  // OF coeff
  void SetOFC(const d_float c) { ofc=c; }
  void SetObj(const d_float c) { ofc=c; }
  d_float GetObj() const { return ofc; }

  void * addi_info;
  void * & GetAddiInfo() { return addi_info; }

  explicit Column(const int m=0) :ofc(0), c(0), d(0),
    addi_info(0)
  { id.reserve(m); } // Now [] goes but values lost then
  Column (d_vec &v):ofc(0), c(0), d(0),
      addi_info(0)
  {
    id.reserve(v.size());
    for (int i=0;i<v.size();++i)
      if (v[i]) PushID(i,(colentry)v[i]);
  }
  void clear() { id.clear(); }
  void PushID(const int i,const colentry d)
  { id.push_back(ID(i,d)); }
//  d_float GetRedCost(const d_vec &d) const;

  void MakeSlack(const int i, const colentry v)
  { clear(); PushID(i,v); d=v; }
  void MakeSlack(LPCut * const p, const colentry v)
    { c = p; d = v; }
  LPCut * GetCutSlackCut() const { return c; }
  int GetCutSlackCoef() { return d; }
  bool IsSlack() const { return (d!=0); }
  colentry GetSlackCoef() const { return d; }

template <class vec>

  d_float VecProd(vec & v) {
    double sum=0;
    for_each_in(id,iid,iterator) sum += v[iid->i] * iid->d;
    return sum;
  }

  bool operator < (const Column & c) const
  { return id.size() ? (c.id.size() ?
    lexicographical_compare(id.begin(),id.end(),
    c.id.begin(),c.id.end()) : true) :
     (c.id.size() ? false
      : GetCutSlackCut() < c.GetCutSlackCut()); }
  
  ostream& print(ostream& os) {
    for_each_in(id,it,iterator)
      os<<(it->i)<<':'<<(it->d)<<' ';
    os<<GetObj();
    return os;
  }

  void Sort() const { // before comparison
    sort(id.begin(),id.end());
  }
};
/* // MSVC has errors with this.
inline ostream& operator<< (ostream& os,Column& c) {
  for_each_in(c.id,it,Column::iterator)
    os<<(it->i)<<':'<<(it->d)<<' ';
  os<<c.GetObj();
  return os;
} */

class LPCut;

struct ColumnSet {
  set<Column> cs;
  typedef set<Column>::iterator iterator;
  typedef set<Column>::const_iterator const_iterator;

  void Add(const Column& c) { cs.insert(c); }
};
typedef ColumnSet ColSet;

struct ColumnList {
  list<Column> cl;
  typedef list<Column>::iterator iterator;
  typedef list<Column>::const_iterator const_iterator;

  void Add(const Column& c) { cl.push_back(c); }
  Column& Add()
  { cl.push_back(Column()); return cl.back(); }
};
typedef ColumnList ColList;

//  enum LPStatusValues
//  { none, ok, infeas, unbounded, LPCYCLE, error };
const char* const LPStatusName[] =
  { "none", "ok", "infeas", "unbnd", "cycle", "error" };

}                                                // SS_END_NAMESPACE__

#endif // __LPDEFS_H__32
