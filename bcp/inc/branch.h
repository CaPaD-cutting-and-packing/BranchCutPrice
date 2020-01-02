/*****************************************************
                          branch.h  -  description
                             -------------------
    begin                : Mon Dec 9 2002
    copyright            : (C) 2002 by Gleb Belov
    email                : belov@math.tu-dresden.de
 *****************************************************/

#ifndef __BRANCH_H__32
#define __BRANCH_H__32

//#include "lpcol.h"
#include "cuts.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

class VarBnd //:public LPCut
{
public:
  /*virtzual*/ int Type() const { return 0; }
  int j; // main var index // we keep a main index
//    Column * c; // the column
  bool upper; // or lower
  int bnd;
  explicit VarBnd(int jCol=0,bool u=false,int b=0)
    :j(jCol), upper(u), bnd(b) { }
  bool operator<(const VarBnd &bd) const
  {
    if (Type()==bd.Type()) {
      int a[2],b[2];
      ProduceVec(a);
      bd.ProduceVec(b);
      return lexicographical_compare(a,a+2,b,b+2);
    }
    return Type()<bd.Type();
  }
  template <class A>
  void ProduceVec(A *v) const
  { v[0]=j;v[1]=(A)upper;/*v[2]=bnd;*/ }
};

}                                                // SS_END_NAMESPACE__

#endif // __BRANCH_H__32

