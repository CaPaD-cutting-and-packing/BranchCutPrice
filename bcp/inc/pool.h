/*******************************************************
                          pool.h  -  description
                             -------------------
    begin                : Sat Dec 7 2002
    copyright            : (C) 2002 by Gleb Belov
    email                : belov@math.tu-dresden.de
 ********************************************************/

#ifndef __POOL_H__32
#define __POOL_H__32

namespace bgn {                     // SS_BEGIN_NAMESPACE__

// insert SORTED cols/cuts
template <class D> //, class Cmp = less<D> >
class Pool : public set<D> {
public:
  D* Add(const D &d) { // returns the addr of the added copy
//    d->Sort(); // for lexico. ordering
    pair<typename set<D>::iterator,bool> res =
      this->insert(d);
    if (!res.second) return NULL; // if already
    return (D*)&*res.first;
  }
  D* Find(const D &d) {
//    d->Sort(); // do it yourself
    typename set<D>::iterator it=set<D>::find(d);
    if (it != this->end()) return (D*)&*it;
    return NULL; // if not found
  }
  D* FindVirtual(const D &d) {
    D* pd = Find(d);
    if (pd) if (not pd->Hidden()) return pd;
    return NULL; // if not found
  }
};

}                                                // SS_END_NAMESPACE__

#endif // __POOL_H__32
