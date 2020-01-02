#ifndef __SSHDRS_H
#define __SSHDRS_H


//#define namespace bgn {                     // SS_BEGIN_NAMESPACE__ namespace ss {
//#define }                                                // SS_END_NAMESPACE__ }


// BEGIN_COMMON_NAMESPACE__
// END_COMMON_NAMESPACE__


namespace bgn {                     // SS_BEGIN_NAMESPACE__
  using namespace COMMON_NAMESPACE__;
/*#ifdef DBG_ON // !!!!!!!!!!!
template<class D>
class Vector : public std::Vector<D> {
  public:
    Vector() : std::Vector<D>() { }
    explicit Vector(int m) : std::Vector<D>(m) { }
    D & operator [](int i)
    { return std::Vector<D>::operator[](i); }
    const D & operator [](int i) const
    { return std::Vector<D>::operator[](i); }
}; // also with resize() ??
#endif*/

}                                                // SS_END_NAMESPACE__


#define CFG_NAME__ "bcp.cfg"


#include "alg.h"
//#include "sstypes.h"
//#include "ssopt.h"


#endif // __SSHDRS_H
