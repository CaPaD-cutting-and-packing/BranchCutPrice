#ifndef __BB_MOS_H
#define __BB_MOS_H

#include "bb.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

class BBMos: public BB { // for min of open stacks
public:
//  BBMos() : BB() { }
  /////////// INPUT: ////////////////
  Vector<int> fOpen;
  int nOpenMax;
  bool fDominance; // whether check dom.
protected:
  int nOpenIn;
  int nStillOpen;
  void Optimize();
}; // BBMos
  //____________________________________________________

}                                                // SS_END_NAMESPACE__

#endif // __BB_MOS_H
