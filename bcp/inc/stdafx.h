// stdafx.h : include file for standard system include files,
//  or project specific include files that are used frequently, but
//      are changed infrequently
//

#ifndef _BITS_ATOMICITY_H
#define _BITS_ATOMICITY_H       1

typedef int _Atomic_word;

/*

  _Atomic_word
  __attribute__ ((__unused__))
  __exchange_and_add(volatile _Atomic_word* __mem, int __val)
  {
    register _Atomic_word __result;
    __asm__ __volatile__ ("lock; xadd{l} {%0,%1|%1,%0}"
			  : "=r" (__result), "=m" (*__mem)
			  : "0" (__val), "m" (*__mem));
    return __result;
  }

  void
  __attribute__ ((__unused__))
  __atomic_add(volatile _Atomic_word* __mem, int __val)
  {
    __asm__ __volatile__ ("lock; add{l} {%1,%0|%0,%1}"
			  : "=m" (*__mem) : "ir" (__val), "m" (*__mem));
  }

  */

static inline _Atomic_word
__attribute__ ((__unused__))
__exchange_and_add (volatile _Atomic_word *__mem, int __val)
{
  register _Atomic_word __result;
    __asm__ __volatile__ ("lock; xadd{l} {%0,%1|%1,%0}"
			  : "=r" (__result), "=m" (*__mem)
			  : "0" (__val), "m" (*__mem));
/*  __asm__ __volatile__ ("lock; xadd{l} {%0,%1|%1,%0}"
                        : "=r" (__result), "+m" (*__mem)
                        : "0" (__val)
                        : "memory");
*/  return __result;
}

static inline void
__attribute__ ((__unused__))
__atomic_add (volatile _Atomic_word* __mem, int __val)
{
  __asm__ __volatile__ ("lock; add{l} {%1,%0|%0,%1}"
                        : "+m" (*__mem) : "ir" (__val) : "memory");
}

#endif /* atomicity.h */



#include <sstream>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <memory>
#include <cstring>
#include <list>
//#include <valarray>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <numeric>
#include <functional>
//#include <utility>
#include <stdexcept>
#include <csignal>
#include <cmath>
#include <climits>
#include <ctime>
#include <cerrno>
#include <cstdlib>
#include <cstdio>
#include <sys/time.h>

//#include <assert>


using namespace std;


#include "firsthdr.h"

