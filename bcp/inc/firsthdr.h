#ifdef  __FIRSTHDR_H
#error should be included only once
#endif
#define __FIRSTHDR_H


/* Author: Gleb Belov <bg37@gmx.net> */


#define COMMON_NAMESPACE__ bgn
//#define BEGIN_COMMON_NAMESPACE__ namespace bgn {
//#define END_COMMON_NAMESPACE__   }


namespace bgn {
//using namespace std;
}


//#define DBG_ON


#include "mydebug.h"
#include "mytools.h"
#include "myopt.h"
#include "mystat.h"


//#include "spxdefines.h"
#include "timer.h"
#include "random.h"
namespace bgn {
typedef soplex::Timer Timer;
typedef soplex::Random Random;
}


#include "sshdrs.h"
#include "mymath.h"


