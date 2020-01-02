#ifndef __MYMATH_H__32
#define __MYMATH_H__32

namespace bgn {                     // SS_BEGIN_NAMESPACE__

// VECTOR PRODUCT: 2 arrays by iterators
template <class iter1,class iter2>
double VectorProduct(iter1 a1b,iter1 a1e,iter2 a2b) {
  double sum=0;
  while (a1b < a1e)
    sum += (*(a1b++)) * (*(a2b++));
  return sum;
} //____________________________________________________
#define NOT_ZERO(v) (0!=(v)) // should be tested always
  // care when float=double usw.
const double INFINITY__ = 1e+100;

template <class LPint>
LPint __gcd(LPint a,LPint b) {
  a=fabs(a); b=fabs(b);
  while(a*b!=0) if(a>b) a=fmod(a,b); else b=fmod(b,a);
  return a+b;
}

template <class D,class E>
inline D ceilEps(const D d,const E eps)
{ return ceil(d - (D)eps); }

template <class iter>
inline double Norm2(iter i1,iter i2)
{ return VectorProduct(i1,i2,i1); }

template <class iter>
inline double Norm(iter i1,iter i2)
{ return sqrt(Norm2(i1,i2)); }

template <class D>
inline D frac(const D d) { return d-floor(d); }
template <class D>
inline D Frac(const D d) { return d-floor(d); }

template <class D>
inline D Round(const D d) { return (D)floor((d+0.5)); }
template <class D>
inline D round(const D d) { return (D)floor((d+0.5)); }
template <class D>
inline D smfrac(const D d) { return d-round(d); }

#define round(x) floor((x)+0.5)
template <class N>
inline int sgn(N v) { return v<0 ? -1: v>0 ? 1:0; }
template <class N> inline int sgn_eps(N& v,N eps)
{ return v<-eps?-1:v>eps?1:0; }
template <class N>
inline N abs(const N&v) {return v>=0 ? v: -v; }

}                                                // SS_END_NAMESPACE__

#endif // __MYMATH_H__32

