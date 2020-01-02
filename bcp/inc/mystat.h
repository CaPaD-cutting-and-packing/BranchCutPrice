#ifndef __MYSTAT_H
#define __MYSTAT_H

namespace bgn {

namespace mystat {

class Accumulator {
  int flags;
  double nTimes;
  double accu;
  double last;
  const char * name;
  static int ww;  // Output width
public:
  enum { bPrintRatio=1, bPrintAbs=2, bPrintTotal=4 };
  Accumulator(const char*n,const int fl=7)
    : name(n), flags(fl) { reset(); }
  void Do(int what, ostream &);
  void AddTry() { ++nTimes; last = 0; }
  void AddAmount(double a) { accu += a; last = a; }
  void AddPositive() { AddAmount(1); }
  void reset() { nTimes=accu= 0.0; }
  double GetLast() { return last; }
  bool GetTrigger() { return last!=0; }
};//____________________________________________________


template <class D>
class VectorStat { // each time fixes the next accumulated
  // value and tracks the differences. their ave. and dev.
  int n;
  D div;
  int addingAt;
  Vector<D> prev;
  Vector<D> accu;
  Vector<D> accu_sqr, accu_log, accu_log_sqr;
public:
  VectorStat() : n(0), div(1), addingAt(0) { }
  int getN() const { return n; }
  void FirstTimeSetDivisorForEach(const D d) {
    if (!n) {div=d; assert(d);}
  }
  VectorStat<D>& NextRun() { ++n; addingAt=0; return *this; }
  VectorStat<D>& operator<<(D d); // d is the accumulated value!!!
  void print_ave(ostream&,const char) const;
  void print_ave_log(ostream&,const char) const;
  void print_dev(ostream&,const char) const;
  void print_dev_log(ostream&,const char) const;
  void print_diff(ostream&,const char) const;
};

template <class D>  //inline // !!!
 VectorStat<D>&
  VectorStat<D>::operator<<(D d) {
    d /= div;
    if (addingAt==accu.size()) {
      if (n>1) { assert(0); }
      else {
        // assert(
 //       cout << " Ad1: "<<d<<" at:"<<addingAt;
        ++ addingAt;
        prev.push_back(0); accu.push_back(d); accu_sqr.push_back(d*d);
        accu_log.push_back(log(d));
        accu_log_sqr.push_back(log(FMax(0.01, d*d)));
      }
    }
    else {
//      cout << " Ad2: "<<d<<" at:"<<addingAt;
      assert(addingAt < accu.size());
      prev[addingAt] = accu[addingAt];
      accu[addingAt] = d;
      accu_sqr[addingAt] += pow((d-prev[addingAt]),2);
      accu_log[addingAt] += log(FMax(0, d-prev[addingAt])+1);
      accu_log_sqr[addingAt] += pow(log(FMax(0, d-prev[addingAt]))+1,2);
      ++ addingAt;
    }
    return *this;
  }

//// LOGARITMIC VALUES GIVEN

template <class D>  //inline // !!!
  void VectorStat<D>::print_ave(ostream& os, const char dv) const {
  assert(n);
  for (int i=0;i<accu.size();++i)
    os << accu[i]/n << dv;
}
template <class D>  //inline // !!!
  void VectorStat<D>::print_dev(ostream& os, const char dv) const {
  assert(n>=2);
  for (int i=0;i<accu.size();++i)
    os << sqrt((accu_sqr[i]-pow(accu[i],2)/n) / (n-1)) << dv;
}
template <class D>  //inline // !!!
  void VectorStat<D>::print_diff(ostream& os, const char dv) const {
  for (int i=0;i<accu.size();++i)
    os << (accu[i] - prev[i]) << dv;
}
template <class D>  //inline // !!!
  void VectorStat<D>::print_ave_log(ostream& os, const char dv) const {
  assert(n);
  for (int i=0;i<accu.size();++i)
    os << exp(accu_log[i]/n)-1 << dv;
}
template <class D>  //inline // !!!
  void VectorStat<D>::print_dev_log(ostream& os, const char dv) const {
  assert(n>=2);
  for (int i=0;i<accu.size();++i)
    os << exp(sqrt((accu_log_sqr[i]-pow(accu_log[i],2)/n) / (n-1)))-1 << dv;
}


} // namespace mystat

} // END_COMMON_NAMESPACE__

#endif // MYSTAT_H
