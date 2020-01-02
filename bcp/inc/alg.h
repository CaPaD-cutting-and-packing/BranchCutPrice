#ifndef __ALG_H__32xx
#define __ALG_H__32xx

namespace bgn { // BEGIN_COMMON_NAMESPACE__

class Env {
  ofstream ofs2,
    ofs3;
  double outputLevel;
public:
  char outname1[1024];
  ofstream & GetOFS2() { return ofs2; }
  ofstream & GetOFS3() { return ofs3; }
  ofstream & GetLog2() { return ofs2; }
  ofstream & GetLog3() { return ofs3; }
  double GetOutputLevel() { return outputLevel; }
};

class Alg {
  Env * env;
public:
  Env & GetEnv() const { return *env; }
  Alg(Env & e) : env(&e) { }
  Alg(Env * pe) : env(pe) { }
  Alg(Alg & a) : env(&(a.GetEnv())) { }
  Alg(Alg * pa) : env(&(pa->GetEnv())) { }
//  Alg & operator=(const Alg& a) { 

  double GetOutputLevel() { return env->GetOutputLevel(); }
};

} // END_COMMON_NAMESPACE__

#endif // __ALG_H__32xx
