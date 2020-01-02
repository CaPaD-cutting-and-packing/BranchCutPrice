#ifndef __BGN__MYTOOLS_H
#define __BGN__MYTOOLS_H


// #include "myautoptr.h" use simple delete

//typedef double Real; // for SoPlex

extern string __glb_file;
extern int __glb_inumber;


namespace bgn {


//using namespace std;

//#ifdef DBG_ON // !!!!!!!!!!!
template<class D>
class Vector : public std::vector<D> {
  public:
    Vector() : std::vector<D>() { }
    explicit Vector(int m) : std::vector<D>(m) { }
    Vector(int m, const D& d) : std::vector<D>(m,d) { }
#ifdef DBG_ON
    D & operator [](int i)
    {
      assert(i>=0&&i<this->size());
      return std::vector<D>::operator[](i);
    }
    const D & operator [](int i) const
    {
      assert(i>=0&&i<this->size());
      return std::vector<D>::operator[](i);
    }
#endif
}; // also with resize() ??
//#endif

//#define FOR(v,s,e) for(v=s;v<=e;++v)
//#define FORDN(v,s,e) for(v=s;v>=e;--v)
#define FOR_EACH_IN(cont,it,it_type) \
 for(it_type it=(cont).begin();it!=(cont).end();++(it))
#define for_each_in FOR_EACH_IN
#define for_each_r_in(cnt,rit,rit_type) \
 for(rit_type rit=(cnt).rbegin();rit!=(cnt).rend();++(rit))

// rammelvoll

//typedef Vector<double> d_vector;
typedef Vector<double> d_vec;
//typedef Vector<int> i_vector;
typedef Vector<int> i_vec;


/*
template <class M> inline
void MemSet(M*p,int v,int n) { memset(p,v,n*sizeof(M)); }
template <class M> inline
void MemCpy(M*d,M*s,int n) { memcpy(d,s,n*sizeof(M)); }
*/
/* // THIS IS DOOF!!:
template <class D,class E>// inline
D IMax(D d,E e) { return (d>(D)e) ? d : (D)e; }
template <class D,class E>// inline
D Min(D d,E e) { return (d<(D)e) ? d : (D)e; }
template <class D,class E>// inline
void Swap(D &d,E &e) { D t=d; d=(D)e; e=(E)t; }
*/
template <class D,class E>// inline
D Max(D d,E e) { return (d>(D)e) ? d : (D)e; }
template <class D,class E>// inline
D Min(D d,E e) { return (d<(D)e) ? d : (D)e; }
inline int IMax(const int a,const int b) { return a>b ? a: b; }
inline long LMax(const long a,const long b) { return a>b ? a: b; }
inline double FMax(const double a,const double b) { return a>b ? a: b; }
inline int IMin(const int a,const int b) { return a<b ? a: b; }
inline long LMin(const long a,const long b) { return a<b ? a: b; }
inline double FMin(const double a,const double b) { return a<b ? a: b; }
template <class D>// inline
void Swap(D &d,D &e) { D t=d; d=e; e=t; }
template <class D,class E>// inline // care!!
void Swap2Types(D &d,E &e) { D t=d; d=(D)e; e=(E)t; }
#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif
template <class D,class E> //inline
D Del0(D a,E b) { return b==0 ? 0 : D(a/(D)b); }

template <class D>
void Reconstruct(D& d) {
  d.~D();
  new (&d) D();
}

//______________________________________________________
// MY STL
////////////////////////////////////////////////////////


template <class vec,class vec2>
void erase_whichnot(vec &v,vec2 &whichNot) {
  int i,j;
  assert(v.size() >= whichNot.size());
  for (i=j=0;i<whichNot.size();++i)
    if (whichNot[i]) v[j++] = v[i];
  v.erase(v.begin()+j,v.end());
}


template <class C>
void resize(C & c,int m) {
  c.resize(m);
} //____________________________________________________


template <class D>
class cycle : public list<D> {
  private: int n;
  public: void push_back(const D&d) {
            if (this->size()>=n) this->pop_front();
            list<D>::push_back(d);
          }
  public: void push_front(const D&d) {


            if (this->size()>=n) this->pop_back();
            list<D>::push_front(d);
          }
          explicit cycle(int n_=1) :n(n_) { }
};
template <class V1,class V2>
inline void CopyVec(V1&d,V2&s) {
  if (d.size() < s.size()) // only then: m.b. reserving
    d.resize(s.size());
  int i; for(i=0;i<s.size();++i) d[i] = s[i];
}


template <class D>
class CmpByArray {
  D* d; public: CmpByArray(D* const d_) :d(d_) { }
  bool operator()(int i1,int i2) const
  { return d[i1] < d[i2]; }
};


template <class vec>
class CmpByVec {
  vec v; public: CmpByVec(vec & v_) :v(v_) { }
  bool operator()(int i1,int i2) const
  { return v[i1] < v[i2]; }
};


template <class D>
class ReverseCmpByArray {
  D* d; public: ReverseCmpByArray(D* const d_) :d(d_) {}
  bool operator()(int i1,int i2) const
  { return d[i1] > d[i2]; }
};


template <class Ptr>
class CmpPtrByVal {
  public:
  bool operator()(const Ptr p1, const Ptr p2) const
  { assert( p1 && p2 ); return *p1 < *p2; }
};


#if (defined(__GNUC__))
  #define bool int
  #define true 1
  #define false 0
#endif
#if (defined(_MSC_VER)/* || defined(__GNUC__)*/)
  #define and &&
  #define not !
  #define or ||
#endif


// IN/OUTPUT:
ofstream & GetMyLog__();

void error2_handler__();
extern int nErr2__;

#define assert2(expr) if (!(expr)) { PRINT_ERROR("Error2: "<<#expr/*<<" __FILE__ " :" __LINE__   __DATE__*/);  error2_handler__(); }

// errors in GNU when having the macros!!!
// to use in before else only in brackets

#define __asErtm(x,emsg) if (!(x)) { PRINT_ERROR("Error2: "<<emsg);  error2_handler__(); }

#define assert2m __asErtm

#define PRINT__(e) { std::cout << e; }
#define PRINTLN(e) PRINT__(e << '\n')
#define PRINT(e) PRINTLN(e)


#define log_out__ PRINT_LOG__
#define log_ln(e) PRINTLN(e)
#define log_out PRINT__


#define log__ log_out__

#define log_n(n,expr) { if (OUTP_LEV__ >= n) log_ln(expr); }
#define log_n_(n,expr) { if (OUTP_LEV__ >= n) log__(expr); }


template <class vec>
inline void PrintVec(ostream &os,const vec &v) {
  int i;
  for (i=0;i<v.size();++i)
    os << v[i] << ' ';
}


#define OUTP_LEV__ \
  FMin(outputLevel,opt::GlobalOutputLevel())
#define DEF_OUTP_LEVEL 5
extern double outputLevel;


template <class D>
inline void PrintVar(ostream&os,int w,const D& d)
{ os << ' ' << setw(w) << d; }


extern bool fShouldExit;
class UserBreak__ { };
inline bool CheckForUserBreak__() {
  if (fShouldExit)
    throw UserBreak__();
  // Maybe check if keypressed .. and some menu
  return false;
} //____________________________________________________


template <class A>
class Incrementor {
  A & a;
public:
  Incrementor(A& aa) :a(aa) { ++a; }
  ~Incrementor() { --a; }
};//____________________________________________________


template <class T>
inline auto_ptr<T> make_auto_ptr(T* pt)
{ return auto_ptr<T>(pt); }


//______________________________________________________
// MY CONTAINERS
////////////////////////////////////////////////////////
  class StringComparator {
  public:
    bool operator()
      (const char * s1, const char * s2) const
    { return (0 > strcmp(s1,s2)); }
  };//__________________________________________________


  template <class IE, class E>
    inline  IE  Ptr2Iter (const E * pe) {
      assertm(sizeof(IE) == sizeof (E*),
        "Cannot convert pointer to iterator !!! Recode");
    /// Assuming iterator contains a shifted ptr:
    static IE it; // maybe 'volatile' would help to correctly compute the difference...
    *(char**)(&it) = (char*)pe; // but 'volatile' needs -fpermissive to compile!
    IE it2; // so now, at least 'it' points to some valid address.
    *(char**)(&it2) = (char*)pe + (*(char**)(&it) - (char*)(&*it));
    assert(&*it2 == pe); // actually this should be enough to check everything...
    return it2;
  }

//////////////// My application ////////////////////////
class MyApp {
public:
  MyApp();
  ~MyApp();
  void Tests();
  void InitHandlers();
};//____________________________________________________


struct MyTime {
  double t;
  MyTime &operator = (const double _) { t= _; return *this; }
  operator double() const { return t; }
  MyTime(const double _=0) : t(_) { }
};
inline ostream & operator <<
  (ostream &ofs,const MyTime & m) {
/* //ios::fmtflags f=os.flags();
  char buf[20];
  ostrstream os(buf,sizeof(buf));
//os << right;
  if (m.t>3600) os << long(m.t)/3600 << 'h';
  if (m.t>60) os << int(m.t)%3600/60 << 'm';
  os << int(m.t)%60;
  if (m.t < 60) os << '.' << int((m.t-floor(m.t))*100);
  os << ends;
  //os.flags(f);


  ofs << buf;*/
  double tm = m;
  if (tm>=100) tm=ceil(tm);
  else if (tm>=10) tm= ceil(tm*10.0)/10.0;
  else tm = ceil(tm*100.0)/100.0;
  ofs << tm;
  return ofs;
}

inline const char * GetDateAndTime() {
    time_t time0=time(0);
    return ctime(&time0);
}

double Alarm (double seconds);


} // END_COMMON_NAMESPACE__
#endif // __BGN__MYTOOLS_H
