#ifndef __MYOPT_H
#define __MYOPT_H

namespace bgn {

//______________________________________________________
// MY OPTIONS
////////////////////////////////////////////////////////
namespace opt {
  class VirtualOption {
  public:
    const char *name;
    const char *descr;
    VirtualOption(const char*n,const char *d)
      : name(n), descr(d) { }
    virtual void read(istream &)=0;
    virtual void write(ostream &)=0;
    virtual void SetDefault()=0;
  };//____________________________________________________
  typedef list<VirtualOption*> OptContainer;
  typedef
    map<const char*,VirtualOption*,
      StringComparator>
    OptAssocCont;
}
// In the above namespace:
inline opt::OptContainer & operator <<
(opt::OptContainer& oc,opt::VirtualOption * vo) {
  oc.push_back(vo);
  return oc;
}
namespace opt {
  template <class type>
  class Option :public VirtualOption {
  public:
    type * pVar;
    const type defVal;
    void read(istream& is) { is >> (*pVar); }
    void write(ostream& os) { os << *pVar; }
    void SetDefault() { *pVar = defVal; }
    Option
    (type *pv,const type def,const char*n,const char *d)
    :pVar(pv), defVal(def), VirtualOption(n,d) { }
  };//____________________________________________________
  template <class type,class t2>
  VirtualOption * MakeOpt
  (type *pv,const t2 def,const char*n,const char *d) {
    return (VirtualOption*)
      (new Option<type>(pv,(type)def,n,d));
  } //____________________________________________________
  class Config;
  class OptSection {
  public:
    const char * name;
    const char * descr;
    OptContainer options;
    OptAssocCont optionsAssoc;
    Config * pCfg;
    int priority;
    OptSection
      (const char *n,const char *d,const OptContainer & c,
       Config * pc, const int pr);
  };//____________________________________________________
  typedef map<const char*,OptSection*,StringComparator>
    SectionContainer;
  inline SectionContainer& operator <<
    (SectionContainer&sc,OptSection* os) {
    sc.insert(make_pair(os->name,os));
    return sc;
  }
  class Config {
  public:
    const char * fileName;
    SectionContainer sections;
    Config(const char * n) :sections()
    { fileName = n; }
    void ReadOptions(), WriteOptions(),
      WriteOptionsShort(ostream&);
    void SetDefault();
    void ReadOptionName(istream&,char*);
    void ReadSectionName(istream&,char*);
    void SkipComments(istream&);
  };//____________________________________________________

  double & GlobalOutputLevel();
  Config * SolverCfg();
//  Config * OutputCfg();
  extern bool writeParams;

  void ReadOptions(), WriteOptions();

} //_namespace opt________________________________________

namespace myopt = opt;

} // END_COMMON_NAMESPACE__

#endif // __MYOPT_H

