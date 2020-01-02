#ifndef __BGN__MYDEBUG_H
#define __BGN__MYDEBUG_H

void error_handler__();

namespace bgn {

#ifdef assert
#undef assert
#endif
#ifdef assertm
#undef assertm
#endif

#define ASSERT_STR(s) \
 "Error. " s  " :" __FILE__ ":" /*__LINE__ ":  " #__DATE__*/
#define assert(x) \
{ if (!(x)) \
  its_error("Error. " #x " :" __FILE__ " :", __LINE__, \
  __DATE__); }
#define assertm(x,emsg) \
{ if (!(x)) its_error(emsg,"",""); }

#define its_error(msg,m2,m3) \
 { PRINT_ERROR(msg<<m2<<':'<<m3<<flush); error_handler__(); }

#if (defined(DBG_OUTPUT) || defined(DBG_ON))
#define dbg_out(str) \
{ mylog << str << flush; }  // Redirection ?
//    setw(2/*int(FunctionEnterer::rdep/sp_fct____)*/) << " "
#define dbg_outn(level,str) \
{ if (OUTP_LEV__>=level) dbg_out(str<<'\n'); } // no else after ?
#define dbg_outn_(level,str) \
{ if (OUTP_LEV__>=level) dbg_out(str); }
#define ENTER_FN(fn,dep) FunctionEnterer fefefefe(fn,dep)
#else
#define dbg_out(str)
#define dbg_outn(level,str)
#define dbg_outn_(level,str)
#define ENTER_FN(fn,dep)
#endif

#ifdef DBG_ON
#define dbg_assert assert
#define dbgpr(n,expr) { if (OUTP_LEV__>=n) log__(expr); }
#define dbgcout(n,e) { if (OUTP_LEV__>=n) cout << e << flush; }
#else
  #define dbg_assert(e)
  #define dbgpr(n,expr)
  #define dbgcout(n,e)
#endif


#define LOGSTREAM__ cout //bgn::GetMyLog__()


#define PRINT_LOG__(e) { LOGSTREAM__ << e << flush; }
#define PRINT_LOGLN(e) PRINT_LOG__(e << '\n')
#define PRINT_LOG(e) PRINT_LOGLN(e)


extern ofstream thelog__;
extern ofstream * pmylog__; // = &thelog__;
#define mylog cout //bgn::GetMyLog__()
extern int floatWidth;



#ifdef DBG_ON
  #define PRINT_ERROR(str) \
{ dbg_out( str<<endl ); PRINT_LOGLN(str); }
#else
  #define PRINT_ERROR(str) { PRINT_LOGLN(str); \
    ofstream ofs("errata.txt",ios::out|ios::app); \
    ofs << __glb_file << " #" << __glb_inumber << ": " << str << \
    '\n'; }
#endif


extern ofstream dbg_os__;
extern int dbg_level__;
extern double fne_level__;
extern double sp_fct____;

/*
struct FunctionEnterer {
  char *nm;
  double dep;
  static int rdep; //Function stack depth. Single process!
  FunctionEnterer(char*fn, double dp) :nm(fn), dep(dp) {
    if (fne_level__ >= dp)
      dbg_out("--> Fn="<<fn<<" depth="<<dp);
    ++rdep;
  }
  ~FunctionEnterer() {
    --rdep;
    if (fne_level__ >= dep)
      dbg_out("<-- Fn="<<nm<<" depth="<<dep);
  }
}; // --------------------- 07.08.01 -------------------
*/

}
#endif // __BGN__MYDEBUG_H

