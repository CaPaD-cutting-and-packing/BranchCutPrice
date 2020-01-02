// FILE: mytools.cpp
// AUTHOR: Gleb <Belov@web.de>

#include "stdafx.h"
//#include "mytools.h"
#include "lasthdr.h"

string __glb_file;
int __glb_inumber;

namespace bgn {          //  BEGIN_COMMON_NAMESPACE__

int nErr2__=0;

void error2_handler__() {
  ++ nErr2__;
//  if (signal_error__)
//    throw logic_error("Some error.");
}


double outputLevel=2;
ofstream *pmylog__ = &thelog__;
//ostream & mylog = thelog__;
ofstream & GetMyLog__()
 { return *(COMMON_NAMESPACE__::pmylog__); }

MyApp myApp;

MyApp::MyApp() {
//  opt::ReadOptions(); // done in the algorithm
  InitHandlers();
  Tests();
} //____________________________________________________
MyApp::~MyApp() {
//  opt::WriteOptions();// could be changed when executing
} //____________________________________________________
void MyApp::Tests() {
  double v=0;
  double v1=v;
  assertm(not NOT_ZERO(v1),
    "0!=0 !!! Modify NOT_ZERO macro.");
} //____________________________________________________

bool fShouldExit=false; // Can be modified by the appl.
 double  Alarm (double seconds) 
     { 
       struct itimerval old, new_; 
//       new.it_interval.tv_sec = 0; 
       new_.it_value.tv_usec = (long int)((seconds-floor(seconds))*(1000000.0));
       new_.it_value.tv_sec = (long int) seconds;
       new_.it_interval = new_.it_value;
       if (setitimer (ITIMER_VIRTUAL, &new_, &old) < 0)
         return -1.34; 
       else
         return (double)old.it_value.tv_sec;
     } 
void SigTermHandler(int signum) {
  throw UserBreak__();
}
void SigErrHandler(int signum) {
  char *s = "Unknown error";
  switch (signum) {
  case SIGFPE: s="Floating point exception";
    break;
  case SIGSEGV: s="Segment violation";
    break;
  case SIGILL: s="Illegal instruction";
    break;
  }
  throw runtime_error(s);
} //____________________________________________________
void MyApp::InitHandlers() {
//#ifndef DBG_ON
  signal(SIGTERM,SigTermHandler);
  signal(SIGINT,SigTermHandler);
#ifdef WIN32
  signal(SIGBREAK,SigTermHandler); // Ctrl-C
#else // UNIX
  signal(SIGHUP,SigTermHandler);
  signal(SIGQUIT,SigTermHandler);
#endif
  signal(SIGFPE,SigErrHandler);
  signal(SIGSEGV,SigErrHandler);
  signal(SIGILL,SigErrHandler);
//#endif
} //____________________________________________________


opt::OptSection::OptSection
  (const char *n,const char *d,const OptContainer & c,
  Config * pc, const int pr)
    :name(n),descr(d),options(c),pCfg(pc),priority(pr) {
  for_each_in(c,io,OptContainer::const_iterator)
    optionsAssoc.insert(make_pair((*io)->name,*io));
  pCfg->sections << this;
} //____________________________________________________

opt::Config * opt::SolverCfg() {
  static auto_ptr<opt::Config> pC
    = make_auto_ptr(new opt::Config(CFG_NAME__));
  return pC.get();
}
void opt::ReadOptions() {
  SolverCfg()->ReadOptions();
}
void opt::WriteOptions() {
  if (writeParams)
    SolverCfg()->WriteOptions();
}

void opt::Config::ReadOptions() {
  SetDefault();
  ifstream is(fileName);
  if (!is) {
    PRINT_LOGLN(fileName << " will be created after work.");
    writeParams = true;
    return;
  }
//  is.setf(ios::skipws,ios::skipws);
  char buf1[1024], buf2[1024];
  while (is) {
    ReadSectionName(is,buf1); // omitting comments
    if (!is)
      break;
    ReadOptionName(is,buf2);
    if (!is)
      break;
    SectionContainer::iterator iSec
      = sections.find(buf1);
    if (sections.end() != iSec) {
      OptAssocCont::iterator iOpt
        = iSec->second->optionsAssoc.find(buf2);
      if (iSec->second->optionsAssoc.end() != iOpt) {
        iOpt->second->read(is);
        if (!is) {
          PRINT_ERROR(fileName<<": error reading option "
            <<buf1<<'.'<<buf2);
          is.ignore(INT_MAX,'\n');
        }
      } else {
        PRINT_ERROR(fileName<<": option not found: \""
          <<buf2<<"\"");
        is.ignore(INT_MAX,'\n');
      }
    } else {
      PRINT_ERROR(fileName
        <<": option section not found: \""<<buf1<<"\"");
      is.ignore(INT_MAX,'\n');
    }

    is.ignore(INT_MAX,'\n');
  }
  if (!is.eof())
    PRINT_ERROR(fileName<<": could not read")
} //____________________________________________________
void opt::Config::SetDefault() {
  for_each_in(sections,iSec,
    SectionContainer::iterator)
    for_each_in(iSec->second->options,iOpt,
      OptContainer::iterator)
      (*iOpt)->SetDefault();
} //____________________________________________________
void opt::Config::WriteOptions() {
 ofstream os(fileName);
#ifdef _WIN32
  os << left;
#endif
  os << "'Options. Delete a line and run BCP without params to restore default\n\n";
  for_each_in(sections,iSec,
    SectionContainer::iterator)
  if (os) {
    os << '[' << iSec->second->name
      << "]: " << iSec->second->descr << '\n';
    OptContainer::iterator iOpt;
    int wMax=0;
    for_each_in(iSec->second->options,iOpt,)
			wMax = IMax(wMax,(int)strlen((*iOpt)->name));
    for_each_in(iSec->second->options,iOpt,)
    if (os) {
      os << iSec->second->name << '.';
      os.width(wMax);
      os << (*iOpt)->name << "  ";
      os.width(6);
      (*iOpt)->write(os);
      if ((*iOpt)->descr)
        os << "  \'" << (*iOpt)->descr;
      os << '\n';
      // here could be deleted.

    }
    os << '\n';
  }
//  os.close();
  if (!os)
     PRINT_ERROR(fileName<<": could not write")
} //____________________________________________________
void opt::Config::WriteOptionsShort(ostream & os) {
  for_each_in(sections,iSec,
    SectionContainer::iterator)
  if (os) {
    os << '[' << iSec->second->name
      << "]: ";
    OptContainer::iterator iOpt;
    int wMax=0;
    for_each_in(iSec->second->options,iOpt,)
			wMax = IMax(wMax,(int)strlen((*iOpt)->name));
    for_each_in(iSec->second->options,iOpt,)
    if (os) {
//      os << iSec->second->name << '.';
  //    os.width(wMax);
      os << (*iOpt)->name << ':';
    //  os.width(6);
      (*iOpt)->write(os);
      os << ' ';
      //if ((*iOpt)->descr)
        //os << "  \'" << (*iOpt)->descr;
      //os << '\n';
      // here could be deleted.

    }
//    os << '\n';
  }
//  os.close();
  if (!os)
     PRINT_ERROR(fileName<<": could not write")
} //____________________________________________________
void opt::Config::SkipComments(istream& is) {
  for (;;) {
    char c;
    for (;;) {
      c = (char)is.peek();
      if (isspace(c)) {
        is.get();
        continue;
      } else
        break;
    };    
    if (('\''==c) or ('#'==c) or ('['==c)) {
      is.ignore(INT_MAX,'\n');
      continue;
    } else
      break;
  };
}

void opt::Config::ReadSectionName(istream &is, char*buf) {
  SkipComments(is);
  is.get(buf,1023,'.');
  is.get();
}
void opt::Config::ReadOptionName(istream &is, char*buf) {
  SkipComments(is);
  is.width(1023);
  is >> buf; // until space
} //____________________________________________________

double glbOutputLevel__=5;
double & opt::GlobalOutputLevel() {
  return glbOutputLevel__;
}
bool opt::writeParams=true;;
static opt::OptContainer MyToolsOptions() {
  opt::OptContainer oc;
//  using namespace opt;
  oc 
    << opt::MakeOpt(&opt::writeParams,true,
    "writeParams","Whether this params file will be re-written "
    "before execution")
    << opt::MakeOpt(&glbOutputLevel__, 1,
    "glbOutputLevel","local will be not higher");
  return oc;
} //____________________________________________________

static opt::OptSection mytoolsopt
  ("~sys","system parameters",MyToolsOptions(),
  opt::SolverCfg(), 10000);

/*
int mystat::Accumulator::ww=20;
void mystat::Accumulator::Do(int wh,ostream &os) {
  switch (wh) {
  case 0: reset();
    break;
  case 1:
    if (flags|bPrintAbs) {
      if (flags|bPrintTotal)
        PrintVar(os,ww,nTimes);
      PrintVar(os,ww,accu);
    }
    if (flags|bPrintRatio)
      PrintVar(os,ww,Del0(accu,nTimes));
    break;
  case 2:
    char buf[1024];
    strncpy(buf,name,1020);
    int L=strlen(buf);
    if (flags|bPrintAbs) {
      if (flags|bPrintTotal) {
        strcpy(buf+L,"_N");
        PrintVar(os,ww,buf);
      }
      strcpy(buf+L,"_A");
      PrintVar(os,ww,buf);
    }
    if (flags|bPrintRatio) {
      strcpy(buf+L,"_R");
      PrintVar(os,ww,buf);
    }
    break;
  }
}
*/

ofstream thelog__("log.txt");
//ofstream 
int floatWidth = 10;

} //  END_COMMON_NAMESPACE__
int FUCK;
//____________________________________________________

using namespace COMMON_NAMESPACE__;

#ifdef _TST_MYTOOLS
double xx,xp;
static bgn::opt::OptContainer MyToolsOptions2() {
  using namespace bgn::opt;
  OptContainer oc;
  oc 
    << MakeOpt(&xx,2.55,
    "xx","test1")
    << MakeOpt(&xp,3.28,
    "xp","test2");
  return oc;
} //____________________________________________________

static bgn::opt::OptSection mytoolsopt2
  ("tst","test options section",MyToolsOptions2(),
  bgn::opt::SolverCfg(), 10000);

int main() {
  try {
    COMMON_NAMESPACE__::MyApp myApp;
    for (char i=32;i<=126;i++)
      cout << i;
    cin.get();
//    cin.ignore(INT_MAX,'\n');
    return 0;
  }
  catch (const UserBreak__&)
  { PRINT_LOGLN("User break."); }
  catch (const exception & e)
  { PRINT_LOGLN(e.what()); }
  catch (...)
  { PRINT_LOGLN("Unknown exception."); }
  return 1;
}
#endif
