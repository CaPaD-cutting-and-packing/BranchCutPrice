// main.cpp : Defines the entry point for the console application.
// Author : Gleb Belov <bg32767@yahoo.com>

#include "stdafx.h"
#include "solver.h"
#include "bcp.h"  // uses problem.h
#include "bbcuts.h"
#include "probl_csp1.h"
//#include "probl_csp2.h"
#include "probl_cp22.h"
#include "probl_pmp1.h"
#include "bcp.h"
#include "bcp2.h"
#include "raster.h"
#include "lasthdr.h"

////////////////////////////////////
#include <ctype.h>

using namespace COMMON_NAMESPACE__;

// Option vars: initialize!

const char * Logo =
  "Branch-Cut-Price Algorithm for the 1d cutting stock problem,\n"
  "   the 2d 2-staged (=>guillotine) cutting problem\n"
  "   and Integer Rounding for the 2d cutting stock problem.\n"
  "\nGleb <Belov@math.tu-dresden.de>\n";
//  << bgn::BCP::Version() <<
const char * manual =
  "\nPlease indicate a list of data files,\n"
  "each file containing a list of problem instances.\n"
  "See solver.cfg for parameters, including PROBLEM TYPE.\n"
  "Integer data only.\n"
  "CmdLine options: -in, n=starting instance from a file\n"
  "-ln last instance n; -In single instance n; -cn cut type n; "
  "-pn problem type n (see *.cfg for values); -s save these temp. values to *.cfg; "
  "-h0: first cut along the second dimension (CP22); -mn: max. n cuts in LP "
  "-on: global output level n; -tn: max time per instance n sec. "
  "-^n: MSVC.pow_l=n\n; -dn: PMP1::deltaKpercent=n"
  "1D CSP input format:\n"
  "[<instance name/text line(s) comment> \\newline]\n"
  "m L l1 b1 ... lm bm \n"
  "2D input format:\n"
  "[<instance name/text line(s) comment> \\newline]\n"
  "L W m\nl1 w1 b1\n ...\nlm wm bm (upper-level cut in L-direction)\n"
  "  or:  L W m\nl1 w1 c1 b1\n ...\nlm wm cm bm (weighted)\n"
  //"Warning: do not add level cut if not gcd(c_i)\\in\\Z\n"
  ;

int iBr; // for testing

int main(int argc, char** argv)
{

int ret=0;
#ifdef DBG_ON
cout << "Warning: DBG_ON defined!!!"<<endl;
cout << "Capacity of size variables: sizeof(size)=" << sizeof(bgn::size)
 << ", SIZE_MAX__="<<bgn::SIZE_MAX__ << endl;
#endif

bgn::Solver::RWOptions(); // even if no files
if (bgn::BCP::fCheckCuts)
  cout << "Warning: BCP::fCheckCuts set in options!!!"<<endl;
if (bgn::BBCuts::fCheckBnd)
  cout << "Warning: BBCuts::fCheckBnd set in options!!!" << endl;
if (bgn::PMP1::pF * 1e6 != round(bgn::PMP1::pF * 1e6)) {
  cout << "bgn::PMP1::pF was read as "<<bgn::PMP1::pF;
  bgn::PMP1::pF = round(bgn::PMP1::pF * 1e6) / 1e6;
  cout << " Setting to "<<bgn::PMP1::pF << endl;
}
if (bgn::PMP1::pV * 1e6 != round(bgn::PMP1::pV * 1e6)) {
  cout << "bgn::PMP1::pV was read as "<<bgn::PMP1::pV;
  bgn::PMP1::pV = round(bgn::PMP1::pV * 1e6) / 1e6;
  cout << " Setting to "<<bgn::PMP1::pV << endl;
}

try
{
  if (argc<2) {
    PRINT(Logo);
    PRINT(manual);
    cout << "Capacity of size variables: SIZE_MAX__=" << bgn::SIZE_MAX__
      <<", DEMAND_MAX__="<<bgn::DEMAND_MAX__<< endl;

    return 0;
  }
  bgn::Solver slv;
  int i1;
  double d1;
  char * ts;
  for (int i=1;i<argc;i++)
  try
  {
    if (('-'!=argv[i][0]) and ('/'!=argv[i][0])) {
      for (iBr=0;iBr<IMin(3,1+2*bgn::Solver::fTestParams);++iBr) {
        if (bgn::Solver::fTestParams) {
        log_ln("PERFORMING 3*8 test-param runs .");
        switch (iBr) {
          case 0:
            bgn::BCP::nBranch = 0;
            { ofstream os1("allfiles.res", ios::app);
            os1 << "\nBranching on G&G vars: only b&b col.gen., Ax>=b formulation\n"; }
            break;
          case 1:
            bgn::BCP::nBranch = 1;
            bgn::CSP1::nHP = 0;
            { ofstream os2("allfiles.res", ios::app);
            os2 << "\nBranching on CVRP vars: only Ax>=b formulation\n"; }
            break;
          case 2:
            bgn::BCP::nBranch = 1;
            bgn::CSP1::nHP = 1;
            { ofstream os3("allfiles.res", ios::app);
            os3 << "\nBranching on AFF vars: only Ax=b formulation\n"; }
            break;
        }
        }
        slv.ProcessFile(argv[i]);
      }
      bgn::Solver::RWOptions(); // after each file
    }
    else {
      switch (argv[i][1]) {
      case 's': // Save options
        bgn::Solver::WriteOptions();
        break;
      case 'i': //  Starting instance
        ts = argv[i]; ++++ts;
        istringstream(ts) >> i1;
        log_ln("Starting from instance "<<i1);
        bgn::Solver::instFirst = i1;
        bgn::Solver::instLast = IMax(bgn::Solver::instLast,i1);
        break;
      case 'I': // A single instance
        ts = argv[i]; ++++ts;
        istringstream(ts) >> i1;
        log_ln("Processing only instance "<<i1);
        bgn::Solver::instFirst = i1;
        bgn::Solver::instLast = i1;
        break;
      case 'l': //  Last instance
        ts = argv[i]; ++++ts;
        istringstream(ts) >> i1;
        log_ln("Last instance "<<i1);
        bgn::Solver::instFirst = IMin(bgn::Solver::instFirst,i1);
        bgn::Solver::instLast = i1;
        break;
      case 'c':
        ts = argv[i]; ++++ts;
        istringstream(ts) >> i1;
        log_ln("Setting cut type "<<i1);
        bgn::BCP::cutType = i1;
        break;
      case 'p':
        ts = argv[i]; ++++ts;
        istringstream(ts) >> i1;
        log_ln("Setting problem type "<<i1);
        bgn::Solver::problemType = i1;
        break;
      case 'h':
        ts = argv[i]; ++++ts;
        istringstream(ts) >> i1;
        log_ln("Setting first cut (CP22) along first dimension: "<<i1);
        bgn::CP22::fFirstCut1stD = i1;
        break;
      case 'm':
        ts = argv[i]; ++++ts;
        istringstream(ts) >> i1;
        log_ln("max Cuts: "<<i1);
        bgn::BCP::maxCuts__ = i1;
        break;
      case 'o':
        ts = argv[i]; ++++ts;
        istringstream(ts) >> d1;
        log_ln("Glb out lev: "<<d1);
        opt::GlobalOutputLevel() = d1;
        break;
      case 't':
        ts = argv[i]; ++++ts;
        istringstream(ts) >> i1;
        log_ln("Time limit per instance, sec: "<<i1);
        bgn::BCP::TimeLimit__ = i1;
        break;
      case '^':
        ts = argv[i]; ++++ts;
        istringstream(ts) >> d1;
        log_ln("MSVC.pow_l: "<<d1);
        bgn::CSP1::MSVC::pow_l = d1;
        break;
      case 'd':
        ts = argv[i]; ++++ts;
        istringstream(ts) >> d1;
        log_ln("PMP1.deltaKpercent = "<<d1);
        bgn::PMP1::deltaKpercent = d1;
        break;
      }
    }
  }
  catch (const exception & e) {
    PRINT_LOGLN(e.what()<<flush); ret += 10;
  }
  catch (...)
  { PRINT_LOGLN("Unknown exception.");
    ret += 100;
    return ret; }

  return ret;
}
  catch (const UserBreak__&) {
    PRINT_LOGLN("User break.");
    ret +=2;
  }
  catch (const exception & e)
  { PRINT_LOGLN(e.what()<<flush); ret += 1; }
  catch (...)
  { PRINT_LOGLN("Unknown exception."); ret += 3; }
  return ret;
}


namespace bgn {                     // SS_BEGIN_NAMESPACE__


void Solver::RWOptions()
{
  opt::SolverCfg()->ReadOptions();
  if (opt::writeParams)
    opt::SolverCfg()->WriteOptions();
}

void Solver::WriteOptions()
{
  opt::SolverCfg()->WriteOptions();
}

/// probl_others must be also included if used in the proc
int ProblemLoader::Read(istream &ifs) {
  long long n1; int i;
  prNm[0] = '\0';
  for (i=0;i<10;i++) { // Omitting beginning text
    ifs >> n1;
    if (ifs) break; // success
    if (ifs.eof()) return 3; // eof
    ifs.clear();
//    if (0==i)
      ifs.get(prNm, sizeof(prNm)); // not read eol
    ifs.ignore(INT_MAX,'\n');
  }
// Now not only CSP1 but also CSP2
  if (Solver::problemType==1 || Solver::problemType==7)
    pr = new CSP1(penv, infile,inst,prNm);
  else
//  if (Solver::problemType==2)
//    pr = new CSP2(penv, infile,inst,prNm);
//  else
  if (Solver::problemType==3)
    pr = new CP22(penv, infile,inst,prNm);
  else
  if (Solver::problemType==4
    or Solver::problemType==6)
    pr = new PMP1(penv, infile,inst,prNm);
  else
    assertm(0, Solver::problemType<<": unknown problem type");

  return pr->Read(ifs,n1);
}

BCP * ProblemLoader::GetSolver() {
  if (Solver::problemType != 4)
    return new BCP(pr);
  return new BCP2(pr);
}


void Solver::ProcessFile(const char *fln) {
  int inst;
// SIMPLE STATISTICS:
  int nInst = 0;
  int nOpt = 0, nFeas=0, nErr = 0, nInfeas = 0, nLPErr=0;
  double LP0=0, LP=0, IP=0;
  double tLP0 = 0, tIP = 0, tIPBest=0, tLPBest=0, tOpt=0;
  double tCG1=0, tCG2=0, tRnd=0;
  double nIter = 0, nLPCol=0, nIPCol =0, nIterAll=0;
  double nNodes = 0;
  double nMIPSolBest = 0;
  int nnIRUP = 0;
  int nTooLongCG = 0;
  nErr2__ = 0;
  // double lOpt=0, lInst=0, lErr=0, lErr2=0;

  double s_m=0, IPmin=1e100, IPmax=-1e100;
  mystat::VectorStat<double> stat, stat_ave;

      for (int iBas=0;iBas<IMin(4,1+3*bgn::Solver::fTestParams);++iBas) {
        if (bgn::Solver::fTestParams)
          bgn::CSP1::nStartBasis = (iBas<2)? iBas: iBas+1;
        for (int iCG=0;iCG<IMin(2,1+bgn::Solver::fTestParams);++iCG) {
          if (0==iBr) { // G&G
            if (bgn::Solver::fTestParams)
//              if (bgn::BCP::fLocalUB) //continue;
                bgn::BCP::BrVarFrac = 0.5f + 1.5f*(float)iCG;

          }
          else
          if (bgn::Solver::fTestParams)
            bgn::CSP1::fCompareBySPPRC_DP = iCG*2; // but not for local UB!!
//          for (int iEq=0;iEq<IMin(2,1+bgn::Solver::fTestParams);++iEq) {
            if (bgn::Solver::fTestParams) {
              if (iBr<2)
                bgn::CSP1::fModelEquality = 0;//iEq*2;
              else
                bgn::CSP1::fModelEquality = 1;//iEq*2;
              log_ln("\n\n  TEST-PARAM RUN: "<<iBr<<iBas<<iCG//<<iEq
                <<'\n');
            }


  ifstream ifs(fln);
  if (!ifs) goto IOErr;
//   strncpy(__glb_file,fln,1024); // for err reporting
	__glb_file = fln;
//  InitFile(); // Done in BCP
  for (inst=1; ifs && (inst<=instLast); ++inst) {
    Env env;
    ProblemLoader pl(env,fln,inst);
    __glb_inumber = inst; // for error reporting
    int readRes=pl.Read(ifs);
    if (3==readRes) goto Finish; // eof
    if (readRes>0)
//      if (2==readRes) // I/O err
//        goto IOErr; else
      goto Finish; // bad format=>next file
    if (inst < instFirst) continue;
    try {
      if (6 == Solver::problemType) {
        ProcessPr6(&pl);
        continue;
      }
      if (7 == Solver::problemType) {
        ProcessPr7(&pl);
        continue;
      }
      auto_ptr<BCP> bcp(pl.GetSolver());
      bcp->firstinstance = (inst==instFirst);
      bcp->Run();
      ++ nInst;
      if (BCP::opt==bcp->status)
      { ++ nOpt; tOpt += bcp->timeIP; }
      else if (BCP::error==bcp->status) ++nErr;
      else if (BCP::infeas==bcp->status) ++ nInfeas;
      else if (BCP::LPErr==bcp->status) ++nLPErr;
      if (BCP::opt == bcp->status or BCP::feas==bcp->status) {
        ++nFeas;
        LP0 += bcp->glv0;
        LP += bcp->glv;
        IP += bcp->gub;
        if (bcp->gub < IPmin)
          IPmin=bcp->gub;
        if (bcp->gub > IPmax)
          IPmax=bcp->gub;
      }
      s_m += bcp->pr->the_m();
      tLP0 += bcp->timeLP0;
      tIP += bcp->timeIP;
      tIPBest += bcp->timeIPBest;
      tLPBest += bcp->timeLPBest;
      tCG1 += bcp->timeCG1;
      tCG2 += bcp->timeCG2;
      tRnd += bcp->timeRnd;
      nIter += bcp->iter; // +1?
      nLPCol += bcp->lpCols;
      nIPCol += bcp->ipCols;
      nIterAll += bcp->iterAll;
      nNodes += bcp->cntInitNode;
      nTooLongCG += bcp->nTooLongColGen;
      if (bcp->fMIPSolBest) ++ nMIPSolBest;
      if (bcp->glb0 < bcp->glb) ++nnIRUP;

//      stat.FirstTimeSetDivisorForEach(nInst); // ????
      stat.NextRun() << (nOpt) << (nInst) << (nErr) << nErr2__ << tLP0 << tIP << nNodes;


      ofstream ofs ("stattemp.txt");

	  if (nInst > 0)
  {
//  ofs.precision(4);
  ofs <<'\n' << fln;
  if (instFirst>1 || instLast<10000)
    ofs << " Inst:"<<instFirst<<'-'<<instLast
	<<" (all "<<nInst<<") ";
  ofs
    << " Feas:" << nFeas
    << " Opt:" << double(nOpt)/nInst*100 <<"% ("
    << nOpt<<'/'<<nInst <<')';
  if (nErr) ofs << " nErr:" << nErr;
  if (nErr2__) ofs << " nErr2__:" << nErr2__;
  if (nInfeas) ofs << " nInfeas:" << nInfeas;
  if (nLPErr) ofs << " nLPErr:" << nLPErr;

  ofs << " AVE: m:"<<s_m/nInst;

  if (nFeas)
    ofs << " LP0: "<<LP0/nFeas<<" LP: "<<LP/nFeas
      <<" IP: "<<IP/nFeas
      << " IPmin:"<<IPmin<<" IPmax:"<<IPmax;

  ofs << " tLP: " << (tLP0/nInst) << " tIP: " << (tIP/nInst)
    << " tLPBest: " << (tLPBest/nInst)
    << " tIPBest: " << (tIPBest/nInst);
  ofs << " tOpt: " << ((nOpt) ? tOpt/nOpt: 1e100)
    << " tCGRoot: " << (tCG1/nInst)
    << " tCGBranch: " << (tCG2/nInst)
    << " tRnd: " << (tRnd/nInst);
  ofs << " nnode: " << (nNodes/nInst);
  ofs <<" Iter: "<<(nIter/nInst)
    << " LPC: "<<nLPCol/nInst << " IPC: "<<nIPCol/nInst
    <<" IterAll: "<<(nIterAll/nInst)
    <<" nMIPBest: "<<nMIPSolBest
    <<" nBetterLP: "<<nnIRUP
    <<" nTooLongColGen: "<<nTooLongCG
   <<'\n';

      bcp->PrintStatistics(ofs);
  ofs << '\n';

  time_t time0=time(0);
  ofs << " Date "<<ctime(&time0);

//  BCP::PrintOptions(ofs);
  opt::SolverCfg()->WriteOptionsShort(ofs);
  ofs << '\n';
	  }
    } catch (const exception & e) {
      PRINT_ERROR(e.what()); // How to note errors
    }
  }

Finish:

      stat_ave.FirstTimeSetDivisorForEach(nInst); // ????
      stat_ave.NextRun() << (nOpt) << (nInst) << (nErr) << nErr2__ << tLP0 << tIP << nNodes;

  {
    ofstream os("allfiles.res", ios::app|ios::out);
//    if (1== stat.getN())
  //    os >>
    os <<iBr<<iBas<<iCG//<<iEq
      <<'\t';
    os.precision(4);
    stat.print_diff(os,'\t');
    os << '\n';
  }

  IOErr:
  if (not (ifs or ifs.eof()))
//  DoneFile();
    PRINT_ERROR(fln<<": "<<strerror(errno));

  }}//} // the 3 test for's
  if (nInst < 1) goto SkipStat;
  {
  ofstream ofs("allfiles.res", ios::out|ios::app);
//  ofs.precision(4);

  ifstream its("stattemp.txt");
  //char c;
  while (its) ofs.put((char)its.get()); // to clean input stack
  ofs << '\n';

  if (fTestParams) {
    ofs << "\nNumber of instances: " << stat.getN()<< '\n';
    ofs << "\tnOpt\t" << "nInst\t" << "nErr\t" << "nErr2\t"
     << "tLP\t" << "tIP\t" << "nNodes\n";

    ofs.precision(4);
    ofs<<"AVE:\t"; stat.print_ave(ofs,'\t');
    ofs<<"\nDEV:\t"; stat.print_dev(ofs,'\t');
    ofs<<"\naveLOG:\t"; stat.print_ave_log(ofs,'\t');
    ofs<<"\ndevLOG:\t"; stat.print_dev_log(ofs,'\t');
    ofs<<'\n';
  }

  }

SkipStat:
  return;
} //____________________________________________________



struct ID1 {
  int i; double d;
  bool operator <(const ID1 & id) const
  { return d < id.d; }
};

void Solver::ProcessPr6(ProblemLoader *pl) {
  if (nTries < nRounds) {
    nTries = nRounds;
    log_n(1,"Setting nTries=nRounds");
  }
  int iRound,i,iGroup,iTry,
    n1; // - current group size
  PMP1 * pr = dynamic_cast<PMP1*>(pl->pr);
    // - should be really that
  i_vec iTypes;
  Vector<i_vec> iPart(nTries);
  CSP1::MSVC::SolContainer sAll; // solutions for the whole
  CSP1::MSVC::SolContainer::iterator
    its, its1, itsOS, itsND, itsEQ;
  rndGrp.setSeed(0.5); // for each 1st round !!!
//  char fln[1024];
  d_vec NZ(nTries);
  Vector<ID1> vID(nTries),
    spl(pr->m);
  // STAT:
  static double s_no=0, s_nd=0, s_nz=0,
    s_noAll=0, s_ndAll=0, s_nzMin=0, s_n=0, tmTries=0;
  static int nCspNotOpt=0, grszMin=INT_MAX, grszMax = 0;


  // LOOKING FOR GOOD SPLITTINGS: // + note max/min group size
  Timer tryTimer;
  tryTimer.start();
  for (iTry=0;iTry<nTries;++iTry) {
    // CREATING A SPLITTING:
    iPart[iTry].resize(pr->m);
    // THE FIRST METHOD:
    for (i=0;i<pr->m;++i)
      iPart[iTry][i] = int (rndGrp * nParts);
    // The SECOND METHOD: (smth. worse)
    /// SPLITTING OVER

    NZ[iTry] = 0;
    for (iGroup=0;iGroup<nParts;++iGroup) {
      string fln;
      {
      ostringstream(fln)
        << "part" << iGroup << '_' << iTry << ends;
      ofstream os(fln.c_str());
      os << pr->L << ' '
        << (n1=(int)count(iPart[iTry].begin(), iPart[iTry].end(), iGroup)) << '\n';
      if (n1 < grszMin) grszMin = n1;
      if (n1 > grszMax) grszMax = n1; // / AND PRINT
      for (i=0;i<pr->m;++i) {
        if (iGroup == iPart[iTry][i])
          os << pr->pc[i].l << ' ' << pr->pc[i].b << '\n';
      }
      } // closing os
      log_n(1, "Try "<<iTry<<", PARTIAL PROBLEM No. " << iGroup
        <<" N PRODUCT TYPES: "<<n1);
      string prNm;
      ostringstream (prNm) << pl->prNm
        << ':' << iGroup << '_' << iTry << ends;
      CSP1 *csp1 = new CSP1(pl->penv, fln.c_str(), pl->inst, prNm.c_str());
      ifstream ifs(fln.c_str());
      ifs >> csp1->L;
      csp1->Read(ifs, csp1->L);
      BCP bcp(csp1);
      bcp.firstinstance = //(0==ind) &&
        (pl->inst==instFirst);
      bcp.Run();
      // SAVING RESULTS OF THE PARTIAL PROBLEM:
      NZ[iTry] += bcp.gub;
      nCspNotOpt += not (BCP::opt == bcp.status);
    }
    log_n(1,"Sum IP: "<<NZ[iTry]);
    vID[iTry].i = iTry;
    vID[iTry].d = NZ[iTry];
  }
  tryTimer.stop();
  tmTries += tryTimer.userTime();
  // SORTING OUT THE BEST SPLITTINGS:
  sort(vID.begin(), vID.end());
  for (iRound=0;iRound<nRounds;++iRound) {
    list<CSP1::MSVC::SolContainer> sols;
    for (iGroup=0;iGroup<nParts;++iGroup) {
      string fln;
      ostringstream(fln)
        << "part" << iGroup << '_' << vID[iRound].i << ends;
      log_n(1, " ROUND "<<iRound<<", PARTIAL PROBLEM No. " << iGroup
        //<<" N PRODUCT TYPES: "<<n1
        );
      string prNm;
      ostringstream (prNm) << pl->prNm
        << ':' << iGroup << '_' << vID[iRound].i << ends;
      PMP1 *pmp1 = new PMP1(pl->penv, fln.c_str(), pl->inst, prNm.c_str());
      ifstream ifs(fln.c_str());
      ifs >> pmp1->L;
      pmp1->Read(ifs, pmp1->L);
      BCP2 bcp2(pmp1);
      bcp2.firstinstance = //(0==iRound) &&
        (pl->inst==instFirst);
      double timeLimit = bcp2.TimeLimit__;
      bcp2.TimeLimit__ = splitTimeLimit;
      minZ = 1e100; // solving partial
      bcp2.Run();
      bcp2.TimeLimit__ = timeLimit;
      // SAVING RESULTS OF THE PARTIAL PROBLEM:
      sols.push_back(pmp1->ssvc.sols);
      assertm(not pmp1->ssvc.sols.empty(),
        "Turn on sequencing in PMP1, no sequenced solutions found");
        // + save times?
    }
    log_n(1,"Sum IP: "<< vID[iRound].d);
    ////////////////////////////////////////////////////
    // Constructing Pareto-best sols for the whole problem:
    ////////////////////////////////////////////////////
    // Calculate min/max nOS:
    double minOS = 1e100, maxOS=0;
    double nos;
    list<PMP1::SSVC::SolContainer>::iterator igs;
    PMP1::SSVC::SolContainer::iterator iss, issBest;
    for_each_in(sols,igs,)
      for_each_in(*igs,iss,) {
      if (iss->no > maxOS) maxOS = iss->no;
      if (iss->no < minOS) minOS = iss->no;
    }
    log_n(2,"All parts: nOS between "<<minOS<<" and "<<maxOS);
    // Explicitly consider 2 solutions with min os / min nd ?
    double nosLast=0;
    for (nos=minOS;nos<=maxOS;++nos) {
      double os1=0; // actual nos, ie max {os <= nos}
      double nDD=0; // sum of N patterns
      double nZ=0; // material input
      for_each_in(sols,igs,) { // for all subproblems
        double nD1 = 1e100;
        issBest = igs->begin();
        for_each_in(*igs,iss,) { // for all pareto-best in subpr.
          if (iss->no <= nos) { // if a valid N open stacks
            if (iss->nd < nD1
              or (iss->nd == nD1 and iss->nz < issBest->nz)) { // if a better N patterns
              nD1 = iss->nd;
              issBest = iss;
            }
          }
        }
        if (nD1 >= 1e50) // for this subproblem,
        { os1=-1; break; } //  no sol. with no <= nos found
        if (os1 < issBest->no) os1 = issBest->no; // max n open
        nDD += nD1;
        nZ += issBest->nz;
      }
      if (os1 != nosLast and os1 > 0) {// only if a new NOS value
        bool fBetter = true;
        for_each_in(sAll,its,) // SEARCH NOT WORSE:
          if (its->no <= os1
            and its->nz <= nZ and its->nd <= nDD) {
            fBetter = false;
          }
        // KILLING PARETO-WORSE:
        its = sAll.begin();
        if (fBetter) // ONLY THEN DELETE WORSE
        while (its != sAll.end()) {
          its1= its;
          ++ its1;
          if (its->no >= os1
            and its->nz >= nZ and its->nd >= nDD) {
            sAll.erase(its);
//            fBetter = 1; // - ERROR
          }
          its = its1;
        }
        // ADDING CURRENT
        if (fBetter) {
          PMP1::SSVC::Solution s;
          s.nz = nZ; s.no = os1; s.nd = nDD;
          sAll.insert(s);
          log_n(2,"Combined sol. found: ("<<os1<<' '
            <<nDD<<' '<<nZ<<')');
        }
      }
      nosLast = os1;
    }
  }
  // Investigate the set of pareto-best combined solutions:
  double minOS=1e100, minND=1e100;
  minZ=1e100; // static
  for_each_in(sAll,its,) {
    if (its->no < minOS) {
      minOS = its->no;
      itsOS = its;
    }
    if (its->nd < minND) {
      minND = its->nd;
      itsND = its;
    }
    if (its->nz < minZ) minZ = its->nz;
  }
  assertm(minOS<1e100 and minND<1e100,
    "min OS and min ND must be findable. solsAll.size()=="
    <<sAll.size());
  // Find the "equivocal" Pareto-best:
  double minFF = 1e100;
  for_each_in(sAll,its,) {
    if (its->no * minND * minZ
      + its->nd * minOS * minZ
      + its->nz * 20 * minND * minOS < minFF) {
      minFF
        = its->no * minND * minZ
        + its->nd * minOS * minZ
        + its->nz * 20 * minND * minOS;
      itsEQ = its;
    }
  }
  if (OUTP_LEV__ >= 2) {
    log__("The 3 solutions: ");
    log__(" ("<<itsOS->no<<' '<<itsOS->nd<<' '<<itsOS->nz<<')');
    log__(" ("<<itsND->no<<' '<<itsND->nd<<' '<<itsND->nz<<')');
    log__(" ("<<itsEQ->no<<' '<<itsEQ->nd<<' '<<itsEQ->nz<<')');
    log_ln("");
  }
  // SOLVE THE WHOLE PROBLEM
  log_n(1,"SOLVING THE WHOLE PROBLEM.");
  BCP2 bcp2(pr);
  bcp2.firstinstance = (pl->inst == instFirst);
  bcp2.Run();
  // Printing results for the problem:
  char outfile1[1024];
    strncpy(outfile1,pr->infile,sizeof(outfile1)-5);
    strcat(outfile1,".txt");
    ofstream ofs(outfile1,ios::out|ios::app); // Here app
    ofs
      << " The 3: ("<<itsOS->no<<' '<<itsOS->nd<<' '<<itsOS->nz<<')'
      << " ("<<itsND->no<<' '<<itsND->nd<<' '<<itsND->nz<<')'
      << " ("<<itsEQ->no<<' '<<itsEQ->nd<<' '<<itsEQ->nz<<')'
      ;
    ofs << " ALL COMBINED:";
  for_each_in(sAll,its,)
    ofs
      << " ("<<its->no<<' '<<its->nd<<' '<<its->nz<<')'
      ;
  ofs << '\n';
  // Statistics:
  if (pl->inst==instFirst) {
    s_no= s_nd= s_nz= s_noAll= s_ndAll= s_nzMin= s_n = tmTries =0;
    nCspNotOpt=0; grszMin=INT_MAX; grszMax = 0;
  }
  s_no += itsEQ->no;
  s_nd += itsEQ->nd;
  s_nz += itsEQ->nz;
  s_noAll += pr->ssvc.nOpenMaxMin;
  s_ndAll += pr->ssvc.nDiffPatMin;
  s_nzMin += pr->cspIP;
  ++ s_n;
  char outstat[1024];
    strncpy(outstat,pr->infile,sizeof(outfile1)-5);
    strcat(outstat,".sta");
    ofstream osta(outstat); // Here app
    osta
       << s_n << ' '
       << s_no/s_n << ' ' << s_nd/s_n << ' ' << s_nz/s_n
       << ' ' << s_no/s_noAll << ' ' << s_nd/s_ndAll << ' ' << s_nz/s_nzMin
       << ' ' << tmTries/s_n
       << ' ' << nCspNotOpt/s_n << ' ' <<  grszMin << ' ' << grszMax
       <<"\nN instances; AVE: no, nd, nz for EQUIVOCAL solution;\n"
       "and relations to that (each independently best) for the whole problem;\n"
       "ave time for all tries, ave nCSPNotOpt, group size min/max";
}


void Solver::ProcessPr7(ProblemLoader *pl) {
    int i;
    CSP1* pr1 = dynamic_cast<CSP1*>(pl->pr);
    int m = pr1->m;
    Vector<int> sz(m), bi(m), rp1;
    for (i=0;i<m;++i) {
      assertm(pr1->pc[i].l <= INT_MAX,
              "For this problem type, only 'int' sizes are possible now.");
      sz[i] = (int)pr1->pc[i].l; // IMPORTANT: w's. 64->32 bits? TODO
    }
    for (i=0;i<m;++i) bi[i] = pr1->pc[i].b;
    int L1 = minL;
    int L2 = maxL;
    double slb=0;
    for (i=0;i<m;++i)  slb += sz[i]*bi[i];

    if (L1<=0)
      for (i=0;i<m;++i) if (L1 < sz[i]) L1 = sz[i];
    if (L2<=0)
      for (i=0;i<m;++i)  L2 += sz[i];
    cout << "CSP1 with unknown L. Only the 1st instance of the file will be processed. Constructing raster points..." << flush;
    ConstructRP(sz,bi,L2,rp1);
    int p, p1 = FindRPUnder(L1,rp1);
    cout << " N raster points="<<rp1.size()
      << "\nPLEASE NOTE: minimal unit=1, so divide b[i] in the input by min.partial size!"<< endl;
//    auto_ptr<BCP> bcp(pl->GetSolver());
//    pr1 = dynamic_cast<CSP1*>(bcp->pr); remains the same
    for (p=p1;p<rp1.size();++p) {
      CSP1 *csp1 = new CSP1(pl->penv, pl->infile, pl->inst, pl->prNm);
      ifstream ifs(pl->infile);
      ifs >> csp1->L;
      csp1->Read(ifs, csp1->L);
      auto_ptr<BCP> bcp(new BCP(csp1));
      bcp->firstinstance = (p==p1); // ...
      csp1->L = csp1->L0 = rp1[p];
      bcp->Run();
      cout << "Solved with L="<<rp1[p];
      if (BCP::opt == bcp->status)
        cout << ": optimal. ";
      if (BCP::error == bcp->status)
        cout << ": error!! ";
      cout << bcp->gub << ", " << slb/(rp1[p]*bcp->gub)*100 << "%, LPB="<<bcp->glb << ", time="<<bcp->timeIP;//<<endl;
//      bcp->lp->Close();
    }
    cout << "If option BCP::fPrintIntSol active, find solutions in the file ...__.txt" << endl;
    delete pr1;
}

//namespace bgn {                     // SS_BEGIN_NAMESPACE__



void Solver::CalculateTimeUnit() {
} //____________________________________________________


// options: instFirst, instLast
opt::OptContainer Solver::Options() {
  opt::OptContainer oc;
  oc
    << opt::MakeOpt(&problemType, 1,
      "problemType", "int: 1: CSP1, 2: CSP2, 3: CP22, 4: PMP1, 6: PMP1+MOSP1, 7: CSP1 with search for opt. L")
    << opt::MakeOpt(&instFirst, 1,
      "instFirst", "The 1st instance in the file is indexed by 1")
    << opt::MakeOpt(&instLast, INT_MAX,
      "instLast",  "")
    << opt::MakeOpt(&fTestParams, 0,
      "fTestParams",  "For each of 3 branching methods, do 16 runs ...")
    << opt::MakeOpt(&nTries, 30,
      "nTries",  "PMP+MOSP: how many different problem separations to try")
    << opt::MakeOpt(&nRounds, 3,
      "nRounds",  "PMP+MOSP: how many best problem separations to solve")
    << opt::MakeOpt(&nParts, 3,
    "nParts",  "PMP+MOSP: how many parts of a problem")
    << opt::MakeOpt(&splitTimeLimit, 50,
    "splitTimeLimit",  "PMP+MOSP: time limit for a partial problem")
    << opt::MakeOpt(&maxL, 0,
    "maxL",  "CSP1 with unknown L, ==0 means sum{l[i]} without b[i]")
    << opt::MakeOpt(&minL, 0,
    "minL",  "CSP1 with unknown L, ==0 means max{l[i]}")
    << opt::MakeOpt(&outputLevel, DEF_OUTP_LEVEL,
      "outputLevel", "");
  return oc;
} //____________________________________________________
opt::OptSection Solver::opt
  ("Solver", "The solver framework",
  Solver::Options(), opt::SolverCfg(), 5000);

int Solver::problemType = 1;
int Solver::instFirst=1;
int Solver::instLast=INT_MAX;
int Solver::fTestParams;

int Solver::nTries, Solver::nRounds, Solver::nParts;
int Solver::minL, Solver::maxL;
double Solver::splitTimeLimit;
double Solver::minZ;

double Solver::outputLevel=DEF_OUTP_LEVEL;

}                                                // SS_END_NAMESPACE__


/////////////////////////////////////////////////
//////// ...................

//int __ctype_b[256];
//int __ctype_tolower[256];
//int __ctype_toupper[256];


//    for (char i=32;i<=126;i++)
//      cout << i;
//    cin.get();
