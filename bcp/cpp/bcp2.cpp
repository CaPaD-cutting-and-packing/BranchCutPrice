// FILE: bcp2.cpp, branch&cut&price for 1D pattern minimization
// Author: Gleb <Belov@math.tu-dresden.de>

#include "stdafx.h"
#include "bcp2.h"
#include "solver.h"
#include "lasthdr.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

char * BCP2::Version() {
  return
    "BCP for PMP1: Version 1.0.\n"
    "\\Delta K also < \\infty; branching2"
    "\nBuild "__DATE__"  "__TIME__;
} //____________________________________________________

void BCP2::InitRun() {
  try {
  lp->Close();

  if (!CheckData())
    return;
  SolveCSP();
  SolveBPP(); // Lower bnd. Use it .
//  AdjustOptions();
//  int nBr = nBranch;
  nBranch = 0;
  lp->Open();
  double gub__ = gub;
  BCP::InitRun();
  gub = pr->raster_ceil(gub__);
  pr->SetHeurBnd(gub);
  dbg_outn(4,"gub0 = "<<gub);
  //  nBranch = nBr;
  if (firstinstance) {
    NCSP = ND0 = ND = 0; nFeas = 0;
  }
  } catch (const exception& e) {
    PRINT_ERROR(e.what());
    PRINT_LOG("Error while init BCP_PMP. Continuing...");
  }
}

bool BCP2::CheckData() {
  // if the objective equivalent to CSP, solve CSP and ret
  // if a BPP, solves it as such and returns 0
  return 1;
}

void BCP2::SolveCSP() {
  CSP1 *csp1 = new CSP1(*(dynamic_cast<CSP1*>(pr)));
    // creating a CSP1 problem;
  dbg_outn_(1, "Solving CSP (m="<<csp1->Dim()<<")...");
  BCP bcp(csp1);
  bcp.firstinstance = 0; // firstinstance;
  bcp.fSubproblem = 1;
//  firstinstance = 0; // no more log init
  double oL = opt::GlobalOutputLevel();
    opt::GlobalOutputLevel()=pr->CSPOutpLev;
  double TL = BCP::TimeLimit__;
    BCP::TimeLimit__=pr->CSPTiLim;
  bcp.Run();
  fCSPOpt = (BCP::opt == bcp.status);
  pr->cspIP = bcp.gub; // what if err?
  pr->cspLP = bcp.glb; // what if err?
  pr->deltaK = 0; // for the time
  // COUNTING N DIFFERENT PATTERNS IN THE CSP SOL: // save it ?
  // need it ?
  pr->zi = 1e100;
  pr->b = csp1->b; // for SSVC
  assertm(not csp1->patBest.empty(),
    "UUUPS: No CSP solution found. ??? dK="
    <<pr->deltaK);
  pr->ExtractBestSolution(csp1);
  assertm(not pr->patBest.empty(),
    "Need the best CSP solution for start. dK="
    <<pr->deltaK);
  pr->ubD = (double)pr->xiBest.size();
  gub = //pr->raster_ceil // not yet initialized
    (pr->GetHeurBnd());
//  pr->SetHeurBnd(gub);
  BCP::TimeLimit__ = TL;
  opt::GlobalOutputLevel() = oL;
  dbg_outn_(1, " cspLP="<<pr->cspLP<<" cspIP="<<pr->cspIP
    <<" CSPnDiff="<<pr->ubD<<" objVal="<<gub);
}

void BCP2::SolveBPP() {
  dbg_outn(1, "    Solving BPP...");
  int i;
  // CSP1::fMergePieces works.
  CSP1 *bpp1 = new CSP1(*(dynamic_cast<CSP1*>(pr)));
    // creating a BPP1 problem;
  for (i=0;i<bpp1->m;++i) bpp1->pc[i].b = 1;

  BCP bcp(bpp1);
  bcp.firstinstance = 0; // here not: don't
  bcp.fSubproblem = 1;
    // rewrite logs
  double oL = opt::GlobalOutputLevel();
    opt::GlobalOutputLevel()=pr->CSPOutpLev;
  double TL = BCP::TimeLimit__;
    BCP::TimeLimit__=pr->BPPTiLim;
  bcp.Run();
  fBPPOpt = (BCP::opt == bcp.status);
  pr->lbD = bcp.glb; // GLB!!!
  BCP::TimeLimit__ = TL;
  opt::GlobalOutputLevel() = oL;
  dbg_outn(1, " ... minDiff="<<pr->lbD);
}

void AdjustOptions() {
/*
Check options: integer bounding (?) + select: most infeas
CSP:: rUseX - own option, =1
Model equality, then start basis: only diag.
or try greedy, if not feas => diag.
Artificial column instead of starting basis? (deltaK < inf)
Then not branch on it. Mark as slack
Adding col: consider K-_D+1
Where LPZF const terms?
  + change LP obj coef to =pV for bounded from below
LP: Ax = b .
*/
}

void BCP2::UpdateVarBounds() {
  BCP::UpdateVarBounds();
  int i;
  lpzfAdd = 0;
  for (i=0;i<lp->NCols();++i) // now the actual cols in f.
    if (not IsSlackCol(i)) {
      int j = GetColIndex(i);
      if (lb[j] > 0) {
        lpzfAdd += pr->pF;
        lp->SetObjCoef(i, pr->pV);
        GetCol(i)->SetObj(pr->pV); // For Phase I
      } else {
        int lub = LocalUpperBound(GetCol(i));
        int newUB = IMin(lub, ub[j]);
        if (newUB) { // For PhaseI
          lp->SetObjCoef(i, pr->pV + pr->pF/newUB);
          GetCol(i)->SetObj(pr->pV + pr->pF/newUB);
        }
      }
    }
  dbg_outn(4,"Init node: lpzfConst="<<lpzfAdd);
}

/// Returns 1 if a better solution:
bool BCP2::Rounding(int fFast) {
  // fFast=2: e.g. after dual, feas != opt
  // maybe separate: =2 only feas. test (prim not opt)
  // some counter: later in the CPA of a single node
  // => less intensive, more randomized
  int i;
  double z_old = gub;
  double obj = 0;
  for (i=0;i<lpx.size();++i) {
    if (lpx[i] - floor(lpx[i]+pr->GetXEps())
      > FMax(pr->GetXEps(), pr->GetXEps() * lpx[i])
      && not IsSlackCol(i))
        goto LPNonInt;
    if (lpx[i] > pr->GetXEps() and not cols[i].slackCut)
      if (GetCol(i)->fInfeasible())
      goto LPNonInt;
  }

  for (i=0;i<lpx.size();++i)
    if (not IsSlackCol(i) and lpx[i] > 0.5)
        obj += pr->pF + pr->pV * round(lpx[i]); ///////////////
  obj = pr->raster_ceil(obj); // against rnd errors

  if ((OUTP_LEV__ >=2 && !fFast) or OUTP_LEV__ >= 4) {
    log_ln(" ALREADY INTEGRAL: LP Bnd " << clb << " Orig ObjVal: "<<obj);
  }
  if (obj < pr->GetHeurBnd()) {
    timeIPBest = timer.userTime();
    ConvertLPCols(lpx,pr->patBest,pr->xiBest); // no 0-intens!!
    for (i=0;i<pr->xiBest.size();++i) {
      pr->xiBest[i] = round(pr->xiBest[i]);
      pr->patBest[i].SetObj(1);
    }
    pr->ExtractSolution(pr->patBest, pr->xiBest);
      // -- before changing bound
    pr->SetHeurBnd(gub = obj); // ALSO gub
    fMIPSolBest = false;
  }
  if (LocalOptimum()) {
    if (!fFast) dbg_outn_(3," lo_");
//    return true;
  }
//  return false;
    return (gub != z_old);

LPNonInt:
  if (fFast) return false; // for the first version !!!
  if (RNDNodeInterval < 1) goto MIP; 
     //RNDNodeInterval=10;
  if (RNDInterval < 1) goto MIP; //RNDInterval=1000;
  if (RNDColInterval < 1) goto MIP;
//  if (depth>4) // only
  if (colpool.size() >= nextRNDNCols)
    nextRNDNCols += RNDColInterval;
  else
  if (cntNode%RNDNodeInterval
    or iter%RNDInterval) goto MIP;
  ConvertLPCols(lpx,pr->pat,pr->lpx);
  pr->lpd = lpd;
  pr->SetFastRounding(fFast);
  z_old = pr->GetHeurBnd();
  dbg_outn_(2, " RND.");
  pr->ConstructIntSol();
  gub = pr->raster_ceil(pr->GetHeurBnd());
  pr->SetHeurBnd(gub); // ... assure gub \in \RP
  if (z_old > pr->GetHeurBnd()) {
    timeIPBest = timer.userTime();
    fMIPSolBest = false;
    dbg_outn_(2, '='<<gub);
  }
  pr->SetFastRounding(0);
MIP: // lpb rounding: -> numerics
  // PMP: another obj => no MIP
  if (gub < 1e50) status = feas;
  if (z_old >= 1e50 && gub < 1e50) {
    // -- THIS is true only once
  // -- to setup raster in 2D/1MCSP
    SetCurrentLPV(clv);
    GetLocalLagrBnd();
    SetLocalLPV(llv);
    glb_others = pr->raster_ceil(glv_others);
    UpdateGlobalLPV(); // FMin(glb_others,llb);
    // Because new raster was set up
    // OR SMTH LIKE IF (pr->Need...) .....
  }
  if (LocalOptimum()) {
    if (!fFast) dbg_outn_(3," lo_");
//    return true;
  }
//  return false;
    return (gub != z_old);
}

pair<int, double> BCP2::SelectBrVar_MostInfeas() { // + least frac
  double level1Min=1e100, infMin=1e100, t;
  int i, iLevel1Min = -1, iInfMin=-1;
  for (i=0;i<lpx.size();++i)
    if (LP::basic == cstat[i] and not IsSlackCol(i)) {
      double lbi, ubi;
      lp->GetVarBnds(i, lbi, ubi);
      ubi = FMin (ubi, LocalUpperBound(GetCol(i))); // !!! if no ub
      // why not on slacks?
      if (lbi < 1e-6 and lpx[i] > 1e-6
        and lpx[i] < ubi-1e-6) { // the first level
		  if (4 == nVarSel)
        t = 1.0 - lpx[i]/ubi; // second best LocalUpperBound(GetCol(i));
		  else if (3==nVarSel)
        t =  - lpx[i]; //LAST BEST LocalUpperBound(GetCol(i));
		  else //if (2==nVarSel)
	      t = fabs(lpx[i] - 1);
 //       t = fabs(frac(lpx[i]) - 0.5); // -0.5) // bloed
        if (t < level1Min)
        { level1Min = t; iLevel1Min = i; }
      }
    }
  dbg_outn(4,"BrVarSelect: iLevel1="<<iLevel1Min
    <<", infMin="<<level1Min);
  if (iLevel1Min != -1) { // BLOCK !!!!!!
    if (2 == nVarSel)
      return pair<int,double>(iLevel1Min, 0);
    return pair<int,double>(iLevel1Min, FMax(floor(lpx[iLevel1Min]),0));
  }

  // LEVEL 2:
  for (i=0;i<lpx.size();++i)
    if (LP::basic == cstat[i] and not IsSlackCol(i))
      // why not on slacks?
     { // the fractional level
      double lbi, ubi;
      lp->GetVarBnds(i, lbi, ubi);
      ubi = FMin (ubi, lbi + // !!!
        LocalUpperBound(GetCol(i))); // !!! if no ub
      if (fabs(lbi - ubi) > 1e-6)
        if (fabs(frac(lpx[i]) - 0.5)
          + wTillUB * (ubi - lpx[i])/(ubi-lbi) // very good
          < infMin) { // but better somehow measure the abs value of lpxi
          infMin = fabs(frac(lpx[i]) - 0.5) + wTillUB * (ubi - lpx[i])/(ubi-lbi);
          iInfMin=i;
        }
      }
  dbg_outn(4,"BrVarSelect: iLevel2="<<iInfMin
    <<", infMin="<<infMin);
  assertm(iInfMin >= 0, "No branching var ???");
  double newub = floor(lpx[iInfMin]);
      double lbi, ubi;
      lp->GetVarBnds(iInfMin, lbi, ubi);
      assertm(fabs(lbi - ubi) > 1e-6, " fixed var is basic and no better chosen ??? ");
  if (newub > ubi-1e-6) newub = ubi - 1; // >=
  if (newub < lbi-1e-6) newub = lbi; // <
  return pair<int,double>(iInfMin, newub);
}


void BCP2::DoneRun() {
  BCP::DoneRun();
      if (BCP::opt == status or BCP::feas==status) {
        ++nFeas;
        NCSP += accumulate(pr->xiBest.begin(), pr->xiBest.end(), 0.0);
        ND0 += pr->ubD;
        ND += (double)pr->xiBest.size();
      }
}

void BCP2::InitOptimize() {
  BCP::InitOptimize();
  lpzfAdd = 0;
  dbg_outn(4,"Starting optimize: gub="<<gub<<" glb="<<glb);
}

void BCP2::InitBounding() {
  BCP::InitBounding();
  pr->q0max = 0; // max bi
  int i;
  for (i=0;i<pr->Dim();++i)
    if (pr->q0max < pr->b_cg[i]
      and pr->b_cg[i] < double(INT_MAX))
      pr->q0max = pr->b_cg[i];
}

void BCP2::PrintIter()
{
//  if (not SmthChanged() and (iter>0))
//    return;
  if (OUTP_LEV__ > 1.99
    or (OUTP_LEV__ >0.49 && fOutputTimer)) {
    mylog 
      << "\n\"" << pr->infile
      << "\" #" << pr->inst
      << " U" << pr->GetHeurBnd()
      << " L" << glv//pr->GetLPBnd()
      << " LL" << llb
      << ' '
        << //double(int(1000.0*
        setprecision(4) <<
        Del0(double(pr->GetHeurBnd() - glv),fabs(glv)) * 100.0
          //))/10.0
            << '%'
            << setprecision(6)
//      << " (" << pr->GetHeurBnd() - glv << ')'
      <<" n" << colpool.size()
      << " m" << cuts.size()
      << " nd"<<theNode->no
      ;
    if (theNode->parent)
      mylog<<"<-"<<theNode->parent->no;//cntNode
    mylog
      << " Nn"<<nodes.GetSize() // this info only sometimes
      << " d"<<depth
      << " lpC" << lpzfAdd
//      << " It "<<iter
      ;
    if (iter)
      mylog << " i"<<iter;
    fOutputTimer = 0;
  }
} //____________________________________________________

void BCP2::PrintLog()
{
  if (OUTP_LEV__ <= 0.00001)
    return; // for batch usage as a tool
  int prec=7;
  if (status != error && status != infeas && status != LPErr)
    if (pr->Optimum()) status = opt;
  // time_t time0=time(0); // in the log
  // printing ...
  cout << setprecision(prec) << setw(3);
    strncpy(outfile1,pr->infile,sizeof(outfile1)-5);
    strcat(outfile1,".txt");
  ofstream ofs(outfile1,ios::out|ios::app); // Here app
  ofs
#ifdef _MSC_VER
      << left
#endif
      << setprecision(prec)
      << setw(3) << (pr->inst) << ' '
      << setw(5) << statusName[status] << ' '
      << setw(prec+1) << pr->GetHeurBnd()
      << (fMIPSolBest ? "* " : " ")
      << setw(prec+1) << pr->xiBest.size() << ' '
      << setw(prec+1)
        << accumulate(pr->xiBest.begin(),
          pr->xiBest.end(), 0.0) << ' '
      << setw(prec+1)
        << //pr->GetHeurBnd() - 
        glv << ' ' // best
      << setw(prec+1)
        << glv0 << ' '
      << setw(prec+1) << pr->cspLP << ' '
      << setw(prec+1) << pr->cspIP << ' '
      << setw(prec+1) << pr->lbD << ' '
      << setw(prec+1) << timeIP << ' '
      << setw(prec+1) << timeLP0 << ' '
      << setw(prec+1) << timeLPBest << ' '
      << setw(prec+1) << timeIPBest << ' '
//      << ctime(time0) <<
      << setprecision(3) << setw(5)
//        << nRuns // - 1 + FMin(tLastRun/tRunMax,0.9999)
//        << ' ' // float
      << setw(5) << cntInitNode << ' ' // only really slv
      << setw(5) << sum_iter << ' '
      << setw(5) << lpCols << ' '
      << setw(5) << ipCols << ' '
// Info : IRUP, ..., errors ? // see log
//      "lMin   lMax   L      "
      << pr->prName;
  if (glb0 < glb
    && (LPBnd==status || opt==status || feas==status)
      && Solver::problemType==1) // CSP1
    ofs << " nIRUP";
//  else
  if (iter) ofs << " I";
  if (timeIP > TimeLimit__) ofs << " T";
  if (nErr2 != nErr2__) ofs << " e2";
  pr->PrintLog(ofs);
  ofs << '\n';
  PrintSolutions(); // with corr. params
}

void BCP2::InitLog1()
{
  if (OUTP_LEV__ <= 0.00001)
    return; // for batch usage as a tool
  int prec = 7;// The compact output for each problem:
//  if (instFirst==inst) {
    strncpy(outfile1,pr->infile,sizeof(outfile1)-5);
    strcat(outfile1,".txt");
    time_t time0=time(0);
    ofstream ofs(outfile1,ios::out|ios::app); // Here app
    ofs 
#ifdef _MSC_VER
      << left
#endif
      << //setprecision(prec) <<
      "\n\n" << Version() << '\n'
      << ctime(&time0);
//    PrintOptions(ofs);
    opt::SolverCfg()->WriteOptionsShort(ofs);
    ofs << "\npF="<<pr->pF << " pV="<<pr->pV/*<<" gcd="<<pr->gcdPP*/<<'\n';
    ofs << "\n\n#  status "
      << setw(prec+2) << "IP_best"
      << setw(prec+2) << "nDiff"
      << setw(prec+2) << "sum_xj"
      << setw(prec+2) << "LPBnd"
      << setw(prec+2) << "LPBnd0 "
      << setw(prec+2) << "cspLP "
      << setw(prec+2) << "cspIP "
      << setw(prec+2) << "bppLP "
      << setw(prec+2) << "time"
      << setw(prec+2) << "timeLP0"
      << setw(prec+2) << "tLPbest"
      << setw(prec+2) << "tIPbest"
//      << setw(6) << " NRuns" // float

      << " NNode NIter LPCol IPCol "
//      "lMin   lMax   L      "
      "[remark]\n"
      "*) IPBest*: found with MIP optimizer"
      "\n\n";
//    ofs << "\n\n";
//  }
}

void BCP2::PrintStatistics(ostream &os) {
  BCP::PrintStatistics(os);
  if (!nFeas) return;
  os << " PMP: NCSPresult(<=cspIP+dK): "<<NCSP/nFeas
    << " ND0: "<<ND0/nFeas
    << " ND: "<<ND/nFeas
    << " redu%: "<<((ND0-ND)/ND0)*100;
}

double BCP2::NCSP, BCP2::ND0, BCP2::ND;
int BCP2::nFeas;

double BCP2::outputLevel;
int BCP2::nVarSel;

opt::OptSection BCP2::opt__
  ("BCP2", "Branch Cut Price for PMP1",
  BCP2::Options(), opt::SolverCfg(), 2500);
opt::OptContainer BCP2::Options() {
  opt::OptContainer oc;
  oc
    << opt::MakeOpt(&nVarSel, 3,
      "nVarSel", "1-4. 1,2: closest to 1 (2: ub=0), 3,4: largest")
    << opt::MakeOpt(&outputLevel, DEF_OUTP_LEVEL,
      "outputLevel", "0-5")
     ;
  return oc;
} //____________________________________________________


}                                                // SS_END_NAMESPACE__

