// FILE: bcpstuff.cpp, branch&cut&price stuff (i/o...)
// Author: Gleb <Belov@math.tu-dresden.de>

#include "stdafx.h"
#include "bcp.h"
#include "solver.h"
#include "bbcuts.h" // for options
#include "probl_cp22.h" // also
#include "lasthdr.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

void BCP::ReadGivenSolution() {
  if (fSubproblem or not fCheckGivenSol) return;
  ifstream is("solution.dat");
  if (!is) {
    perror("Could not read solution.dat: "); return;
  }
  char line[1024];
  int a,b;
  is >> givsolobj;
  log_n(2,"Reading help solution: obval = " << givsolobj);
  while (is) {
    is.getline(line, sizeof(line));
    Column col;
    istringstream iss(line);
    iss >> col.index; // as x
    iss.get(); // :
    while (iss) {
      iss >> a;
      iss.get();
      iss >> b;
      if (iss)
        col.PushID(a,b);
    }
    if (OUTP_LEV__ >= 2) {
      log__("col: ");
      int i;
      for (i=0;i<col.id.size();++i)
        log__(col.id[i].i<<':'<<col.id[i].d<<' ');
      log_ln("");
    }
    if (not col.id.empty())
      givsol.insert(col);
  }
}

void BCP::CheckGivenSolution() {
  // ADD: considering ALL columns of the givsol ___
  // Still incorrect
  if (fSubproblem or not fCheckGivenSol
    or givsol.empty()) return;
  int i;
//  bool feas=1;
  double lbi, ubi;
  Column *pc;
  Pool<Column>::iterator ic;
//  if (llv < givsolobj)
    log_n_(3,"Local LP val "<<llv<<", GIVEN sol. "<<givsolobj);
  // FIRST: Existent columns. What if hidden ?
  // Thus, checking only the current lp.
  for_each_in(givsol,ic,) ic->nHidden = 0;
  for (i=0;i<lp->NCols();++i) {
    pc = givsol.Find(*GetCol(i));
    if (pc) {
      pc->nHidden = 1;
      lp->GetVarBnds(i,lbi,ubi);
      if (pc->index < lbi-1e-6 or pc->index > ubi+1e-6) {
        log_n_(2," GivenINFEAS0");
        return;
      }
    }
  }
  // SECOND: the other columns
  for_each_in(givsol,ic,)
  if (not ic->nHidden) {
    int ub = LocalUpperBound(&*ic);
    if (ic->index > ub) {
        log_n_(2," GivenINFEAS1");
        return;
    }
  }
  log_n_(2," GivenFeas.");
  if (llv < givsolobj+pr->GetVEps())
  { log_n_(2," Obj OK."); }
  else
    assertm(0,"Given feas BUT Local WORSE ____ ");
  return;
}

volatile bool BCP::fOutputTimer=0;

void BCP::PrintIter()
{
//  if (not SmthChanged() and (iter>0))
//    return;
	/// Saving LP:
	if (fOutputTimer)
		lp->WriteModel("LPModelNow.lp");
		
  if (OUTP_LEV__ > 1.99
    or (OUTP_LEV__ >0.49 && fOutputTimer)) {
    mylog 
      << "\n\"" << pr->infile
      << "\" #" << pr->inst
      << " U" << pr->GetHeurBnd()
      <<setprecision(16)
      << " L" << glv//pr->GetLPBnd()
      <<setprecision(6)
//       << " LL" << llb
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
      << " Nn"<<nodes.GetSearchSize()
        <<'/'<<nodes.GetSize() // this info only sometimes
      << " d"<<depth
//      << " It "<<iter
      << " LR "
        << nLeft << '/' << nRight
//      << " Toff"
//        << ( nCG ? int (double(nCG_TOff)*100 / nCG) : 0 ) << '%'
      ;
    if (iter)
      mylog << " i"<<iter;
    fOutputTimer = 0;
 }
} //____________________________________________________

bool BCP::TimeLimit() {
	if ( timer.userTime() > TimeLimit__ ) {
		log_ln("  TIME LIMIT REACHED: " << TimeLimit__ );
		GetEnv().GetLog2() << "  TIME LIMIT REACHED: " << TimeLimit__ << '\n';
		return 1;
	}
	return 0;
}

void BCP::AlarmHandler(int sig) {
  fOutputTimer=1;
  signal (sig, AlarmHandler);
}

void BCP::PrintLog()
{
  if (OUTP_LEV__ <= 0.00001 or fSubproblem)
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
      << ' ' << setw(2) << (pr->inst) << ' '
      << setw(5) << statusName[status] << ' '
      << setw(prec+1) << pr->GetHeurBnd()
      << (fMIPSolBest ? "* " : " ")
			<< setprecision(15)
			<< setw(prec+1)
        << //pr->GetHeurBnd() - 
        glv << ' ' // best
      << setw(prec+1)
        << glv0 << ' '
       << setprecision(prec)
			<< setw(prec+1) << timeIP << ' '
      << setw(prec+1) << timeLP0 << ' '
      << setw(prec+1) << timeLPBest << ' '
      << setw(prec+1) << timeIPBest << ' '
//      << ctime(time0) <<
      << setprecision(3) << setw(5)
//        << nRuns // - 1 + FMin(tLastRun/tRunMax,0.9999)
//        << ' ' // float
      << setw(5) << cntNode << ' ' // only really slv
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
  if (sum_iter) ofs << " I";
  if (timeIP > TimeLimit__) ofs << " T";
  if (nErr2 != nErr2__) ofs << " e2";
  pr->PrintLog(ofs);
  ofs << '\n';
  PrintSolutions(); // with corr. params
}

void BCP::InitLog1()
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
    ofs << "\n\n # status "
      << setw(prec+2) << "IP_best"
      << setw(prec+2) << "LPBnd"
      << setw(prec+2) << "LPBnd0 "
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

void BCP::PrintContSol() {
  if (OUTP_LEV__ <= 0.00001)
    return; // for batch usage as a tool
  GetEnv().GetLog2()
    << "\nFile: "<<pr->infile
    <<"  Instance No. "<<pr->inst<<endl;
  GetEnv().GetLog2()
    << GetDateAndTime() <<endl;

      GetEnv().GetLog2()   
      << "LP solution: lpval="
      << setprecision(16)
      << clv << " lpbnd=" << glb <<"\n\n";
    int i;
    for (i=0;i<lpx.size();++i) {
      if (!lpx[i]) continue;
      GetEnv().GetLog2() << lpx[i] << ": ";
      if (not cols[i].slackCut) // ?????
        pr->PrintColumn(GetEnv().GetLog2(),GetCol(i));
      else // and slacks ???
        GetEnv().GetLog2() << " cut slack: " << cols[i].slackCut;
      GetEnv().GetLog2() << "\\\\\n";
    }
    GetEnv().GetLog2()  << setprecision(6) << "\n\n" << endl;
}

void BCP::PrintSolutions() {
  if (OUTP_LEV__ <= 0.00001)
    return; // for batch usage as a tool
  // fPrintContSol0, fPrintContSolBest,
  if (!(fPrintContSolBest || fPrintIntSolBest))
    goto PrUnsolved;

  GetEnv().GetLog2()
    << "\n\n\nFile: "<<pr->infile
    <<"  Instance No. "<<pr->inst<<endl;
  GetEnv().GetLog2()
    << GetDateAndTime() <<endl;

  if (fPrintIntSolBest) {
    GetEnv().GetLog2() 
      << "LP0 Time: " << timeLP0 << " IPTime: " << timeIP
      << " The best integer solution ("
      <<pr->GetHeurBnd()<<"):\n\n";
    int i;
    for (i=0;i<pr->xiBest.size();++i) {
      Column c;
      if (!pr->xiBest[i]) continue;
      GetEnv().GetLog2() << pr->xiBest[i] << ": ";
      pr->MakeColumn(&c,&pr->patBest[i]);
      pr->PrintColumn(GetEnv().GetLog2(),&c);
      GetEnv().GetLog2() << "\\\\\n";
    }
    pr->PrintBestSolution(GetEnv().GetLog2());
  }
  // + Smt fPrintCuts, fPrintProblem
  if (fPrintContSolBest) {
    GetEnv().GetLog2()
      << "\n\n\nThe last LP solution \n\n";
    PrintContSol();
  }

PrUnsolved:
  if (fPrintUnsolved
    and (opt != status or (sum_iter>0 or cntNode>0)))
    ; else
    return;
  char uns[1024];
    strncpy(uns,pr->infile,sizeof(uns)-13);
    strcat(uns,"_compl.txt");
    ofstream ofs(uns, ios::out|ios::app);
    pr->PrintProblem(ofs);
}

void BCP::InitLog2() { // None ?
  if (OUTP_LEV__ <= 0.00001)
    return; // for batch usage as a tool
//  if (!(fPr
    strncpy(outfile2,pr->infile,sizeof(outfile2)-7);
    strcat(outfile2,"__.txt");
  if (firstinstance) {
    GetEnv().GetLog2().open(outfile2);
    GetEnv().GetLog2().close();
  }
//  GetEnv().GetLog2().close();
  if (Solver::problemType!=7 or firstinstance) // because the same input file
    GetEnv().GetLog2().open(outfile2, ios::app);
  if (!GetEnv().GetLog2()) {
    cerr << "Opening Log2 for append: "<<outfile2 << endl;
    perror("Error: ");
  }

  if (not fPrintUnsolved)
    return;
  char uns[1024];
    strncpy(uns,pr->infile,sizeof(uns)-13);
    strcat(uns,"_compl.txt");
  if (firstinstance) {
    ofstream ofs(uns);
  }
}

void BCP::InitLog3() {
// The detailed output: MAYBE should be able to switch off
  if (OUTP_LEV__ >= 2) {
/*    strncpy(outfile3,"__verbose.txt",sizeof(outfile2)-7);
    GetMyLog__().close();
    GetMyLog__().open(outfile2);
    GetMyLog__().close(); // where closed before ?
    GetMyLog__().open(outfile2,ios::out|ios::app); // Here not app

  */}
}

BCP::BCP(Problem *p_) 
	: Alg(p_), pr(p_), lp(LP::CreateLP()), mip(LP::CreateLP()),
  subgr(new Subgr()), fSubproblem(0), nColsPrintLPVal(0),
  glb0(-1e100), glv0(-1e100)
{ }

void BCP::PrintOptions(ostream &ofs) { // obsolete
/*
  ofs
      <<" TimeLimit:"<<BCP::TimeLimit__
    <<" nRunsMax:"<<BCP::nRunsMax
      <<" maxDepth:"<<BCP::maxDepth
      <<" nBFS:"<<BCP::nBFS
      <<" nDFS:"<<BCP::nDFS
      <<" BrStratRL:"<<BCP::BrStratRL
      <<" nBranch:"<<BCP::nBranch
      <<" dInfeasFracPart:"<<BCP::dInfeasFracPart
      <<" Alfa1:"<<BCP::Alfa1
      <<" wTillUB:"<<BCP::wTillUB

      <<" fRedCostBnd:"<<BCP::fRedCostBnd
      <<" fLocalUB:"<<BCP::fLocalUB
      <<" fLocalReduce:"<<BCP::fLocalReduce
     << " MIPTimeLimit:"<<BCP::MIPTiLim
      <<" RNDNodeInterval:"<<BCP::RNDNodeInterval
      <<" RNDInterval:"<<BCP::RNDInterval
      <<" RNDColInterval:"<<BCP::RNDColInterval
      <<" MIPNodeInterval:"<<BCP::MIPNodeInterval
      <<" MIPInterval:"<<BCP::MIPInterval
      <<" fOriginalRedCost:"<<BCP::fOriginalRedCost
      <<" fUseLagrange:"<<BCP::fUseLagrange
    <<" _CUTS: iterRootMax:"<<BCP::iterRootMax
     << " iterMax:"<<BCP::iterMax
     << " kNodeNewCuts:"<<BCP::kNodeNewCuts
     << " maxCuts:"<<BCP::maxCuts__
     << " cutType:"<<BCP::cutType
     << " nPoolIterMax:"<<BCP::nPoolIterMax
     << " nCutsAddedMax:"<<BCP::nIterkMax
     << " nIterTailOff:" <<  BCP:: iterTailOff
     << " fSkipColGenAtLPBound:"<<BCP::fSkipColGenAtLPBound
     << " CGNormMax:"<<CGNormMax
     << " CGFracParts:"<<CGFracParts
     << " cntDel0:"<<cntDel0
     << " cntDel0inc:"<<cntDel0inc
     << " delReset0:"<< delReset0
     << " delResetInc:"<< delResetInc
     << " fCheckCuts:"<< fCheckCuts
    <<" BBCuts::nStepsTooMuch:"<<BBCuts::nStepsTooMuch
    <<" Subgr::iterMaxRatio:"<<Subgr::iterMaxRatio
    <<" Subgr::weight1:"<<Subgr::weight1
    <<" Subgr::rhoDef:"<<Subgr::rhoDef
    <<" CP22::fFirstCut1stD:"<<CP22::fFirstCut1stD
    ;*/
}

opt::OptContainer BCP::Options() {
  opt::OptContainer oc;
  oc
//    << opt::MakeOpt(&fLPOnly, DEF_LP_ONLY,
//      "fLPOnly", "bool: if LP bound only")
    << opt::MakeOpt(&TimeLimit__, 3600,
      "TimeLimit", "Overall time limit")
    << opt::MakeOpt(&fNoOpt, 0,
      "fNoOpt", "Optimality is not recognized -- to test cuts etc.")
     << opt::MakeOpt(&nRunsMax, N_RUNS_MAX_DEF,
      "nRunsMax", "N trial runs with diff params, "
      "if not fLPOnly")
    << opt::MakeOpt(&maxDepth, 999999,
      "maxDepth", "of the BCP tree")

    << opt::MakeOpt(&nBFSDFS, 1,
      "nBFSDFS", "Number of nodes selected according to the standard"
      " strategy, i.e. sorted by LB levels, then by depth")
    << opt::MakeOpt(&nDive, 2,
      "nDive", "m*nDive is the number of nodes selected according to DFS in the subtree"
      " of the last node, then again BFSDFS")
    << opt::MakeOpt(&BrStratRL, 0,
      "BrStratRL", "prefer nodes: 0: less left branches (xj<=), "
      "1: more left, "
      "between: probability of '=1'.  ")
    << opt::MakeOpt(&fCGTailOff, 0,
      "fCGTailOff", "Quit optimizing subnode if local LPval <= globalLPval")

    << opt::MakeOpt(&fBranchPsCosts, 0, // default: no
      "fBranchPsCosts", "Branch on Vars, var selection: 1: pseudo-costs, 0: most fractional")
    << opt::MakeOpt(&nBranch, 0, // default: no
      "nBranch", "1: branching on hyperplanes (defined in Problem), 0: branching on variables")
    << opt::MakeOpt(&BrVarFrac, 0.5,
      "BrVarFrac", "between 0 and 1: is the attracting fractional part for branching choice with BrOnVar; greater than 1: random choice")
    << opt::MakeOpt(&dInfeasFracPart, 0.001,
      "dInfeasFracPart", "variable infeasibility distance from int")
    << opt::MakeOpt(&Alfa1, 0,
      "Alfa1", "pseudo-cost (var) = Alfa1*Min(psDn,psUp) + Max(psDn,psUp)")
    << opt::MakeOpt(&wTillUB, 0.0,
      "wTillUB", "Selecting most infeas var: minimize "
      "fabs(frac(lpx[i]) - 0.5) + wTillUB * (ubi - lpx[i])/(ubi-lbi)"  )

    << opt::MakeOpt(&fRedCostBnd, false,
      "fRedCostBnd", "reduced cost bounding")

    << opt::MakeOpt(&fLocalReduce, true,
      "fLocalReduce", "Reduce RHS for col gen in nodes")
    << opt::MakeOpt(&fLocalUB, true,
      "fLocalUB", "Set local upper bounds for all vars. Set 1 for CP22 "
      /*"BUT only with identity start matrix for CSP1"*/)
    << opt::MakeOpt(&nReducePool, 0,
      "nReducePool", "!=0: In each node, master is initialized only with parent basis")
    << opt::MakeOpt(&fSkipColGenAtLPBound, false,
      "fSkipColGenAtLPBound",
      "Skip CG if LP Bnd == that in the last CPA iteration")
    << opt::MakeOpt(&fUseLagrange, false,
      "fUseLagrange",
      "Use Lagrange bound to stop col gen")
    << opt::MakeOpt(&RNDNodeInterval, 60,
      "RNDNodeInterval", "Every ...th [processed] node to apply rounding heuristics."
      "For CSP smth like 30, for PMP 10")
    << opt::MakeOpt(&RNDInterval, 60,
      "RNDInterval", "Every ...th CPA iteration to apply")
    << opt::MakeOpt(&RNDColInterval, 30,
      "RNDColInterval", "Apply when N cols increases by it. Def: 30 like Vanderbeck")
    << opt::MakeOpt(&MIPTiLim, 200,
      "MIPTimeLimit", "... in each run. SET THIS 0 FOR NO MIP ROUNDING")
    << opt::MakeOpt(&MIPNodeInterval, 500,
      "MIPNodeInterval", "Every ...th [processed] node to apply MIP Solver")
    << opt::MakeOpt(&MIPInterval, 40,
      "MIPInterval", "Every ...th CPA iteration to apply")
    << opt::MakeOpt(&fOriginalRedCost, true,
      "fOriginalRedCost",
      "Subgradient: take only new columns with negative _original_"
      " reduced cost")
    << opt::MakeOpt(&M__, 1e5,
      "M__", "initial obj. coef of infeasible slacks")
//    << opt::MakeOpt(&kNodeOutput, 100,
 //     "kNodeOutput", "each ... node output allowed")
    << opt::MakeOpt(&fPrintProblem, true,
      "fPrintProblem", "bool: print what we solve")
    << opt::MakeOpt(&fPrintContSol0, true,
      "fPrintContSol0", "bool: print initial LP solution")
    << opt::MakeOpt(&fPrintContSolBest, true,

      "fPrintContSolBest", "bool: print the best (better: last) LP solution")
    << opt::MakeOpt(&fPrintIntSolBest, true,
      "fPrintIntSolBest", "bool: print int solution")
    << opt::MakeOpt(&fPrintUnsolved, true,
      "fPrintComplicated", "bool: print problems unsolved after initial LP")
    << opt::MakeOpt(&fCheckGivenSol, 0,
      "fCheckGivenSol", "Solution from solution.dat will be checked for "
      "feasibility in every node. Format: objvalue {x: i1:n1 .. ik:nk[newline]} "
      "like mine output in ...__.txt BUT NO slacks. Not with cuts yet. "
      "THE PROC is currently wrong.")
    << opt::MakeOpt(&outputLevel, DEF_OUTP_LEVEL,
      "outputLevel", "0-5")
    << opt::MakeOpt(&outputInterval, 1.0,
      "outputInterval", "If outputLevel==1, you will see only the main info each .. sec.")
     ;
  return oc;
} //____________________________________________________
opt::OptSection BCP::opt__
  ("BCP", "Branch Cut Price",
  BCP::Options(), opt::SolverCfg(), 2500);
opt::OptContainer BCP::OptCuts() {
  opt::OptContainer oc;
  oc
    << opt::MakeOpt(&firstCutAfterNode, 0,
      "firstCutAfterNode", "")
    << opt::MakeOpt(&iterRootMax, 10,
      "iterRootMax", "Root node: IMax main CPA iterations (with col generation)")
    << opt::MakeOpt(&iterMax, 3,
      "iterMax", "Subnodes: IMax main CPA iterations (with col generation)")
    << opt::MakeOpt(&kNodeNewCuts, 1,
      "kNodeNewCuts", "each k-th node new cuts are constructed")
    << opt::MakeOpt(&maxCuts__, 1,
      "maxCuts", "IMax N cuts present. SET =0 IF DON'T NEED CUTS. Most effective: 0 or 1")
    << opt::MakeOpt(&cutType, 1,
      "cutType", "0: Mixed-Integer Gomory cuts by superadditive funcs\n"
      "'1: Lifted CHVATAL-GOMORY cuts")
//    << opt::MakeOpt(&fSE, 0,  "fSE", "if slack elimination")
    << opt::MakeOpt(&nPoolIterMax, 30,
      "nPoolIterMax", "N CPA iterations w/o adding columns")
//    << opt::MakeOpt(&fPoolSearch, 1,
 //     "fPoolSearch", "Bool: try to find violated cuts in pool")
    << opt::MakeOpt(&nIterkMax, 10,
      "nCutsAddedMax",
      "N cuts added each time in each CPA subiter."
      " Set as large as poss but so that col gen is not too long")
    << opt::MakeOpt(&iterTailOff, 15,
      "nIterTailOff",
      "N CPA iter w/o bound change => 'tail off', break")
//    << opt::MakeOpt(&cutSelect, 1,
    //  "cutSelect", "Choosing cuts acc. to violation: 0-near 0.5, -1-min, 1-max, 2-all")
    << opt::MakeOpt(&rndViolDev, 0.001,
      "rndViolDev",
      "adding cuts according to randomized violation grade,\n"
      "'sort.wgt = viol * (1+(rnd-0.5)*rndViolDev)")
//    << opt::MakeOpt(&nMarkedMax, 100,
    //  "nMarkedMax", "IMax N cuts protected from deletion (more => reset) NOT  IMPL")
    << opt::MakeOpt(&CGNormMax, 0,
      "CGNormMax", "Strengthening CG cuts: -1: mult those u with fr(ub)<1/2 by -1; "
      "0: mult by int(1/fr(ub); k>0: mult by int n: fr(nub) = max, n<=k")
    << opt::MakeOpt(&CGFracParts, 1,
      "CGFracParts", "Bool: take only frac parts of u's keeping the sign")
    << opt::MakeOpt(&cntDel0, 1,
      "cntDel0",
      "initial n of iterations till inactive cut del.")
    << opt::MakeOpt(&cntDel0inc, 1.2,
      "cntDel0inc",
      "coef. to mult. cntDel0 after cut deletion")
    << opt::MakeOpt(&delReset0, 30,
      "delReset0",
      "After how many iterations all inactive cuts are del"
      " in spite of the non-del counter (initially). Set as large as"
      " possible but so that col gen is not too long")
    << opt::MakeOpt(&delResetInc, 1.05,
      "delResetInc",
      "increase ratio for delReset")
//    << opt::MakeOpt(&fAddLevelCut, false,
    //  "fAddLevelCut", "BETTER NOT FOR CP22")
    << opt::MakeOpt(&fCheckCuts, false,
      "fCheckCuts", "bool: for debugging")
  ;
  return oc;
} //____________________________________________________
opt::OptSection BCP::optCuts
  ("BCP_Cuts", "The Cutting Plane part of BCP",
  BCP::OptCuts(), opt::SolverCfg(), 2500);

bool BCP::fLocalReduce=0;
bool BCP::fLocalUB=0;
int BCP::nReducePool=0;
int BCP::kNodeNewCuts=8;
int BCP::maxDepth=32000;
bool BCP::fNoOpt=0;

int BCP::nBFSDFS=10;
double BCP::nDive=2;
double BCP::BrStratRL=0;
bool BCP::fCGTailOff=0;

int BCP::nBranch;
float BCP::BrVarFrac=0.5;
bool BCP::fBranchPsCosts;
double BCP::dInfeasFracPart;
double BCP::Alfa1=2.0;
double BCP::wTillUB=1;

bool BCP::fRedCostBnd;
bool BCP::fCheckGivenSol=0;

double BCP::outputLevel;
double BCP::outputInterval;
int BCP::kNodeOutput=100;
double BCP::TimeLimit__ = 900;
int BCP::nRunsMax=N_RUNS_MAX_DEF;

int BCP::firstCutAfterNode=200;
int BCP::iterMax=0;
int BCP::iterRootMax=35;
int BCP::maxCuts__=20;
int BCP::iterTailOff=15;
int BCP::cutType = 0;
int BCP::fSE = 0;
int BCP::nPoolIterMax = 5;
int BCP::fPoolSearch = 1;
int BCP::nIterkMax = 7;
int BCP::cutSelect = 1;
//int BCP::nMarkedMax = 10;
int BCP::CGNormMax = -1;
int BCP::CGFracParts = 1;
double BCP::rndViolDev = 1;
int BCP::cntDel0 = 5;
double BCP::cntDel0inc = 1.5;
double BCP::delReset0 = 50; // all ... iter all inactive wb dl
double BCP::delResetInc = 1.04; // all ... iter all inactive wb dl
bool BCP::fLPOnly=DEF_LP_ONLY;
bool BCP::fSkipColGenAtLPBound=true;

bool BCP::fOriginalRedCost=true;
bool BCP::fUseLagrange=false;
double BCP::M__ = 1e8; // mult by 1.5
double BCP::MIPTiLim=30; // - time limit?
int BCP::MIPInterval=40; // how often to apply
int BCP::MIPNodeInterval=20; // how often to apply
int BCP::RNDInterval=40; // how often to apply
int BCP::RNDNodeInterval=20; // how often to apply
int BCP::RNDColInterval=10; // how often to apply
bool BCP::fAddLevelCut = 0;
bool BCP::fPrintContSol0=true,

  BCP::fPrintContSolBest=true, BCP::fPrintIntSolBest=true,
  BCP::fPrintProblem=true, BCP::fPrintUnsolved=true;
bool BCP::fCheckCuts=false;
const char * BCP::statusName[] =
  { "none", "OPT", "LPBnd", "feas", "LPErr", "infeas", "error"};

}                                                // SS_END_NAMESPACE__
