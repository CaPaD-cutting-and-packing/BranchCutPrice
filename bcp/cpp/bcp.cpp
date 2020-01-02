// FILE: bcp.cpp, branch&cut&price for CSP and 2CP
// Author: Gleb <Belov@math.tu-dresden.de>

#include "stdafx.h"
#include "bcp.h"
#include "solver.h"
#include "lasthdr.h"

#include "probl_csp1.h"             // to change params
#include "../inc/bcp.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

// + see other authors
// + some printing
// cut pool s. contain all
// where branching constr?
// RANDOM: in node choice, var choice
// ABACUS: also separation if pricing too long but
//   NO specific separation cycle w/o pricing
// BUT FOR US: If col gen too long,
// this means already that too many cuts.
// pruning/combinatorial bnds
// RERUN with cuts: root LP solution is always known
// Tailing Off on CG: no complete CG?
// Nodes' ordering! (all by bfs/dfs, diving specially,
//  to fathom specially)

// to get feasible need Phase 1 ?
// compare bb & spprc for cg (no integer bounding)

char * BCP::Version() {
  return (char*)
    "BCP: Version 2.0.  Build "__DATE__"  "__TIME__ 
    ;
} //____________________________________________________

void BCP::Run() {
try {
  lp->Open(); // Closing in the destructor
	mip->Open();
  InitRun(); // Randomization
  while (NextOptimize()) {
		int nTry=0;
Retry:
		++nTry;
    try {
      Optimize();
    } catch (const exception & e) {
      PRINT_ERROR(e.what());
      status = error;
    }
    catch (const double dErr) {
			assert(8.1 == dErr);
			assert(1==nTry);
			CSP1::fCompareBySPPRC_DP=1;
			fLocalUB=0;
			nBranch=1;
			CSP1::nHP=1;
			CSP1::fModelEquality=1;              // ?????????
			fLocalReduce=0;
			Subgr::iterMaxRatio=0;
			maxCuts = maxCuts__ = 0;
			lp->Close();
			lp->Open();
			mip->Close();
			mip->Open();
			allcols.clear();
			colpool.clear();
			cols.clear();
			patMIP.clear();
			cutpool.clear();
			cuts.clear();
			nodes.clear();
			backtrack.clear();
			pr->Init();
			log_ln(" RESTARTING INSTANCE... ");
			goto Retry;
		}
  }
} catch (const exception & e) {
  PRINT_ERROR(e.what());
  PRINT_LOG("Error while init BCP.");
      status = error;
}
try {
  DoneRun();
} catch (const exception & e) {
  PRINT_ERROR(e.what());
  PRINT_LOG("Error while closing BCP.");
}
} //____________________________________________________

void BCP::Optimize() {
  InitOptimize();
  do {
    if (not Bounding())  // if llb < gub and nodeLP feas
      if (not Stop()) Branch();
  } while (not Stop() && Select());
  DoneOptimize();
}

bool BCP::Stop()
{ return TimeLimit() or Optimum() or maxDepth<1
    or fLPOnly or not pr->BranchingPossible(); }

// bound comparison: no eps 'cause pr->raster_ceil
bool BCP::Bounding() { // true if the node to fathom
  bool ret=false;
  InitBounding();
  for (;;) {
    if (!SolveLP()) goto _Contradict;
    if (fLPIntOpt) goto _Exit;
    if (TimeLimit()) { ret = 0; goto _Exit; }
    if (RoundingMIP())  { // a better UB obtained
			StartReoptTree();
			if (theNode->fLPOpt && llb >= gub && !fNoOpt)
				goto _LocalOpt;
			if (fReoptTree > 0) break;
		}
    if (llb >= gub && !fNoOpt)  goto _Fathom;
    if (Rounding()) { // a better UB obtained
      StartReoptTree();
      if (theNode->fLPOpt && llb >= gub && !fNoOpt)
          goto _LocalOpt;
      if (fReoptTree > 0) break;
    }
    if (!FixAndSet()) goto _Contradict;
    if (NewValues()) {
      if (fLPIntOpt) goto _Exit;
      continue;
    }
    if (CPATailingOff() or StopCPA()) break; // not ABACUS
    Separate(); // no pricing
    if (fLPIntOpt) goto _Exit;
    PrintIter();
  }
//  theNode->llb = llb;
  theNode->MarkProcessed();
  dbg_outn_(2," _P");
  ret = false; goto _Exit;
_Fathom: dbg_outn_(3," f_");
_LocalOpt: Fix();
_Contradict:
  Fathom(theNode);
  dbg_outn_(2," _F");
  ret = true;
_Exit: DoneBounding();
  PrintIter();
  return ret;
}

/// LP lower bound, returns 0 if infeasible or error
bool BCP::SolveLP() {
    InitLP(); // meaning of LP::error ???
    theNode->fLPOpt = true;
    do {  // infeas: can still price
      SolvePrimal(); // + some fast rounding
      if (fLPIntOpt)
        return 1; // "all ok"
      if (LP::error==lp->status) { // +unbounded? not infeas
        if (!depth/* && !iter*/) status = LPErr;
	  // not in root
        log_n_(1," LPerr?_"); return 0;
      }
      if (LP::infeas == lp->status) {
        // sure this cannot be corrected by CG (inf. cols)
        dbg_outn_(3," LPinf!_");
	return 0;
      }
      if (fCGTailOff && clv <= glv // clv <= gub SENSE?
        && depth) { // Tailing off, only in subnodes
        theNode->fLPOpt = false;
	dbg_outn_(2," TO_");
	break;
      }
      if (TimeLimit()) return 1; // "everything o.k."
    } while (Price()/* && not StopLP()*/);
		if (LP::infeas == lp->status) // still after pricing
    { dbg_outn_(3," i_"); return 0;}
    PricingOver(); // everything remains ???????
    return 1;
}

bool BCP::CPATailingOff()
{ return iter - iterLLBLast > iterTailOff; }

bool BCP::StopCPA()
{ return not maxCuts or TimeLimit() // or LocalOpt() // ?
    or (depth ? iter >= iterMax : iter >= iterRootMax)
    or LP::error==lp->status or fLPOnly
      or not pr->SACutsPossible()
      or (cntNode < firstCutAfterNode)
      ;
}

// Attention: add it to the lb to get actual ub!
int BCP::LocalUpperBound(const Column *c) {
  int ub = INT_MAX, ub1;
  double ub1t;
  for_each_in(c->id,iid,Column::const_iterator) {
    assert(iid->d);
    ub1t = pr->b_cg[iid->i]/iid->d;
    if (ub1t > double(INT_MAX)) ub1 = INT_MAX;
    else ub1 = int(ub1t);
    if (ub1<0) if (pr->validsign[iid->i] > 0) ub1=0; //
    if (ub1<ub) ub=ub1;
  }
  return ub;
}

/// returns 1 if a better solution found:
bool BCP::Rounding(int fFast) {
  // fFast=2: e.g. after dual, feas != opt
  // maybe separate: =2 only feas. test (prim not opt)
  // some counter: later in the CPA of a single node
  // => less intensive, more randomized
//  if (fFast) return false; // for testing pure enum.
  tmRnd.start();
  int i;
  double z_old = -1e200;
  for (i=0;i<lpx.size();++i) {
    if (lpx[i] - floor(lpx[i]+pr->GetXEps())
      > FMax(pr->GetXEps(), pr->GetXEps() * lpx[i])
      && not IsSlackCol(i))
        goto LPNonInt;
    if (lpx[i] > pr->GetXEps() and not cols[i].slackCut)
      if (GetCol(i)->fInfeasible())
      goto LPNonInt; // but infeas? needed during fast computation
    if (lpx[i] > pr->GetXEps() and -2 == cols[i].j)
      goto LPNonInt;
  }
  if (clb < gub) // only when a better solution
  if ((!depth && OUTP_LEV__ >= 1.3)
    or (OUTP_LEV__ >=2 && !fFast) or OUTP_LEV__ >= 4) {
    log_ln(" LP SOLUTION ALREADY INTEGRAL: " << clb);
  }
  if (!fFast) dbg_outn_(3," lo_");
  if ((clb) < pr->GetHeurBnd()) {
    timeIPBest = timer.userTime();
    ConvertLPCols(lpx,pr->patBest,pr->xiBest);
    for (i=0;i<pr->xiBest.size();++i)
      pr->xiBest[i] = round(pr->xiBest[i]);
    pr->SetHeurBnd(gub = (clb)); // ALSO gub
    fMIPSolBest = false;
    if (clb == glb)
      fLPIntOpt = true;
  }
//  return true;
  tmRnd.stop();
  // assert(z_old != -1e200); // ???
  return (gub != z_old);

LPNonInt:
  if (fFast) {tmRnd.stop(); return false;} // for the first version !!!
  if (RNDNodeInterval < 1) goto MIP;
     //RNDNodeInterval=10;
  if (RNDInterval < 1) goto MIP; //RNDInterval=1000;
//  if (depth>4) // only
  if (RNDColInterval < 1) goto MIP;
//  if (depth>4) // only
  if (colpool.size() >= nextRNDNCols)
		nextRNDNCols = (int)colpool.size() + RNDColInterval;
  else
  if (cntNode%RNDNodeInterval
    or iter%RNDInterval) goto MIP;
  ConvertLPCols(lpx,pr->pat,pr->lpx);
  pr->lpd = lpd;
  pr->SetFastRounding(fFast);
  z_old = pr->GetHeurBnd();
  log_n_(2, " RND.");
  pr->ConstructIntSol();
  gub = pr->raster_ceil(pr->GetHeurBnd());
  pr->SetHeurBnd(gub); // ... assure gub \in \RP
  if (z_old > pr->GetHeurBnd()) {
    timeIPBest = timer.userTime();
    fMIPSolBest = false;
    log_n_(2, '='<<gub);
  }
  pr->SetFastRounding(0);
MIP: // lpb rounding: -> numerics
  /// DOING NOTHING, EXTRA FUNC
// FINISH:
  if (gub < 1e50) status = feas;
  if (z_old >= 1e50 && gub < 1e50) {
    // -- THIS is true only once
  // -- to setup raster in 2D/1MCSP
    SetCurrentLPV(clv);
    GetLocalLagrBnd();
    SetLocalLPV(llv);
    glb_others = pr->raster_ceil(glv_others);
    UpdateGlobalLPV(); // FMin(glb_others,llb);
    // because new raster was set up
    // OR SMTH LIKE IF (pr->Need...) .....
  }
  if (LocalOptimum()) {
    if (!fFast) log_n_(3," lo_");
//    return true;
  }
//  return false;
    tmRnd.stop();
  return (gub != z_old);

}

bool bgn::BCP::RoundingMIP()
{
	tmRnd.start();
	double z_old = pr->GetHeurBnd();
	
	if (MIPInterval && MIPNodeInterval)
		if (     // !LocalOptimum()  and         // local !! ?? TODO
			MIPTiLim > 0
			and (!(iter%MIPInterval))
			and (!(cntNode%MIPNodeInterval))
			and (cntNode)                           // to skip MIP at root ????????? Better for hard28? TODO
// 			and !fFast
			and (nMIPColsLast + 30 <= patMIP.size())           // TODO
		) {
			/// Deleting cuts:
			// 		CutList c00 = cuts;
			// 		i_vec l00 = lb, u00 = ub;
			// 		for (int i=cuts.size()-1; i>=0; --i) {
			// 			DelCut(i, true);
			// 		}
			// 		/// Updating var bounds:
			// 		cout << " lb.size() = " << lb.size() << endl;
			// 		lb.clear();
			// 		lb.resize(lp->NCols(), 0);
			// 		ub.clear();
			// 		ub.resize(lp->NCols(), INT_MAX);
			// 		UpdateVarBounds();
			// 		SolutionTree::node_iterator ndLast = theNode;
			// 		theNode = nodes.GetRoot();
			// 		InitNewNode();
			nMIPColsLast = (int)patMIP.size();                        // to prevetn solving the same
			Vector<double> xMIP(mip->NCols());
			assert(xMIP.size() == patMIP.size());
			for (int i=0;i<mip->NCols();++i)
				xMIP[i] = (true);
			double z_own = pr->GetHeurBnd();
			log_n_(2, " MIP.");
			double z_mip = mip->SolveMIP(xMIP, MIPTiLim, z_own-1);  // ????? TODO
			// Use own solution for start ? Would help cutting
			log_n_(2, "="<<z_mip);
			//    z_mip = pr->raster_ceil(z_mip); // ?????????
			z_mip = pr->raster_ceil(z_mip);
			pr->UpdateHeurBnd(z_mip);
			//      cout << setprecision(20) << z_mip << ' ' << z_own << endl;
			if (z_mip < z_own) {
				timeIPBest = timer.userTime();
				fMIPSolBest = true;
				StoreMIPPatterns(xMIP,pr->patBest,pr->xiBest);
				gub = z_mip; // rounding errors!!!
			}
			/// Restoring:
			//     theNode = ndLast;
			// 		InitNewNode();
			// 		lb = l00;
			// 		ub = u00;
			// 		UpdateVarBounds();
			// 		for (int i=0; i<c00.size(); ++i) {
			// 			AddCut(c00[i]);
			// 		}
		}
		if (gub < 1e50) status = feas;
		if (z_old >= 1e50 && gub < 1e50) {
			// -- THIS is true only once
			// -- to setup raster in 2D/1MCSP
			SetCurrentLPV(clv);
			GetLocalLagrBnd();
			SetLocalLPV(llv);
			glb_others = pr->raster_ceil(glv_others);
			UpdateGlobalLPV(); // FMin(glb_others,llb);
			// because new raster was set up
			// OR SMTH LIKE IF (pr->Need...) .....
		}
// 		if (LocalOptimum()) {
// 			if (!fFast) log_n_(3," lo_");
// 			//    return true;
// 		}
		//  return false;
		tmRnd.stop();
		return (gub != z_old);
}


void BCP::ConvertLPCols(Vector<double> &lpx,
  Vector<Pattern> &pat,Vector<double> &patx) {
// Converting non-0 cols to patterns:
  pat.clear(); patx.clear(); patx.reserve(Dim());
  int i; // does reserving prevent alloc. problems ???
  for (i=0;i<lpx.size();++i)
    if (fabs(lpx[i]) > pr->GetXEps() && not IsSlackCol(i))
      patx.push_back(
      //round
      (lpx[i]));
  pat.resize(patx.size());
  int j=0;
  for (i=0;i<lpx.size();++i)
    if (fabs(lpx[i]) > pr->GetXEps() && not IsSlackCol(i)) {
      pr->MakePattern(&pat[j++], GetCol(i));
      if (OUTP_LEV__ >=6) {
        mylog<<"Taking col: ";
        pr->PrintColumn(mylog,GetCol(i));
        mylog<<endl;
      }
    }
  assert(pat.size() == j);
}


void bgn::BCP::StoreMIPPatterns(Vector< double >& xMIP, Vector< Pattern >& pat, Vector< double >& patx)
{
	assert(xMIP.size() == patMIP.size());
// Converting non-0 cols to patterns:
  pat.clear(); patx.clear(); patx.reserve(Dim());
  int i; // does reserving prevent alloc. problems ???
  for (i=0;i<xMIP.size();++i)
		if (fabs(xMIP[i]) > pr->GetXEps())
      patx.push_back(
      //round
				(xMIP[i]));
	for (i=0;i<xMIP.size();++i)
		if (fabs(xMIP[i]) > pr->GetXEps()) {
			pat.push_back(patMIP[i]);
    }
}

bool BCP::FixAndSet() {
  Nnv = 0;
  if (!fRedCostBnd
    or fSkipColGenAtLPBound or fUseLagrange)
    return true; // nothing happened
  if (LocalOptimum())
    return true;
  lp->GetRedCosts(redCosts);
  int i;
  double gub_low = pr->raster_below(gub);
  for (i=0;i<lpx.size();++i)
  if (not IsSlackCol(i)
    and fabs(redCosts[i])>pr->GetRCEps()) {
    double ddx=(int)floor((gub_low - llb)/fabs(redCosts[i]));
    if (ddx> INT_MAX)
      continue;
    assert(ddx >= 0);
    int dx = (int)ddx;
    double lb,ub;
    VarBnd * vb;
    lp->GetVarBnds(i,lb,ub);
    if (LP::atLower==cstat[i]) {
      NnvTst ++;
      if (lb+dx < ub) {
        Nnv ++;
        dbg_outn_(5," x["<<i<<":("<<lb<<','<<ub<<")->("
          <<lb<<','<<lb+dx<<')');
        int ub1=int(lb)+dx;
        vb = theNode->rcBnds.Find
          (VarBnd(GetColIndex(i),1,int(lb)+dx));
        if (vb) ub1=vb->bnd = IMin(vb->bnd,int(lb)+dx);
        else
          theNode->rcBnds.Add
          (VarBnd(GetColIndex(i),1,int(lb)+dx));
        lp->ChangeVarBnds(i,lb,ub1);
      }
    } else if (LP::atUpper==cstat[i]) {
      NnvTst ++;
      if (ub-dx > lb) {
        Nnv ++;
        dbg_outn_(5," x["<<i<<":("<<lb<<','<<ub<<")->("
          <<ub-dx<<','<<ub<<')');
        int lb1=int(ub)-dx;
        vb = theNode->rcBnds.Find
          (VarBnd(GetColIndex(i),0,int(ub)-dx));
        if (vb) lb1=vb->bnd = IMax(vb->bnd,int(ub)-dx);
        else
          theNode->rcBnds.Add
          (VarBnd(GetColIndex(i),0,int(ub)-dx));
        lp->ChangeVarBnds(i,lb1,ub);
      }
    }
  }

  NnvAll += Nnv;
  return true;
  dbg_outn_(3," fs_");
  return false;
}
bool BCP::NewValues() {
  if (0==Nnv) return false;
  dbg_outn_(3," nv_"); // Frequency ? Print at end
  SolveDual();
  return true;
}
void BCP::Fix() { // if glb or gub changed
}

void BCP::Separate() { // we do this cycle always
  if (maxCuts < 1) return; // even if no new cuts in node
  ++ iter; // number of iterations with col.gen.
  int iterPool=0;
  int nPoolIter = int(rndCuts * nPoolIterMax );
  do { // at least once
    MarkTightCuts(); ////////////////////////
    DeleteSomeCuts(0);
    if (!SolvePrimal()) // why?
      assertm(false,"Could not solve primal after del cuts!!!");
    if (fLPIntOpt)
      return;
    if (cuts.size() >= MaxAllCuts()) return;

    if (not FindViolCutsInPool()) {
      dbg_outn_(2," no viol"); 
      if (depth && !kNodeNewCuts)
        return;
      if (kNodeNewCuts)
      if ((cntNode % kNodeNewCuts))
        return;
      ConstructNewCuts(); // assuming violated cuts can
    }  // always be constructed
    AddSomeCuts();
    PrintNCuts();

    if (!SolveDual()) { // like in InitNewNode()
      if (!pr->CanBeInfeasOnRestrPool()) // ?????
        return; // + SolvePrimal ? Globally infeas.
      else {
        if (!StartPhaseI())
          return;
      }
    }
    if (fLPIntOpt)
      return;
    ++ iterAll;
    dbg_outn_(2," SubIter="<<iterAll);
    RecalcCutsRHS();
  } while (not StopCPA() and (++iterPool<nPoolIter));
  PrintCuts();
}

bool BCP::Optimum() {
  if (fNoOpt)
    return false;
// 	if (fabs(132.0000140922412-llv)<1e-11) {
// 		cout << " llv="<<setprecision(16)<<llv;
// 	}
  if (glb <= gub + pr->GetVEps()) //// !!!!!
    return fabs(glb-gub) < pr->GetVEps();
//  outputLevel=5;
//  opt::GlobalOutputLevel() = 5;
  mylog << setprecision(16) << "glb="<<glb<<" > gub="<<gub<<"!! llv="<<llv<<"\n";
//  PrintCuts();
  int i, j;
  Vector<int> iNZ; // indices of non-0 vars
  Vector<double> xNZ; // the values
/*  iNZ.push_back(19); xNZ.push_back(3.0);
  iNZ.push_back(24); xNZ.push_back(1.0);
  iNZ.push_back(44); xNZ.push_back(2.0);
  iNZ.push_back(49); xNZ.push_back(1.0);
  iNZ.push_back(50); xNZ.push_back(1.0);
  CutList::iterator ic; // iterator(Cut*) ?
  CutRefSet::iterator ipc;
  cutViol.clear();
  for_each_in(cutpool,ipc,) (*ipc)->fPresent = false;
  for_each_in(cuts,ic,) (*ic)->fPresent = true;

  map<LPCut*,double> slackVals;*/
  // Filling with known slacks:
/*  for (i=0;i<cols.size();++i)
    if (cols[i].slackCut)
      slackVals.insert(
        make_pair(cols[i].slackCut,lpx[i]));
    else if (lpx[i]) {
      iNZ.push_back(GetColIndex(i)); // orig. index!
      xNZ.push_back(lpx[i]);
    }*/
/*  double sv;
  for (i=0;i<cuts.size();++i) {
      cout << "\n\nCut "<<cuts[i]<<": ";
      cuts[i]->Print(cout);
      ClearCoefs();
      sv = cuts[i]->CalcSlackValue(iNZ,xNZ,slackVals);
//      log_ln("Cut "<<cuts[i]<<": slack="<<sv);
  }*/
  // Checking the cols of the best solution:
  log_ln("Checking the cols of the best solution:");
  double rc;
  for (j=0;j<pr->patBest.size();++j) {
    Column col;
    pr->MakeColumn(&col,&pr->patBest[j]);
    log__("Col "<<j<<". x="<<pr->xiBest[j]<<": ");
    pr->PrintColumn(mylog,&col);
    if (allcols.Find(col)) log__(" IS IN. ") else log__(" IS NOT IN. ");
    rc = pr->CalcRedCost(&col,&cuts,lpd); // what if bounded down?
    log_ln("rc="<<rc);
  }
  log_ln('\n');
  SolvePrimal();
  log_ln("Resolved LP value = "<<setprecision(16)<<GetLPValue());
  // CHECKING CUT COEFS:
  vector<d_vec> cf(cuts.size()); //(colpool.size());
  for (i=0;i<cuts.size();++i) {
    cf[i].resize(colpool.size());
    for (j=0;j<colpool.size();++j) {
      cf[i][j] = cuts[i]->GetCoef(j);
    }
    ((SACutSE*)(&cuts[i]))->rawCoef.clear();
  }
  for (i=0;i<cuts.size();++i) {
    for (j=0;j<colpool.size();++j) {
      if (fabs(cf[i][j] - cuts[i]->Calc__(GetMainCol(j),j)) >1e-6)
        log__("!!");
    }
    ClearCoefs();
    double rh = cuts[i]->GetRHS();
    if (fabs(cuts[i]->CalcRHS__(pr->b) - rh) > 1e-6)
        log__("!!");
  }
  log_ln("Checking that all invloved cuts are in...");
  int nn=0;
  for (i=0;i<cuts.size();++i) cuts[i]->no = -1;
  for (i=0;i<cuts.size();++i) cuts[i]->Number(nn);
  assert(cuts.size() == nn);

  CheckCuts();

  log_ln("Writing LP to glbGgub.lp.");
  lp->WriteModel("glbGgub.lp");

  assert(glb <= gub); return glb == gub;
}

bool BCP::LocalOptimum()
{ /*assert(llb <= gub); No !!*/
  if (llb >= gub) return 1;
  assert( gub - llb > pr->GetVEps() );
  return 0;
}

void BCP::InitRun()
{
  assert(pr->Dim());
  //assert(delReset); // no div-0
  nRuns=0;
  gub = +1e100;
  llb = glb = llv = glv = llbLast=llvLast = -1e100;
  fLevelCut = !fAddLevelCut;
  // PARAMS FOR THIS INSTANCE WHICH CAN BE MODIFIED:
  maxCuts = maxCuts__;

  if (firstinstance) {
    InitLog1();
    pr->InitStat();
    if (nBranch) // branching on hyperplanes
      if (pr->NoForbiddenColsWithHyperplanes())
        if (fLocalUB or fLocalReduce) {
	  log_n(1, "Warning: INTEGER BOUNDING _COMPLETELY_ turned off because branching on hyperplanes");
	  fLocalUB = fLocalReduce = false;
	}
    if (nBranch) // branching on hyperplanes
      if (pr->NoCutsWithHyperplanes())
        if (maxCuts>0) {
	  log_n(1, "Warning: NO CUTTING PLANES (yet) when branching on hyperplanes");
	  maxCuts = 0;
	}
    if (nBranch) // branching on hyperplanes
      if (fBranchPsCosts)
      {
	  log_n(1, "Warning: NO PsCosts when branching on hyperplanes");
	  fBranchPsCosts = 0;
      }
  }
  pr->fLocalUB = fLocalUB;
  /*if (firstinstance) */InitLog2(); // all stat +++
  InitLog3(); // too detailed
  pr->cuts = &cuts; // enough?
  pr->forbidden = &allcols; // before Init();
  pr->Init();
//  lp->Init();
  timer.reset();  timer.start();
  tmRnd.reset(); // not start
  tmCG1.reset();
  tmCG2.reset();
  status = none;
  if (fPrintProblem) pr->PrintProblem(GetEnv().GetLog2());

  // +PrintInitRun(); // s^tart time
  // fLPOnly = ...
  // OUTPUT TIMER:
  AlarmHandler();
  if (outputInterval > 0)
    if (-1.34 == Alarm(outputInterval))
     log_ln("Output timer probably not installed");
  // INIT STATISTICS:
  NnvAll = NnvTst = 0;
  nTooLongColGen = 0;
  nErr2 = nErr2__;
} //____________________________________________________

void BCP::InitOptimize() {
//  gub = +1e100; // here not
// Here: mixing code & calls in 1 level!!
  llb = glb = llv = glv = glb0 = glv0 = llbLast=llvLast = -1e100;
  glb_others = glv_others = 1e100;
  depth = nHP = 0;
  fReoptTree = 0;
  InitRootNode();
  subgr->SetRHS(pr->b,pr->validsign);
  cntNode = cntInitNode = 0;
  maxObjCoef = 1; // or -1e100 ??? Before InitBasis() !

  InitBasis();
  fLPOpt = fMIPSolBest = false;
  fLevelCut = !fAddLevelCut;
  sum_iter = sum_iterAll = 0;
  iter = iterAll = 0; // BCP
  nextRNDNCols = 0;
	nMIPColsLast = -1;
  fInfeasColsIn = 0;
   // not 0 for the case of empty basis
  nCG = nCG_TOff = 0;
  nLeft = nRight = 0;
	
	nColsPrintLPVal = 0;

  ReadGivenSolution();
//  srand(12345); // for CPLEX -- no reaction
//  kNodeNewCuts = IMax(1,kNodeNewCuts); // 0: no local
} // init all vars properly

void BCP::InitRootNode() {
  nodes.AddRoot(Node());
  theNode = nodes.GetRoot();
  theNode->depth = depth;
  theNode->state = Node::current;
  pr->b_cg = pr->b; // no rhs reductions
//  list<Column*> empty;
//  pr->CopyForbCols(empty);
//  InitNewNode(); // not here. Only pr->b_cg like there
}

void BCP::InitBounding() { // init CPA
  iter = iterAll = 0; // BCP
  delReset = delReset0; delResetIterNext = int(delReset);
  llvLast = llv; //llbLast = llb; //-1e100;
  iterLLBLast = 0;
  assert(!nBranch || nHP == depth);
  // cuts' cntDel?
//  pr->InitBounding(); // e.g. copy forbidden cols !!! -- done in InitNewNode()
}

void BCP::DoneBounding() {
  sum_iter+=iter; sum_iterAll+=iterAll;
  theNode->llv = llv;
  theNode->llb = llb; // for DoneRun()
  UpdatePseudoCost();

  if (OUTP_LEV__ > 5.1) lp->WriteModel("donenode.lp");
    // makes problems!!!!

/// OUTPUTTING LP Solution:
  if (/*theNode->no == 612 or theNode->no == 1091 or
    theNode->no == 260.5 or theNode->no == 258*/0) {
    char fln [1924];
//    ostringstream(fln, sizeof(fln)) << "node" << theNode->no<<".lp" << ends;
    lp->WriteModel(fln);
    log_ln("LP solution in LP indexation:");
        for (int j=0;j<lpx.size();++j)
          if (lpx[j]) log__('['<<j<<"]:"<<lpx[j]<<' ');
    log_ln("LP solution in pool indexation:");
        for (int j=0;j<lpx.size();++j)
          if (lpx[j]) log__('['<<GetColIndex(j)<<"]:"<<lpx[j]<<' ');
  }

/// CHECK THAT THE BOUND IS NOT WORSE THAN ABOVE:
  bool BndMonot=true;
  if (not nTooLongColGen && not fSkipColGenAtLPBound
    && not fUseLagrange) {
    for (Node* pnode = theNode->parent;
      pnode;pnode=pnode->parent) {

    if (pnode ->fLPOpt && pnode->llv > llv + pr->GetVEps()) {
//      if (nErr2__ < 7)
        log__("_BND_WRS! To parent "<<pnode->no
          <<", diff="<<pnode->llv - llv);
      BndMonot = false;
    }
    }
    if (!BndMonot) { // extra for debug
//    if (OUTP_LEV__>=3) {
      if (nErr2__ < 7) {
        log__("\nWRS BND. LP SOL: ");
        for (int j=0;j<lpx.size();++j)
          if (lpx[j]) log__('['<<GetColIndex(j)<<"]:"<<lpx[j]<<' ');
      }
//    }
      CheckCuts();
      if (nErr2__ < 7)
      { __asErtm(false," Bound worse in a subproblem!!!"); }
      else
        ++ nErr2__;
    }
  }
}

void BCP::DoneOptimize() {
//   glv = glb = gub;                         WHY WAS THIS HERE?????????? TODO
  if (not nodes.SearchEmpty())
    glv = FMin(glv,(*nodes.GetBeginSearch())->llv);
  glb = pr->raster_ceil(glv);
  dbg_outn(1.23," N non-fathomed nodes w/o children (still open): "<<nodes.GetSearchSize());
}

void BCP::DoneRun() {
  if (//none==status or
   (error!=status && LPErr!=status)) {
    status = infeas;
    if (glb > -1e50) status = LPBnd;
    if (gub < 1e50) {
      if (Optimum()) status = opt; else status = feas;
    }
  }
  timeIP = timer.stop();
  timeRnd = tmRnd.userTime();
  timeCG1 = tmCG1.userTime();
  timeCG2 = tmCG2.userTime();
	ipCols = (int)colpool.size();
  if (OUTP_LEV__ > .99)
    lp->WriteModel("last.lp");
  fOutputTimer = 1;  PrintIter();
  log_n(1," STAT " << statusName[status]);
  if (nErr2 != nErr2__) mylog << " errors2: "<< (nErr2__-nErr2)<<' ';
  if (OUTP_LEV__ > 0.99) {
    if (nTooLongColGen)
      mylog << "nTooLongColGen:" << nTooLongColGen << endl;
    if (fRedCostBnd)
      mylog << "RedCostBounding:" << NnvAll <<" cases of "<<NnvTst << endl;
  }
  if (OUTP_LEV__ >= 2)
    mylog << '\n';
  pr->Done();
  PrintLog();
//  cntNode = 0; // ???
}

// =====================================================
// SERVICE UTILS

bool BCP::NextOptimize()
// Between runs: keep bounds, columns, else ?
{
  if ((nRuns >= nRunsMax)
    or ((nRuns==1) and fLPOnly
    // and not (LPErr == status) // no rnd yet
    )
    or Optimum())
    return false;
  // also:if time break or optimal or err or infeas or ?
  ++nRuns;
  // RANDOMIZATION: ensure effective changes
  return true;
} //____________________________________________________


}                                                // SS_END_NAMESPACE__
