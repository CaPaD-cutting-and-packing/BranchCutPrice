// FILE: bcp_cuts.cpp, cut management for branch&cut&price
// Author: Gleb <Belov@math.tu-dresden.de>

#include "stdafx.h"
#include "bcp.h"
#include "lasthdr.h"

#include "bbcuts.h"  // TODO: actually should use Problem's methods to access

namespace bgn {                     // SS_BEGIN_NAMESPACE__

//#define CHECK_CUTS

// How maxNCuts if del after add cuts?
// ClearCoefs called properly when needed..?

void BCP::RecalcCutsRHS() {
// NEEDED AFTER LEVEL CHANGE
// AND BEFORE COL GEN TO KNOW THE rhs
// USED TO CALC THE LAGRANGE BND.
//  Column col0(pr->b);
//  col0.SetObj(pr->GetLPBnd()); // !!!!
//  for (i=0;I < pr->b.size();++i) col0.PushID(i,pr->b[i]);
  CutList::iterator ic;
//  for_each_in(cuts,ic,) // also adding cut!
//    (*ic)->Clear(); // - this is too long
  ClearCoefs();
  pr->ba.resize(pr->Dim() + cuts.size());
  pr->GetLevelCut(); // Resets the level
  int i=pr->Dim();
  for_each_in(cuts,ic,) {
    pr->ba[i] = (*ic)->CalcRHS__(pr->b); // not b_cg
    lp->ChangeRHS(i,pr->ba[i]);
    ++ i;
  }
}
void BCP::ClearCoefs() {
    for_each_in(cutpool,ic,CutRefSet::iterator)
      (*ic)->ClearNonRec(); // to avoid too much recusrion in 2D,
//    for_each_in(bnds,ic,)  //  ???????
  //    (*ic)->ClearNonRec();
    pr->GetLevelCut()->ClearNonRec();
}
/////////////////// CUTS:
void  BCP::MarkTightCuts() { ////////////////////////
  if (cuts.empty()) return;
  // marking inactive:
  assert(lpd.size() >= pr->Dim() + cuts.size());
  log_n_(2," cuts");
	int nact=(int)cuts.size(); //, ndel=0;
  CutList::iterator ic;
  fNTD.resize(cuts.size());
  bool fDelAllInact = (iterAll == delResetIterNext);
  if (fDelAllInact) {
    delReset *= delResetInc;
    delResetIterNext += int(delReset);
    dbg_outn_(2,'%');
  }
  int i = pr->Dim();
  for_each_in(cuts, ic, ) {
    if (fabs(lpd[i]) < 1e-8 and (*ic)->CanBeDeleted()) {
      -- nact; // only now:
      fNTD[i-pr->Dim()]
        = not (0 >= --((*ic)->cntDel) || fDelAllInact);
    } else
      fNTD[i-pr->Dim()] = true; // not del active
    ++ i;
  }
  log_n_(2, "!"<<nact);
//  lp->BeforeModifCuts(); ////////////////////////
}

void BCP::DeleteSomeCuts(int mode) {
  if (cuts.empty()) return;
//  assert(lp->d.size() >= pr->Dim() + cuts.size());
  int i = pr->Dim();
  int ndel=0;
  CutList::iterator
    // ic1=cuts.begin(), 
    ic;
  if (1 <= mode) { // too long col gen, emergency del half
		int nc = (int)cuts.size();
    if (1 == mode) { // must reduce maxCuts
      log_ln("Too long col gen, trying to reduce maxCuts by 1<=n<=20%...");
      maxCuts = IMin(maxCuts,int(0.8*IMin(nc,maxCuts)));
      if (maxCuts < 1) {
				// THIS DOESN't HELP:
// 				maxCuts=1;
// 				BBCuts::nStepsTooMuch *= 4;
//         log_ln("WARNING: maxCuts=1, increasing iteration limit for BBCuts to "
// 					<< BBCuts::nStepsTooMuch);
				log_ln( " \n\n          WARNING: NO CUTS USED ANY MORE. \n\n\n");
			}
    }
    int nn=0;
    for_each_in(cuts,ic,) (*ic)->no = -1;
    for_each_in(cuts,ic,) (*ic)->Number(nn);
    assertm(nn==nc, // as long as branches are no cuts
      "n involved cuts="<<nn<<" n present cuts="<<nc);
    fNTD.resize(nc);
    fill_n(fNTD.begin(),nc,true);
    int i=0;
    int ndel = IMax(1+int(0.2 * rndCuts * nc), nc-maxCuts);
    for_each_in(cuts,ic,) {
      if ((*ic)->no >= nc-ndel) // >=: not forever
        fNTD[i] = false; // they _WILL_ be deleted
      ++i;
    }
  }

  int nn=0;
  for_each_in(cuts,ic,) (*ic)->no = -1;
  i=0;
  for_each_in(cuts,ic,)
    if (fNTD[i++])
      (*ic)->Number(nn);
  i=0;
  for_each_in(cuts,ic,)
     fNTD[i++] = ((*ic)->no != -1);

	for (i=(int)cuts.size()-1;i>=0;--i)
    if (not fNTD[i]) {
        DelCut(i); // to pool, mark as deleted
        ++ ndel;
    }
  log_n_(2,"-"<<ndel);
}

void BCP::DelCut(int iCut, bool fFakeDel) {   // fake del: only deleting to solve MIP
  int i, nd;
  CutList::iterator ic = cuts.begin()+iCut;
  if ((*ic)->Type() < 0 or (*ic)->Type() == 345)
    -- nHP;
	if (!fFakeDel) {
		(*ic)->cntDel0 *= cntDel0inc; // set to smth at creation
		(*ic)->cntDel = (*ic)->cntDel0;
	}
// ORDERING CONSISTENCY BETWEEN LP AND CUT LIST
  lp->DelRow(pr->Dim()+iCut);
  nd = 0;
	for (i=(int)cols.size()-1;i>=0;--i) // BACKWARDS !!!
    if ((*ic) == cols[i].slackCut) {
      ++nd; DelCol(i);
    }//goto SlackFound; NO! 2 slacks
  assertm(1+(pr->CanBeInfeasOnRestrPool()>0) == nd,
    "DelCut: not all cut slacks found in the LP");
// SlackFound:
  cuts.erase(ic);
}
bool BCP::FindViolCutsInPool() {
/* Choose some subset of cuts from the pool
  for them call ->CalcSlackValue()
*/
  int i;
  Vector<int> iNZ; // indices of non-0 vars
  Vector<double> xNZ; // the values
  iNZ.reserve(Dim()+depth);
  xNZ.reserve(Dim()+depth);
  CutList::iterator ic; // iterator(Cut*) ?
  CutRefSet::iterator ipc;
  cutViol.clear();
  for_each_in(cutpool,ipc,) (*ipc)->fPresent = false;
  for_each_in(cuts,ic,) (*ic)->fPresent = true;

  map<LPCut*,double> slackVals;
  // Filling with known slacks:
  for (i=0;i<cols.size();++i)
    if (cols[i].slackCut) {
     if (-1 == cols[i].j) // normal cut slack
      slackVals.insert(
        make_pair(cols[i].slackCut,lpx[i]));
    }
    else if (lpx[i]) {
      iNZ.push_back(GetColIndex(i)); // orig. index!
      xNZ.push_back(lpx[i]);
    }
  double sv;
  for_each_in(cutpool,ipc,)
    if (!(*ipc)->fPresent) {
      ClearCoefs();
      sv = (*ipc)->CalcSlackValue(iNZ,xNZ,slackVals);
      if (sv < -1e-9) // also '>=' cuts
        cutViol.push_back(CGIV(*ipc,sv));
    }

  return not cutViol.empty();
}

void BCP::ConstructNewCuts() {      ////// + level cuts
  if (not fLevelCut) return;
  if (nBranch) {
//    while (ModifyAndResolveLP());
    Node::CutCont1 cts;
    pr->SeparateHP(colpool, cols, lpx, lpd, cts);
    if (cts.empty())
      { }
    RestoreLP(); // here pr should provide original rhs etc.
    for_each_in(cts, ic, Node::CutCont1::iterator) {
      cutpool.insert(*ic);
      cutViol.push_back(CGIV(*ic, 1 )); // because we expect only 1 cut...
    }
    theNode->cutcont1.splice(theNode->cutcont1.end(),cts);
    return;
  }
  int i;
  Vector<Vector<double> > bi;
  Vector<double> x;
//  lp->GetSolution(&cstat,&lpx); // done in Solve ?
  try {
    lp->GetBasisInverse(bi,x); // x: only bas values
  }
  catch (...) {
    CheckCuts();
    throw;
  }
  // THE SOURCE OF DIFFS!!
  for (i=0;i<x.size();++i) {
    double xfi = frac(x[i]);
    int iCol = 0;
    // iCol = the col index of this bas var
    if (cols[iCol].slackCut ?
      cols[iCol].slackCut->IntegerSlack() : true and
      xfi>0.001 and xfi<0.999) // not too small fractionality
      ConstructSACut(bi[i],xfi,cstat,lpx); // marking as violated
  }
}

bool BCP::ModifyAndResolveLP() {
  if (not pr->TempModifyLP
    (colpool, cols, lpx, lpd, theNode->cutcont1))
    return false;
  // Modifying rhs:
  RestoreLP(); // just copy it from pr
    // care with older capacity cuts!
  /// Resolving:
  SolveDual_GetMultsOnly();
    do {  // infeas: can still price
      SolvePrimal(); // + some fast rounding
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
//      if (TimeLimit()) return 1; // "everything o.k."
    } while (Price()/* && not StopLP()*/);
  return true;
}

/// now just change rhs:
void BCP::RestoreLP() {
//  CutList::iterator ic;
//  ClearCoefs();
  int i=0;
//  for_each_in(cuts,ic,) {
//  for ( ; i < pr->Dim(); ++i)
//    lp->ChangeRHS(i,pr->b_cg[i]);
  for ( ; i < Dim(); ++i) {
//    pr->ba[i] = (*ic)->CalcRHS__(pr->b); // not b_cg
    lp->ChangeRHS(i,pr->ba[i]);
  }
}

void BCP::ConstructSACut
  (Vector<double> &bii, double xfi,
    Vector<LP::ColStatus> &cstat, d_vec &lpx) {
  int i=0;
  CutList::iterator ic;
  SACutSE cut;
  cut.cntDel0 = cut.cntDel = cntDel0; // ...
  cut.fSE = fSE;
// THE ALFA:
  if (!cutType) // MI cuts
    cut.alfa = xfi;
  else { // CG cuts
    cut.alfa = 1;
    // SOME HEURISTIC ENSMALLING OF Us?
    cut.cutType = cutType; // LIFTING CG cuts:
    if (0==CGNormMax) {
      int k = int(1/xfi);
      if (xfi * k > 0.999999)  -- k;
      for (i=0;i<bii.size();++i)  bii[i] *= k;
      xfi *= k;
    }
    else if (-1==CGNormMax) {
      if (xfi < 0.5) {
        for (i=0;i<bii.size();++i)  bii[i] = -bii[i];
        xfi = 1-xfi;
      }
    } else if (0 < CGNormMax) {
      int k = int(1/xfi);
      if (xfi * k > 0.999999)  -- k;
      int kMax = k;
      double frMax = xfi * k;
      for (;k<=CGNormMax;++k) {
        double fr = xfi * k;
        fr -= floor(fr);
        if (fr > frMax) { frMax = fr; kMax = k; }
      }
      for (i=0;i<bii.size();++i)  bii[i] *= kMax;
      xfi = frMax;
    }
    // ENSMALLING Us:
    if (CGFracParts)
    for (i=0;i<bii.size();++i)
      bii[i] = (bii[i] > 0) ? (bii[i]-floor(bii[i])) : (bii[i]-ceil(bii[i]));
  }
// u's for orig constr:
  cut.resize(pr->Dim());
  for (i=0;i<pr->Dim();++i) {
    cut.u[i] = bii[i];
/*    if (fSE)
      cut.uSE[i] = (0==pr->validsign[i]) ? 0:
        (pr->validsign[i] > 0) ? // negative slack coefs
          cut.F_Alfa( -cut.u[i])
          : - cut.F_Alfa( cut.u[i]);*/
  }
// u's for the cuts:
//  if (fLevelCut) ////// ???????
  for (i=pr->Dim(), ic=cuts.begin(); ic!=cuts.end();
    ++ ic, ++ i)
      if /*(0!=bii[i]) {*/(fabs(bii[i]) > 1e-15) {
        double uSE = 0;   // ??? 0-epsilon
/*        if (fSE) {
          uSE = (*ic)->LHS() > -1e50 ?
            cut.F_Bar( - bii[i]):
            - cut.F_Bar( bii[i]); // negated!
          // BUT
          if ((*ic)->IntegerSlack())
            uSE = (*ic)->LHS() > - 1e+50 ?
            cut.F_Alfa( -bii[i]) :
            - cut.F_Alfa( bii[i]);
        }*/
        cut.dep.push_back
          (SACutSE::Dependence(bii[i],uSE,*ic));
    }
  cut.dep.sort(); // FOR COMPARISON OF CUTS
  // Filling involved with their slack coefs:
  ClearCoefs();
  invCuts.clear();
  cut.ProduceListOfInvolved(invCuts);
  assert(invCuts.back() == &cut); // _cut_ disappears then
  invCuts.erase(invCuts.end()-1);
  for_each_in(invCuts,ic,)
    cut.invCuts.Add(SACutSE::Dependence(-1e100,*ic));
//////////////////// Dependence on var bounds: /////////
/// (upper + non0 lower)
  for (i=0;i<cstat.size();++i)
    if (LP::atUpper == cstat[i]) {
      assert(GetColIndex(i) >= 0);
      cut.bnds.Add
        (VarBnd(GetColIndex(i),true,(int)Round(lpx[i])));
    } else
    if (LP::atLower == cstat[i])
    if (fabs(lpx[i])>1e-6 && not IsSlackCol(i)) {
      assert(GetColIndex(i) >= 0); // no slacks
      cut.bnds.Add
        (VarBnd(GetColIndex(i),false,(int)Round(lpx[i])));
    }
  CutRefSet::iterator icut = cutpool.find(&cut); //no bnd
  if (cutpool.end() == icut) {
    SACutSE * pcut = theNode->cutcont.Add(cut);              // IS it modified here??????????
    if (!pcut)                       // HOW ??? TODO  Happens with hard28/22 on CentOS
        cerr << "\n\n\n    **************** CUT NOT IN THE POOL BUT IN THE NODE'S POOL?  ************** \n\n\n";
    else {
      cutpool.insert(pcut);
      cutViol.push_back(CGIV(pcut,xfi));
// CALC COEFS OF ALL CUTS FOR ALL COLS (needed for viol)
      pcut->rawCoef.reserve(int(1.3*(double)colpool.size()));
      for (i=0;i<colpool.size();++i)
	pcut->Calc__(GetMainCol(i),i);
      for_each_in(pcut->invCuts,iic,
	Pool<SACutSE::Dependence>::iterator) {
	ClearCoefs(); // before each!
	pcut->CalcCutSlackCoef__(iic->c);
      }
    }
  }
}
// ADDING: ALSO IN cuts
void BCP::AddLevelCut() { // always ?
  assert(false);
/*
  cuts.push_back(pr->GetLevelCut());
  pr->GetLevelCut()->cntDel = (long)1e+9;
  // - never to delete before a "cleanup"
//  lp->AddCut(pr->GetLevelCut()); // the niveau set alr.
  // lp ? AddFromPool ?
  // Then modify its rhs when changes. where ?
*/}
void BCP::AddSomeCuts() { // FROM THE POOL
  int i, nAddCuts;
  list<CGIV>::iterator icg;
  CutList::iterator ic;
  CutRefSet::iterator ipc;
  int nadd = 0;
//  fNTD.resize(cuts.size());
  if (not fLevelCut) {
    fLevelCut = true;
    AddLevelCut();
    dbg_outn_(2," +1");
//    fNTD.push_back(true);
    goto AfterAdd;
  }
//  Sorting cuts by randomized violation:
  if (cutViol.empty()) dbg_outn(2,"UUPS: No violated cuts!!!");
  for_each_in(cutViol,icg,)
    if (1==cutSelect) icg->wgt = float(icg->v * (1+(rndCuts - 0.5)*rndViolDev));
    else if (-1==cutSelect) icg->wgt = float(- icg->v * (1+(rndCuts - 0.5)*rndViolDev));
    else if (0==cutSelect) icg->wgt = float(- fabs(icg->v-0.5) * (1+(rndCuts - 0.5)*rndViolDev));
    else icg->wgt = float((1+(rndCuts - 0.5)*rndViolDev)); // all
  cutViol.sort();
	nAddCuts = IMin(MaxAllCuts() - (int)cuts.size(),
    1 + int(rndCuts * nIterkMax));
  if (nAddCuts<1) { log_ln("no cuts???"); return; }

// as cutpool contains also added cuts,
// care not to add again:
  for_each_in(cutpool,ipc,) (*ipc)->fPresent = false;
  for_each_in(cuts,ic,) (*ic)->fPresent = true;
  for(i=0;i<nAddCuts && not cutViol.empty()
    && nadd<nAddCuts;++i) {
    ClearCoefs();
    invCuts.clear();
    cutViol.front().c->ProduceListOfInvolved(invCuts);
    for (int j=0;j<invCuts.size();++j)
      if (not invCuts[j]->fPresent) {
        invCuts[j]->fPresent = true;
        AddCut(invCuts[j]);
        ++ nadd;
      } // adding involved first for max. feasibility
    cutViol.pop_front();
  }
  dbg_outn_(2," +"<<nadd);

AfterAdd:;
  if (OUTP_LEV__ > .99) lp->WriteModel("AddCuts.lp");
  // Checking that all 
  if (!fCheckCuts) return;
  CheckCuts();
  for (i=0;i<cuts.size();++i) {
    int j,k;
    invCuts.clear();
    ClearCoefs();
    cuts[i]->ProduceListOfInvolved(invCuts);
    for (j=0;j<invCuts.size();++j) {
      bool here=false;
      for (k=0;k<cuts.size();++k)
        if (invCuts[j] == cuts[k])
        { here = 1; break; }
      assertm(here, "Some involved cuts are not in the formulation!!!");
    }
  }
}
void BCP::AddCut(LPCut *pcut) {
// ORDERING CONSISTENCY BETWEEN LP AND CUT LIST
  int i;
  matind.clear(); matval.clear();
  double v;

  log_n_(3, "ADD CUT "<<pcut<<" alf="
    <<((SACutSE*)pcut)->alfa );
  cuts.push_back(pcut);
  if (pcut->Type() < 0 or pcut->Type() == 345)
    ++ nHP;
  CutList::iterator icc;
  ClearCoefs();
  pcut->CalcRHS__(pr->b); // value will be used
// Does it use the correct rhs of the level cut ??????
  for (i=0;i<cols.size();++i) {
    if (cols[i].slackCut) {
     if (-1 == cols[i].j) // only normal slack
     {
      ClearCoefs();
      v = pcut->GetCutSlackCoef__(cols[i].slackCut);
//      cout << "coef["<<cols[i].slackCut<<"]="<<v<<' ';
     }
    }
    else
      v = pcut->Calc__(GetCol(i),GetColIndex(i));
    if (v) { matind.push_back(i); matval.push_back(v); }
  }
  lp->AddRow
  ((int)matind.size(),&matind.front(),&matval.front(),pcut->GetRHS());
// + SLACK !!!
  matind.clear(); matval.clear(); ClearCoefs();
  for (i=0;i<cuts.size();++i) {
    v = cuts[i]->GetCutSlackCoef__(pcut); // also own
    if (v) { matind.push_back(pr->Dim() + i); matval.push_back(v); }
  }
  lp->AddCol(0,(int)matind.size(),&matind.front(),&matval.front(),0,1e100);
  cols.push_back(ColId(-1,pcut));
  if (pr->CanBeInfeasOnRestrPool()) { // Ax=b
    matind.clear(); matval.clear();
		matind.push_back(pr->Dim() + (int)cuts.size()-1);
    /// The coefficient of the infeasible column is the negative of that
    /// of the slack:
    matval.push_back(-(pcut->GetCutSlackCoef__(pcut))); // identity
    lp->AddCol(maxObjCoef*M__,
       (int)matind.size(),&matind.front(),&matval.front(),0,1e100);
    cols.push_back(ColId(-2,pcut)); // -2: for the "second" slack
      // -1 or -2 is distinguished?
      // deletion: both?
  }
  if (OUTP_LEV__ >= 5)
    lp->WriteModel("AddCut.lp");
}
void BCP::PrintCuts() {
//  cout<<"Number of cuts: " << lp->Dim()-pr->Dim()<<endl;

  if (OUTP_LEV__ < 5) return;
  cout<<"Cuts:\n";
  for_each_in(cuts,ic,CutList::iterator) {
    cout<<"\n\nCut "<<*ic<<": "<<flush;
    (*ic)->Print(cout); cout << endl;
  }
}
void BCP::PrintNCuts() {
  if (OUTP_LEV__ < 2) return;
  dbg_outn_(2," ="<<cuts.size());
}

struct CutInfo {
  SACutSE *pcut;
  d_vec u;
//  double alfa;
  d_vec lb, ub;
  void realloc(int n) {
    lb.resize(n); ub.resize(n);
    fill(ub.begin(), ub.end(), 1e+75);
  }
};
  class CmpPCut {
  public:
    bool operator() (const LPCut* c1,const LPCut* c2) const
    { return c1->no < c2->no; }
  };


#define Log__
#define Log_ln
void BCP::CheckCuts() {
  int i,j,r;
  int nn=0;
  if (nBranch)
    return;
	Vector<d_vec> A(pr->Dim() + (int)cuts.size());
  for (i=0;i<pr->Dim()+cuts.size();++i)
    A[i].resize(cols.size());
	Vector<d_vec> Araw((int)cuts.size());
  for (i=0;i<cuts.size();++i)
    Araw[i].resize(cols.size());

  for (i=0;i<cuts.size();++i) cuts[i]->no = -1;
  for (i=0;i<cuts.size();++i) cuts[i]->Number(nn);
  assert(cuts.size() == nn);
  CutList pc = cuts;
  sort(pc.begin(), pc.end(), CmpPCut());

	d_vec lb((int)cols.size()), ub((int)cols.size());
  for (j=0;j<cols.size();++j) {
    lp->GetVarBnds(j,lb[j],ub[j]);
  }

  Vector<CutInfo> cu(nn);
  for (i=0;i<nn;++i) {
    cu[i].pcut = (SACutSE*)(pc[i]);
		cu[i].realloc((int)cols.size());
    cu[i].u = cu[i].pcut->u;
    cu[i].u.resize(pr->Dim() + nn);
//    log_ln("Fill cut u's of cut"<<i<<" and checking var bounds...");
    for_each_in(cu[i].pcut->dep, idep, SACutSE::DepContainer::iterator)
      cu[i].u[pr->Dim()+idep->c->no] = idep->u;
    Pool<VarBnd>::iterator ib;
		d_vec lb_((int)allcols.size()), ub_((int)allcols.size(),1e+75);
    for_each_in(cu[i].pcut->bnds,ib,)
      if (ib->upper) ub_[ib->j] = ib->bnd;
      else lb_[ib->j] = ib->bnd;
    // CONVERT TO LP INDEXING:
    for (j=0;j<cols.size();++j)
      if (not IsSlackCol(j)) {
        cu[i].lb[j] = lb_[GetColIndex(j)];
        assert(cu[i].lb[j] <= lb[j]+1e-6);
        cu[i].ub[j] = ub_[GetColIndex(j)];
        assert(cu[i].lb[j] > -1e-6);
        if (cu[i].lb[j] > 1e-6) assert(cu[i].ub[j] > 1e+50);
        if (cu[i].ub[j] < 1e+50) {
          assert(fabs(cu[i].lb[j]) < 1e-6);
          assert(cu[i].ub[j] >= ub[j]-1e-6);
        }
      }
  }
  // EXTRACT LP COLS, ALSO CUT SLACKS:
//  log_ln("Filling LP cols from own memory. N Cols="<<cols.size());
  for (j=0;j<cols.size();++j) {
//    log__(j);
    if (-1 == cols[j].j) {
      A[pr->Dim() + cols[j].slackCut->no][j] = 1;
//      log__('s');
    }
    else {
      int i;
      for (i=0;i<GetCol(j)->id.size();++i)
        A[GetCol(j)->id[i].i][j] = GetCol(j)->id[i].d;
//      log__(' ');
    }
  }
  // CALCULATE CUT COEFS:
  //for (r=0;r<nn;++r) log__("u["<<r<<"].size() == "<<cu[r].u.size()<<' ');
//  log__("\ncut coefs.");
  for (j=0;j<cols.size();++j) {
//    log__("\nCol "<<j<<": ");
    for (r=0;r<nn;++r) {
      double sum=0;
      for (i=0;i<pr->Dim() + r;++i)
        sum += cu[r].u[i] * A[i][j];
      Araw[r][j] = sum;
      if (-1 == cols[j].j) {
        if (r > cols[j].slackCut->no)
          A[pr->Dim() + r][j] = SACutSE::F_Bar(cu[r].pcut->alfa,sum,cutType);
      }
      else
      if (cu[r].ub[j] < 1e+50)
        A[pr->Dim() + r][j] = -SACutSE::F_Alfa(cu[r].pcut->alfa, - sum,cutType);
      else
        A[pr->Dim() + r][j] = SACutSE::F_Alfa(cu[r].pcut->alfa, sum,cutType);
//      log__(A[pr->Dim() + r][j]<<' ');
    }
  }
//  log__('\n');

   // COMPARING LP coefs: 1st basic constraints
//      log__("\n cmp orig coefs:\n");
  for (j=0;j<cols.size();++j)
    for (i=0;i<pr->Dim();++i)
      if (fabs(lp->GetLPCoef(i,j) - A[i][j]) > 1e-6) {
        log__(" ("<<i<<','<<j<<"): "<<lp->GetLPCoef(i,j) <<','<<A[i][j]);
        if (NULL != cols[j].slackCut) log__(" -slk ");
        assert(0);
      }
   // cut coefs:
//      log__("\n cmp cut coefs:\n");
  for (j=0;j<cols.size();++j)
    for (r=0;r<nn;++r)
      if (fabs(lp->GetLPCoef(pr->Dim()+r,j) - A[pr->Dim()+cuts[r]->no][j]) > 1e-6) {
        log__(" ("<<pr->Dim()+r<<','<<j<<"): "
          <<lp->GetLPCoef(pr->Dim()+r,j) <<','<<A[pr->Dim()+cuts[r]->no][j]);
        if (NULL != cols[j].slackCut) log__(" -slk ");
        assert(0);
      }
  // RHS:
  d_vec b = pr->b;
  b.resize(pr->Dim() + nn);
//    log__("\nRHS. ");
    for (r=0;r<nn;++r) {
      double sum=0;
      for (i=0;i<pr->Dim() + r;++i)
        sum += cu[r].u[i] * b[i];
      for (j=0;j<cols.size();++j) {
        if (cu[r].lb[j] > 1e-6) sum -= Araw[r][j] *cu[r].lb[j];
        if (cu[r].ub[j] < 1e+50) sum -= Araw[r][j] *cu[r].ub[j];
      }
      sum = SACutSE::F_Alfa(cu[r].pcut->alfa, sum,cutType);
      for (j=0;j<cols.size();++j) {
        if (cu[r].lb[j] > 1e-6) sum += A[pr->Dim()+r][j] *cu[r].lb[j];
        if (cu[r].ub[j] < 1e+50) sum += A[pr->Dim()+r][j] *cu[r].ub[j];
      }
      b[pr->Dim() + r] = sum;
//      log__(b[pr->Dim() + r]<<' ');
      if (fabs(b[pr->Dim() + r] - lp->GetLPCoef(pr->Dim()+r,-1)) > 1e-6) {
        log__(" but in LP: "<<lp->GetLPCoef(pr->Dim()+r,-1)<<' ');
        assert(0);
      }
    }
}

}                                                // SS_END_NAMESPACE__
