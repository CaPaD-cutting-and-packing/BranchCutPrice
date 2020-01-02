// FILE: bcp_barnch.cpp, branch&cut&price for CSP and 2CP
// Author: Gleb <Belov@math.tu-dresden.de>

#include "stdafx.h"
#include "bcp.h"
#include "solver.h"
#include "lasthdr.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

// Write level-wise

void BCP::StartReoptTree() {
  if (fCGTailOff) { // only if the corresp. option
    // and the number of not node->fLPOpt ??
    fReoptTree = 1;
    log_n_(1.3, " REOPT ");
  }
}

void BCP::Branch() {
/* storing info _IN_ the current node + cleaning WHAT?,
  getting LP sol,
  saving basis +ub indexes but care: del slacks?
  +cuts in the node
  +add sons knowing the father AND father is to know them
*/
  SaveBasis(); // with array reserving
   // even if theNode->fLPOpt? Delete after that (Reconstruct) ?

  if (depth >= maxDepth)
    return; // No branching deeper
  if (fReoptTree)
    return; // exit here!

  Create2Sons();
  CopyBasTo(theNode->sonLeft);
  CopyBasTo(theNode->sonRight);

  if (0 == nBranch) // attention: meaning of nBranch
    BranchOnVar();
  else
    BranchOnHyperplanes();
}

void BCP::CopyBasTo(Node * pn) {
  pn->iUpper = theNode->iUpper;
  pn->slkBas = theNode->slkBas;
  pn->iBas = theNode->iBas;
  pn->cuts = theNode->cuts; // STORE IN  PARENT!!! But when del?
}

void BCP::SaveBasis() {
  int nU=0, nCS=0, nB=0; unsigned i;
  for (i=0;i<lpx.size();++i) {
    if (LP::atUpper == cstat[i])
      ++ nU; // may be quite different dep. on formulation
    if (LP::basic == cstat[i]) {
      if (cols[i].slackCut)  ++ nCS;
      else  ++ nB;
    }
  }
  theNode->iUpper.clear(); theNode->iUpper.reserve(nU);
  theNode->slkBas.clear(); theNode->slkBas.reserve(nCS);
  theNode->iBas.clear(); theNode->iBas.reserve(nB);

  /// Saving basis & vars on upper bounds
  dbg_outn_(3," Savebas. Non0Cols not at LB: ");
  for (i=0;i<lpx.size();++i) {
    if (LP::basic == cstat[i])
      if (cols[i].slackCut) {
        theNode->slkBas.push_back(cols[i].slackCut);
        dbg_outn_(3,cols[i].slackCut<<"c="<<lpx[i]<<' ');
      }
      else {
        theNode->iBas.push_back(cols[i].j);
        dbg_outn_(3,cols[i].j<<':'<<lpx[i]<<' ');
      }
    else
      if (LP::atUpper == cstat[i]) {
        theNode->iUpper.push_back(cols[i].j);
        dbg_outn_(3,cols[i].j<<"u="<<lpx[i]<<' ');
      }
  }
  theNode->cuts = cuts; // STORE IN  PARENT!!! But when del?
    // and when storing branching hyperplanes in cuts, need it
}

void BCP::Create2Sons() {
  /// Creating one node
  Node n1; // right, >=

  n1.parent = &*theNode;
  n1.no = cntNode+1;
  n1.depth = depth+1;
  n1.llv = llv;
  n1.llb = llb;
  n1.state = Node::open;
  n1.nDesc = 0;
  n1.nLeft = theNode->nLeft;

/// Copying the node:
  Node n2 = n1; // left, <=
  n2.nLeft ++;
  n2.no += 0.5;
  theNode->nDesc = 2; // decrease when del.

/// Adding to the tree:
  theNode->sonRight = &*(sonRight=nodes.AddOpen(n1));
  theNode->sonLeft = &*(sonLeft=nodes.AddOpen(n2));

/// Preparing backtracking for diving
  if (BrStratRL >= 0.5) { // more left sons
    nextNodeDFS = sonLeft;
    backtrack.push_back(sonRight);
  } else { // less left sons
    nextNodeDFS = sonRight;
    backtrack.push_back(sonLeft);
  }
  fLastNodeSplit = true;
}

void BCP::BranchOnVar() {
  int iMin, brUB;
  pair<int, double> varsel;
  Node *pnode;

  /// SELECT VARIABLE:
  if (fBranchPsCosts) varsel=SelectBrVar_PsCosts();
  else varsel=SelectBrVar_MostInfeas();
  iMin = varsel.first;
  brUB = (int)(varsel.second);
  // CHECKING THAT THIS VAR DOES NOT VIOLATE A BOUND:
  const int jMin = GetColIndex(iMin);
  log_n(2,"VarSel: x"<<jMin<<", val="<<lpx[iMin]<<", ub="<<brUB);

  for (pnode = &*theNode;pnode;pnode=pnode->parent)
    for_each_in(pnode->bnds,ib,
      Node::BndCont::iterator)
      if (jMin == ib->j) {
        if (ib->upper)
          assertm(lpx[iMin] - ib->bnd < 1e-6,
          "\nUpper bnd on "<<iMin<<": "
            <<ib->bnd<<", val:"<<lpx[iMin])
        else
          assertm(lpx[iMin] - ib->bnd > -1e-6,
          "\nLower bnd on "<<iMin<<": "
            <<ib->bnd<<", val:"<<lpx[iMin])
  }
  sonLeft->xValue
    = sonRight->xValue = lpx[iMin];
  sonRight->bnds.reserve(1);
  sonRight->bnds.push_back(
    VarBnd(GetColIndex(iMin),false,brUB+1));
  sonLeft->bnds.reserve(1);
  sonLeft->bnds.push_back(
    VarBnd(GetColIndex(iMin),true,brUB));
}

void BCP::BranchOnHyperplanes() {
  LPCut *pcL, *pcR;
  pr->BrOnHP(colpool, cols, lpx, pcL, pcR);
  /// Store the pointers (which are allocated by new()):
  sonLeft->cutcont1.push_back(pcL); // not my_auto_ptr
//  sonLeft->cutcont1.back().reset(pcL);
  sonRight->cutcont1.push_back(pcR);
//  sonRight->cutcont1.back().reset(pcR);
  /// Now say that this new "cut" belongs to the formulation:
  sonLeft->cuts.push_back(pcL);
  sonRight->cuts.push_back(pcR);
  // how should we then know that these are branching cuts??
}


pair<int, double> BCP::SelectBrVar_MostInfeas() {
  double dMin=1e100; int i, iMin=-1;
  for (i=0;i<lpx.size();++i)
    if (LP::basic == cstat[i]
      && not IsSlackCol(i)) { // why not on slacks?
      double lbi, ubi;
      lp->GetVarBnds(i, lbi, ubi);
      ubi = FMin (ubi, lbi + // !!!
        LocalUpperBound(GetCol(i))); // !!! if no ub
      double the_frac = frac(lpx[i]);
      if (BrVarFrac>1)
        the_frac = RndBrFrac;
      if (fabs(lbi - ubi) > 1e-6)
        if (fabs(lpx[i]-floor(lpx[i]+1e-3))>1e-3 // fractional enough; 1e-6 too small when BrVarFrac != 0.5
         and 
         fabs(the_frac - BrVarFrac)
          + wTillUB * (ubi - lpx[i])/(ubi-lbi) // not so good
// 				the_frac * (the_frac-1.0)          // SCIP !!
          < dMin) {
          dMin = fabs(the_frac - BrVarFrac) + wTillUB * (ubi - lpx[i])/(ubi-lbi);
// 					dMin = the_frac * (1.0-the_frac);
          iMin=i;
        }
    }
  assert(iMin>=0);
  double newub = floor(lpx[iMin]);
      double lbi, ubi;
      lp->GetVarBnds(iMin, lbi, ubi);
      assert2m(fabs(lbi - ubi) > 1e-6, " fixed var is basic ??? ");
  if (newub > ubi-1e-6) newub = ubi - 1; // >=
  if (newub < lbi-1e-6) newub = lbi; // <
  return pair<int,double>(iMin, newub);
}

pair<int, double> BCP::SelectBrVar_PsCosts() {
  int i, iMax=-1, nInit=0;
  double dFrac = dInfeasFracPart;
  if (dFrac>=0.5 or dFrac<0) dFrac=0.001;
  int nf=0; // N candidates
  for (i=0;i<lpx.size();++i)
    if (LP::basic == cstat[i])
    if (fabs(frac(lpx[i]) - 0.5) < (0.5 - dFrac)
      && not IsSlackCol(i)) // why not on slacks?
      nf++;
  if (nf<3) dFrac = pr->GetXEps();
// ALL CANDIDATES: new cols when testing ??
// Some informations will be changed, incl.
// those needed later in node splitting. SAVING:
  lpx__ = lpx; cstat__ = cstat;
//  double clv__ = llv, clb__ = llb; // !!!! local to current
    // + global too
  double PsCMax = -1e100;
  double dInfeasMax = 0;
  for (i=0;i<lpx__.size();++i)
    if (LP::basic == cstat__[i])
    if (fabs(frac(lpx__[i]) - 0.5) < (0.5 - dFrac)
      && not IsSlackCol(i) // why not on slacks?
      && not TimeLimit())
    {
      if (0 == GetCol(i)->nLB) {
        InitPsCost(i,0,lpx__[i]);
        ++ nInit;
      }
      if (0 == GetCol(i)->nUB) {
        InitPsCost(i,1,lpx__[i]);
        ++ nInit;
      }
      double f=frac(lpx__[i]);
      double t=PseudoCost(GetCol(i), f);
      if (t>PsCMax or
        (t>PsCMax - pr->GetRCEps() and // equal psc. VEps ?
          FMin(f,1-f)>dInfeasMax)) {
        iMax=i; PsCMax=t; dInfeasMax=FMin(f,1-f);
      }
    }
  lpx = lpx__; cstat = cstat__; // used later? Yes, in node creation
//  clv = clv__; clb = clb__;
//  UpdateLocalLPV(); // return to the pre-test values
  // Reoptimize if InitPsC called ?
  // Reason: provide starting point
  // but CPLEX should be able to start from any node...
//  if (nInit)
  //  SolvePrimal();

  double newub = floor(lpx__[iMax]);
      double lbi, ubi;
      lp->GetVarBnds(iMax, lbi, ubi);
      assert2m(fabs(lbi - ubi) > 1e-6, " fixed var is basic ??? ");
  if (newub > ubi-1e-6) newub = ubi - 1; // >=. MIST
  if (newub < lbi-1e-6) newub = lbi; // <

  return pair<int,double>(iMax, newub);
}

void BCP::InitPsCost(int iCol, bool upper, double x) {
  int i;
	int n__ = (int)allcols.size();
	int nc__ = (int)cols.size();
	d_vec lb0((int)allcols.size()),
	ub0((int)allcols.size());
  for (i=0;i<cols.size();++i)
    if (!IsSlackCol(i))
      lp->GetVarBnds(i,lb0[GetColIndex(i)],ub0[GetColIndex(i)]);
  lb.resize(allcols.size());
  ub.resize(allcols.size());
  for (i=0;i<lb0.size();++i) {
    lb[i]=(int)lb0[i];
    ub[i]=(ub0[i]>(double)INT_MAX) ? INT_MAX: (int)ub0[i];
  }
  if (upper) ub[GetColIndex(iCol)]=int(floor(x));
  else lb[GetColIndex(iCol)]=int(ceil(x));
//  lp->ChangeVarBnds(iCol,lb,ub);

  bool lpOpt = 1; //SolveDual_GetMultsOnly(); // for col gen
  // PREPROC:
  if (!TightenVarBounds())
  { lpOpt = 0; goto CGdone; }
  UpdateVarBounds();

  if (!SolveDual()) {
    if (!pr->CanBeInfeasOnRestrPool()) // ?????
      lpOpt = 0; // + SolvePrimal ?
    else {
      if (!StartPhaseI())
        lpOpt = 0;
    }
  }

  /////////// NOW GENERATE NEW COLS:
  if (lpOpt) {
    InitLP(); // meaning of LP::error ???
    do {  // infeas: can still price
      SolvePrimal(); // + some fast rounding
      if (LP::error==lp->status) { // +unbounded? not infeas
        lpOpt = 0;
        log_n_(3," e!!?_"); goto CGdone;
      }
    } while (Price()/* && not TimeLimit()*/);
    if (LP::infeas == lp->status)
      lpOpt = 0;
  }

CGdone:
  log_n_(2," +"<<allcols.size()-n__<<'c');
  double z=gub;
  if (lpOpt) z=GetLPValue();
  double f=frac(x);
  if (upper) {
    GetCol(iCol)->nUB=1;
    GetCol(iCol)->sumPsCUB=(z-llv)/f;
  } else {
    GetCol(iCol)->nLB=1;
    GetCol(iCol)->sumPsCLB=(z-llv)/(1.0-f);
  }
  dbg_outn(2," tst psc["<<GetColIndex(iCol)
    <<"]: fUpper="<<upper<<", ZDual="<<z);
  // NO COL GEN

//  lp->ChangeVarBnds(iCol,lb0,ub0);
  //////////// RESTORE BOUNDS:
  // FOR OLD COLS:
  for (i=0;i<nc__;++i)
    if (!IsSlackCol(i))
      lp->ChangeVarBnds(i,lb0[GetColIndex(i)],ub0[GetColIndex(i)]);
  // FOR NEW COLS:
  // ? In UpdateBounds ?
}

double BCP::PseudoCost(Column *c, double f) {
  assert(c->nLB && c->nUB);
  double psu=c->sumPsCUB/c->nUB*f;
  double psl=c->sumPsCLB/c->nLB*(1.0-f);
  return Alfa1*FMin(psu,psl) + FMax(psu,psl);
//  return Alfa1*psl + psu;
}

void BCP::UpdatePseudoCost() { // After solving node
  if (!depth || !fBranchPsCosts) return;
  int iC=theNode->bnds.front().j;
  double f = frac(theNode->xValue);
  Column *c=GetMainCol(iC);
  double dZ=llv - theNode->parent->llv;
  if (LP::infeas == lp->status)
    dZ = (gub-theNode->parent->llv);
  if (theNode->bnds.front().upper) {
    c->nUB ++;
    c->sumPsCUB += dZ/f;
  } else {
    c->nLB ++;
    c->sumPsCLB += dZ/(1.0-f);
  }
}

void BCP::Fathom(SolutionTree::node_iterator itF) {

/// Delete the node physically only later because it
/// stores its cuts which may be in the formulation
    Node * pnode = &*itF,  * parent;
    if (pnode->nDesc) { // GOING DOWN
      assert(fReoptTree);
      for (;;) {
        if (pnode->sonLeft) {
          assert(pnode == pnode->sonLeft->parent);
          -- (pnode->nDesc);
          assert((pnode->nDesc) >= 0);
          pnode = pnode->sonLeft;
          pnode->parent->sonLeft = NULL; // to forbid going again
          nodes.MarkFathomed(pnode);
        }
        else 
        if (pnode->sonRight) {
          assert(pnode == pnode->sonRight->parent);
          -- (pnode->nDesc);
          assert((pnode->nDesc) >= 0);
          pnode = pnode->sonRight;
          pnode->parent->sonRight = NULL; // to forbid going again
          nodes.MarkFathomed(pnode);
        }
        else {
          if (pnode == &*itF) break;
          assert(pnode->parent);
          pnode = pnode->parent;                              
        }
      }
    }
    /// After that,  GOING UP:
    for ( ; ; pnode=parent) {
//      if (pnode->nDesc)  break;
      nodes.MarkFathomed(pnode);
      parent = pnode->parent;
      if (!parent) break;
      -- (parent->nDesc);
      assert((parent->nDesc) >= 0);
      if (parent->sonLeft == pnode)
        parent->sonLeft = NULL;
      else {
          assert(parent->sonRight == pnode);
          parent->sonRight = NULL;
      }
      if (0 < parent->nDesc) break;
    }
}

void BCP::DeleteFathomed() {
  nodes.DeleteFathomed();
}

bool BCP::Select() {
Loop:
//  if (ListEmpty()) return false;
  if (not SelectNextNode()) return false;  // sets llb
//  PrintConstr(); /////////-------------------------
  if (theNode->llb >= gub)
    { Fathom(theNode); goto Loop; }
  if (ContradictionsInNewNode())
    { Fathom(theNode); goto Loop; }
//  SetByLogImp();
  if (not InitNewNode())
  { Fathom(theNode); goto Loop; }
  return true; // dual simplex
}

bool BCP::SelectNextNode() {
  /* dfs: deepest node & best bound.
  WHAT IS BEST BOUND ? FMax or FMin? (Balas/ABACUS) Option.
  */
  // all nodes are left in the tree at first

  /// FATHOMING or REOPTIMIZING suboptimal nodes:
  // if (fReoptTree) { // starting from the highest
  //   bn = nodes.GetRoot();
  //   assert(!bn->depth);
  //   ....
  // } else
    SolutionTree::INode bn;
    while (nodes.GetLastOpen(bn)) {
     // actually quite the opposite direction as when reopt
      if (bn->llb > gub) {
        Fathom(bn);
	nodes.DeleteSearch(bn);
      }
      else break;
    }

  if (nodes.SearchEmpty()) return false;
  ++ cntNode;
	
	/// Just testing Emergency AFF switch:
// 	static int nSwitch=0;
// 	if (100 == cntNode)
// 		if (!nSwitch++)
// 			throw double(8.1);

// 		throw double(8.0);
	

  int nStr // Strategy. 0=bfsdfs, 1=dive
    = (nBFSDFS + nDive > 0) ?
    (cntNode % (nBFSDFS+int(nDive*pr->Dim())) < nBFSDFS) ? 0: 1 : 0;
/*  if (nBFS + nDFS__ > 0)
    if (cntNode % (nBFS+nDFS__) == nBFS)
      nStr = 2; // just the beginning of BFS in the cycle
      */
  // When mixed, at first DFS ????
  if (nStr) {// when subtree over before end of cycle,
    // do exactly nBFSDFS selections?
    if (fLastNodeSplit)
    {
      bn = nextNodeDFS;
      nodes.DeleteSearch(bn);
      log_n_(3," NextNd");
      goto TakeNode;
    }
    if (not backtrack.empty()) // no backtracking
    {
      bn = backtrack.back();
      backtrack.pop_back();
      nodes.DeleteSearch(bn);
      log_n_(3," BTrk");
      goto TakeNode;
    }
    // when NOTHING FOUND from the last node:
//    if (2==nStr) {// DFS just started
    nStr = 0; // then still BFS hoping that
    log_n_(1.3," Bfs -- subtree over ");
//    }
  }
  else
    backtrack.clear(); // when bfs

  // BFSDFS:
  bn = nodes.ExtractBestOpen();

//  log_n_(2.3, '|');

//  oldNodes.splice(theNode .. where to ?
// somehow mark old nodes + the last one
TakeNode:
  fLastNodeSplit = false; // reset;
  llb = bn->llb;  // ???
  llv = bn->llv;  // ???
  depth = bn->depth; // ........
// erase last node if no sons ??? No, it is marked...
// and if erase then not here...
  lastNode = theNode;
  theNode = bn;
  theNode->state = Node::current;
  return true;
}

bool BCP::ContradictionsInNewNode() { return false; }
// + clearance? or in Branch() or Select()?

bool BCP::InitNewNode() {
  // when reconsider active cuts
  /* Fill the set of forbidden cols
  update the cuts & branches uncommon with the last node
  also the slacks of the cuts
  compile the cut pool from all parents' local cuts
  reset cuts' cntDel and cntDel0 (for a new cutting cycle)
  reduce problem's RHS (INTEGER BOUNDING),
  recalc cut coefs/rhs, set starting basis
  */ // CUTLIST !!!!
  llv = theNode->llv;  // ???
  llb = theNode->llb;  // ???
  depth = theNode->depth; // ........
  nLeft += theNode->nLeft;
  nRight += depth - theNode->nLeft;

  log_n_(2,"\n INInode " << cntNode << " (" <<theNode->no);
	if (theNode->depth)
		log_n_(2, "<-"<<theNode->parent->no);
	log_n_(2, ") ");

  CompileColumns();
  CompileVarBounds();
  if (!TightenVarBounds()) return false;

  CompileCuts(); // including "branching" cuts
  UpdateVarBounds();
  UpdateCuts();

  RestoreBasis();

    // + var fixing

  if (!SolveDual()) {
    if (!pr->CanBeInfeasOnRestrPool()) // ?????
      return false; // + SolvePrimal ?
    else {
      if (!StartPhaseI())
        return false;
    }
  }
  // COL GEN NEEDED HERE ?
  ++ cntInitNode;

  RecalcGLB();
  DeleteFathomed();
// + set of forbidden ? Automatically

  if (OUTP_LEV__ > 5.1) lp->WriteModel("node.lp");
  return true;
}

void BCP::CompileColumns() {
  if (not nReducePool)
    return;
  Node * pnode;
  pnode = &*theNode;
//  ubnds.clear();
//  dbg_outn_(3.5,"\nBNDS:");
  int i,j;
//  for (i=0;i<cols.size();++i)
  //  if (not IsSlackCol(i))
    //  GetCol(i)->nHidden = 
  dbg_outn_(3," ColsMainBeforeModif*: ");
  for (i=0;i<colpool.size();++i) {
        dbg_outn_(3,i<<(char)('a'+GetMainCol(i)->Hidden()));
  }
  for (;pnode;pnode = pnode->parent) {
//    dbg_outn_(3.5," depth="<<pnode->depth);
    for_each_in(pnode->bnds,ib,Node::BndCont::iterator) {
//      dbg_outn_(3.5," bnd"<<(ib->upper?'U':'L')<<'['<<ib->j<<"]="<<ib->bnd);
      if (1 == GetMainCol(ib->j)->Hidden())
        GetMainCol(ib->j)->nHidden = 2;
      else // can be double bounds
      if (0 == GetMainCol(ib->j)->Hidden())
        GetMainCol(ib->j)->nHidden = 3;
    }
  // AND NOT FOR RED COST BOUNDS: (added if only possible)
  }
   // ALSO BASIC: (can be in the bounds => not modify any more)
  for (i=0;i<theNode->iBas.size();++i)
    if (not GetMainCol(theNode->iBas[i])->IsSlack()) { // ???
      if (1 == GetMainCol(theNode->iBas[i])->Hidden())
        GetMainCol(theNode->iBas[i])->nHidden = 2;
      else
      if (0 == GetMainCol(theNode->iBas[i])->Hidden())
        GetMainCol(theNode->iBas[i])->nHidden = 3;
    }
   // ALSO At Upper Bound:
  for (i=0;i<theNode->iUpper.size();++i)
//    if (not GetMainCol(theNode->iUpper[i])->IsSlack()) // ???
      if (1 == GetMainCol(theNode->iUpper[i])->Hidden())
        GetMainCol(theNode->iUpper[i])->nHidden = 2;
      else
      if (0 == GetMainCol(theNode->iUpper[i])->Hidden())
        GetMainCol(theNode->iUpper[i])->nHidden = 3;
  dbg_outn_(3," ColsMain*: ");
  for (i=0;i<colpool.size();++i)
/*      if (cols[i].slackCut) {
        dbg_outn_(3,cols[i].slackCut<<"c ");
      }
      else */{
        dbg_outn_(3,i<<(char)('a'+GetMainCol(i)->Hidden()));
      }
  // ADD MISSING:
  for (j=0;j<colpool.size();++j)
    if (2 == GetMainCol(j)->nHidden)
      AddColToLP(lp.get(), GetMainCol(j), j);
  // DELETE UNNEEDED:
  for (i=(int)cols.size()-1;i>=0;--i) // backwards
    if (not IsSlackCol(i)) // NOT SLACKS !!
      if (0 == GetCol(i)->nHidden) {
        GetCol(i)->nHidden = 1;
        DelCol(i);
      }
  for (i=0;i<cols.size();++i)
    if (not IsSlackCol(i))
      GetCol(i)->nHidden = 0;
  dbg_outn_(3," ColsInLP*: ");
  for (i=0;i<cols.size();++i)
      if (cols[i].slackCut) {
        dbg_outn_(3,cols[i].slackCut<<"c ");
      }
      else {
        dbg_outn_(3,cols[i].j<<' ');
      }
}

void BCP::CompileVarBounds() {
  // COMPILING Var bounds, also by reduced cost bounding
  lb.resize(allcols.size()); ub.resize(allcols.size());
  fill(lb.begin(),lb.end(),0);
  fill(ub.begin(),ub.end(),INT_MAX);
  Node * pnode;
  pnode = &*theNode;
//  ubnds.clear();
  dbg_outn_(3.5,"\nBNDS:");
  for (;pnode;pnode = pnode->parent) {
    dbg_outn_(3.5," depth="<<pnode->depth);
    for_each_in(pnode->bnds,ib,Node::BndCont::iterator) {
      dbg_outn_(3.5," bnd"<<(ib->upper?'U':'L')<<'['<<ib->j<<"]="<<ib->bnd);
      if (ib->upper) {
        ub[ib->j] = IMin(ub[ib->j],ib->bnd);
//        ubnds.push_back(*ib);
      } // IMin: as long as ub, lb int
      else lb[ib->j] = IMax(lb[ib->j],ib->bnd); // lower
    }
  // AND ALSO RED COST BOUNDS: (added if only possible)
    for_each_in(pnode->rcBnds,irc,Pool<VarBnd>::iterator) {
      dbg_outn_(3.5," rcBnd"<<(irc->upper?'U':'L')<<'['<<irc->j<<"]="<<irc->bnd);
      if (irc->upper) {
        ub[irc->j] = IMin(ub[irc->j],irc->bnd);
//        ubnds.push_back(*ib);
      } // IMin: as long as ub, lb int
      else lb[irc->j] = IMax(lb[irc->j],irc->bnd); // lower
    }
  }
}

char BCP::TightenVarBounds() {
  int i;
  // NOW removing lower-bounded cols from b_cg:
  // (r-h side for col gen):
  if (fLocalReduce or fLocalUB) {
    pr->b_cg = pr->b;
    for (i=0;i<lb.size();++i)
      if (lb[i]) {
        assert( not GetMainCol(i)->Hidden());
        for_each_in(GetMainCol(i)->id,iid,Column::iterator) {
          pr->b_cg[iid->i] -= iid->d * lb[i];
					if (pr->b_cg[iid->i]<0)
						pr->b_cg[iid->i]=0;                        // !!!!!!!!!!!
				}
      }
  }
  // RESETTING upper bounds:
  if (OUTP_LEV__ >= 4)
  { mylog << ' '; PrintVec(mylog,pr->b_cg); }
  if (fLocalUB) { // INTEGER BOUNDING
    for (i=0;i<ub.size();++i) {
    if (not GetMainCol(i)->IsSlack()) {
      int lub = LocalUpperBound(GetMainCol(i));
      dbg_outn_(4," lub["<<i<<"]="<<lub);
      if (lub<0) { dbg_outn_(4," _F!!"); return 0; }
      int ub1 = lb[i] + lub;
      if (ub1<ub[i]) ub[i] = ub1;
    }
    }
  }
  else if (fLocalReduce) {
    for (i=0;i<ub.size();++i) {
    if (not GetMainCol(i)->IsSlack()) {
      int lub = LocalUpperBound(GetMainCol(i));
      dbg_outn_(4," lub["<<i<<"]="<<lub);
      if (lub<0) { dbg_outn_(4," _F!!"); return 0; }
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (0 == lub) ub[i] = lb[i];
    }
    }
  }
  return 1;
}

void BCP::CompileCuts() {
  int i;
  CutList::iterator ic;
  CutRefSet::iterator ipc;
  CutContainer::iterator icc;
  Node *pnode;
// GATHER all cuts (ptrs) from above nodes to cutpool:
  cutpool.clear();
  pnode = &*theNode;
  for (;pnode;pnode = pnode->parent) {
    /// Inserting regular cuts:
    for_each_in(pnode->cutcont,icc,)
      cutpool.insert((LPCut*)&*icc);
    /// Inserting branching cuts:
    /// because we need pooling for marking
    log_n_(2.5," cuts1: "<<cutpool.size()<<'+'<<pnode->cutcont1.size());
    for_each_in(pnode->cutcont1,ic1,Node::CutIter1) {
      if (OUTP_LEV__ >= 2.5)
        (*ic1)->Print();
      cutpool.insert((LPCut*)(*ic1));
    }
  }
// + cntDel
  for_each_in(cutpool,ipc,)
    (*ipc)->cntDel0 = (*ipc)->cntDel = cntDel0;
// CALC COEFS OF ALL CUTS FOR ALL COLS (needed for viol)
  for_each_in(cutpool,ipc,)
    for (i=0;i<colpool.size();++i)
      (*ipc)->Calc__(GetMainCol(i),i);
// UPDATING THE CUT SET OF THE FORMULATION:
  // -- last node must not have been deleted
  for_each_in(cutpool,ipc,) (*ipc)->Mark(0); // to add
  for_each_in(cuts,ic,) (*ic)->Mark(1); // old
  for_each_in(theNode->cuts,ic,)
    if (1==(*ic)->Mark()) (*ic)->Mark(2); // already here
  /* Del superfl., not changing marks inside:*/
  for (i=(int)cuts.size()-1;i>=0;--i)
    if (1==cuts[i]->Mark())
      DelCut(i); // +slack; not from newNode->cuts!
      // maybe not del them from last node 'cause needed?
}

void BCP::UpdateVarBounds() {
  int i;
// ACTUAL UPDATING BRANCHES:
	for (i=0;i<lp->NCols()         // WRONG: min(lp->NCols(), int(lb.size()))
		;++i) // now the actual cols in f.
  if (not IsSlackCol(i)) {
    int j = GetColIndex(i);
    assert2(lb[j]<=ub[j]);
    if (lb[j]==ub[j])
    { dbg_outn_(4," fx"<<j); } // do it BY A VECTOR COMMAND?
    lp->ChangeVarBnds(i,lb[j],ub[j]<INT_MAX?ub[j]:1e100);
  }
}

void BCP::UpdateCuts() {
  CutList::iterator ic;
  CutRefSet::iterator ipc;
  CutContainer::iterator icc;
  /* NOW ADD new: (MARKING not lost???) */
  for_each_in(theNode->cuts,ic,)
    if (0==(*ic)->Mark())
      AddCut(*ic); // +slack; not from newNode->cuts!
  assert(!nBranch or (cuts.size() == nHP and cuts.size() == depth))
  if (!nBranch && cuts.size() > MaxAllCuts()
      && nTooLongColGen) {
    DeleteSomeCuts(2);
    assert(cuts.size() <= MaxAllCuts());
  }

  if (fCheckCuts)
    CheckCuts();

//  RecalcCutsRHS(); // -> to InitNewNode() ?
}

void BCP::RestoreBasis() {
  int i;
  CutList::iterator ic;

  Node * pnode = theNode->iBas.empty() ? theNode->parent
    : &*theNode; // choosing parent's basis if no own data
  assert(pnode->iBas.size());

  // BASIS STATES OF CUT SLACKS:
  for_each_in(cuts,ic,) (*ic)->Mark(0);
  for (i=0;i<pnode->slkBas.size();++i)
    pnode->slkBas[i]->Mark(1); // NOW OF MAIN VARS:
  for (i=0;i<lb.size();++i) lb[i] = (int)LP::atLower;
  for (i=0;i<pnode->iBas.size();++i)
    lb[pnode->iBas[i]] = (int)LP::basic;
  for (i=0;i<pnode->iUpper.size();++i)
    lb[pnode->iUpper[i]] = (int)LP::atUpper;
  cstat.resize(cols.size());
  for (i=0;i<cstat.size();++i)
    if (cols[i].slackCut)
      if (1==cols[i].slackCut->Mark()) cstat[i] = LP::basic;
      else cstat[i] = LP::atLower;
    else cstat[i] = (LP::ColStatus)lb[GetColIndex(i)];
  lp->LoadBasis(cstat);

  if (pnode == &*theNode
    // || (pnode != theNode && pnode->fLPOpt) no, the 2nd child may need
    ) {
    Reconstruct(pnode->iBas);
    Reconstruct(pnode->iUpper);
    Reconstruct(pnode->slkBas);
    Reconstruct(pnode->cuts); // + clear fathomed subtrees
  }
}

/// GLV only approximately:
void BCP::RecalcGLB() {
// BUT glb cannot change unless reopt tree
  SolutionTree::ISrch is1 = nodes.GetBeginSearch();
  glb_others = glv_others = gub;
  // looking for other nodes in the BFSDFS order
  while (nodes.GetEndSearch() != is1 &&
    ( *is1== theNode // normally yes at first
      || 0 < (*is1)->nDesc // has children
      || Node::fathomed == (*is1)->state ))
     nodes.Inc(is1);
  if (nodes.GetEndSearch() != is1) {
    glv_others = (*is1)->llv; // must not be minimal
    glb_others = (*is1)->llb;
  }
  double glb_old = glb;
  UpdateGlobalLPV();
  if (glb < glb_old) // But it cannot change unless reopt tree
    log__(" GLB got smaller !! ");
}


}                                                // SS_END_NAMESPACE__
