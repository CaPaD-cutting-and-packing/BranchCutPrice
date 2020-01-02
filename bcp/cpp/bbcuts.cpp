// FILE: bbcuts.cpp
// Branch-and-bound to solve column generation
// with general-purpose (superadditive) cutting planes
// AUTHOR: Gleb Belov  <belov@math.tu-dresden.de>
// Arrays from 0.




// How data storage ?
// Print problems/results optionally


// EVERY DETAIL IS IMPORTANT




#include "stdafx.h"
#include "bbcuts.h"
#include "lasthdr.h"




namespace bgn {                     // SS_BEGIN_NAMESPACE__




#undef dbgpr
#define dbgpr(n,e)  // if (OUTP_LEV__ >=n) cout<<e;




//#define CHECK // No bound usage
//#define CHECKCALC // Check objective


int cnt=0; // for debug


void BBCuts::Optimize() {
  cnt ++;
  dbg_outn_(5,"BBCuts.cnt="<<cnt<<' ');
  InitEnumeration();
  ConsiderCalcObj1(); // for the ini. solution, see CP22
  if (ms <= 0) return; // if no more pieces
StepAdd:             //:: The adding step
  k++;             dbgpr(5,"\nS: k="<<k<<" lb="<<lb);  
                   dbgpr(5," bnd="<<ub1+UB2(k));
//  if (ms == k) goto StepB;
//#ifndef CHECK
  if (ub1+UB2(k) < lb)    // Check upper bound(+eps)
    goto StepB;           // for the whole length L
//#endif
             dbgpr(5," li="<<pieces[k].l<<" dL="<<dL);
  if (dL < pieces[k].l)   // no space for current piece
    goto StepAdd;         // -- take next
  AddMaxItemsOfThePiece();
  ConsiderCalcObj1();
  if (dL >= pieces[k].lma) // if space enough for further
    goto StepAdd;
////////////////////////////////////////////////////////
//  CalcObjective();      // This is done after every
                          //piece type.
                        // Though don't use it for bound
////////////////////////////////////////////////////////
StepB:               //:: The backtracking step
  if (EnumerationOver())
    return;
  RemoveLastItem();       // What when full patts only ?
                          // Some addi removals ? (opt)
  ConsiderCalcObj2();
  if (dL < pieces[k].lma) // if no space for any further
    goto StepB;




  if (0 == fmod( (++nSteps), 2048))
    if (CheckForEarlyTermination())
      return;
  goto StepAdd;           // and so on.
} //____________________________________________________
/** Remarks.
+ Several calls with almost the same data for M>1 and 2d
 => Flexible! Hierarchy (+ of solver classes)
+ eg. level cuts
+ if sorting not acc. to bound then next pieces not nec.
 the best!
+ asserts that data are prepared?
+ Perturbations: only actual?
+ Use DSU to exclude useless piecess
 - Use DSL to omit PI's calcul
+ Theoretically, PIs don't grow. Practice?
+ Some inline, all static.
+ What when last types fully in? Only for monotone c().
 - But this can be softened by the bound.
  */
void BBCuts::InitEnumeration() {
  // asserts that prepared
  dL = L;
  pk=k = -1;
  nSteps=0;
  fETerm = found = false;
//  trgETerm.AddTry();
  if (fUpdateBnd)
    lb=FMax(zMin, zLowerInitial) + eps; // !
  z = - INFINITY__; // !!! to indicate. But use (found).
  RecalcSums();
} //____________________________________________________
void BBCuts::RecalcSums() {
  // Because of rounding errors we recalc accumulators
  // after a while. + In the initialization.
  // Using pk as the current end of the compact list.
  ub1 = apprError + bndConstTerm; // + apprError !!!
  sad = zfConstTerm;
  // Set constant terms for cut sums accumulators:
  cuts.AssignConstTerms(); // must be calc. somewhere
  int pk1; // Index for piece type in the compact solution
  for (pk1=0;pk1<=pk;pk1++) {
    int k=ix[pk1].i;
    int k0=pieces[k].i0;
    cuts.AddToSums(k0,ix[pk1].x);
    ub1 += pieces[k].dAppr * ix[pk1].x;
    sad += pieces[k].d     * ix[pk1].x;
  }
} //____________________________________________________
void BBCuts::AddMaxItemsOfThePiece() {
  Piece * pc = &pieces[k];
  int x = int(dL / pc->l);      dbgpr(5," x="<<x);
  if (x > pc->b)
    x = pc->b;           // Proper relaxation
        dbgpr(5," x="<<x);
  ++ pk;                   dbgpr(5," pk="<<pk);
  ix[pk].set(k,x);
  dL -= pc->l * x;
  ub1 += pc->dAppr *x;        // the bound
  sad += pc->d * x;           // the linear portion of c()
                           dbgpr(5," sad="<<sad);
                           dbgpr(5," ub1="<<ub1);
               dbgpr(5," dL="<<dL<<" lma="<<pc->lma);
  cuts.AddToSums(pc->i0,x);
} //____________________________________________________
void BBCuts::RemoveLastItem() {
  k = ix[pk].i;        dbgpr(5,"\n      B: k="<<k);
  Piece * pc = &(pieces[k]);
  sad -= pc->d;
  ub1 -= pc->dAppr;
  dL  += pc->l;
                           dbgpr(5," sad="<<sad);
                           dbgpr(5," ub1="<<ub1);
               dbgpr(5," dL="<<dL<<" lma="<<pc->lma);
  cuts.Sub1FromSums(pc->i0);
  if (0==(--ix[pk].x)) {
    -- pk;                    // Taking the prev. pc. type
    fPrevPieceType=true;
  } else {
    fPrevPieceType=false;
  }   dbgpr(5," pk="<<pk);
} //____________________________________________________
void BBCuts::CalcObjective() {
  int i;
//#ifndef CHECK
  if (ub1 < lb)      // Checking the bound
    return;          // for the occupied length
//#endif
  double sDPi = cuts.CalcCoefsUsingIntermSums();
  double zCurrent = sad+sDPi;
     dbgpr(5," zC="<<zCurrent);
//  if (zCurrent<lb) return; // it was calc. with eps
#ifdef CHECKCALC
     Column cc = *colBest; // keeping its obj & some elems
//     int i;
     for (i=0;i<=pk;++i)
       cc.PushID(pieces[ix[i].i].i0,ix[i].x);
     double zz=0;
     for_each_in(cc.id,iid,Column::iterator)
       zz += d__[iid->i] * iid->d;
     double zzz = zz;
     for_each_in(cuts.invCuts,ic,CutList::iterator)
       (*ic)->ClearNonRec();
     for_each_in(cuts.dep,idep,
       CutSet::DepContainer::iterator)
       zzz += idep->u * idep->c->Calc__(&cc);
     dbgpr(5," zC="<<zzz);
     assertm(fabs(zzz-zCurrent) < eps,
       zzz << " -actual cg zf diff from fast calc-d: "
       << zCurrent);
     if (27550 == cnt or 27551== cnt or 27733==cnt or 28969==cnt) {
       dbg_outn_(4," COL WITH BETTER BOUND u1b="
         <<ub1<<" (z="<<zCurrent
	 <<", zRecalc="<<zzz<<", lb="<<lb<<"): ");
       for (i=0;i<=pk;++i)
         dbg_outn_(4,pieces[ix[i].i].i0<<':'<<ix[i].x<<' ');
       dbg_outn(4,"");
     }
#endif
//     dbg_outn_(5," zCurrent="<<zCurrent);
  if (zCurrent > lb) { // A better total value
         dbgpr(5," **** SAVING **\n");
//    dbg_outn_(4,'^'); // will be shown 'f'
    SaveSolutionToTemp();
    if (NULL == forbidden__->FindVirtual(colTmp)) {
       dbg_outn_(4," BETTER COL (z="<<zCurrent<<"): ");
       for (i=0;i<=pk;++i)
         dbg_outn_(4,pieces[ix[i].i].i0<<':'<<ix[i].x<<' ');
       dbg_outn(4,"");
      if (fUpdateBnd) {
        found=true;
        z = zCurrent;
        lb = z+eps;
      } else if (zCurrent > z) {
        z = zCurrent;
      }
      SaveBetterSolutionFromTemp();
      SignalBetterSolution(); // no nsteps=0;
    } else {
       dbg_outn_(3,'f');
       dbg_outn_(4," THE COL IS FB (z="<<zCurrent<<"): ");
       for (i=0;i<=pk;++i)
         dbg_outn_(4,pieces[ix[i].i].i0<<':'<<ix[i].x<<' ');
       dbg_outn(4,"");
//       goto StepB3; // to check immed. after del
     }


  }
} //____________________________________________________
void BBCuts::SignalBetterSolution() {
  if (FMin(outputLevel,opt::GlobalOutputLevel())>3)
    PRINT__('^');
  if (FMin(outputLevel,opt::GlobalOutputLevel())>4) {
    // Printing some info
    return;
  }
} //____________________________________________________
bool BBCuts::CheckForEarlyTermination() {
// Increase nstepsmax, nsavesmax steadily till INT_MAX/2
// Doubling (?) interval ~m
  if (nSteps > nStepsMin) {
    if (found || foundInit) {
      fETerm = true;// Signal early term. above
      dbg_outn_(4,'x');
      return true;
    }
    else if (nSteps > nStepsTooMuch)
      throw int(1);
  }
  if (CheckForUserBreak__()) // signal Ctrl-C or timer
    return true; // if still nece. (not thrown)
  RecalcSums();
  return false;
} //____________________________________________________
// Precondition: Init() must have been called
void BBCuts::Run() {
  // Auto: if no good found then with NotOnlyFull=true
  // ++:: For my rnd proc, not full pts can be better
  // modify it: consider 0-multipliers.
  // cleanup: internals
  Print();
  assert(pieces__.size()==m);
  fETerm = found=false;
  z = - INFINITY__; // !!! to indicate. But use (found).
//  if (!ms) { dbg_outn_(5," bbc:m=0"); return; }
  // EMPTY SOLUTION IS HANDLED ALSO. FOR CP22
/*  if (fTryAlsoNotOnlyFullPatterns)
    fConsiderNotOnlyFullPatterns = false;*/
  fUpdateBnd = true;
  Optimize();
  if (!found)
    if (fCheckBnd) {
      fUpdateBnd = false;
      lb = -1e100;
      Optimize();
      if (z > FMax(zMin, zLowerInitial) + eps) {
        if (outputLevel<5) 
          outputLevel=5;
        dbg_outn(5,"RUNNING BOTH ONCE MORE:");
        Print();
        fUpdateBnd = true;
        Optimize();
        fUpdateBnd = false;
        lb = -1e100;
        Optimize();
        assertm(false, "BBCuts: only run w/o bound found a pattern!!!");
      }
    }
/*fdef BBCuts__NotOnlyFullPatterns
  if (fTryAlsoNotOnlyFullPatterns) {
    if (!found) { // No proper solution found
      trg2ndRunFinds.AddTry();
      fConsiderNotOnlyFullPatterns = true;
      Optimize();
      if (found) {
        trg2ndRunFinds.AddPositive();
        cout << "2nd run finds!" << endl;
        its_error("...","","");
      }
    }
  }
#endif*/
//  dbg_outn_(3,'.');
  if (found) CopyResult();
  dbgcout(3.5,'.');
} //____________________________________________________
// Precondition: All input vars set!
void BBCuts::Init() {
  assert(pieces__.size()==m);
  assert(mc>=0);
  assert(d__.size()>=m+mc);
  assert(forbidden__);
//  assert(addi__.size()==mc);
  CopyInput();
  ReduceNCuts();
  CalcBounds();
  CalcConstTerms();  // =UpdateAddiConstr(). But 2d!
     // for the level cut: before CalcApprError()! ???
  SortPieces();
  InitServiceData();
} //____________________________________________________
void BBCuts::CopyInput() {
  // Copying info from input:
  int i,j;
  for (j=i=0;i<m;i++)
  if (pieces__[i].b > 0 && pieces__[i].l <= L) {
    pieces[j].l = pieces__[i].l;
    pieces[j].b = pieces__[i].b;
    pieces[j].i0= i;
    pieces[j].d = pieces__[i].d;
    ++ j;
  }
  ms = j;
  cuts = cuts__;
} //____________________________________________________
void BBCuts::ReduceNCuts() {
  // Leaving only somehow needed cuts: done at constr.
  // cuts.ReduceWithZeroWeights(); // Now mu invalid
} //____________________________________________________
void BBCuts::CalcBounds() {
  //////////////////////////////////////////////////////
  //// Calc the linear approx of multipliers
  // assert (muS>0) ? Hpw reduced u-containers
  //////////////////////////////////////////////////////
  int i;
  dAppr.clear(); dAppr.resize(m+mc);
  for (i=0;i<ms;++i) {
    double sDUAppr = cuts.CalcApprCoefs(pieces[i].i0);
      // and mult them with d's at once to return
    dAppr[pieces[i].i0]
      = pieces[i].dAppr = pieces[i].d + sDUAppr;
  } // for each pieces
  // Also for already added pieces:
  for_each_in(colBest->id,iid,Column::iterator)
    if (!dAppr[iid->i]) {
      double sDUAppr
        = cuts.CalcApprCoefs(iid->i); // !!!
      dAppr[iid->i] = d__[iid->i] + sDUAppr;
    }
  // + for each addi constr (2D):
  for (i=m;i<m+mc;i++) {
    double sDUAppr
      = cuts.CalcApprCoefs(i); // !!!
    dAppr[i] = d__[i] + sDUAppr;
  }
  //////////////////////////////////////////////////////
  //// Calc the abs error of the appr when using it as
  //// an upper bound
  //////////////////////////////////////////////////////
  apprError = cuts.CalcApprError();
/*  double aErr = cuts.GetAprrErrorAlt();
  assert(fabs(aErr - ApprError) < 1e-6);*/
} //____________________________________________________
void BBCuts::CalcConstTerms() {
  // How do addi const influence the bound ?
// FOR 2D: NEED THE WHOLE d-Vector, then vec. prod.
  zfConstTerm  = colBest->VecProd(d__);
  bndConstTerm = colBest->VecProd(dAppr);
  cuts.CalcConstTerms(colBest);
} //____________________________________________________
void BBCuts::SortPieces() {
  int i;
  for (i=0;i<ms;++i) {
    pieces[i].q = pieces[i].dAppr / double(pieces[i].l);
    pieces[i].w = pieces[i].q;
    if (fRandomizeWeights)
      RandomizePieceWeightRatio(&pieces[i]);
  }
//  list<Piece> pl;
//  pl.assign(pieces.begin(),pieces.end());




  sort(pieces.begin(),pieces.begin()+ms,greater<Piece>());
//  pieces.assign(pl);
}
void BBCuts::InitServiceData() {
  size lma=SIZE_MAX__;
  double qma = - INFINITY__;
  double dama = qma;
  double lsa=0;
  double dSApprA=0;




  int i;
  for (i=ms-1;i>=0;i--) {
    Piece * pc = &pieces[i];
    pc->lma=lma;
    if (pc->l < lma)
      lma=pc->l;
//    pc->lma1=lma;
    if (pc->q > qma)
      qma=pc->q;
    pc->qma=qma;
    if (pc->dAppr > dama)
      dama=pc->dAppr;
    pc->dApprMA=dama;




// IMPORTANT: add before
    lsa += double(pc->l) * pc->b;
    dSApprA += FMax(0, pc->dAppr * pc->b);
    if (lsa > SIZE_MAX__)
      pc->lsa=SIZE_MAX__;
    else {
      pc->lsa=(size) lsa;
      pc->dSApprA=dSApprA;
    }
  }
  if (OUTP_LEV__ >= 4) {
    dbg_outn_(4,"\nBBCuts. Pieces order: ");
    for (i=0;i<ms;++i) dbg_outn_(4,pieces[i].i0<<' ');
    dbg_outn(4,"");
  }
/*
  /// FORBIDDEN
  forbidden.clear();
  i_vec ai1(m); // eliminated pieces stay i1=0
  fill(ai1.begin(),ai1.end(),INT_MAX); // all eliminated not!
  for (i=0;i<ms;++i) ai1[pieces[i].i0] = i;
  pair<set<Vector<IX> >::iterator,bool> res;


  dbg_outn_(4,"\nBB: FORBIDDEN COLS (all "<<forbidden__.size()
    <<"). Pieces order: ");
  for (i=0;i<ms;++i) dbg_outn_(4,pieces[i].i0<<' ');
  dbg_outn(4,"");


  for_each_in (forbidden__,iff,
      list<Vector<IX> >::const_iterator) {
    Vector<IX> vix = (*iff);
    for (int i=0;i<vix.size();++i) {
      if (vix[i].i >= m) {
        mylog << "forbidden vecs:\n";
        for_each_in (forbidden__,iff,
          list<Vector<IX> >::const_iterator) {
            for (int i=0;i<iff->size();++i)
              mylog << ' '<<(*iff)[i].i<<':'<<(*iff)[i].x;
            mylog << '\n';
          }
        assertm(false," Memory fleak");
      }
      vix[i].i = ai1[vix[i].i];
      if (INT_MAX == vix[i].i) goto NextVec;
    }
    for (int i=0;i<vix.size();++i)
      dbg_outn_(4,pieces[vix[i].i].i0<<':'<<vix[i].x<<' ');
    dbg_outn(4,"");
    sort(vix.begin(), vix.end());
    res = forbidden.insert(vix);
//    assert(res.second); // no: when ms small...
NextVec: ;
  }
  */
} //____________________________________________________
void BBCuts::RandomizePieceWeightRatio
(BBCuts::Piece *pc) {
} //____________________________________________________
void BBCuts::SaveSolutionToTemp() { // SORTING!!
  // TO SEARCH FORBIDDEN
  colTmp = *colBest;
  colTmp.id.reserve(colBest->id.size()+pk+1);
  // Restoring the original numeration (before sorting):
  for (int i=0;i<=pk;++i)
    colTmp.PushID(pieces[ix[i].i].i0,ix[i].x);
  colTmp.SortWithMerging1st(); // if pre-set piece in the beginning...
} //____________________________________________________
void BBCuts::SaveBetterSolutionFromTemp() {
  colTmpBetter = colTmp; // colBest stays untouched
// NONE
} //____________________________________________________
void BBCuts::CopyResult() { // Only once
  *colBest = colTmpBetter;
//  colBest->Sort(); // !!!!! before
  colsetRes->cs.clear();
  colsetRes->cs.insert(*colBest);
} //____________________________________________________
void BBCuts::Reallocate(int m/*,int mc*/) {
// INPUT:
  resize(pieces__,m);
//  resize(addi__,mc);
//  resize(dAddi__,mc);
// VARS:
//  resize(colRes.ix,m);
  resize(colTmp.id,m);
  resize(ix,m);
  resize(pieces,m);
//  resize(dApprAddi,mc);
} //____________________________________________________




// INPUT
int BBCuts::m/*,BBCuts::mc ,mu*/;
int BBCuts::mc;
size BBCuts::L; // Stock length
BBCuts::PieceContainer BBCuts::pieces__;
// Each knows l,b,d
/// The additional constraints part of the column, fixed:
//valarray<addi_float> BBCuts::addi__;
//valarray<double> BBCuts::dAddi__; // their weights
/// Cutting planes:
CutSet BBCuts::cuts__;    // with weights incaps.
Pool<Column> * BBCuts::forbidden__;
d_vec BBCuts::d__;
/// Initial, eg incumbent, lower bound:
double BBCuts::zLowerInitial;
/// A priori lower bound, eg for primal simplex or
/// obtained from Lagrange relaxation.
/// E.g. CSP1: 1+eps
double BBCuts::zMin;
bool BBCuts::foundInit;
/// The tolerance by which a better solution must be >=:
double BBCuts::eps; // =smth*bbpi_eps
/// Early termination controls:
double BBCuts::nStepsMin=0; // =16384*8=nStepsMax0
// Use global prms nStepsMaxIncrRatio + interval
bool BBCuts::fRandomizeWeights;
// + How -- params
bool
    BBCuts::fConsiderNotOnlyFullPatterns,
    BBCuts::fTryAlsoNotOnlyFullPatterns;
double BBCuts::outputLevel;
  //////////////////////////////////////////////////////
  /// OUTPUT::
  //____________________________________________________
/// Result flag. =true e.g. if (z > zMin).
bool BBCuts::found;
/// The value of the best found solution
double BBCuts::z;
bool BBCuts::fETerm;
/// A corresponding solution in compact form
Column * BBCuts::colBest;
ColSet * BBCuts::colsetRes;
  //////////////////////////////////////////////////////
  /// STATISTICS::
// Whether the 2nd run (not only full...) better:
/*mystat::Accumulator
    BBCuts::trg2ndRunFinds("BBCuts.2ndRunFinds"),
    BBCuts::trgETerm("BBCuts.ETerm");*/




  //////////////////////////////////////////////////////
  /// VARIABLES::
  /// (now all static, no parallelization)
  /// (even the latter can use static with templates)
  /// or use separate executables
  //____________________________________________________
/// A lower bound, initally max(zMin,zLowerInitial)+eps,
/// thereafter (bestSolution + eps):
double BBCuts::lb;
/// An upper bound for the partial solution:
double BBCuts::ub1; // = dAppr*a - alfL;
double BBCuts::sad;
BBCuts::PieceContainer BBCuts::pieces;
size BBCuts::dL;
int
    BBCuts::ms,
    BBCuts::k,
    BBCuts::pk;
Vector<IX> BBCuts::ix;
double BBCuts::nSteps;
CutSet BBCuts::cuts;
//set<Vector<IX> > BBCuts::forbidden;
//valarray<double> BBCuts::dApprAddi;
Column BBCuts::colTmp, BBCuts::colTmpBetter;
double
    BBCuts::apprError,
    BBCuts::bndConstTerm,
    BBCuts::zfConstTerm;
Vector<double> BBCuts::dAppr;
bool BBCuts::fPrevPieceType;
// _____________________________________________________
double BBCuts::nStepsTooMuch=1e+7;
bool BBCuts::fCheckBnd, BBCuts::fUpdateBnd;


opt::OptContainer BBCuts::Options() {
  opt::OptContainer oc;
  oc
/*    << opt::MakeOpt(&fRandomizeWeights, true,
      "fRandWghts",
      "Bool: randomizing the approximated weights "
      "of pieces for sorting")*/
    << opt::MakeOpt(&nStepsTooMuch, 1e+8,
      "nStepsTooMuch",
      "So many steps means there are too many cuts")
    << opt::MakeOpt(&fCheckBnd, false,
      "fCheckBnd",
      "whether the procedure should be run again w/o bound in the case no sol is found (very slow)")
/*    << opt::MakeOpt(&fConsiderNotOnlyFullPatterns, false,
      "fNotOnlyFullPatt",
      "Bool: Consider not only full patterns")
    << opt::MakeOpt(&fTryAlsoNotOnlyFullPatterns, false,
      "fTryAlsoNotOnlyFull",
      "Bool: Try also with not only full patterns "
      "if failed without")*/
    << opt::MakeOpt(&outputLevel, 4,
      "outputLevel", "");
  return oc;
} //____________________________________________________
opt::OptSection BBCuts::opt
  ("BBCuts", "Br&Bound for Col Gen with cuts",
  BBCuts::Options(), opt::SolverCfg(), 5000);




// Print some solution-relevant data.
// Control: outputLevel
void BBCuts::Print() {
  if (FMin(outputLevel,opt::GlobalOutputLevel())<=5.99)
    return;
  log_ln(" BBCuts. Piece lengths:");
  PieceContainer::iterator ipc;
//  int i;
  for_each_in(pieces,ipc,)
    PrintVar(mylog,floatWidth,ipc->l);
  log_ln(" BBCuts. Piece orig. numbers:\n");
  for_each_in(pieces,ipc,)
    PrintVar(mylog,floatWidth,ipc->i0);
  log_ln(" BBCuts. Piece weights:\n");
  for_each_in(pieces,ipc,)
    PrintVar(mylog,floatWidth,ipc->d);
/*  log__("\nAddi weights: ");
  for (i=0;i<mc;i++)
    PrintVar(mylog,floatWidth,dAddi__[i]);*/
  log__("\nCut weights: ");
  for_each_in(cuts.dep,idep,
    CutSet::DepContainer::iterator) {
      mylog << ' '<<idep->c << ':';
      PrintVar(mylog,floatWidth,idep->u);
  }
  log__("apprError: " << apprError);
  log_ln(" Piece approx weights:");
  for_each_in(pieces,ipc,)
    PrintVar(mylog,floatWidth,ipc->dAppr);
/*  log__("\nAddi approx weights: ");
  for (i=0;i<mc;i++)
    PrintVar(mylog,floatWidth,dApprAddi[i]);
  log__('\n'); */
/*  if (trgETerm.GetTrigger())
    log__("EarlyTerm. ");
  if (trg2ndRunFinds.GetTrigger())
    log__("2ndRunFinds.");*/
  log__('\n');
} //____________________________________________________




}                                                // SS_END_NAMESPACE__
