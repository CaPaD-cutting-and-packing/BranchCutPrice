// FILE: bb.cpp

#include "stdafx.h"
#include "bb.h"
#include "lasthdr.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

// HOW: take care in B2 that you can start from the 
// beginning of the sparse list (pk=-1)

#undef dbgpr
#define dbgpr(n,e)

void BB::Optimize() {
  assert((m>0) and (L>=0));
  /// VARIABLES::
/// A lower bound, initally max(zMin,zLowerInitial)+eps,
/// thereafter (bestSolution + eps):
  double lb=FMax(zMin, zLowerInitial) + eps; // !
  double sad=0; // the accumulated o.f. value
  PieceContainer pieces(m); // = pieces__;
  assert(forbidden__);
  Column colTmp, colTmpBetter;
//  set<Vector<IX> > forbidden;
  size dL = L;
  int
    k=-1,
    pk=-1,
    x;
  Vector<IX> ix(m);
  Vector<double> zfBuf(m+1);
  double * zf = &zfBuf[1]; // for zf[-1]
    zf[-1] = 0;
  Vector<size> l2(m);
  Vector<int> bb2(m);
  found=false;
  Piece *pc=NULL;
  assert(colsetRes && colBest);
//  colsetRes->cs.clear();
//  cycle<Column> colcclRes(1);

// Reading test data:
//  ifstream is("hard28_hardcg.txt");
//  {for (int i=0;i<m;++i) is >> pieces__[i].d;}
  // For testing:
  // colRes.ix.resize(0); // ???

// INIT:
  assert(pieces__.size()==m);
// Copy input:
  int i,j;
  for (i=0;i<m;i++) {
    pieces__[i].i0= i;
    pieces__[i].w = pieces__[i].l; // sorting in 2D
  }
  PieceContainer pc__ = pieces__;
  sort(pc__.begin(),pc__.end(),greater<Piece>());
// DOMINANCE:
  size dh;
// l2 - the undominated size:
// the max. size left by dominating items
  for (i=0;i<m; ++i)  bb2[i]=pc__[i].b; // copy demands
  for (i=0;i<m; ++i)  l2[i]=L;
  for (i=1; i<m; ++i) assert(pc__[i-1].l >= pc__[i].l);
// Item Dominance (items sorted non-incr.! !!!):
// - only if no forbidden cols -it is also possible but ??
  bool fNoForbidden = (forbidden__->empty() or fHeur);
  if (fNoForbidden)
  for (i=m-1; i>=0; --i)
  if(pc__[i].l <= l2[i]) {
   bb2[i]=IMin(bb2[i],l2[i]/pc__[i].l);
   dh=pc__[i].l*bb2[i];
   for (j=i-1; j>=0; --j)
     if ((pc__[i].d > pc__[j].d-eps) 
       && (l2[j] > 0))
       l2[j] -= dh;
  } // Compressing references:
  j=0;
  for (i=0;i<m;++i)
    if ((l2[i] >= pc__[i].l)
      && ((pc__[i].d>eps) || not fNoForbidden)
      && (bb2[i]>0)) {
        pieces[j]=pc__[i];
        pieces[j].b = bb2[i];
        ++j;
    }
  int ms=j;
  pieces.resize(ms);
  const int ms_1 = ms - 1;

// SORT:
  for (i=0;i<ms;++i) {
    pieces[i].q = pieces[i].d / pieces[i].l;
    pieces[i].w = pieces[i].q;
/*    if (fRandomizeWeights)
      RandomizePieceWeightRatio(&pieces[i]);*/
  }
  sort(pieces.begin(),pieces.end(),greater<Piece>());

// Init service data:
  size lma=SIZE_MAX__;
  double qma = - INFINITY__;
  double lsa=0;
  double dsa=0;
  for (i = pieces.size() - 1;i>=0;i--) {
    Piece * pc = &pieces[i];
    pc->lma=lma;
    if (pc->l < lma)
      lma=pc->l;
    if (pc->q > qma)
      qma=pc->q;
    pc->qma=qma;

  }
  if (OUTP_LEV__ >= 4) {
    dbg_outn_(4,"\nBB. Pieces order: ");
    for (i=0;i<ms;++i) dbg_outn_(4,pieces[i].i0<<' ');
    dbg_outn_(5,"lower bounds: "<<zMin<<' '<<zLowerInitial<<' '<<lb);
    dbg_outn(4,"");
  }
  if (OUTP_LEV__ >= 4) { // so is faster
    log_ln("BB: ms="<<ms<<" L="<<L);
    for (i=0;i<ms;++i)
      log__(pieces[i].l<<':'<<pieces[i].d
      <<':'<<pieces[i].b<<' ');
  }
  dbgcout(4,'.'<<flush);
  double ub;
  double nSteps = 0;
  fETerm = 0;

  // CHECK THE INITIAL SOLUTION, NEEDED FOR CP22.
  // COMPARE lb TO ITSELF ? NO
  // ALSO FORBIDDEN!!
  if (sad > lb) {  // <=> lb<0
      // SAVING _initial_ SOLUTION TO TEMP:
      colTmp = *colBest;
      colTmp.SortWithMerging1st(); // if pre-set piece in the beginning...
     if (fHeur or NULL == forbidden__->FindVirtual(colTmp))
     {
       colTmpBetter = colTmp;
       dbg_outn_(4," BETTER INITIAL COL (z="<<sad<<"): ");
       for (i=0;i<=pk;++i)
         dbg_outn_(4,pieces[ix[i].i].i0<<':'<<ix[i].x<<' ');
       dbg_outn(4,"");
//      colcclRes.push_back(*colBest);
      z=sad;  // here always, also w/o dominance!
      lb=z+eps;         dbgpr(5," **** SAVING **\n");
      found=true;       dbg_outn_(3,'^'<<flush);
     } else
     {
       dbg_outn_(3,'f');
       dbg_outn_(4," THE INITIAL COL IS FB (z="<<sad<<"): ");
       for (i=0;i<=pk;++i)
         dbg_outn_(4,pieces[ix[i].i].i0<<':'<<ix[i].x<<' ');
       dbg_outn(4,"");
     }
  }
  if (!ms) { // Only after checking the initial solution (CP22, PMP)
    dbg_outn(4,"No pieces for the knapsack problem w/o cuts");
    return;
  }

StepS:                //:: The adding step
  ++ k;             dbgpr(5,"\nS: k="<<k<<" lb="<<lb);
  pc=&pieces[k];  dbgpr(5," bnd="<<sad+pc->qma*dL);
  ub = (pc->qma)*dL; if (ub<0) ub=0;
  if (sad + ub > lb) { // Check ub(+eps)
    if (dL < pc->l)        // no space for current piece
      goto StepS;      dbgpr(5," li="<<pc->l<<" dL="<<dL);
    x = dL/pc->l;      dbgpr(5," x="<<x);
    if (x>pc->b) {         // PROPER pattern
      x = pc->b;       dbgpr(5," x="<<x); }
    ++ pk;             dbgpr(5," pk="<<pk);
    ix[pk].set(k,x);
    zf[pk] = (sad += pc->d*x);    dbgpr(5," sad="<<sad);
    dL -= pc->l*x;  dbgpr(5," dL="<<dL<<" lma="<<pc->lma);
    if (sad > lb) {
CheckBetter:
      // SAVING SOLUTION TO TEMP:
      colTmp = *colBest;
      colTmp.id.reserve(colBest->id.size()+pk+1);
      // Restoring the original numeration (before sorting):
      for (i=0;i<=pk;++i)
        colTmp.PushID(pieces[ix[i].i].i0,ix[i].x);
      colTmp.SortWithMerging1st(); // if pre-set piece in the beginning...
     if (fHeur or NULL == forbidden__->FindVirtual(colTmp))
     {
       colTmpBetter = colTmp;
       dbg_outn_(4," BETTER COL (z="<<sad<<"): ");
       for (i=0;i<=pk;++i)
         dbg_outn_(4,pieces[ix[i].i].i0<<':'<<ix[i].x<<' ');
       dbg_outn(4,"");
      z=sad;  // here always, also w/o dominance!
      lb=z+eps;         dbgpr(5," **** SAVING **\n");
      found=true;       dbg_outn_(3,'^'<<flush);
     } else
     {
       dbg_outn_(3,'f');
       dbg_outn_(4," THE COL IS FB (z="<<sad<<"): ");
       for (i=0;i<=pk;++i)
         dbg_outn_(4,pieces[ix[i].i].i0<<':'<<ix[i].x<<' ');
       dbg_outn(4,"");
     }
    }
    if (dL >= pieces[k].lma) goto StepS;
    goto StepB2;
  }
StepB2:                 dbgpr(5,"\n  B2:");
    if (0 > pk) goto Exit; // when coming directly from..
    k = ix[pk].i;       dbgpr(5," k="<<k<<" lma="<<pieces[k].lma);
    if (0 == -- ix[pk].x) {
      sad = zf[--pk];   dbgpr(5," pk="<<pk);
    }
    else
    {zf[pk] = (sad-=pieces[k].d); dbgpr(5," x[pk]="<<ix[pk].x);}
    dL += pieces[k].l; dbgpr(5," sad="<<sad<<" dL="<<dL);

    if (0 == fmod( (++nSteps), 1024)) {
     if (fEvenIfNotFound or found) {
      if (nSteps > nStepsTooMuch) { // even if not found ?
        dbg_outn_(4," et");
        fETerm = 1;
        goto Exit;
      }
     }
    }
    if (dL >= pieces[k].lma) goto StepS; // to get the best
    if (sad > lb) // still
      goto CheckBetter;
    goto StepB2;

Exit:                             dbgpr(4,flush);
  *colBest = colTmpBetter;
//  colBest->Sort();
  colsetRes->cs.clear();
  if (found) colsetRes->cs.insert(*colBest);
} //____________________________________________________
/** Remarks.
+ Several calls with almost the same data for M>1 and 2d
 => Flexible! Hierarchy (+ of solver classes)
+ eg. level cuts
+ if sorting not acc. to bound then next pieces not nec.
 the best!
+ asserts that data are prepared?
+ Perturbations: only actual?

+ Use lower bounds to exclude useless piecess
 - Use DSL to omit PI's calcul
+ Theoretically, PIs don't grow. Practice?
+ Some inline, all static.
+ What when last types fully in? Only for monotone c().
 - But this can be softened by the bound.
  */
void BB::SignalBetterSolution() {
  if (FMin(outputLevel,opt::GlobalOutputLevel())>3) {
    // Printing some info
    return;
  }
  if (FMin(outputLevel,opt::GlobalOutputLevel())>2)
    PRINT__('^');
} //____________________________________________________
// Precondition: Init() must have been called
void BB::Run() {
  // ++:: For my rnd proc, not full pts can be better
  // modify it: consider 0-multipliers.
  // cleanup: internals
  assert(pieces__.size()==m);
  // assert(addi__.size()==mc);
  Print();
  Optimize();
  dbgcout(3.5,'.');
} //____________________________________________________
void BB::RandomizePieceWeightRatio
(BB::Piece *pc) {
} //____________________________________________________
void BB::Reallocate(int m) {
// INPUT:
  resize(pieces__,m);
//  resize(addi__,mc);
//  resize(dAddi__,mc);
// VARS:
  //resize(colRes.ix,m);
} //____________________________________________________

// Print some solution-relevant data.
// Control: outputLevel
void BB::Print() {
  if (FMin(outputLevel,opt::GlobalOutputLevel())<3)
    return;
  /*g__("\nBB. Weights: ");
  for_each_in(pieces__,ipc,PieceContainer::iterator)
    log__(ipc->l<<':'<<ipc->b<<'='<<ipc->d<<' ');
  log_ln("");*/
//  if (found) log_ln("  Found: z="<<z);
} //____________________________________________________

double BB::outputLevel;
double BB::nStepsTooMuch;
bool BB::fEvenIfNotFound;
//bool BB::fRandomizeWeights;

// _____________________________________________________

opt::OptContainer BB::Options() {
  opt::OptContainer oc;
  oc
/*    << opt::MakeOpt(&fRandomizeWeights, true,
      "fRandWghts",
      "Bool: randomizing the weights "
      "of pieces for sorting")  // + aRndStdDev*/
    << opt::MakeOpt(&nStepsTooMuch, 1e+6,
      "nStepsTooMuch",
      "So many backtrack steps => termination")
    << opt::MakeOpt(&fEvenIfNotFound, 0,
      "fEvenIfNotFound",
      "nStepsTooMuch backtrack steps => termination even if no good sol. found")
      
    << opt::MakeOpt(&outputLevel, 4,
      "outputLevel", "");
  return oc;
} //____________________________________________________
opt::OptSection BB::opt
  ("BB", "Br&Bound for Col Gen (w/o cuts)",
  BB::Options(), opt::SolverCfg(), 4500);

void BB::DoStatistics(int wh,ostream &os) {
} //____________________________________________________

}                                                // SS_END_NAMESPACE__
