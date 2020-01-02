// FILE: bb_mos.cpp

#include "stdafx.h"
#include "bb_mos.h"
#include "lasthdr.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

// HOW: take care in B2 that you can start from the 
// beginning of the sparse list (pk=-1)

#undef dbgpr
#define dbgpr(n,e)

void BBMos::Optimize() {
  assert((m>0) and (L>=0));
  /// VARIABLES::
/// A lower bound, initally max(zMin,zLowerInitial)+eps,
/// thereafter (bestSolution + eps):
  double lb=FMax(zMin, zLowerInitial) + eps; // !
  dbg_outn(5,zMin<<' '<<zLowerInitial<<' '<<lb<<' '<<eps);
  double sad=0; // the accumulated o.f. value
  PieceContainer pieces(m); // = pieces__;
//  assert(forbidden__);
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
    pieces__[i].w = (double)pieces__[i].l; // sorting in 2D
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
  bool fNoForbidden = (//forbidden__->empty() or
    fHeur);

  // better no dominance with nOpen
  // but yes without

  if (fDominance &&
    fNoForbidden)
  for (i=m-1; i>=0; --i)
  if(pc__[i].l <= l2[i]) {
   bb2[i]=Min(bb2[i],l2[i]/pc__[i].l);
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
  if (!ms) {
    dbg_outn(4,"No pieces for the knapsack problem w/o cuts");
    return;
  }
  pieces.resize(ms);
  const int ms_1 = ms - 1;

// SORT:
  for (i=0;i<ms;++i) {
    pieces[i].q = pieces[i].d / (double)pieces[i].l;
    pieces[i].w = pieces[i].q;
/*    if (fRandomizeWeights)
      RandomizePieceWeightRatio(&pieces[i]);*/
  }
  sort(pieces.begin(),pieces.end(),greater<Piece>());

// Init service data:
  size lma=SIZE_MAX__;
  double qma = - INFINITY__;
  // double lsa=0;
  //double dsa=0;
  for (i = (int)pieces.size() - 1;i>=0;i--) {
    Piece * pc = &pieces[i];
    pc->lma=lma;
    if (pc->l < lma)
      lma=pc->l;
    if (pc->q > qma)
      qma=pc->q;
    pc->qma=qma;

/*    if (lsa > SIZE_MAX__)
      pc->lsa=SIZE_MAX__;
    else {
      pc->lsa=(size) lsa;
      pc->dsa=dsa;
      lsa += double(pc->l) * pc->b;
      dsa += pc->d * pc->b;
    }*/
  }
  if (OUTP_LEV__ >= 4) {
    dbg_outn_(4,"\nBB. Pieces order: ");
    for (i=0;i<ms;++i) dbg_outn_(4,pieces[i].i0<<' ');
    dbg_outn(4,"");
  }
  /*
  /// FORBIDDEN
  forbidden.clear();
  i_vec ai1(m);
  fill(ai1.begin(),ai1.end(), INT_MAX);
  for (i=0;i<ms;++i) ai1[pieces[i].i0] = i; // others 0?
  pair<set<Vector<IX> >::iterator,bool> res;
  if (!fHeur)
  for_each_in (forbidden__,iff,
      list<Vector<IX> >::const_iterator) {
    Vector<IX> vix = *iff;
    for (int i=0;i<vix.size();++i) {
      vix[i].i = ai1[vix[i].i];
      if (INT_MAX == vix[i].i)
        goto NextVec;
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
  if (OUTP_LEV__ >= 4) { // so is faster
    log_ln("\nBB: ms="<<ms);
    for (i=0;i<ms;++i)
      log__(pieces[i].l<<':'<<pieces[i].d
      <<':'<<pieces[i].b<<' ');
  }
  dbgcout(4,'.'<<flush);
  double ub;

  // CHECK THE INITIAL SOLUTION, NEEDED FOR CP22.
  // ALSO FORBIDDEN!!
  if (sad > lb) {  // <=> lb<0
      // SAVING _initial_ SOLUTION TO TEMP:
      colTmp = *colBest;
//      colTmp.id.reserve(colBest->id.size()+pk+1);
      // Restoring the original numeration (before sorting):
//      for (i=0;i<=pk;++i)
 //       colTmp.PushID(pieces[ix[i].i].i0,ix[i].x);
      colTmp.SortWithMerging1st(); // if pre-set piece in the beginning...
     if (fHeur// or NULL == forbidden__->Find(colTmp)
       )
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
     }
  }
  // BUT NOW QUIT IF NO MORE PIECES
  if (ms <= 0) goto Exit;

  /// OPEN:
  nStillOpen = accumulate(fOpen.begin(), fOpen.end(), int(0));
  nOpenIn = 0;

StepS:                //:: The adding step
  if (++ k > ms_1)
/*    if (sad > lb) goto CheckBetter;
    else*/ goto StepB2;             dbgpr(5,"\nS: k="<<k<<" lb="<<lb);
  pc=&pieces[k];  dbgpr(5," bnd="<<sad+pc->qma*dL);
  ub = (pc->qma)*(double)dL; if (ub<0) ub=0;
  if (sad + ub > lb) { // Check ub(+eps)
    if (dL < pc->l)        // no space for current piece
      goto StepS;      dbgpr(5," li="<<pc->l<<" dL="<<dL);
    if (nStillOpen - nOpenIn + pk +1 >= nOpenMax) // open full
      if (not fOpen[pc->i0])
//        if (k<ms_1)
          goto StepS;
//        else goto StepB2; // ??? or del all last ? Prob. not
    x = int(dL/pc->l);      dbgpr(5," x="<<x);
    if (x>pc->b) {         // PROPER pattern
      x = pc->b;       dbgpr(5," x="<<x); }
    ++ pk;             dbgpr(5," pk="<<pk);
    /// OPEN:
    if (fOpen[pc->i0]) ++ nOpenIn;
//    assert(nStillOpen - nOpenIn + pk + 1 <= nOpenMax);
//    assert(nStillOpen >= nOpenIn);
    ix[pk].set(k,x);
    zf[pk] = (sad += pc->d*x);    dbgpr(5," sad="<<sad);
    dL -= pc->l*x;  dbgpr(5," dL="<<dL<<" lma="<<pc->lma);
//    if (dL >= pc->lma
      /*and nStillOpen - nOpenIn + pk < nOpenMax*/
    //) // 1 more fits
  //    goto StepS;
    // save always because of N Open:
    if (sad > lb) { CheckBetter:
      // SAVING SOLUTION TO TEMP:
      colTmp = *colBest;
      colTmp.id.reserve(colBest->id.size()+pk+1);
      // Restoring the original numeration (before sorting):
      for (i=0;i<=pk;++i)
        colTmp.PushID(pieces[ix[i].i].i0,ix[i].x);
      colTmp.SortWithMerging1st(); // if pre-set piece in the beginning...
       colTmpBetter = colTmp;
       dbg_outn_(4," BETTER COL (z="<<sad<<"): ");
       for (i=0;i<=pk;++i)
         dbg_outn_(4,pieces[ix[i].i].i0<<':'<<ix[i].x<<' ');
       dbg_outn(4,"");
//      colcclRes.push_back(*colBest);
      z=sad;  // here always, also w/o dominance!
      lb=z+eps;         dbgpr(5," **** SAVING **\n");
      found=true;       dbg_outn_(3,'^'<<flush);
    }
    if (dL >= pieces[k].lma
      /*and nStillOpen - nOpenIn + pk < nOpenMax*/)
      goto StepS;
    if (ms_1==k)
      goto StepB1; // what with forbidden?
    goto StepB2;
  }
  if (0 > pk)  goto Exit;
StepB:
  while (k-1 == ix[pk].i) {
    -- k;               dbgpr(5,"\n      B: k="<<k);
StepB1:                 dbgpr(5,"\n    B1:");
    dL += pieces[k].l*ix[pk].x; dbgpr(5," dL="<<dL);
    if (0 == pk) goto Exit;
    sad = zf[--pk];     dbgpr(5," sad="<<sad<<" pk="<<pk);
    dbg_assert(0 != ix[pk].x);  dbgpr(5," bi="<<pieces[k].b);
    if (fOpen[pieces[k].i0])  -- nOpenIn;
    if (pieces[k].b == ix[pk+1].x)
      goto StepB;
    goto StepB2;
  }
StepB2:                 dbgpr(5,"\n  B2:");
//  do {
    if (0 > pk) goto Exit; // when coming directly from..
    k = ix[pk].i;       dbgpr(5," k="<<k<<" lma="<<pieces[k].lma);
    if (0 == -- ix[pk].x) {
//      dbg_assert(0<=pk); // goto Exit;
      sad = zf[--pk];   dbgpr(5," pk="<<pk);
      /// OPEN:
      if (fOpen[pieces[k].i0])  -- nOpenIn;
    }
    else
    {zf[pk] = (sad-=pieces[k].d); dbgpr(5," x[pk]="<<ix[pk].x);}
    dL += pieces[k].l; dbgpr(5," sad="<<sad<<" dL="<<dL);
//  } while (dL < pieces[k].lma);
//  goto StepS;
    if (dL >= pieces[k].lma
      /*and nStillOpen - nOpenIn + pk + 1 < nOpenMax*/)
      goto StepS; // to get the best
    if (sad > lb) // still
      goto CheckBetter;
    goto StepB2;
/*StepB3:                 dbgpr(5,"\n  B3:"); // after forb.
//    if (0 > pk) goto Exit; // when coming directly from..
    k = ix[pk].i;       dbgpr(5," k="<<k<<" lma="<<pieces[k].lma);
    if (0 == -- ix[pk].x) {
      dbg_assert(0<=pk); // goto Exit;
      sad = zf[--pk];   dbgpr(5," pk="<<pk);
    }
    else
    {zf[pk] = (sad-=pieces[k].d); dbgpr(5," x[pk]="<<ix[pk].x);}
    dL += pieces[k].l; dbgpr(5," sad="<<sad<<" dL="<<dL);
    if (dL >= pieces[k].lma) goto StepS; // to get the best
    if (sad > lb) // still
      goto CheckBetter;
    goto StepB2;*/
Exit:                             dbgpr(4,flush);
  *colBest = colTmpBetter;
//  colBest->Sort();
  colsetRes->cs.clear();
  if (found) colsetRes->cs.insert(*colBest);
//  for_each_in(colcclRes,iix,cycle<Column>::iterator)
//    colsetRes->Add(*iix);
} //____________________________________________________

}                                                // SS_END_NAMESPACE__
