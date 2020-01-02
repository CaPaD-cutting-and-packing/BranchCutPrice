// 2D 2-staged (=>guillotine) cutting (=>knapsack) problem
// transposition -> rotation
//
//////////////////////////////////////////////////////////////////////

// ALL VARS IN THE BEGINNING -- LIKE pascal.
// Output

// ATTENTION: with cuts, dimension = pr->Dim() = ma != m
// MINIMIZATION PROBLEM

#include "stdafx.h"
#include "probl_cp22.h"
#include "solver.h"
#include "raster.h"
#include "lasthdr.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

void CP22::Init() {
  int i;
  opttest = 0; // not initialized

  bb->Reallocate(m);
  BBCuts::Reallocate(m);
  //veps = deps;
  zi = +INFINITY__;
  lpb=llrv=lrv1= -INFINITY__; // lpb1 in FindPrimal?

  FillSigns();
  // + filling b:
  b.resize(m+1); // reserved FMax(100,m*2) or so (reall)
  ba.resize(m+1);
  b_cg.resize(m+1);
  for (i=0;i<m;++i) b_cg[i] = ba[i] = b[i] = pc[i].b;
  b_cg[m] = ba[m] = b[m] = (double)W; // the width constr.

// Algorithms:
  bb->m = m; bb->L = L;
  BBCuts::m = m;
  BBCuts::mc = 1;
  BBCuts::d__.resize(m+1);
//  BBCuts::mc = mc;
  BBCuts::L = L; // effective ?
  cMax =0; cMin=1e100;
  wMin = SIZE_MAX__;
  areaPriceMin = 1e+100;
  for (i=0;i<m;++i) {
    bb->pieces__[i].l = pc[i].l;
    bb->pieces__[i].b = pc[i].b;
    BBCuts::pieces__[i].l = pc[i].l;
    BBCuts::pieces__[i].b = pc[i].b;
    if (pc[i].c > cMax) cMax = pc[i].c;
    if (pc[i].c < cMin) cMin = pc[i].c;
    wMin = Min(wMin,pc[i].w);
    areaPriceMin = FMin(areaPriceMin,
      pc[i].c / double(pc[i].l) / double(pc[i].w));
  }
  bb->zMin = 1;
  bb->zLowerInitial = -INFINITY__;
  bb->eps = GetBBEps();

// FORBIDDEN:
  BBCuts::forbidden__ = bb->forbidden__ = forbidden;

  // LEVEL CUT: is a specialized superadditive cut
  // (alfa = 0)
/*  levelCut.resize(Dim());
//  levelCut.uSE.resize(Dim());
  for (i=0;i<m;++i)
    levelCut.u[i] = - pc[i].c;
*/
  gcdC = pc[0].c;
  for (i=1;i<m;++i)
    gcdC = __gcd(gcdC, pc[i].c);
  dbg_outn_(4," gcd(Ci)="<<gcdC);
  if (!gcdC) {
    gcdC = cMax / 100000;
    log_ln("WARNING: GCD(Ci) == 0, setting to "
      <<gcdC<<". Trying...");
  }
}

void CP22::FillSigns() {
  validsign.resize(m+1);
  fill_n(&validsign[0],m+1,-1); // Ax <= b
}

// returns 2 if I/O err, 1 if bad format, 0 else
int CP22::Read(istream &ifs,long long L_)
{
  char buf[1024];
  L0 = (size)L_;
  ifs >> W0 >> m0;
  if (!ifs) return 2;
  m=m0;
  if ((m0<2) or ((L0<3) or W0<3))
  { PRINT_ERROR("Bad data."); return 1; }
  pc0.resize(m0);
  int i;
  double d1,d2;
  for (i=0;i<m0;i++) {
    ifs >> pc0[i].l >> pc0[i].w;
    ifs.getline(buf,sizeof(buf));
    istringstream iss(buf); iss >> d1;// >> d2;
    if (!iss) {
     pc0[i].b = 1;
     pc0[i].c = double(pc0[i].l) * double(pc0[i].w);
    }
    else {
      iss >> d2;
      if (!iss) {
       pc0[i].b = (int)d1;
       pc0[i].c = double(pc0[i].l) * double(pc0[i].w);
      }
      else { pc0[i].b = (int)d2; pc0[i].c = d1; }
    }
    if ((pc0[i].l<=0) or (pc0[i].w<=0)
      or (pc0[i].b<0)) { // b < 0 !!
      PRINT_ERROR(infile<<", instance "<<inst
        <<": bad item data");
      return 1;
    }
    if ((pc0[i].l > L0) or (pc0[i].w > W0)) {
      PRINT_ERROR(infile<<", instance "<<inst
        <<": too large piece");
      return 1;
    }
  }
  if ((!ifs) && (!ifs.eof())) return 2;
  InitProblem();
  if (m<2) return 1;
  return 0;
}

void CP22::InitProblem()
{
  int i;
  // Sorting & compressing pieces' list (optionally)
  // + effective rod length ?
  pc = pc0;
  L = L0; // maybe calc. effective
  W = W0;
  if (not fFirstCut1stD) {
    Swap(L,W);
    for (int i=0;i<pc.size();++i)
      Swap(pc[i].l, pc[i].w);
  }
  //for (int i=0;i<m0;++i)
  //  pc[i].i_ = i;   // original numbers
  if (true/*fSortPieces*/) {
    for (i=0;i<m0;++i)
      pc[i].wgt = float(pc[i].w);    // weights for sorting
    sort(pc.begin(),pc.end(),greater<Piece>());
    if (fMergePieces) {
      int i,j;
      for (j=0,i=1; i<m0; ++i) {
        if (((pc[i].l != pc[j].l) or (pc[i].w != pc[j].w)
            or (pc[i].c != pc[j].c)) // FUCK !!!
          and not (0 == pc[i].b)) {
          ++j;
          pc[j] = pc[i];
        }
        else
          pc[j].b += pc[i].b;
      }
      m = j+1;
      pc.resize(m);
      if (m<2) PRINT_ERROR("m<2.");
      for (i=0;i<m;++i)
        pc[i].i0 = i;   // to know after sorting
    }
  }

  // EFFECTIVE LENGTHS
  Vector<int> sz(m), bi(m), rp1;
  for (i=0;i<m;++i) bi[i] = pc[i].b;
  if (W <= INT_MAX) {
		for (i=0;i<m;++i) sz[i] = (int)pc[i].w; // IMPORTANT: w's
    ConstructRP(sz,bi,(double)W,rp1);
		size W1 = rp1[FindRPUnder((int)W,rp1)];
    if (W1 != W) {
      dbg_outn(2," Old W="<<W<<", effective W="<<W1);
      W=W1;
    }
  } else {
    cout << "W > INT_MAX => no effective W computation" << endl;
  }
  if (L <= INT_MAX) {
		for (i=0;i<m;++i) sz[i] = (int)pc[i].l;
    rp1.clear();
    ConstructRP(sz,bi,(double)L,rp1);
		size L1 = rp1[FindRPUnder((int)L,rp1)];
    if (L1 != L) {
      dbg_outn(2," Old L="<<L<<", effective L="<<L1);
      L=L1;
    }
  } else {
    cout << "L > INT_MAX => no effective L computation" << endl;
  }
}

// ALL VARS in the beginning like PASCAL

void CP22::FFSBasis
  (const PieceContainer &pc__,ColumnList &bas)
{
//#warning Don't decrease bi's at FFSBasis
  
  // maybe sort pc with randomized weights
  // then use ->i0 to restore indices
  if (nStartBasis) { // GreedyBasis
    if (2==nStartBasis)
      return; // no one at all
    d_vec d(m), b(m);
    int i;
    for (i=0;i<m;++i) d[i] = pc__[i].c;
    for (i=0;i<m;++i) b[i] = pc__[i].b;
    ColSet css;
    GenColForHeur(&css,d,b,W);
    for_each_in(css.cs,ic,ColSet::iterator) {
      bas.Add(*ic);
      SetWConstr(&bas.cl.back());
    }
    return;
  }
  int i;
  PieceContainer pc=pc__;
	int m=(int)pc.size();
  for (i=0;i<m;++i)
    pc[i].i0 = i;
  // if (fRndFFSBas) { ... sort(pc...); }
  Vector<int> b0(m);
  for (i=0;i<m; ++i)
    b0[i] = pc[i].b;
  double sx = 0; // Solution value
  int j=0;  // column
  for ( ;(j<m); ++j ) {
    size L1=L;
    int x1=INT_MAX;
    if (0==pc[j].b) x1 = 0;
		Pattern a1((int)pc.size()); // stack<IX> ?
    int i;
    for (i=j; i<m; ++i) {
      int a=int(L1/pc[i].l);
      if (x1) {if (a > pc[i].b) a=pc[i].b;}
      else {if (a > b0[i]) a=b0[i];} // PROPER pattern
      L1 -= pc[i].l * a;
      if (a>0) {
        a1.PushIX(i,a);
        int x=pc[i].b/a;
        if (x<x1)
          x1 = x;
      }
    } // next i
/*    {for_each_in (a1.ix,iix,Pattern::iterator)
      pc[iix->i].b -= x1*iix->x;} //Need sorted numeration
      */
    sx += x1;
    // Restore original numeration if sorted:
    for_each_in (a1.ix,iix,Pattern::iterator)
      iix->i = pc[iix->i].i0;
    Column * c = &bas.Add();
    MakeColumn(c,&a1);
    SetObj(c);
//    SetWConstr(c); // -in MakeCol
  } // next j
//  UpdateHeurBnd(sx); // -- not in this form of FFD. */
} //____________________________________________________

void CP22::SetObj(Column * c) {
  Column::iterator iid;
  double ofc = 0;
  for_each_in(c->id,iid,)
    if (iid->i < m) ofc -= pc[iid->i].c * iid->d;
  c->SetObj(ofc); // -= : MINIMIZATION PROBLEM
}
void CP22::SetWConstr(Column *c) {
  Column::iterator iid;
  size w = 0;
  for_each_in(c->id,iid,) {
    assertm(iid->i < m, "Double adding W-constr");
    w = Max(pc[iid->i].w, w);
  }
  c->PushID(m, (int)w);
}
double CP22::GetW(Pattern *p) {
  Pattern::iterator iix;
  size ww = 0;
  for_each_in(p->ix,iix,) {
    assert(iix->i < m);
    ww = Max(ww,pc[iix->i].w);
  }
  return (double)ww;
}

void CP22::GenCol(ColSet &cs,const d_vec &d) {
  assert(cuts);
  assert(forbidden&&bb->forbidden__&&BBCuts::forbidden__);
  if (cuts->size())
    if (cuts->size()==1
      and cuts->front()==GetLevelCut())
        GenColWithCuts(cs,d); // Could be a spec proc
    else
      GenColWithCuts(cs,d); // CMP with multiwidth col gen
  else  // when one reports rdc=0 then the other also
    GenColPure(cs,d);
}

// ColGen: assumption: pieces sorted with non-incr w's
// Use all cols for the LP

// For Heur: ColSet? Then we may choose the most valuable
void CP22::GenColForHeur
  (ColSet *css,const d_vec &d,
  Vector<double> &b, size WLeft) {
 
  int i, im;
  int iMaxNZ = 0; // the last i with b_cg[i] > 0
  // double val;
  bb->fHeur = 1;

  css->clear();
  assert(d.size() >= m);
// Cleaning info from the prev generation: (+assert)
  for (i=0;i<m;++i) {
    bb->pieces__[i].d = d[i];
    bb->pieces__[i].b = int(b[i]);
    if (b[i] > 0) iMaxNZ = i;
  }
  dbg_outn_(5,"\nColGen for Heur. MULTS: ");
  if (OUTP_LEV__ >= 5) {
    for (i=0;i<m;++i) log__(d[i]<<' '); dbg_outn(4,"");
    log__("iMaxNZ="<<iMaxNZ);
  }

for (im=0;im<=iMaxNZ;++im)
if (b[im] > 0)
if (pc[im].w <= WLeft) // the reduced width
{
  for (i=0;i<im;++i) bb->pieces__[i].b = 0;
  bb->pieces__[im].b = int(b[im]) - 1;
  bb->L = L - pc[im].l;
  // THE BOUND USING PREVIOUS SOLUTIONS:
  bb->zMin = 
    bb->zLowerInitial = -INFINITY__;
  bb->eps = GetBBEps();

  ColSet cs; Column col;
  col.PushID(im,1); // the beginning piece is already there
//  col.PushID(m,pc[im].w); // the width constraint -- here not.
  bb->colsetRes = &cs;
  bb->colBest = &col; // not forget clBest

  bb -> Run();
//  if (not (bb->found)) continue;

  SetObj(&col);
  if (not col.id.empty())
    css->cs.insert(col);
}

} // GenColForHeur()

/*// To be called in each B&B node
void CP22::CopyForbCols(list<Column*> & fc) {

  list<Column*>::iterator ifc;
  bb->forbidden__.clear();
  BBCuts::forbidden__.clear();
  list<Vector<IX> > * pfc=
//    (cuts->size()) ?      &(BBCuts::forbidden__) :
      &(bb->forbidden__);
  Vector<IX> vix;
  for_each_in(fc,ifc,) {
    vix.clear();
//    pix->reserve((*ifc)->id.size());
    for_each_in((*ifc)->id,iid,Column::iterator)
      if (iid->i < m) // CP22 !!!
        vix.push_back(IX(iid->i,iid->d));
    pfc->push_back(vix);
  }
  BBCuts::forbidden__ = *pfc;
}*/

// To be called only by LP
void CP22::GenColPure
  (ColSet &css,const d_vec &d) {
  int i, im;
  int iMaxNZ = 0; // the last i with b_cg[i] > 0
  double val;
  bb->fHeur = false;

  css.clear();
  assert(d.size() >= m+1);
// Cleaning info from the prev generation: (+assert)
  redCostBest = 0;  // SORTING: PIECES' W's non-increasing
  for (i=0;i<m;++i) {
    bb->pieces__[i].d = d[i] + pc[i].c;
    bb->pieces__[i].b = int(b_cg[i]);
    if (b_cg[i] > 0) iMaxNZ = i;
  }
  dbg_outn_(4,"\nMULTS: ");
  if (OUTP_LEV__ >= 4) {
    for (i=0;i<=m;++i) log__(d[i]<<' '); dbg_outn(4,"");
    log__("iMaxNZ="<<iMaxNZ);
  }

for (im=0;im<=iMaxNZ;++im)
if (b_cg[im] > 0)
if (pc[im].w <= b_cg[m]) // the reduced width
{
  for (i=0;i<im;++i) bb->pieces__[i].b = 0;
  bb->pieces__[im].b = int(b_cg[im]) - 1;
  bb->L = L - pc[im].l;
  // THE BOUND USING PREVIOUS SOLUTIONS:
  bb->zMin =  - redCostBest
    - bb->pieces__[im].d - (double)pc[im].w * d[m];
  bb->zLowerInitial = -INFINITY__;
  bb->eps = GetBBEps();

  ColSet cs; Column col; // to search forbidden:
  col.PushID(im,1); // the beginning piece is already there
  col.PushID(m,(int)pc[im].w); // the width constraint
  bb->colsetRes = &cs;
  bb->colBest = &col; // not forget clBest

  bb -> Run();
  if (not (bb->found)) continue;

  SetObj(&col);
  // FIX THE BEST COL & ITS REDUCED COST:
  val =  - bb->z
    - bb->pieces__[im].d - (double)pc[im].w * d[m];
  if (val < redCostBest) {
    redCostBest =  val;
    clBest = col;
  // Now leave only this best col:
    css.clear(); css.Add(col);
  }
}
  UpdateLagrBnd(d);  //  !!!!!!!!!!
}

// NEW PROC: COMAPRE with the old
//  -- if not found, try old
void CP22::GenColWithCuts (ColSet &css,const d_vec &d)
{
  int i, im;
  int iMaxNZ=0;
  css.clear();
  assert(d.size() >= Dim() + cuts->size());
  redCostBest = -1e100;
  BBCuts::nStepsMin = FMax( // increasing with ... ??
    BBCuts::nStepsMin, nStepsMin0)
    * nStepsMinInc; // ... each col gen!

// + multipliers: eliminating 0's:
  BBCuts::cuts__.dep.clear();
  i=Dim();
  CutList::iterator ic;
  for_each_in(*cuts,ic,) {
//    cout << d[i]<<' '<<fabs(d[i])<<' '<<GetDEps()
//      <<(fabs(d[i]) > GetDEps())<<endl;;
    if (fabs(d[i]) > 1e-15) {
//      cout << '+';
      BBCuts::cuts__.dep.push_back
        (SACutSE::Dependence(d[i],*ic));
    }
    ++ i;
  }
  BBCuts::cuts__.AddInvolved(cuts);

  for (i=0;i<m;++i) {
    BBCuts::d__[i] =
      BBCuts::pieces__[i].d = d[i] + pc[i].c;
    BBCuts::pieces__[i].b = int(b_cg[i]); // again
    if (b_cg[i] > 0) iMaxNZ = i;
  }
  for (i=m;i<Dim();++i)
    BBCuts::d__[i] = d[i]; // for w-constr
  if (OUTP_LEV__ >= 4) {
    dbg_outn_(4,"\nMULTS: ");
    for (i=0;i<=m;++i) log__(d[i]<<' '); dbg_outn(4,"");
    log__("iMaxNZ="<<iMaxNZ);
  }
  if (OUTP_LEV__ >= 4) {
    dbg_outn_(4,"\nRHS: ");
    for (i=0;i<=m;++i) log__(b_cg[i]<<' '); dbg_outn(4,"");
    log__("iMaxNZ="<<iMaxNZ);
  }

  double rdcBest = 0;
  BBCuts::foundInit=false;

for (im=0;im<=iMaxNZ;++im)
if (b_cg[im] > 0)
if (pc[im].w <= b_cg[m]) // the reduced width
{ // <= !!!
  for (i=0;i<im;++i)
    BBCuts::pieces__[i].b = 0;
  BBCuts::pieces__[im].b = int(b_cg[im]) - 1;
  BBCuts::L = L - pc[im].l;
  // THE BOUND USING PREVIOUS SOLUTIONS:
  BBCuts::zMin =  - rdcBest;
//    - pc[im].w * d[m]; // this comes later
  BBCuts::zLowerInitial = -INFINITY__;
  BBCuts::eps = GetBBEps();

  ColSet cs;
  Column col;
  BBCuts::colsetRes = &cs;
  BBCuts::colBest = &col; // not forget clBest
  col.PushID(im,1); // the beginning item
  col.PushID(m,(int)pc[im].w); // adding the width constr coef
  // THIS MAY BE A BETTER COL! BBCuts must handle appr.

  // THE BOUND USING PREVIOUS SOLUTIONS:
  BBCuts::zMin =  - rdcBest;
  BBCuts::Init();
  BBCuts::Run();
  if (not (BBCuts::found)) continue;

  BBCuts::foundInit=true;
// in the following runs early term. can be made
  // if a good solution found in previous runs!!!

  SetObj(&col);
  double rdcCheck = CalcRedCost(&col,cuts,d);
  assertm(fabs(BBCuts::z + rdcCheck)<1e-6,
    "Recalc red cost: -BBCuts::z="<<-BBCuts::z
      <<", recalc="<<rdcCheck);
  // FIX THE BEST COL & ITS REDUCED COST:
  assert (-BBCuts::z < rdcBest);
    rdcBest =  - BBCuts::z;
    css.clear(); //leave only the best col
    css.Add(col);
  if (!BBCuts::fETerm) {
    redCostBest = rdcBest; // !!!
    clBest = col;
  } else redCostBest = -1e100;
}
  if (redCostBest > -1e50) // after all widths
    UpdateLagrBnd(d);  //  !! here early termination!!!
// TEST:
/*  Column c1(5);
  c1.PushID(24,1); c1.PushID(25,1); c1.PushID(29,1);
  SetWConstr(&c1); SetObj(&c1);
  mylog<<"Col 24 25 29: rc="<<CalcRedCost(&c1,cuts,d)<<endl;
*/  // ?????????????????????????????????
} //____________________________________________________

void CP22::PrintColumn(ostream& os,Column *c) {
  int i;
//  if (c->GetCutSlackCut())
//    os << "Cut slack: " << c->GetCutSlackCoef()<< ' ';
      for (i=0; i<c->id.size(); ++i)
        if (c->id[i].i < m) // not the width constr.
//          os << pc[c->id[i].i].l <<'x'<< pc[c->id[i].i].w
          os << c->id[i].i
            <<':'<<c->id[i].d<<' ';
        else os << 'w' << c->id[i].d << ' ';
      os << "obj " << c->GetObj();
}

void CP22::MakePattern(Pattern *p, Column *c) {
    p->clear();
    for_each_in(c->id,iid,Column::iterator)
      if (iid->i < m)
        p->PushIX(iid->i, (int)iid->d);
    p->SetObj(c->GetObj());
} ////////////
/*LP::Cut * CP22::GetLevelCut() {
  levelCut.lhs = GetLocalLPBnd();  //  ??????????????
  return &levelCut;
}
double CP22LevelCut::Calc__(Column * c) {
  return c->IsSlack() ?
    (this == c->GetCutSlackCut() ? c->GetCutSlackCoef() : 0)
    : SACutSE::Calc__(c);
}*/


////////////////////////////////////////////////////////
/////////// Integer Rounding ///////////////////////////
////////////////////////////////////////////////////////

double CP22::raster_ceil(const double v) {
  if (0 == opttest) { // no ratser points yet
    if (GetLPBnd() > -1e50 && zi < 1e50) {
      InitOptTest(); // not calling GetHeurBnd() 'cause..
      goto TestRaster;
    }
    return gcdC * ceilEps(v/gcdC, deps/gcdC);
  }
TestRaster:
  if (v < rpM*rpc[0]) return gcdC * ceilEps(v/gcdC, deps/gcdC);
  if (v > rpM*rpc.back()) return gcdC * ceilEps(v/gcdC, deps/gcdC);
  while (rpM*rpc[lastPos] > v - GetVEps()) -- lastPos;
  while (rpM*rpc[lastPos] < v - GetVEps()) ++ lastPos;
  return rpM * rpc[lastPos];
}
double CP22::raster_below(const double v) {
  if (0 == opttest) { // no ratser points yet
    return gcdC * (ceilEps(v/gcdC, deps/gcdC) - 1);
  }
  if (v < rpM*rpc[0]) return gcdC * (ceilEps(v/gcdC, deps/gcdC)-1);
  if (v > rpM*rpc.back()) return gcdC * (ceilEps(v/gcdC, deps/gcdC)-1);
  while (rpM*rpc[lastPos] > v - GetVEps()) -- lastPos;
  while (rpM*rpc[lastPos] < v - GetVEps()) ++ lastPos;
  return rpM * rpc[IMax(0,lastPos-1)];
}
void CP22::InitOptTest() {
  opttest = 1;
  double nnn = - GetLPBnd() / gcdC;
  rpM = gcdC;
  if (nnn > ((unsigned)(1))<<28) {
    // we allow only 64 MB for the raster point generation
    rpM = gcdC * nnn / (((unsigned)(1))<<28);
  }
  Vector<int> ci(m), bi(m), rpc1;
  int i;
  for (i=0;i<m;++i) ci[i] = int(pc[i].c / rpM);
  for (i=0;i<m;++i) bi[i] = int(pc[i].b);
  ConstructRP(ci,bi,-GetLPBnd()/rpM+1,rpc1);
  int sp = int(-GetHeurBnd()/rpM)-1;
  if (sp<0) sp=0;
  int p = FindRPUnder(sp,rpc1);
//  rpc.assign(rpc1.begin()+p,rpc1.end());
  // rotating and negating:
  rpc.reserve(rpc1.size()-p+1);
  rpc.push_back(-INT_MAX);
	for (i=(int)rpc1.size()-1;i>=p;--i) rpc.push_back(-rpc1[i]);
  lastPos = 0;
}

bool CP22::ConstructIntSol() {
  assert(pat.size() && lpx.size() && lpd.size()
    ); // Input
  Round();
  do {
    CompleteIntSol();
    if (Optimum()) break;
  } while (VaryResProblem()); // normally extending
  lpx.clear();
  lpd.clear();
  if (OUTP_LEV__ >=3) log__("=" << zi);
  return Optimum();
}

double CP22::GetXEps() { return deps; }
void CP22::SaveRoundedPart() {
  patBest.resize(pat.size()); patBest = pat;
  xiBest.resize(xi.size()); xiBest = xi;
}

  class J0Cmp { Vector<double> &xf; public:
    J0Cmp(Vector<double> &x_)  : xf(x_) { }
    bool operator () (int j1, int j2) const
    { return xf[j1] > xf[j2]; }
  };
bool CP22::Round() {
  int i;
  xi.resize(lpx.size());
  xf.resize(lpx.size());    ///////// Rnd down:
	Vector<int> j0((int)xi.size()); // Indices; won't be used outs
  br.resize(m); GetRHS(br);

  Vector<int> dd(3);
  fracused.clear(); fracused.reserve(xi.size());
  zr = 0; WLeft = W;
  dbg_outn(5,"Rnd dn cols:");
  for (i=0;i<xi.size();++i) {
    xi[i] = floor(lpx[i]+GetXEps());
    xf[i] = lpx[i] - xi[i];
    j0[i] = i;   // Updating rhs:
    if (xi[i])
    for_each_in(pat[i].ix,iix,Pattern::iterator) {
      br[iix->i] -= xi[i] * iix->x;
      assertm(br[iix->i] >= 0, "Too many items in rounded sol.");
    }
    zr += pat[i].GetObj() * xi[i];
    WLeft -= size(GetW(&pat[i]) * xi[i]);
    if (OUTP_LEV__ >=5) {
      for_each_in(pat[i].ix,iix,Pattern::iterator)
        log__(pc[iix->i].l<<'x'<<pc[iix->i].w
          <<':'<<iix->x<<' ');
      log__("obj="<<pat[i].GetObj()<<" x="<<xi[i]<<'\n');
    }
  }
  if (OUTP_LEV__ >=3) log__(" Rnd\\" << zr);
  dbg_outn_(5," WLeft"<<WLeft);
  dbg_outn(5,"Rnd up: +cols:");
  sort (&j0[0], &j0[0]+j0.size(), J0Cmp(xf));
//  bool somefit;
//  int i0 = 0;
  i = 0;  nRPE = 0; // Init Variation of Residual
//////////////////////////////////////////////// Rnd up:
//  do {     somefit = false;  // for multiple
    do {
      bool fits = (WLeft >= GetW(&pat[j0[i]]));
      if (!fits) break;
      for_each_in(pat[j0[i]].ix,iix,Pattern::iterator)
        if (/*fabs(lpd[iix->i]) > GetDEps() // signif. demand
          and */ br[iix->i]  <  iix->x)
          { fits = false; break; }
      if (fits) {
        fracused.push_back(j0[i]);
        ++ xi[j0[i]];
        for_each_in(pat[j0[i]].ix,iix,Pattern::iterator) {
          br[iix->i] -= iix->x;
//          if (br[iix->i] < 0) br[iix->i] = 0;
        }


        zr += pat[j0[i]].GetObj();
        WLeft -= size(GetW(&pat[j0[i]]));
    if (OUTP_LEV__ >=5) {
      for_each_in(pat[j0[i]].ix,iix,Pattern::iterator)
        log__(pc[iix->i].l<<'x'<<pc[iix->i].w
          <<':'<<iix->x<<' ');
      log__("obj="<<pat[j0[i]].GetObj()
        <<" x="<<xi[j0[i]]<<'\n');
    }
      }
    } while ( ++i < j0.size() );

  if (OUTP_LEV__ >=3) log__("/" << zr<<flush);
  dbg_outn_(5," WLeft"<<WLeft);
//  } while (somefit)
  return true;
} //////////////////////////////// CP22::Round

void CP22::GetRHS(Vector<double> &br) {

  br.resize(m); int i;
  for (i=0;i<m;++i) br[i] = pc[i].b;
} //////////////////////////////// CP22:: GetRHS

// NEGATE ci to get a minimization problem ?
//namespace {
class CP22::SVC: public Alg {
/////// INPUT via constructor:
  CP22 * pr;
  const Vector<double> & b0;
  const Vector<double> & d0;
  const double lb, ub; // a lower/upper bnd
  const size L, W;
//////// VARS:
  Vector<double> bc;  // the remaining demands

  d_vec d;
  const int m;
  unsigned int kkk; // total iteration number
  int k; // iteration number
  double zz; // current value
  double result; // best value
  int xMin; // intensity of current pattern
  double waste;
  size WLeft;
  size wLastPat;
  double SW;  // randomizer

  typedef Vector<Pattern> PatArray;
  PatArray pat;
  Vector<double> x;
public:
//////////// OUTPUT:
  PatArray patBest;
  Vector<double> xBest;

  SVC(CP22* p_, Vector<double> & b_,
      Vector<double> & d_, const double lb_,
      const double ub_, const size W_)
    : Alg(p_), pr(p_), m(pr->m), b0(b_), d0(d_),
    lb(lb_), ub(ub_), W(W_), L(pr->L), kkk(0)
  { d.resize(m); bc.resize(m); }
  double Run() {
//    try {
      Execute();
      if (result < lb) {
        int i;
        log_ln("SVC CP22: Problem: L="<<L<<" W="<<W);
        for (i=0;i<m;++i)
          log_ln(" pc "<<i<<": "<<pr->pc[i].l<<'x'<<pr->pc[i].w
            <<'='<<pr->pc[i].c<<'.'<<b0[i]);
        log_ln("\n Solution:");
        for (i=0;i<patBest.size();++i) {
          Column col;
          pr->MakeColumn(&col,&patBest[i]);
          mylog << "x=" << xBest[i]<<' ';
          pr->PrintColumn(mylog,&col);
          mylog << '\n';
        }
        assertm(false,"SVC CP22: result="<<result<<" < lb="<<lb);
      }
//      assert(result <= ub); // NO, may be all
/*    } catch (const exception & e) {
      PRINT_ERROR(e.what());
      return 1e+100;
    } catch (...) {
      PRINT_ERROR("Unknown error.");
      return 1e+100;
    }*/
    return result;
  }
private:
  void Execute() {
    result = 1e+100;
    k=0;
    ConvertValues();
    do {
      try {
        zz = 0; bc = b0; pat.clear(); x.clear();
        WLeft = W;
        while (GenPat()) {
          CorrectValues();
        }
        CorrectValuesOfFree();
        ControlSolution();
        if (zz < result) {
          result = zz;
          SaveSolution();
          if (zz <= lb) return;
        }
      } catch (...) {
        if (OUTP_LEV__ >=3) PRINT_ERROR("SVC: Some error.");
      } // Let it be
    } while (++k < iterMax );
  }
  void ConvertValues() { // bec. of scale
    int i; // ATTENTION: when all duals of remaining
          //  pieces very small, our bb procedure doesnt
    for (i=0; i<m; ++i) {
      d[i] = pr->pc[i].c + d0[i];
      if (d[i] < 1) d[i] =1;
    }
  }
  bool GenPat() {
    double sbi = 0;  int i;
    for (i=0; i<m; ++i) sbi += bc[i];
    if (0 == sbi or WLeft < pr->wMin) {
      return false;
    }
    ColSet css;
    pr->GenColForHeur(&css,d,bc,WLeft);
    if (0==(css.cs.size())) return false;
    Column col;
    ChooseCol(&css,&col);
    pr->SetObj(&col);
    pat.push_back(Pattern());
    pr->MakePattern(&pat.back(), &col);
    wLastPat = size(pr->GetW(&pat.back()));
    xMin = int(WLeft/wLastPat); // pattern intensity:
    assert(xMin);
    waste = (double)L;  // original sorting ???
    Pattern::iterator iix;
    for_each_in(pat.back().ix,iix,) {
      assert(iix->x);
      assert(bc[iix->i] > 0);
      if (bc[iix->i] < double(xMin)*iix->x)
        xMin = int(double(bc[iix->i]) / iix->x); // floor
      waste -= double(pr->pc[iix->i].l)
        /* * pr->pc[iix->i].w */ // !!!!!!! double done
        * iix->x;
    }
    if (patternUseRatio > 0 and patternUseRatio <= 1)
      xMin = (int)ceil (double(xMin) * patternUseRatio);
    for_each_in(pat.back().ix,iix,)
       bc[iix->i] -= double(xMin)*iix->x;   // DECR. DMND
    zz += pat.back().GetObj() * xMin;
    x.push_back(xMin);
    WLeft -= wLastPat * xMin;
    waste *= double(wLastPat); // the 2nd dimension
    return true;

  }
  void ChooseCol(ColSet *cs, Column *col) {
    ColSet::iterator icm = cs->cs.begin();
    double valueBest = -1e100;
    for_each_in(cs->cs,ic,ColSet::iterator) {
      size w=0;
      double sumd=0;
      for_each_in(ic->id,iid,Column::iterator) {
        sumd += d[iid->i] * iid->d;
        w = Max(w,pr->pc[iid->i].w);
      }
      double value = sumd - pr->areaPriceMin * double(L) * double(w);
      if (value > valueBest) {
        valueBest = value;
        icm = ic;
      }
    }
    *col = *icm;
//    *col = cs->front(); // THE MOST STUPID
//    cout << " ??? What col ???" << endl;
  }
  void CorrectValues() {
    SW = 1.1 + fabs(double ( ++ kkk % 40 - 15)) / 10;
    assert ( waste < L * wLastPat );
    double LLh = double(L)*double(wLastPat)
      / (double(L)*double(wLastPat)-waste);
    for_each_in(pat.back().ix,iix,Pattern::iterator)
      d[iix->i] = ( LLh * pr->pc[iix->i].c * iix->x
        * double(wLastPat) / double(pr->pc[iix->i].w) // considering
         // the individual waste
        + d[iix->i] * SW * (b0[iix->i] +bc[iix->i]))
        / (SW * (b0[iix->i] +bc[iix->i]) + iix->x);
  }
  void CorrectValuesOfFree() {
    waste = double(L) * double(WLeft);
    double sbi = 0; int i;
    for (i=0;i<pr->m;++i) // not involved types
      if (bc[i] && bc[i] == b0[i]) sbi += pr->pc[i].c * bc[i];
    for (i=0;i<pr->m;++i)
      if (bc[i] && bc[i] == b0[i])
        d[i] += double(waste) * pr->pc[i].c * bc[i] / sbi;

  }
  void ControlSolution() { // the last one
    assert(pat.size() == x.size());
    Vector<double> s_aij(m); // fillled 0 ?
    double zzz = 0;
    int i; size ww = 0;
    fill_n(&s_aij[0],m,0);
    for (i=0;i<pat.size(); ++i) {
      for_each_in(pat[i].ix,iix,Pattern::iterator)
        s_aij[iix->i] += iix->x * x[i];
      zzz += pat[i].GetObj() * x[i];
      ww += size(pr->GetW(&pat[i]) * x[i]);
    }
    assertm(zzz == zz, "zzz="<<zzz<<", zz="<<zz);
    for (i=0;i<m;++i)
      assertm(s_aij[i] <= b0[i],
       "s_aij[i]="<<s_aij[i]<<", b0[i]"<<b0[i]);
    assert(ww <= W);


}

  void SaveSolution()
  {
    patBest = pat; xBest = x;
    if (OUTP_LEV__ >= 5) {
      log_ln("SVC_CP22: Saving better solution:");
      int i;
      for (i=0;i<pat.size();++i) {
        for_each_in(pat[i].ix,iix,Pattern::iterator)
          log__(pr->pc[iix->i].l<<'x'<<pr->pc[iix->i].w
            <<':'<<iix->x<<' ');
        log__("obj="<<pat[i].GetObj()<<" x="<<x[i]<<'\n');
      }
      log_ln("WLeft = "<<WLeft);
    }
  }
  static opt::OptContainer Options();
  static opt::OptSection opt;
public:
  static int iterMax;
  static double outputLevel;
  static double patternUseRatio;
}; ////////////////////// class CP22::SVC

int CP22::SVC::iterMax=20;
double CP22::SVC::patternUseRatio=0.5;
double CP22::SVC::outputLevel=0;

opt::OptContainer CP22::SVC::Options() {
  opt::OptContainer oc;
  oc
    << opt::MakeOpt(&iterMax, 20,
      "iterMax", 0)
    << opt::MakeOpt(&patternUseRatio, 0.5,
      "patternUseRatio",
      "Usage intensity of a new pattern rel. to max")
    << opt::MakeOpt(&outputLevel, DEF_OUTP_LEVEL,
      "outputLevel", "0-5");
  return oc;
} //____________________________________________________
opt::OptSection CP22::SVC::opt
  ("CP22_SVC", "SVC heur for the 2D CP",
  SVC::Options(), opt::SolverCfg(), 5500);
//} ////////////////////// namespace { }

bool CP22::CompleteIntSol() {
  int i;
  double sbi=0;
///// Check if all fit in a single rod:
  for (i=0;i<m;++i) sbi += br[i];
  if (zr < GetHeurBnd() && (0 == sbi || WLeft < wMin)) {

    UpdateHeurBnd(zr);
    SaveRoundedPart();
    return true;
  }
  if (OUTP_LEV__ >=4)
    log__(" Adding to rounded value: " << zr);
  SVC svc(this, br, lpd, GetLPBnd()-zr, GetHeurBnd()-zr, WLeft);
  int iterMax__ = SVC::iterMax;
  if (fFastRounding) SVC::iterMax = IMin(SVC::iterMax,3);
  double SVCres=svc.Run();
  SVC::iterMax = iterMax__;
  if (OUTP_LEV__ >=4) log__(" With SVC: " << zr+SVCres);
  if (zr + SVCres < zi) {
    SaveRoundedPart();
    AddResidualSolution(svc);
    UpdateHeurBnd(zr + SVCres);
    ControlSolution();
    if (0==opttest && GetLPBnd() > -1e+50) {
      InitOptTest();
//      dbg_outn_(4," lpb before RP:"<<GetLPBnd());
//      RecalcLPBnds();
//      dbg_outn_(4," after RP:"<<GetLPBnd());
    }
  }
//  if (Optimum())
  return true;
} //////////////////////////////// CP22::CompleteIntSol

void CP22::AddResidualSolution(CP22::SVC & svc) {
  patBest.insert
    (patBest.end(),svc.patBest.begin(),svc.patBest.end());
  xiBest.insert
    (xiBest.end(),svc.xBest.begin(),svc.xBest.end());
}

void CP22::ControlSolution() { // the best one
    assert(patBest.size() == xiBest.size());
    Vector<double> s_aij(m); // fillled 0 ?
    double zzz = 0;
    int i;
    size ww=0;
    fill_n(&s_aij[0],m,0);
    for (i=0;i<patBest.size(); ++i) {
      for_each_in(patBest[i].ix,iix,Pattern::iterator)
        s_aij[iix->i] += iix->x * xiBest[i];
      zzz += patBest[i].GetObj() * xiBest[i];
      ww += size(GetW(&patBest[i]) * xiBest[i]);
    }
    assert(zzz == zi);

    for (i=0;i<m;++i) assert(s_aij[i] <= pc[i].b);
    assert(ww<=W);

}

bool CP22::VaryResProblem() {
  if (++nRPE > RPEMax) return false;
  if (fFastRounding) return false;
  int i;
  if (fracused.empty()) { // No increased components
/*    for (i=xi.size()-1;i>=0;--i) // normally
      // consider j0 & not fracused...
      if (xi[i]) break;
    if (i<0)*/ return false; 
  } else
  { i = fracused.back(); fracused.pop_back(); }
  -- xi[i] ;  zr -= pat[i].GetObj(); WLeft += (size)GetW(&pat[i]);
  for_each_in(pat[i].ix, iix, Pattern::iterator)
    br[iix->i] += iix->x;
  return true;
} //////////////////////////////// CP22::VaryResProblem

void CP22::PrintProblem(ostream&os) {
  os <<
    "CP22 instance name: "<<prName<<", No. "<<inst<<
    " from file "<<infile<<"\\\\\n"<<
    L<<" "<<W<<" "<<m<<'\n';
  for (int i=0;i<m;++i)
    os << pc[i].l<<' '<<pc[i].w<<' '
    <<pc[i].c<<' '<<pc[i].b<<'\n';
  os << endl;
}
void CP22::MakeColumn(Column * c,Pattern * p)
  { 
    c->clear();
    for_each_in(p->ix,iix,Pattern::iterator)
      c->PushID(iix->i,iix->x);
//    c.id.insert(c.id.end(),id.begin(),id.end());
    c->ofc=p->ofc;
    c->GetAddiInfo() = p->GetAddiInfo();
    SetWConstr(c); // !!!!
    c->Sort(); // !!!
  }


bool CP22::fFirstCut1stD=true;
bool CP22::fSortPieces=true;
bool CP22::fMergePieces=true;
int CP22::nStartBasis=0;
double CP22::nStepsMin0=65536;
double CP22::nStepsMinInc=1.05;
double CP22::deps=1e-6, CP22::bb_eps=1e-6;
int CP22::RPEMax = 10;
double CP22::outputLevel=DEF_OUTP_LEVEL;

opt::OptContainer CP22::Options() {
  opt::OptContainer oc;
  oc
    << opt::MakeOpt(&fFirstCut1stD, true,
      "fFirstCut1stD", "First cut along the first dimension")
    << opt::MakeOpt(&fSortPieces, true,
      "fSortPieces", 0)
    << opt::MakeOpt(&fMergePieces, true,
      "fMergePieces",
      "Bool: Merge equal piece types")
    << opt::MakeOpt(&nStartBasis, 2,
      "nStartBasis",  "0: FFD (best); 1: Greedy, 2: empty.T")
    << opt::MakeOpt(&nStepsMin0, 4096,
      "nStepsMin0",
      "Initial min N steps in B&B with cuts")
    << opt::MakeOpt(&nStepsMinInc, 1.01,
      "nStepsMinInc", "Incr ratio (each generation)")
    << opt::MakeOpt(&deps, 1e-6,
      "deps", "eps for dual multipliers")
    << opt::MakeOpt(&bb_eps, 1e-6,
      "bb_eps", "eps for b&b col gen")
    << opt::MakeOpt(&RPEMax, 10,
      "RPEMax", "IMax residual problem extensions")
    << opt::MakeOpt(&outputLevel, DEF_OUTP_LEVEL,
      "outputLevel", "0-5");
  return oc;
} //____________________________________________________
opt::OptSection CP22::opt
  ("CP22", "The 2D 2-Stage (=>Guillotine) Cutting Problem",
  CP22::Options(), opt::SolverCfg(), 500);

}                                                // SS_END_NAMESPACE__
