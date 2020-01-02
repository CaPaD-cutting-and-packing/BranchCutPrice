// probl_csp2.cpp: Implementierung der Klasse CSP2.
//
//////////////////////////////////////////////////////////////////////

// ALL VARS IN THE BEGINNING -- LIKE pascal.
// Output

// ATTENTION: with cuts, dimension = pr->Dim() = ma != m

#include "stdafx.h"
#include "probl_csp2.h"
#include "lasthdr.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

void CSP2::Init() {
  int i;

  //veps = deps;
  zi = +INFINITY__;
  lpb=llrv=lrv1= -INFINITY__; // lpb1 in FindPrimal?

  FillSigns();
  // + filling b:
  b.resize(m); // reserved IMax(100,m*2) or so (reall)
  ba.resize(m);
  for (i=0;i<m;++i) ba[i] = b[i] = pc[i].b;

// Algorithms:
  ggk2.m = //wang.m    =
    m;
  ggk2.Reallocate();
 // wang.Reallocate();
  ggk2.L = //wang.L    =
    L; // effective sizes ?
  ggk2.W = //wang.W    =
    W; // effective sizes ?
  for (i=0;i<m;++i) {
    ggk2.l[i] = pc[i].l;
    ggk2.w[i] = pc[i].w;
//    wang.b[i] = pc[i].b; // will be set specifically
//    ggk2.c[i] = wang.c[i] = ...;
  }
/*
  if (WangGammaStep<=0 or WangGammaStep>=1)
    WangGammaStep = 0.97;
  wang.gammaStep = WangGammaStep;
  wang.NSol = WangNSol;*/

}

void CSP2::FillSigns() {
  validsign.resize(m);
  fill_n(&validsign[0],m,1); // Ax >= b
}

// returns 2 if I/O err, 1 if bad format, 0 else
int CSP2::Read(istream &ifs,long L_)
{
// Old format for multiple stock:
  bool multi = (L_ < 0);  int M=1;
  L0 = (size)L_;
  if (multi) ifs >> m0 >> M;
  else ifs >> W0 >> m0;
  if (!ifs) return 2;
//  if (!multi && (m0>L0)) Swap(m0,L0);
  m=m0;
  if ((m0<2) or ((L0<3) and !multi))
  { PRINT_ERROR("Bad data."); return 1; }
  pc0.resize(m0);
  int i;
  for (i=0;i<m0;i++) {
    ifs >> pc0[i].l >> pc0[i].w >> pc0[i].b;
    if ((pc0[i].l<=0) or (pc0[i].w<=0)
      or (pc0[i].b<0)) { // b < 0 !!
      PRINT_ERROR(infile<<", instance "<<inst
        <<": bad item data");
      return 1;
    }
    if ((pc0[i].l > L0 or pc0[i].w > W0) && !multi) {
      PRINT_ERROR(infile<<", instance "<<inst
        <<": too large piece");
      return 1;
    }
  }
  if (multi) {
    int dm;
    ifs >> L0  >> dm;
    for (i=1;i<M;++i) ifs >> dm >> dm; // skipping other lengths
  }
  if ((!ifs) && (!ifs.eof())) return 2;
  InitProblem();
  if (m<2) return 1;
  return 0;
}

void CSP2::InitProblem()
{
//  int i;
  // Sorting & compressing pieces' list (optionally)
  // + effective rod length ?
  pc = pc0;
  L = L0; W = W0; // maybe calc. effective sizes
  //for (int i=0;i<m0;++i)
  //  pc[i].i_ = i;   // original numbers
  if (fSortPieces) {
/*    for (i=0;i<m0;++i)
      pc[i].wgt = pc[i].l;    // weights for sorting
    sort(pc.begin(),pc.end(),greater<Piece>());*/
    if (fMergePieces) {
      int i,j;
      for (j=0,i=1; i<m0; ++i) {
        if ((pc[i].l != pc[j].l or pc[i].w != pc[j].w)
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
}

// ALL VARS in the beginning like PASCAL

// NEED at least a diagonal matrix for the beginning!
void CSP2::FFSBasis
  (const PieceContainer &pc__,ColumnList &bas)
{
  int i;
  d_vec d(m);
  Vector<double> bb(m);
  PieceContainer pc=pc__;

  dbg_outn_(3," INIBAS..");

  bas.cl.clear();
  for (i=0;i<m;++i) {
    d[i] = double(pc[i].l) * pc[i].w;
    bb[i] = pc[i].b;
  }

  for (i=0;i<m; d[i++]=0 ) {
    if (bb[i] <= 0) continue; // does d[i++] = 0
    bas.Add(Column());
    bas.cl.back().SetObj(1);
    GenColForHeur(&bas.cl.back(),d,bb);
    int x1 = INT_MAX;
    Column::iterator iid;
    for_each_in(bas.cl.back().id,iid,)
      if (double(iid->d) * x1 > bb[iid->i])
        x1 = (int)ceil(double(bb[iid->i]) / iid->d);
    // NOW subtracting the col from bb x1 times:
    for_each_in(bas.cl.back().id,iid,) {
      bb[iid->i] -= iid->d * x1;
      if (bb[iid->i] <= 0) {
        bb[iid->i] = 0; d[iid->i] = 0;
      }
    }
  }
// AND ADDING A DIAGONAL MATRIX:
  for (i=0;i<m;++i)  { d[i] = 0; }

  for (i=0;i<m; d[i++]=0 ) {
    ColSet cs;
    d[i] = 10;
    GenColPure(cs,d);
    bas.Add(clBest);
  }
  dbg_outn_(3," ENDBAS");

/*  // maybe sort pc with randomized weights
  // then use ->i0 to restore indices
  int i;
  PieceContainer pc=pc__;
  int m=pc.size();
  for (i=0;i<m;++i)
    pc[i].i0 = i;
  // if (fRndFFSBas) { ... sort(pc...); }
  Vector<int> b0(m);
  for (i=0;i<m; ++i)
    b0[i] = pc[i].b;
  double sx = 0; // Solution value
  int j=0;  // column
  for ( ;(j<m); ++j ) {
    int x1=INT_MAX;
//    if (0==pc[j].b) x1 = 0;
    Pattern a1(pc.size()); // stack<IX> ?
    a1.SetObj(1);
    int i;
    for (i=j; i<m; ++i) {
      int a=L1/pc[i].l;
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
    {for_each_in (a1.ix,iix,Pattern::iterator)
      pc[iix->i].b -= x1*iix->x;} //Need sorted numeration
    sx += x1;
    // Restore original numeration:
    for_each_in (a1.ix,iix,Pattern::iterator)
      iix->i = pc[iix->i].i0;
    a1.MakeColumn(bas.Add());
  } // next j
//  UpdateHeurBnd(sx); // -- not in this form of FFD.
*/
} //____________________________________________________

void CSP2::GenCol(ColSet &cs,const d_vec &d) {
  GenColPure(cs,d);
}

class KP2G__FFD {
/////////////// INPUT ///////////////////
public:
  int L, W;
  Vector<int> l,w;
  Vector<int> b0;
  Vector<double> c; // prices
  Vector<int> ib; // The demand/price of piece i: b[ib[i]]
  int maxLevel;
// VARS:
private:
  int m,m2;
  Vector<int> b;
  struct Piece {
    int l,w,
      ii, // the original index as in input
      ib, // ib[i]
      i90; // this item in the rotated list
    void set(int a,int b,int c,int d,int e)
    { l=a; w=b; ii=c; ib=d; i90=e; }
    void rotate() { Swap(l,w); }
    bool operator<(const Piece& p) const
    { return w > p.w ? true:(w==p.w ? l>p.l: false); }
  };
  Vector<Piece> pc[2];
  Vector<int> lMin[2];
  int __maxLevel;
  int ori;
//////////////// OUTPUT //////////////////
public:
  double s,s1; // value for rotated/direct solution
  Vector<int> aa,aa1; // the solution vectors
  Vector<GGK2::IXY> pos,pos1;
private:
  void InitOptimize() {
    b=b0;
    aa.clear(); aa.resize(m2); pos.clear();
  }
  void InitSolve() {
    int i,ori;
    m=l.size();
    assert(m==w.size() && m==b0.size() && m==c.size());
    m2=ib.size(); m2=IMax(m2,m);
    if (ib.empty()) {
      ib.resize(m);
      for (i=0;i<m;++i) ib[i] = i;
    }
    for (i=0;i<m;++i) if (c[i]<1e-6) b0[i] = 0;
    pc[0].resize(m2); pc[1].resize(m2);
    for (ori=0;ori<2;++ori) for (i=0;i<m;++i)
      pc[ori][i].set(l[i],w[i],i,i, // the first m are the same
      0);
    for (ori=0;ori<2;++ori) for (i=m;i<m2;++i)
      pc[ori][i].set(w[ib[i]],l[ib[i]],i,ib[i], // for i>m, ib[i] is the unrotated
      0);
    for (i=0;i<m2;++i) pc[1][i].rotate();
    sort(pc[0].begin(),pc[0].end());
    sort(pc[1].begin(),pc[1].end());
    Vector<int> pci0[2];
    pci0[0].resize(m2); pci0[1].resize(m2);
    for (ori=0;ori<2;++ori) for (i=0;i<m2;++i)
      pci0[ori][pc[ori][i].ii] = i;
    for (ori=0;ori<2;++ori) for (i=0;i<m2;++i) {
      int i90 = pci0[1 - ori][i];
      while (i90)
        if (pc[1-ori][i90-1].w != pc[1-ori][i90].w)
          break; else  -- i90;
      pc[ori][pci0[ori][i]].i90 = i90; // for equal w's.
    }
    lMin[0].resize(m2); lMin[1].resize(m2);
    for (ori=0;ori<2;++ori) {
      int lM=INT_MAX;
      for (i=m2-1;i>=0;--i) {
        if (pc[ori][i].l < lM) lM=pc[ori][i].l;
        lMin[ori][i] = lM;
      }
    }
  }
public:
  KP2G__FFD() : maxLevel(0) { }
  void Reallocate(int m)
  { l.resize(m); w.resize(m); b0.resize(m); c.resize(m); }
  void Run() {
    InitSolve(); if (maxLevel<1) maxLevel = 3200000;
    InitOptimize(); __maxLevel = maxLevel;
    assert(pc[0][0].l<=L && pc[0][0].w <= W);
    PackStrip(0,0,0,L,W,0);
    aa1=aa; pos1=pos;
    InitOptimize(); __maxLevel = maxLevel+1;
    assert(pc[1][0].l<=W && pc[1][0].w <= L);
    PackStrip(1,0,0,W,L,0);
    s=s1=0; int i;
    for (i=0;i<m2;++i) s+=aa[i] * c[ib[i]];
    for (i=0;i<m2;++i) s1+=aa1[i] * c[ib[i]];
  }
private:
  void PackStrip(int level,size Y,size X,size L,size W,int i0) {
    int ori = level % 2;
    if (L < lMin[ori][i0]) return;
    size lLeft = L;
    BeginStrip();
//    AddPiece(Y,X,pc[ori][i0].iUnr);
    do {
      while (b[pc[ori][i0].ib] <= 0
        or pc[ori][i0].l > lLeft)
        if (++ i0 >= m2 or lLeft < lMin[ori][i0]) goto EndS;
      AddPiece(Y,X,pc[ori][i0].ii);
      if (level < __maxLevel) {
        size XS = X, YS = Y;
        (ori ? YS:XS) += pc[ori][i0].w;
        PackStrip(level+1, YS, XS,
        W-pc[ori][i0].w, pc[ori][i0].l, pc[ori][i0].i90);
      }
      lLeft -= pc[ori][i0].l;
      (ori ? X:Y) += pc[ori][i0].l; /// .....
    } while (not lLeft < lMin[ori][i0]);
EndS:    EndStrip();
  }
  void BeginStrip() { cout << " S"; } //AddPiece(-1); }
  void EndStrip() { cout << " E"; } //AddPiece(-2); }
  void AddPiece(size y,size x,int i) {
//    sol.push_back(Item(i));
//    if (i<0) return;
    ++ aa[i]; -- b[ib[i]]; // ORIENTATION ?
    pos.push_back(GGK2::IXY(i,y,x));
    cout << " p"<<i<< "("<<y<<','<<x<<")";
  }
};

void CSP2::GenColForHeur
  (Column *col,const d_vec &d,Vector<double> &b) {
  int i;
/* // WANG:
  for (i=0;i<m;++i)
    ggk2.c[i] = d[i];
// Cleaning info from the prev generation: (+assert)

  ggk2.Init(); ggk2.Run(); ggk2.ProduceColumn();

  col->clear();
  col->SetObj(1);

  Vector<int> bbb(m);
  for (i=0;i<m;++i) bbb[i] = (int)b[i];
if (ggk2.ConstrOK(bbb)) {
  // Filling *col:
  for (i=0;i<m;++i)
    if (ggk2.aa[i]) col->PushID(i,ggk2.aa[i]);
  // + Addi info: coords
  crd.push_back(ggk2.pos); // not too long
    // => not too much memory...
  col->GetAddiInfo() = (void*) &crd.back();

  dbg_outn_(3," EX:"<<ggk2.res);
  return;
}

  for (i=0;i<m;++i) {
    wang.c[i] = d[i];
    wang.b[i] = (int)b[i];
  }
  wang.Run();
  // Filling *col:
  for (i=0;i<wang.ips3->ix.size();++i)
    col->PushID(wang.ips3->ix[i].i,wang.ips3->ix[i].x);

  dbg_outn_(3," HEUR:"<<wang.ips3->InnerValue());

  crd.push_back(wang.pos);
  col->GetAddiInfo() = (void*) &crd.back();
*/
  KP2G__FFD alg;
  alg.Reallocate(m);
  alg.L = L; alg.W = W;
  for (i=0;i<m;++i) {
    alg.l[i] = pc[i].l;
    alg.w[i] = pc[i].w;
    alg.c[i] = d[i];
    alg.b0[i] = (int)b[i];
  }
  alg.Run();
  // Filling *col:
  col->clear();
  col->SetObj(1);
  Vector<int> * paa = alg.s > alg.s1 ? &alg.aa: &alg.aa1;
  for (i=0;i<paa->size();++i)
    if ((*paa)[i])
    col->PushID(i,(*paa)[i]);

  dbg_outn_(3," HEUR:"<<FMax(alg.s, alg.s1));

  crd.push_back(alg.s > alg.s1 ? alg.pos: alg.pos1);
  col->GetAddiInfo() = (void*) &crd.back();
}

// To be called only by LP
void CSP2::GenColPure
  (ColSet &cs,const d_vec &d) {
    if (0) {
      Column c;
      Vector<double> b(m);
      for (int i=0;i<m;++i) b[i] = pc[i].b;
      GenColForHeur(&c,d,b);
      double res = c.VecProd(d);
      if (res > 1+GetRCEps()) {
  clBest.clear();
  clBest = c;
  clBest.SetObj(1);
  dbg_outn_(3," -- Col! LPB:"<<GetCurrentLPBnd());

  cs.cs.insert(clBest); // this is the sign that a col found

  redCostBest = 1 - res;
  UpdateLagrBnd(d);  //  !!!!!!!!!!
  return;
      }
    }
  int i;
  for (i=0;i<m;++i)
    ggk2.c[i] = d[i];
// Cleaning info from the prev generation: (+assert)
  redCostBest = INFINITY__;
  cs.cs.clear();

  ggk2.Init(); ggk2.Run(); ggk2.ProduceColumn();

  // Filling clBest:
  clBest.clear();
  clBest.SetObj(1);
  for (i=0;i<m;++i)
    if (ggk2.aa[i]) clBest.PushID(i,ggk2.aa[i]);
  // + Addi info: coords
  crd.push_back(ggk2.pos);
  clBest.GetAddiInfo() = (void*) &crd.back();

  dbg_outn_(3," LPB:"<<GetCurrentLPBnd()<<" CG:"<<ggk2.res);

  if (ggk2.res < 1+GetRCEps()) return;

  cs.cs.insert(clBest); // this is the sign that a col found

  redCostBest = 1 - ggk2.res;
  UpdateLagrBnd(d);  //  !!!!!!!!!!
}

void CSP2::PrintColumn(ostream& os,Column *c) {
  int i;
//  if (c->GetCutSlackCut())
//    os << "Cut slack: " << c->GetCutSlackCoef()<< ' ';
      for (i=0; i<c->id.size(); ++i)
        os << pc[c->id[i].i].l <<'x'<< pc[c->id[i].i].w
          <<':'<<c->id[i].d<<' ';
      os << "obj " << c->GetObj();
}

LPCut * CSP2::GetLevelCut() {
  levelCut.lhs = GetLPBnd();
  return &levelCut;
}

////////////////////////////////////////////////////////
/////////// Integer Rounding ///////////////////////////
////////////////////////////////////////////////////////

bool CSP2::ConstructIntSol() {
  assert(pat.size() && lpx.size() && lpd.size()
    ); // Input
  Round();
  do {
    CompleteIntSol();
    if (Optimum()) break;
  } while (VaryResProblem()); // normally extending
  lpx.clear();
  lpd.clear();
  if (OUTP_LEV__ >=3) log__(" bnd:" << zi);
  return Optimum();
}

// double CSP2::GetXEps() { return GetBBEps(); }
void CSP2::SaveRoundedPart() {
  patBest.resize(pat.size()); patBest = pat;
  xiBest.resize(xi.size()); xiBest = xi;
}

  class J0Cmp { Vector<double> &xf; public:
    J0Cmp(Vector<double> &x_)  : xf(x_) { }
    bool operator () (int j1, int j2) const
    { return xf[j1] > xf[j2]; }
  };
bool CSP2::Round() {
  int i;
  xi.resize(lpx.size());
  xf.resize(lpx.size());    ///////// Rnd down:
  Vector<int> j0(xi.size()); // Indices; won't be used outs
  br.resize(m); GetRHS(br);
  fracused.clear(); fracused.reserve(xi.size());
  zr = 0;
  for (i=0;i<xi.size();++i) {
    xi[i] = floor(lpx[i]+GetXEps());
    xf[i] = lpx[i] - xi[i];
    j0[i] = i;   // Updating rhs:
    if (xi[i])
    for_each_in(pat[i].ix,iix,Pattern::iterator) {
      br[iix->i] -= xi[i] * iix->x;
      if (br[iix->i] < 0) br[iix->i] = 0;
    }
    zr += pat[i].GetObj() * xi[i];
  }
  if (OUTP_LEV__ >=3) log__(" RND dn:" << zr);
  sort (&j0[0], &j0[0]+j0.size(), J0Cmp(xf));
//  bool somefit;
//  int i0 = 0;
  i = 0;  nRPE = 0; // Init Variation of Residual
//////////////////////////////////////////////// Rnd up:
//  do {     somefit = false;  // for multiple
    do {
      bool fits = true;
      for_each_in(pat[j0[i]].ix,iix,Pattern::iterator)
        if (/*fabs(lpd[iix->i]) > GetDEps() // signif. demand
          and */ br[iix->i]  <  iix->x)
          { fits = false; break; }
      if (fits) {
        fracused.push_back(j0[i]);
        ++ xi[j0[i]];
        for_each_in(pat[j0[i]].ix,iix,Pattern::iterator) {
          br[iix->i] -= iix->x;
          if (br[iix->i] < 0) br[iix->i] = 0;
        }
        zr += pat[j0[i]].GetObj();
      }
    } while ( ++i < j0.size() );
  if (OUTP_LEV__ >=3) log__(" up:" << zr<<flush);
//  } while (somefit)
  return true;
} //////////////////////////////// CSP2::Round

void CSP2::GetRHS(Vector<double> &br) {
  br.resize(m); int i;
  for (i=0;i<m;++i) br[i] = pc[i].b;
} //////////////////////////////// CSP2:: GetRHS

//namespace {
class CSP2::SVC: public Alg {
/////// INPUT via constructor:
  CSP2 * pr;
  const Vector<double> & b0;
  const Vector<double> & d0;
  const double lb, ub; // a lower/upper bnd
//////// VARS:
  Vector<double> bc;  // the remaining demands
  d_vec d;
  const int m;
  unsigned int kkk; // total iteration number
  int k; // iteration number
  double zz; // current value
  double result; // best value
  int xMin; // intensity of current pattern
  double waste; // square can be large-...
  double SW;  // randomizer

  typedef Vector<Pattern> PatArray;
  PatArray pat;
  Vector<double> x;
public:
//////////// OUTPUT:
  PatArray patBest;
  Vector<double> xBest;

  SVC(CSP2* p_, Vector<double> & b_,
      Vector<double> & d_, const double lb_,
      const double ub_)
    : Alg(p_), pr(p_), m(pr->m), b0(b_), d0(d_),
    lb(lb_), ub(ub_), kkk(0) { d.resize(m); bc.resize(m); }
  double Run() {
//    try {
      Execute();
      assert(result >= lb);
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
  void Execute() {   result = 1e+100;
    k=0;
    ConvertValues();
    do {
      try {
        zz = 0; bc = b0; pat.clear(); x.clear();
        while (GenPat()) {
          CorrectValues();
        }
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
      d[i] = d0[i] * pr->L * pr->W;
      if (d[i] < 1) d[i] =1;
    }
  }
  bool GenPat() {
    pat.push_back(Pattern()); Column col;
    pr->GenColForHeur(&col,d,bc);
    pr->MakePattern(&pat.back(), &col);
    xMin = INT_MAX; // pattern intensity:
    waste = double(pr->L)*pr->W;  // original sorting ???
    Pattern::iterator iix;
    for_each_in(pat.back().ix,iix,) {
      assert(iix->x);
      assert(bc[iix->i] > 0);
      if (bc[iix->i] < double(xMin)*iix->x)
        xMin = int(double(bc[iix->i]) / iix->x); // floor
      waste -=
        double(pr->pc[iix->i].l) * pr->pc[iix->i].w
        * iix->x;
    }
    if (patternUseRatio > 0 and patternUseRatio <= 1)
      xMin = (int)ceil (double(xMin) * patternUseRatio);
    for_each_in(pat.back().ix,iix,)
       bc[iix->i] -= double(xMin)*iix->x;   // DECR. DMND
    zz += pat.back().GetObj() * xMin;
    x.push_back(xMin);

    double h = 0;  int i;
    for (i=0; i<m; ++i) h += pr->pc[i].l * bc[i];
    if (0 == h) { return false; } // NOTHING LEFT
    return true;
  }
  void CorrectValues() {
    SW = 1.1 + fabs(double( ++ kkk % 40 - 15)) / 10;
    assert ( waste < pr->L*pr->W );
    double LLh = double(pr->L)*pr->W
      / (double(pr->L)*pr->W - waste);
    for_each_in(pat.back().ix,iix,Pattern::iterator)
      d[iix->i] = ( LLh * iix->x
        * pr->pc[iix->i].l * pr->pc[iix->i].w
        + d[iix->i] * SW * (b0[iix->i] +bc[iix->i]))
        / (SW * (b0[iix->i] +bc[iix->i]) + iix->x);
  }
  void ControlSolution() { // the last one
    assert(pat.size() == x.size());
    Vector<double> s_aij(m); // fillled 0 ?
    double zzz = 0; int i;
    fill_n(&s_aij[0],m,0);
    for (i=0;i<pat.size(); ++i) {
      for_each_in(pat[i].ix,iix,Pattern::iterator)
        s_aij[iix->i] += iix->x * x[i];
      zzz += pat[i].GetObj() * x[i];
    }
    assertm(zzz == zz, "zzz="<<zzz<<", zz="<<zz);
    for (i=0;i<m;++i)
      assertm(s_aij[i] >=b0[i],"s_aij[i]="<<s_aij[i]<<", b0[i]"<<b0[i]);
}
  void SaveSolution()
  { patBest = pat; xBest = x; }
  void FillLastPattern() {
/*    pat.push_back(Pattern());
    int i;
    for (i=0;i<m;++i) if (bc[i] != 0)
      pat.back().PushIX(i,(int)bc[i]);
    pat.back().SetObj(1);
    zz += pat.back().GetObj();
    x.push_back(1);*/
  }
  static opt::OptContainer Options();
  static opt::OptSection opt;
  static int iterMax;
  static double outputLevel;
  static double patternUseRatio;
}; ////////////////////// class CSP2::SVC

int CSP2::SVC::iterMax=1;
double CSP2::SVC::patternUseRatio=0.5;
double CSP2::SVC::outputLevel=0;

opt::OptContainer CSP2::SVC::Options() {
  opt::OptContainer oc;
  oc
    << opt::MakeOpt(&iterMax, 1,
      "iterMax", 0)
    << opt::MakeOpt(&patternUseRatio, 0.5,
      "patternUseRatio",
      "Usage intensity of a new pattern rel. to max")
    << opt::MakeOpt(&outputLevel, DEF_OUTP_LEVEL,
      "outputLevel", "0-5");
  return oc;
} //____________________________________________________
opt::OptSection CSP2::SVC::opt
  ("CSP2_SVC", "SVC heur for the 2D CSP",
  SVC::Options(), opt::SolverCfg(), 5500);
//} ////////////////////// namespace { }

bool CSP2::CompleteIntSol() {
  int i;
  double LL=0;
///// Check if all fit in a single rod:
  for (i=0;i<m;++i) LL += (pc[i].l * br[i]);
  if (0 == LL) {
    UpdateHeurBnd(zr);
    SaveRoundedPart();
    return true;
  }
  if (OUTP_LEV__ >=4)
    log__(" Adding to rounded value: " << zr);
  SVC svc(this, br, lpd, GetLPBnd()-zr, GetHeurBnd()-zr);
  double SVCres=svc.Run();
  if (OUTP_LEV__ >=4) log__(" With SVC: " << zr+SVCres);
  if (zr + SVCres < zi) {
    SaveRoundedPart();
    AddResidualSolution(svc);
    UpdateHeurBnd(zr + SVCres);
    ControlSolution();
  }
//  if (Optimum())
  return true;
} //////////////////////////////// CSP2::CompleteIntSol

void CSP2::AddResidualSolution(CSP2::SVC & svc) {
  patBest.insert
    (patBest.end(),svc.patBest.begin(),svc.patBest.end());
  xiBest.insert
    (xiBest.end(),svc.xBest.begin(),svc.xBest.end());
}

// CSP 1D only:
void CSP2::ControlSolution() { // the best one
    assert(patBest.size() == xiBest.size());
    Vector<double> s_aij(m); // fillled 0 ?
    double zzz = 0; int i;
    fill_n(&s_aij[0],m,0);
    for (i=0;i<patBest.size(); ++i) {
      for_each_in(patBest[i].ix,iix,Pattern::iterator)
        s_aij[iix->i] += iix->x * xiBest[i];
      zzz += patBest[i].GetObj() * xiBest[i];
    }
    assert(zzz == zi);
    for (i=0;i<m;++i) assert(s_aij[i] >= pc[i].b);
}

bool CSP2::VaryResProblem() {
  if (++nRPE > RPEMax) return false;
  int i;
  if (fracused.empty()) { // No increased components
/*    for (i=xi.size()-1;i>=0;--i) // normally
      // consider j0 & not fracused...
      if (xi[i]) break;
    if (i<0)*/ return false; 
  } else
  { i = fracused.back(); fracused.pop_back(); }
  -- xi[i] ;  zr -= pat[i].GetObj();
  for_each_in(pat[i].ix, iix, Pattern::iterator)
    br[iix->i] += iix->x;
  return true;
} //////////////////////////////// CSP2::VaryResProblem

void CSP2::PrintBestSolution(ostream &os) {
  int i,j;
  os <<
    "%------ FROM HERE YOU CAN COPY TO A TeX FILE -------\n\n"
    "CSP2 instance name: "<<prName<<", No. "<<inst<<
    " from file "<<infile<<"\\\\\n"<<
    "L="<<L<<" W="<<W<<" m="<<m<<"\\\\\nPieces (l,w,b): ";
  for (i=0;i<m;++i)
    os << '('<<pc[i].l<<','<<pc[i].w<<','<<pc[i].b<<"), ";
  os<<"\n\nLP Bound: "<<GetLPBnd()
    //<<", LP Value: "<<GetLPValue()
    <<"\\\\\n"<<
  "Best solution value: "<<GetHeurBnd()<<"\n\n";
  os<<"Best solution:\\\\\n";
  os<<"\\unitlength"<<0.49/L<<"\\textwidth\n";
  os<<"%\\newcommand{\\printlabelSIZE}[1]{#1}\n"
    "% REDEFINE EMPTY to omit piece sizes in the pictures\n";

  int np=0;
  for (j=0;j<xiBest.size();++j) {
    Vector<GGK2::IXY> * pos
      = (Vector<GGK2::IXY>*)patBest[j].GetAddiInfo();
    if (!pos) continue;
    if (!xiBest[j]) continue;
    os << "\\begin{picture}("<<L<<","<<W<<")\n"
      "\\put("<<L<<","<<W<<"){\\makebox(0,0)[tr]"
      "{$\\times$"<<xiBest[j]<<"}}\n"
      "\\put(0,0){\\framebox("<<L<<","<<W<<"){}}\n";
    for_each_in(*pos,ip,Vector<GGK2::IXY>::iterator) {
      os << "\\put("<<ip->x<<","<<ip->y<<"){"
        "\\framebox("<<pc[ip->i].l<<","<<pc[ip->i].w
        <<"){\\footnotesize\\printlabelSIZE{"
        <<pc[ip->i].l<<"x"<<pc[ip->i].w
        <<"}}}\n";
    }
    os << "\\end{picture}\n";
    if (!(++np%2)) os << '\n';
  }
  os << '\n';
}

void CSP2::PrintProblem(ostream&os) {
  os <<
    "CSP2 instance name: "<<prName<<", No. "<<inst<<
    " from file "<<infile<<"\\\\\n"<<
    "L="<<L<<" W="<<W<<" m="<<m<<"\\\\\nPieces (l,w,b): ";
  for (int i=0;i<m;++i)
    os << '('<<pc[i].l<<'x'<<pc[i].w<<','<<pc[i].b<<"), ";
  os << endl;
}

bool CSP2::fSortPieces=true;
bool CSP2::fMergePieces=true;
double CSP2::WangGammaStep=0.8;
double CSP2::deps=1e-6;
int CSP2::RPEMax = 0;
double CSP2::WangNSol=1e+6;
double CSP2::outputLevel=DEF_OUTP_LEVEL;

opt::OptContainer CSP2::Options() {
  opt::OptContainer oc;
  oc
    << opt::MakeOpt(&fSortPieces, true,
      "fSortPieces", 0)
    << opt::MakeOpt(&fMergePieces, true,
      "fMergePieces",
      "Bool: Merge equal piece types")
    << opt::MakeOpt(&WangGammaStep, 0.8,
      "WangGammaStep",
      "Multiply Wang-Gamma if it's too high. <1.")
    << opt::MakeOpt(&WangNSol, 1e+6,
      "WangNSol",
      "FMax. number of combinations considered in Wang alg.")
    << opt::MakeOpt(&deps, 1e-6,
      "deps", "eps for dual multipliers")
    << opt::MakeOpt(&RPEMax, 0,
      "RPEMax", "FMax residual problem extensions. 2D: too long")
    << opt::MakeOpt(&outputLevel, DEF_OUTP_LEVEL,
      "outputLevel", "0-5");
  return oc;
} //____________________________________________________
opt::OptSection CSP2::opt
  ("CSP2", "The 2D Cutting Stock Problem",
  CSP2::Options(), opt::SolverCfg(), 500);

}                                                // SS_END_NAMESPACE__
