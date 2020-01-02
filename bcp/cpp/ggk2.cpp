// GENERAL (non-staged) GUILLOTINE KNAPSACK 2D (GGK2.cpp)


// CHECK WHAT IF NO RP REDUCTION


#include "stdafx.h"
#include "ggk2.h"


namespace ss {


double GGK2::recurs(int j, int k) {
  if (valOpt[j][k]) return val[j][k];
  int w=0, // taking the single piece at this pos
    pos=what[j][k].getPos(); // its index
  TryCutsAtLs(j,k,val[j][k],what[j][k]);
  TryCutsAtWs(j,k,val[j][k],what[j][k]);
  valOpt[j][k] = 1;
  return val[j][k];
}


void GGK2::TryCutsAtLs
(int j,int k,double &v,PtInfo &pi) {
  int lmax=rpl[j]/2; double tv;
  int j1 = 1; // no endless recurs
  for (int j2=j;rpl[j1] <=lmax;++ j1) {
    int lr = rpl[j] - rpl[j1];
    while (rpl[j2] > lr)   -- j2;
    if ((tv=recurs(j1,k) + recurs(j2,k)) > v)
    { v=tv; pi.set(1,j1,j2); }
  }
}


void GGK2::TryCutsAtWs
(int j,int k,double &v,PtInfo &pi) {
  int wmax=rpw[k]/2;
  int k1 = 1; double tv;
  for (int k2=k;rpw[k1] <= wmax;++ k1) {
    int wr = rpw[k] - rpw[k1];
    // if (!wr) break; // bec. of wmax
    while (rpw[k2] > wr)   -- k2;
    if ((tv=recurs(j,k1) + recurs(j,k2)) > v)
    { v=tv; pi.set(2,k1,k2); }
  }
}


// CALL EACH TIME prices change:
// Setting table values for all cells fitting a piece:
void GGK2::InitValues() {
  int i,j,k;
  for (j=0;j<nrpl;++j) for (k=0;k<nrpw;++k) { // ZEROS
    val[j][k]=0;
    valOpt[j][k]=1; // For pos fitting no item
    what[j][k].set(3);
  }
// RP indices for l
  Vector<int> il(m); for (i=0;i<m;++i) il[i]=i;
  sort(il.begin(),il.end(),CmpByArray<int>(&l[0]));
  for (i=0,j=0;i<m;++j)
    if (rpl[j]>=l[il[i]]) { // FOR REDUCED RP:
      rplForPiece[il[i]] = j; // SOME L's not in rpl
      while (l[il[i]] <= rpl[j]) {
        rplForPiece[il[i]] = j; // SOME L's not in rpl
        if (++i>=m) break;
      }
    }
// RP indices for w
  Vector<int> iw(m); for (i=0;i<m;++i) iw[i]=i;
  sort(iw.begin(),iw.end(),CmpByArray<int>(&w[0]));
  for (i=0,k=0;i<m;++k)
    if (rpw[k]>=w[iw[i]]) {
      rpwForPiece[iw[i]] = k;
      while (w[iw[i]] <= rpw[k] && i<m) {
        rpwForPiece[iw[i]] = k;
        if (++i>=m) break;
      }
    }
/// SIMPLEST init values for single pieces:
  for (i=0;i<m;++i)
    for (j=rplForPiece[i];j<nrpl;++j)
      for (k=rpwForPiece[i];k<nrpw;++k)
        if (val[j][k] < c[i]) {
          val[j][k] = c[i];
          valOpt[j][k] = 0;
          what[j][k].set(0,i);
        }
}


void GGK2::Reallocate() {
  assert(m);
  l.resize(m); w.resize(m); c.resize(m);
}


void GGK2::ReallocateWorkingData() {
  assert(nrpl&&nrpw); int j;


  val.resize(nrpl);
  valOpt.resize(nrpl);
  what.resize(nrpl);
  for (j=0;j<nrpl;++j) {
    val[j].resize(nrpw);
    valOpt[j].resize(nrpw);
    what[j].resize(nrpw);
  }
  rplForPiece.resize(m);
  rpwForPiece.resize(m);
}


void GGK2::ConstructRP() {
  ss::ConstructReducedRP(l,L,rpl); nrpl = rpl.size();
  ss::ConstructReducedRP(w,W,rpw); nrpw = rpw.size();
}


// CALL EACH TIME DATA CHANGES
void GGK2::Init() {
  assert(l.size() == w.size());
  m = l.size();
  assert(m);
  assert(m==c.size());


  ConstructRP();
  ReallocateWorkingData();
  InitValues();
}


void GGK2::Run() {
  assert(l.size() == w.size());
  res = recurs(nrpl-1,nrpw-1);
}


void GGK2::rec_restore(int j,int k,int x,int y) {
  assert ( (j && k));
  int w = what[j][k].getW();
  if (3==w) return; // noth fits
  if (0==w) { // A real item
    int i;
    ++ aa[i=what[j][k].getPos()];
    pos.push_back(IXY(i,x,y));
    return;
  }
  if (1==w) { // HORIZ
    int j1=what[j][k].getPos();
    rec_restore(j1,k,x,y);
    int j2=what[j][k].getPos2();
    rec_restore(j2,k,x+rpl[j1],y);
  } else { // VERTIC
    int k1=what[j][k].getPos();
    rec_restore(j,k1,x,y);
    int k2=what[j][k].getPos2();
    rec_restore(j,k2,x,y+rpw[k1]);
  }
}
// MOVE this to an ancestor ???????
void GGK2::ProduceColumn() {
  pos.clear();
  aa.clear();
  aa.resize(m);
  rec_restore(nrpl-1,nrpw-1,0,0);
}


bool GGK2::ConstrOK(Vector<int> &b) {
  int i;
  for (i=0;i<m;++i)
    if (aa[i] > b[i]) return false;
  return true;
}


void GGK2::PrintTable(ostream&os) {
  int i,j;
  os << "w/l  ";
  for (i=0; i<nrpl; ++i) os << setw(4)<< rpl[i] << ' ';
  os << endl;
  for (j=0; j<nrpw; ++j) {
    os << setw(4) << rpw[j] << ' ';
    for (i=0;i<nrpl; ++i) os << setw(4)<< val[i][j] << ' ';
    os << endl;
  }
}


void GGK2::PrintBestSolution(ostream &os) {
  os << res << endl;
  int i;
  for (i=0;i<aa.size();++i) os << (aa)[i] << ' ';
  os << endl;
  for (i=0;i<pos.size();++i)
    os<< ' ' << pos[i].i << ':'<<pos[i].x<<','<<pos[i].y;
  os<<endl;
}


} // ss


//#define __TEST_GGK2


#ifdef __TEST_GGK2


int main() {
  int ll[] = {71,31,81};
  int ww[] = {31,41,51};
  double cc[]={8,5,16};
  int m=sizeof(ll)/sizeof(*ll);


  ss::GGK2 ggk2;
  ggk2.l.assign(ll,ll+m);
  ggk2.w.assign(ww,ww+m);
  ggk2.c.assign(cc,cc+m);
  ggk2.m=m; ggk2.L=160;ggk2.W=105;


  ggk2.Init(); ggk2.Run(); ggk2.ProduceColumn();


  ggk2.PrintTable(cout);
  ggk2.PrintBestSolution(cout);


  cin.get();
  return 0;
}


#endif
