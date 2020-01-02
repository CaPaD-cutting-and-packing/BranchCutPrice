#include <cmath>
#include <vector>
#include <iostream>
using namespace std;


#include "stdafx.h"
#include "raster.h"
#include "lasthdr.h"


namespace bgn {
namespace {


#define bpw (sizeof(unsigned)*CHAR_BIT) // bits per word
void CopyBits(unsigned *bf,double s,double SZ) {
 if(!s) return; // ZF-Schw.: beschreiben
 int k, a=int(s/(double)bpw), b=(int)fmod(s,bpw);
 int k2=int(SZ/(double)bpw)-a+1;
 if (b) // ShiftOr:
  for (k=k2;k; -- k)
    bf[k+a] |= ((bf[k]<<b)|(bf[k-1]>>(bpw-b)));
 else
  for (k=k2;k; -- k) bf[k+a] |= bf[k];
}


// No upper bnds
void ConstructRPBitField
(Vector<unsigned> &bf,Vector<int> &sz,double SZ) {
  bf.resize(int(SZ/(double)bpw)+4);
  // MemSet(bf,0,int((zi+128)/bpw));
  bf[0]=0; bf[1]=1; // To begin from bf[1] !!
  int i;
  for (i=0;i<sz.size();++i) {
   int i1=1, i2=int(SZ/sz[i]);
   while(i1*2-1<=i2)
   { CopyBits(&bf[0],double(sz[i])*i1,SZ); i1*=2; }
   CopyBits(&bf[0],double(sz[i])*(i2-i1+1),SZ);
  }
}


// With upper bnds
void ConstructRPBitField
(Vector<unsigned> &bf,
 Vector<int> &sz, Vector<int> &ub, double SZ) {
  bf.resize(int(SZ/(double)bpw)+4);
  // MemSet(bf,0,int((zi+128)/bpw));
  bf[0]=0; bf[1]=1; // To begin from bf[1] !!
  int i;
  for (i=0;i<sz.size();++i) {
   int i1=1, i2=IMin(int(SZ/sz[i]), ub[i]);
   while(i1*2-1<=i2)
   { CopyBits(&bf[0],double(sz[i])*i1,SZ); i1*=2; }
   CopyBits(&bf[0],double(sz[i])*(i2-i1+1),SZ);
  }
}


// what if l[0] == 1?
void ReduceRP(Vector<unsigned> &bf,double SZ) {
  unsigned * const bbff = &bf[1];
  int SZint=(int)SZ; // only smaller tables have sense?
  int j=SZint+1; // begin with pos (S+1)
  int i=0; // begin with the pos 0 to find the 1st j
  for (;i<=SZint; ++i)
    if (bbff[i/bpw]
      & (static_cast<unsigned>(1)<<(i%bpw)))
    if (j>SZint - i) {
      -- j; // PASS IT BY
      for (;j>SZint - i; --j) // CLEAN for j>SZ-i
        bbff[j/bpw]
        &= ~(static_cast<unsigned>(1)<<(j%bpw));
      while (not // SEARCH NEXT SMALLER RP
        (bbff[j/bpw]
        & (static_cast<unsigned>(1)<<(j%bpw))) )
        -- j;
    }
}


int CountBits(unsigned *bf,double SZ) {
  int n=0;
  unsigned * const bbff = &bf[1];
  for (int i=(int)SZ;i>=0; --i)
    if (bbff[i/bpw] & (unsigned(1)<<(i%bpw)))
      ++n;
  return n;
}


void ProduceLengths
(Vector<unsigned> &bf,double SZ,Vector<int> &rps) {
  rps.resize(CountBits(&bf[0],SZ));
  unsigned * const bbff = &bf[1]; int SZint=(int)SZ;
  int j=0;
  for (int i=0;i<=SZint; ++i)
    if (bbff[i/bpw] & (unsigned(1)<<(i%bpw)))
      rps[j++] = i;
}


void PrintBitField(Vector<unsigned> bf,double SZ) {
  unsigned * const bbff = &bf[1]; int SZint=(int)SZ;
  for (int i=0;i<=SZint; ++i)
    if (bbff[i/bpw] & (unsigned(1)<<(i%bpw)))
      cout << i << ' ';
  cout << endl;
}




} // namespace {}


// NO UPPER BNDS:
void ConstructRP
(Vector<int> &sz,double SZ,Vector<int> &rps) {
  Vector<unsigned> bf;
  ConstructRPBitField(bf,sz,SZ);
  ProduceLengths(bf,SZ,rps);
/*  cout << "RP: ";
  for (int i=0;i<rps.size();++i) cout << rps[i] << ' ';
  cout << endl;*/
}


// with UPPER BNDS:
void ConstructRP
(Vector<int> &sz, Vector<int> &ub,
  double SZ,Vector<int> &rps) {
  Vector<unsigned> bf;
  assert(sz.size() == ub.size());
  ConstructRPBitField(bf,sz,ub,SZ);
  ProduceLengths(bf,SZ,rps);
/*  cout << "RP: ";
  for (int i=0;i<rps.size();++i) cout << rps[i] << ' ';
  cout << endl;*/
}
// NO UPPER BNDS:
void ConstructReducedRP
(Vector<int> &sz,double SZ,Vector<int> &rps) {
  Vector<unsigned> bf;
  ConstructRPBitField(bf,sz,SZ);
//  PrintBitField(bf,SZ);
  ReduceRP(bf,SZ);
//  PrintBitField(bf,SZ);
  ProduceLengths(bf,SZ,rps);
/*  cout << "RP: ";
  for (int i=0;i<rps.size();++i) cout << rps[i] << ' ';
  cout << endl;*/
}


int FindRPUnder(int l,Vector<int> &rpl) {
  // LINEAR APPROX
  // BE SURE NO DOUBLED RPs
  // THE FIRST rp = 0?
  // CONVERGENCY? <= j2 yes, j1 yes.
  assert(rpl.size());
  if (rpl.back() == 0) return 0;
	if (l >= rpl.back()) return (int)rpl.size()-1;
	int j1=0, j2= (int)rpl.size()-1; int j; int tl;
  while((tl=rpl[j2]-rpl[j1]) > 0 && (j2-j1 > 1)) {
    j = j1 + int(double(l-rpl[j1]) / tl * (j2-j1));
    if (rpl[j] > l) j2=j;
    else if (rpl[j] < l) j1=j;
    else return j;
    if (rpl[j1+1]>l) return j1;
    ++ j1; // otherwise
  }
  return j1;
}


} // ss


//#define __TEST_RP


#ifdef __TEST_RP


int main() {
  int lll[] = { 6,71,31,81 }; // sorted?
  Vector<int> l(lll,lll+sizeof(lll)/sizeof(lll[0]));
  int L = 160;
  Vector<int> rpl;
  ss::ConstructRP(l,L,rpl);
  int l1=0;
  for (;l1<=160;l1+=10)
  cout << " RP<="<<l1<<": "
    << rpl[ss::FindRPUnder(l1,rpl)];
  cin.get();
  return 0;
}


#endif
