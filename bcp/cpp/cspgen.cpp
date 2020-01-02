// CUTTING STOCK PROBLEM RANDOM GENERATOR
// WITH A POLYNOMIAL SIZE DISTRIBUTION FUNCTION
// Gleb <Belov@math.tu-dresden.de>


#include <fstream>
#include <iostream>
#include <strstream>
#include <iomanip>
#include <cctype>
#include <cmath>
#include <list>
#include <vector>
#include <algorithm>
#include <cassert>
#include "random.h"


using namespace std;


char* GetLogo() {
  return
    "\nCUTTING STOCK PROBLEM RANDOM GENERATOR\n"
    "WITH A POLYNOMIAL SIZE DISTRIBUTION FUNCTION\n"
    "Gleb <Belov@math.tu-dresden.de>" __DATE__
    "\n\nArgument: param file (default: cspgen.cfg)\n"
    "Param file content: param value(s) param ...\n"
    "PARAMS: (default)\n"
    "N: (10) number of instances of each class\n"
    "m: (10) number of piece types\n"
    "d: (1) dimension\n"
    "L (10000), W (2000): stock sizes\n"
    "l1: followed by a few values, lower bounds of pieces lengths\n"
    " (10,100,1500,2500)\n"
    "l2: followed by a few values, upper bnds of l's\n"
    " (2000,3000,...9000)"
    "w1 (100): same as l1 for widths, w2 (1000,1500)\n"
    "b1 (1), b2 (100): lower, upper bnds of order demands\n"
    "f: (0) probability density of a size l is ~ l^(-f),\n"
    "    f=0 means the common uniform distribution\n"
    "r: (0.5) randomizer seed   "
    "p: path for output.";
}


int N=10;
int m=50;
int d=1;
int L=10000, W=2000;
int l1l__[] = {10,100,1500,2500};
list<int> l1l(l1l__,l1l__+sizeof(l1l__)/sizeof(int));
int l2l__[] = {2000,3000,4000,5000,6000,7000,8000,9000};
list<int> l2l(l2l__,l2l__+sizeof(l2l__)/sizeof(int));
int w1w__[] = {100};
list<int> w1w(w1w__,w1w__+sizeof(w1w__)/sizeof(int));
int w2w__[] = {1000,1500};
list<int> w2w(w2w__,w2w__+sizeof(w2w__)/sizeof(int));
int b1=1, b2=100;
double f=0;
double r=0.5;
char path[1024] = "";


void ReadList(istream &is,list<int> &ll) {
  ll.clear();
  int n;
  do {
    is >> n;
    if (!is) break;
    cout << "Putting value: " << n << endl;
    ll.push_back(n);
//    is >> skipws;
  } while (true);
  is.clear();
}


int ReadParams(const char* fln) {
  ifstream ifs(fln);
  if (!ifs) return 1;
  char token[16];
  do {
    ifs.clear();
    ifs >> setw(sizeof(token)) >> token;
    cout << "Read token: " << token << endl;
    if (not strlen(token) or ifs.eof())
      break;
    if (!strcmp(token,"N")) ifs >> N;
    if (!strcmp(token,"m")) ifs >> m;
    if (!strcmp(token,"d")) ifs >> d;
    if (!strcmp(token,"L")) ifs >> L;
    if (!strcmp(token,"W")) ifs >> W;
    if (!strcmp(token,"b1")) ifs >> b1;
    if (!strcmp(token,"b2")) ifs >> b2;
    if (!strcmp(token,"f")) ifs >> f;
    if (!strcmp(token,"r")) ifs >> r;
    if (!strcmp(token,"l1")) ReadList(ifs,l1l);
    if (!strcmp(token,"l2")) ReadList(ifs,l2l);
    if (!strcmp(token,"w1")) ReadList(ifs,w1w);
    if (!strcmp(token,"w2")) ReadList(ifs,w2w);
    if (!strcmp(token,"p")) ifs >> setw(sizeof(path)) >> path;
  } while (ifs);
  return ifs.good();
}


// the integral of s^f
double U(int s,double f) {
  return f==1 ? log(double(s)) : pow(s,1-f);
}
// its reverse
double Ur(double v,double f) {
  return f==1 ? exp(v) : pow(v,1/(1-f));
}


// This function returns the size corresp
// to the quantile p according to the prob. density
// proportional to size^(-f)
int Distr(int s1,int s2,double p,double f) {
  return (int)Ur(U(s1,f) + p*(U(s2,f) - U(s1,f)),f);
}


struct Piece {
  int l,w,b;
  bool operator<(const Piece& pc) const
  { return l > pc.l; } // INVERSE
};
vector<Piece> pc;
soplex::Random rnd, rnd2;


void Generate() {
  assert(m>0 && L>0 && W>0 && b2>0 && b1>0);
  assert(l1l.size() && l2l.size() && w1w.size() && w2w.size());
  int iN,i;
  char filename[1024];
  list<int>::iterator il1,il2,iw1,iw2;
  pc.resize(m);


  for (il1=l1l.begin(),iw1=w1w.begin();
    il1 != l1l.end();  ++ il1 ) {
    for (il2=l2l.begin(),iw2=w2w.begin();
      il2 != l2l.end(); ++ il2 ) {
      rnd.setSeed(r); // BEFORE EACH CLASS !!!
      double seed2
        = r*(*il1)/L*(*il2)/L*b1/b2*m/1000+f; // *d?
      if (d!=1) seed2 *= double(*iw1)/W*(*iw2)/W;
      rnd2.setSeed(seed2); // USE SETTINGS
      // FILENAME:
      ostrstream oss(filename,sizeof(filename));
      oss  << path
        << 'm'<<m
        <<'l'<<*il1<<'l'<<*il2<<'L'<<L;
      if (d!=1)
        oss <<'w'<<*iw1<<'w'<<*iw2<<'W'<<W;
      oss << 'b'<<b1<<'b'<<b2<<'f'<<f<<'r'<<r
        << ends;
      for (char *c=filename;*c;++c)
        if ('.' == *c) *c = '-';
      ofstream ofs(filename);
      if (!ofs) {
        cout << " Could not write: " << filename << endl;
        return;
      }
      for (iN=0;iN<N;++iN) {
        for (i=0;i<m;++i) {
          pc[i].l = Distr(*il1,*il2,rnd,f);
          if (d!=1) pc[i].w = Distr(*iw1,*iw2,rnd,f);
          pc[i].b = b1 + int(rnd * (b2-b1));
        }
        sort(pc.begin(),pc.end());
        // PRINTING
        ofs << "NN=" << iN+1
//          << "CSP-" << int(rnd2*INT_MAX)
          << '\n';
        ofs << L << ' ';
        if (d!=1) ofs << W << ' ';
        ofs << m << '\n';
        for (i=0;i<m;++i) {
          ofs << pc[i].l << ' ';
          if (d!=1) ofs << pc[i].w << ' ';
          ofs << pc[i].b << '\n';
        }
        ofs << '\n';
      }
      if ( ++ iw2==w2w.end() )
        -- iw2 ;
    }
    if (++ iw1==w1w.end())
      -- iw1 ;
  }
}


int main(int argc, char *argv[]) {
  char * cfgname;
  if (2<=argc)
    if ('-'==*argv[1] || '?'==*argv[1]) {
      cout << GetLogo() << endl;
      return 0;
    }
    else cfgname = argv[1];
  else cfgname = "cspgen.cfg";
  ReadParams(cfgname);
  Generate();
  return 0;
}
