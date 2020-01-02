#ifndef __GGK2_H__32
#define __GGK2_H__32


#include "raster.h"


namespace ss {


// GENERAL (:=non-staged...??) GUILLOTINE KNAPSACK 2D (GGK2)
// UNBOUNDED
class GGK2 {
public:
  ////////////// INPUT //////////////
  int m;
  Vector<int> l,w;
  int L,W;
  Vector<double> c;
  ////////////// OUTPUT /////////////
  Vector<Vector<double> > val;
  Vector<Vector<bool> > valOpt;
  class PtInfo {
    int w, // w=0: single piece; w=1/2 cut at l/w
      p, // piece id / cut pos
      p2; // pos of the 2nd part
  public:
    void set(int a,int b=0,int c=0) { w=a;p=b;p2=c; }
    void setPos(int b) { p=b; }
    int getPos() const { return p; }
    int getPos2() const { return p2; }
    int getW() const { return w; }
  };
  Vector<Vector<PtInfo> > what;
  double res;
  Vector<int> aa; // the pattern, non-sparse
  struct IXY
  {int i,x,y; IXY(int a,int b,int c) :i(a),x(b),y(c) {} IXY() {} };
  Vector<IXY> pos; // the positions


  Vector<int> rpl,rpw, rplForPiece, rpwForPiece;
  int nrpl, nrpw;


  ////////////// ACCESS ///////////////
public:
// CALL EACH TIME DATA CHANGES:
  void Init();
// CALL EACH TIME only prices change:
  void InitValues();
  void Run();
  void Reallocate();
  void ProduceColumn();
  bool ConstrOK(Vector<int> &b);
    // -- restore intensity Vector+crd
  void PrintTable(ostream&);
  void PrintBestSolution(ostream&);


  ////////////// PRIVATE //////////////
private:
  double recurs(int j, int k);
  void TryCutsAtLs
  (int j,int k,double &v,PtInfo&);
  void TryCutsAtWs
  (int j,int k,double &v,PtInfo&);
  void ReallocateWorkingData();
  void ConstructRP();
  void rec_restore(int j,int k,int x,int y);
};


} // ss


#endif // __GGK2_H__32
