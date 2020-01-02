#ifndef __RASTER_H__32
#define __RASTER_H__32


//#include <Vector>
using namespace std;


namespace bgn {


void ConstructRP
(Vector<int> &sz,double SZ,Vector<int> &rps);
void ConstructRP
(Vector<int> &sz,Vector<int> &ub, // with upper bounds
 double SZ,Vector<int> &rps);
void ConstructReducedRP
(Vector<int> &sz,double SZ,Vector<int> &rps);


// RETURNS number, not value
int FindRPUnder(int l,Vector<int> &rpl);


} // namespace ss


#endif // __RASTER_H__32
