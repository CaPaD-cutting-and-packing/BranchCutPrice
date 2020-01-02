// probl_csp1.cpp: Implementierung der Klasse CSP1.
//
//////////////////////////////////////////////////////////////////////


// ALL VARS IN THE BEGINNING -- LIKE pascal.
// Output


// ATTENTION: with cuts, dimension = pr->Dim() = ma != m
// + cuts: cannot calc LR if early term.
  // => use DP for ini. solution
// notice ratio how many diff columns generated


#include "stdafx.h"
#include "probl_csp1.h"
#include "bb_mos.h"
#include "raster.h"
#include "solver.h"
#include "HP_ACVRP.h"
#include "HP_AFF.h"
#include "lasthdr.h"


#include "wrpcpx.h"


//ILOSTLBEGIN


namespace bgn {                     // SS_BEGIN_NAMESPACE__

CSP1::CSP1(Env *penv,
    const char *in,const int i,const char *n)
    :Problem(penv,in,i,n), m(0)  //, bb(new BB())
  { }
CSP1::CSP1(int ,Env *penv,
    const char *in,const int i,const char *n)
    :Problem(penv,in,i,n)  //, bb(new BB())
  { }

/// Hyperplane coefs for branching hyperplanes
/// VRP-hyperplanes: 1-based product indexation
double BrVRP::Calc__(Column* c) {
  assert(i<=j);
  assert(i>=0);
  if (c->IsSlack())
    return 0;
  if (0 == i) {
    if (c->id.front().i == j - 1) // -1 !!
      return 1;
    else
      return 0;
  }
  if (i != j) {
    for_each_in(c->id, iid, Column::iterator)
      if (iid+1 != c->id.end()) // not the last item
	if (iid->i == i-1 && (iid+1)->i == j-1)
	  return 1;
  }
  else // i == j
    for_each_in(c->id, iid, Column::iterator)
        if (iid->i == i-1) {
	  assert(iid->d != 0);
	  return (iid->d - 1); // ????? may be 0
	}
  return 0;
}

double Capacity::Calc__(Column* c) {
    double sum = 0; // !!!!!
    int j=0;
    assert(not S.empty());
    if (c->IsSlack())
      return 0;
    /// Considering each arc:
    if (S[c->id[0].i + 1]) // internal indexation 1..m
      ++ sum; // arc from node 0 into S?
    while (++ j < c->id.size())
      if (not S[(c->id)[j-1].i + 1]
        and S[(c->id)[j].i + 1])
	  ++ sum;
    return sum;
}

double BrAFF::Calc__(Column* c) {
  // assume sorted list of items--+++
    size p1=0;
    int i1;
    //assert();
    if (c->IsSlack())
      return 0;
    for (i1=0;i1<c->id.size();++i1)
      if (c->id[i1].i + 1 >= i) // indexation 1..m
        break;
      else p1 += l[c->id[i1].i] * c->id[i1].d;
    if (i1 < c->id.size()) // yes
    if (c->id[i1].i+1 == i && p >= p1) { //space
      int k = int((p-p1) / l[i-1]);
      if (k * l[i-1] == p-p1
        && k < c->id[i1].d)
	  return 1;
    }
    return 0;
}

void CSP1::BrOnHP(const Vector<Column*>& colpool,
    const Vector<ColId>& cols, const Vector<double>& lpx,
    LPCut* &hpU, LPCut* &hpL) {
    if (0 == nHP)
      BrOnHP_VRP(colpool, cols, lpx, hpU, hpL);
    else
      BrOnHP_AFF(colpool, cols, lpx, hpU, hpL);
}

/// VRP-hyperplanes: 1-based product indexation
void CSP1::BrOnHP_VRP(const Vector<Column*>& colpool,
    const Vector<ColId>& cols, const Vector<double>& lpx,
    LPCut* &hpU, LPCut* &hpL)
{ // Branching on Hyperplanes
  // OUTPUT: hpU, hpL
  int i,j,im=-1,jm=-1;
  float minDF = 1e+10f, valueS=-1e35f;
  int i1S=-1, i2S=-1, i3S=-1, iWhatS=0;

  /// INIT: (indexation 1..m)
  if (xVRP.empty()) {
    xVRP.resize(m+1);
    for (i=0;i<=m;++i)
      xVRP[i].resize(m+1);
  }
  /// CLEAR:
  for (i=0;i<=m;++i)
    fill_n(xVRP[i].begin(), m+1, 0);
  /// COMPUTE THE VRP VARIABLES:
  assert(cols.size() == lpx.size());
  for (j=0;j<cols.size();++j) {
    if (cols[j].j >= 0) { // not a cut slack
      if (not colpool[cols[j].j]->IsSlack()) {
        xVRP[0][colpool[cols[j].j]->id.front().i + 1]
	  += (float)lpx[j]; // arc {0,i}
        for_each_in(colpool[cols[j].j]->id,
	  iid,Column::iterator) {
	// the column contains only products, no further entries
	// (cf. multiple stock lengths)
	  if (iid->d > 1) // arc {i,i}
	    xVRP[iid->i + 1][iid->i + 1] += float(lpx[j] * (iid->d - 1));
	  if (iid+1 != (colpool[cols[j].j]->id.end()))
	    xVRP[iid->i + 1][(iid+1)->i + 1] += float(lpx[j]);
	} // arc {i,j}
      }
    }
  }

  /// FIND THE MOST FRACTIONAL VARIABLE:
  for (i=0;i<=m;++i)
  for (j=IMax(i,1);j<=m;++j)
  if (afVRP(xVRP[i][j]) < minDF) /* round(xij) */
  if (fabs(xVRP[i][j] - floor(xVRP[i][j] + 1e-3)) > 1e-3) {
    minDF = afVRP(xVRP[i][j]); // only fractional
    im = i; jm = j;
  }
  assert(im>=0 and jm>=0);
//  assert(frMax >= GetXEps());
/*
  /// Find the most fractional capacity subsets:
  /// 1. Fill "all incoming/outgoing" sums
  xVRPinc.clear(); xVRPout.clear();
  xVRPinc.resize(m+1); xVRPout.resize(m+1);
  for (i=0;i<=m;++i)
  for (j=i+1;j<=m;++j) { // j > i strictly
    xVRPinc[j] += xVRP[i][j];
    xVRPout[i] += xVRP[i][j];
  }
  /// 2. Compute 1-, 2-, and 3-item sets
  /// 2.a) 1-item sets
  /// ????????????? Divide frac by d(S) ????? Yet by |S|
  for (j=1;j<=m;++j)
    if (afVRP(tvS = xVRPinc[j]) < minDF)
    if (fabs(tvS - floor(tvS + 1e-6)) > 1e-6) {
      minDF = afVRP(valueS = tvS);
      iWhatS = 1;
      i1S = j;
    }
  /// 2.b) 2-item sets
  for (i=1;i<=m;++i) // i > 0
  for (j=i+1;j<=m;++j) // j > i strictly
    if (afVRP(tvS = xVRPinc[i] + xVRPinc[j] - xVRP[i][j])
        < minDF)
    if (fabs(tvS - floor(tvS + 1e-6)) > 1e-6) {
      minDF = afVRP(valueS = tvS); // div by 2!
      iWhatS = 2;
      i1S = i; i2S = j;
    }
*/
/*  /// 2.c) 3-item sets
  for (i=1;i<=m;++i) // i > 0
  for (j=i+1;j<=m;++j) // j > i strictly
  for (k=j+1;k<=m;++k) // k > j strictly
    if (afVRP(tvS = xVRPinc[i] + xVRPinc[j] + xVRPinc[k]
          - xVRP[i][j] - xVRP[i][k] - xVRP[j][k]) < minDF) {
      minDF = afVRP(valueS = tvS); // div by 2!
      iWhatS = 3;
      i1S = i; i2S = j; i3S = k;
    }
  /// 3. Compute m-1 and m-2 item sets
  /// 3.a) m-1 item sets

  for (j=1;j<=m;++j)
    if (afVRP(tvS = xVRPout[0] + xVRPout[j] - xVRP[0][j]) < minDF) {
      minDF = afVRP(valueS = tvS);
      iWhatS = 4;
      i1S = j;
    }
  /// 3.b) m-2 item sets
  for (i=1;i<=m;++i) // i > 0
  for (j=i+1;j<=m;++j) // j > i strictly
    if (afVRP(tvS = xVRPout[0] + xVRPout[i] + xVRPout[j]
         - xVRP[0][i] - xVRP[0][j] - xVRP[i][j]
	 ) < minDF) {
      minDF = afVRP(valueS = tvS);
      iWhatS = 5;
      i1S = i; i2S = j;
    }
*/
  /// Choose the best + statistics:
  if (0 == iWhatS) {
    hpU = new BrVRP(im, jm, (int)floor(xVRP[im][jm]), 1);
    hpL = new BrVRP(im, jm, (int)ceil(xVRP[im][jm]), -1);
    log_n_(2," BrOnVRP x["<<im<<','<<jm<<"] = "<<xVRP[im][jm]
      <<" frac: "<<smfrac(xVRP[im][jm]));
  } else { // capacity
    S.clear(); S.resize(m+1);
//    Vector<int> S1(m+1,1);
    assert(i1S>=0 and i2S>=0 and i3S>=0);
    switch( iWhatS ) {
       case 1: S[i1S] = 1; break;
       case 2: S[i1S] = S[i2S] = 1; break;
       case 3: S[i1S] = S[i2S] = S[i3S] = 1; break;
       case 4: S.clear(); S.resize(m+1,1); S[i1S] = 0; break;
       case 5: S.clear(); S.resize(m+1,1); S[i1S] = S[i2S] = 0; break;
       default: assert(0);
    }
    assert(valueS > -1e34);
    hpU = new Capacity(S, (int)floor(valueS), 1);
    hpL = new Capacity(S, (int)ceil(valueS), -1);
    log_n_(2," BrOnCapacity "<<iWhatS<<':'<<valueS);
  } // capacity
  /// Statistics:
  ++ iwstat[iWhatS];

}

void CSP1::BrOnHP_AFF(const Vector<Column*>& colpool,
    const Vector<ColId>& cols, const Vector<double>& lpx,
    LPCut* &hpU, LPCut* &hpL) {
  Vector<map<size,double> > xa(m);
  map<size,double>::iterator it;
  int i,j,k,im=-1;
  size p,pm;
  float minDF = 1e+10;

  /// COMPUTE THE AFF VARIABLES:
  assert(cols.size() == lpx.size());
  for (j=0;j<cols.size();++j) {
    if (cols[j].j >= 0) { // not a cut slack
      if (not colpool[cols[j].j]->IsSlack()) {
        p=0;
        for_each_in(colpool[cols[j].j]->id,
	  iid,Column::iterator) {
	// the column contains only products, no further entries
	// (cf. multiple stock lengths)
          for (k=0;k<iid->d;++k) {
            xa[iid->i][p] += lpx[j];
            p += l[iid->i];
          }
	}
        assert(p <= L);
      }
    }
  }

  for (i=0;i<m;++i)
  for_each_in(xa[i], it, ) {
    double x = it->second;
    if (fabs(x-floor(x+1e-6))>1e-6 and fabs(x-floor(x)-0.5) < minDF) {
      minDF = (float)fabs(x-floor(x)-0.5);
      im = i; pm = it->first;
    }
  }
  assert(im >= 0);
  assertm(pm <= INT_MAX,
          "This part only for 'int' sizes now.");
  hpU = new BrAFF(&l[0], im+1, (int)pm, (int)floor(xa[im][pm]), 1);
  hpL = new BrAFF(&l[0], im+1, (int)pm, (int)ceil(xa[im][pm]), -1);
  log_n_(2," BrOnAFF x["<<im+1<<','<<pm<<"] = "<<xa[im][pm]);

}

//// USER INDEXATION 1..m //////////////
bool CSP1::TempModifyLP(const Vector<Column*>& colpool,
    const Vector<ColId>& cols, const Vector<double>& lpx,
    const d_vec& lpd, list<LPCut*> & cuts1) {
  int i; double sum;
  Vector<int> RE;

  if (fInteractiveCapacity) {
    static int fStartNew = 1;

    cout << "\n\nCapacity cut creation procedure. "
      " Current LP value = "<<GetCurrentLPValue()<<'\n';
    if (fStartNew) goto NewS;
AskYN0: // just after trying a subset
    cout << "Are you happy to add this capacity cut? (y/n) ";
    cin.clear(); cin.ignore(32767, '\n');
    switch (cin.get()) {
    case 'y': fStartNew = 1; return false;
    case 'n': goto NewS;
    default:
      cout << "\nAre you Chinese?? Please answer y or n!_ ";
      goto AskYN0;
    }
Loop: // if not agree with current S
    cout << "\n\nCapacity cut creation procedure. "
      " Current LP value = "<<GetCurrentLPValue()<<'\n';
NewS:
    fStartNew = 0;
    cout << "Listing {i: l[i] b[i] d[i]}:\n";
    for (i=0;i<Dim();++i)
      cout << setw(0) << i+1 <<": "
        << setw(0) << pc[i].l << ' '
        << setw(0) << pc[i].b << ' ' // current rhs
//	<< setw(0) << ba[i] << ' '
	<< lpd[i] << ' ';
    if (not cuts->empty()) {
      cout << "Dual values for cuts:\n";
      CutList::iterator ic;
      i=0;
      for_each_in(*cuts,ic,) {
        cout << "Cut " << i+1 << ": " << lpd[Dim() + i] << '\n';
        ++ i;
      }
    }
    cout << "\nCurrently removed products from LP:\n";
    for (i=0;i<Dim();++i)
      if (pc[i].b != ba[i]) cout << i+1 << ' ';
    cout << "\nNow propose a new removed set: (space-separated, press Ctrl-D or enter 0 & press Enter for end)\n";
    RE.clear();
    cin.clear();
    while (cin) {
      cin >> i;
      if (Dim() < i) {
        cout << "Index " << i << " too big\n";
	goto Loop;
      }
      if (1 > i) break;
      if (cin) RE.push_back(i);
    }
//    if (RE.empty()) {
//      cout << "Nothing entered.\n";
//      goto Loop;
//    }
    sum = 0;
    for (i=0;i<RE.size();++i)
      sum += lpd[RE[i]-1] * pc[RE[i]-1].b;
    cout << "The dual value of the proposed set (";
    for (i=0;i<RE.size();++i)
      cout << RE[i] << ' ';
    cout << ") IN THE CURRENT LP DUAL VALUES is "
      << sum << "\nAre you happy to try this subset? (y/n) ";
AskYN:
    cin.clear();  cin.ignore(32767, '\n');
    switch (cin.get()) {
    case 'y': break;
    case 'n': goto Loop;
    default:
      cout << "\nAre you Chinese?? Please answer y or n!_ ";
      goto AskYN;
    }
  }
  else // non-interactive
  { /*
    static int fStartNew = 1;

    if (fStartNew) goto NewS;
//AskYN0: // just after trying a subset
  //  cout << "Are you happy to add this capacity cut? (y/n) ";
//    case 'y':
/////// TEST VIOLATION
    fStartNew = 1; return false;
//    case 'n': goto NewS;
Loop: // if not agree with current S
NewS:
    fStartNew = 0;
//    cout << "\nNow propose a new removed set: (space-separated, press Ctrl-D or enter 0 & press Enter for end)\n";
    RE.clear();
    sum = 0;
    for (i=0;i<RE.size();++i)
      sum += lpd[RE[i]-1] * pc[RE[i]-1].b;
    cout << "The dual value of the proposed set (";
    for (i=0;i<RE.size();++i)
      cout << RE[i] << ' ';
    cout << ") IN THE CURRENT LP DUAL VALUES is "
      << sum << "\nAre you happy to try this subset? (y/n) ";
//    case 'y': break;
//    case 'n': goto Loop;
*/  }
  // Happy:
    for (i=0;i<Dim();++i)
      ba[i] = pc[i].b;
    for (i=0;i<RE.size();++i)
      ba[RE[i]-1] = 0; // completely removing the type
    for (i=0;i<cuts->size();++i)
      if (((*cuts)[i])->Type() == 345)
        ba[Dim() + i] = 0;
    return true; // to try this subset
} // change b_cg also?


void CSP1::SeparateHP(const Vector<Column*>& colpool,
    const Vector<ColId>& cols, const Vector<double>& lpx,
    const d_vec& lpd, list<LPCut*> & cuts1) {
	
	
	/// SWITCHING OFF:
	return;
	
	
    int i;
    Vector<int> S(1+Dim()); // zeroed. DIM = m+1 !!
    for (i=0;i<Dim();++i)
      if (pc[i].b == ba[i]) // for those kept
        S[i+1] = 1;
      else
        ba[i] = pc[i].b; // restoring ba
    for (i=0;i<cuts->size();++i)
      if ((*cuts)[i]->Type() == 345)
        ba[Dim() + i] = (*cuts)[i] ->GetRHS();
			
		log_ln(" WARNING: CVRP Cuts are rarely efficient and can be slow to generate. "
			"Better disable cutting planes with AFF/HP branching."
		);

/// MIP ======================

  int j,im=-1;
//  float minDF = 1e+10, tvS, valueS;
//  int i1S, i2S, i3S, iWhatS=0;

  /// INIT: (indexation 1..m)
  if (xVRP.empty()) {
    xVRP.resize(m+1);
    for (i=0;i<=m;++i)
      xVRP[i].resize(m+1);
  }
  /// CLEAR:
  for (i=0;i<=m;++i)
    fill_n(xVRP[i].begin(), m+1, 0);
  /// COMPUTE THE VRP VARIABLES:
  assert(cols.size() == lpx.size());
  for (j=0;j<cols.size();++j) {
    if (cols[j].j >= 0) { // not a cut slack
      if (not colpool[cols[j].j]->IsSlack()) {
        xVRP[0][colpool[cols[j].j]->id.front().i + 1]
	  += (float)lpx[j]; // arc {0,i}
        for_each_in(colpool[cols[j].j]->id,
	  iid,Column::iterator) {
	// the column contains only products, no further entries
	// (cf. multiple stock lengths)
	  if (iid->d > 1) // arc {i,i}
	    xVRP[iid->i + 1][iid->i + 1] += float(lpx[j] * (iid->d - 1));
	  if (iid+1 != (colpool[cols[j].j]->id.end()))
	    xVRP[iid->i + 1][(iid+1)->i + 1] += (float)lpx[j];
	} // arc {i,j}
      }
    }
  }
  // INIT MIP::


/* The problem we are optimizing will have 3 rows, 4 columns
   and 9 nonzeros.  */

// int NUMROWS;//    3
// int NUMCOLS;//    4
// int NUMNZ;//      9


//int main (void) {
/* Declare pointers for the variables and arrays that will contain
   the data which define the LP problem.  The setproblemdata() routine
   allocates space for the problem data.  */

   char     probname[] = "cap";
//    int      numcols;
//    int      numrows;
//    int      objsen;
   Vector<double>   obj(m*m*2);
   Vector<double>   rhs(m*m*2);
   Vector<char>     sense(m*m*2);
   Vector<int>      matbeg(m*m*2);
   Vector<int>      matcnt(m*m*2);
   Vector<int>      matind(m*m*2);
   Vector<double>   matval(m*m*2);
   Vector<double>   lb(m*m*2);
   Vector<double>   ub(m*m*2);
   Vector<char>      ctype(m*m*2);
/*   obj      = (double *) malloc (NUMCOLS * sizeof(double));
   rhs      = (double *) malloc (NUMROWS * sizeof(double));
   sense    = (char *) malloc (NUMROWS * sizeof(char));
   matbeg   = (int *) malloc (NUMCOLS * sizeof(int));
   matcnt   = (int *) malloc (NUMCOLS * sizeof(int));
   matind   = (int *) malloc (NUMNZ * sizeof(int));
   matval   = (double *) malloc (NUMNZ * sizeof(double));
   lb       = (double *) malloc (NUMCOLS * sizeof(double));
   ub       = (double *) malloc (NUMCOLS * sizeof(double));
   ctype    = (char *) malloc (NUMCOLS * sizeof(char));
*/
   /* Fill in the data for the problem.  */

/* This function fills in the data structures for the mixed integer program:

      Maximize
       obj: x1 + 2 x2 + 3 x3 + x4
      Subject To
       c1: - x1 + x2 + x3 + 10x4  <= 20
       c2: x1 - 3 x2 + x3         <= 30
       c3:       x2       - 3.5x4  = 0
      Bounds
       0 <= x1 <= 40
       2 <= x4 <= 3
      Integers
        x4
      End
 */

//   Vector<Vector<double> > A(m+1);
   int col=0, pos;
   for (i=0;i<=m;++i) { // filling columns for w[ij]
     for (j=i+1;j<=m;++j) {
       pos = col; // !!!
       matbeg[col]=pos;
       matcnt[col]=1;
       matind[pos]=col; // identity
       matval[pos]=1;

       ++col; //  next col
     }
   }
   pos=col;
   for (im=0;im<=m;++im) {// the y[i]
    int p1=0;
    matbeg[col]=pos;
    matcnt[col]=m+(im>0);
    for (i=0;i<=m;++i) { // filling columns for y[i]
     for (j=i+1; j<=m; ++j) {
     if (i==im or j==im) {
			 if (pos>matind.size()) {
				 log_ln("\n\n BETTER SET maxCuts=0 with nBranch=0, CVRP cuts seem to have a bug\n\n");
				 return;
			 }
       matind[pos]=p1;
       matval[pos]=(i==im)?1:-1;

       ++ pos;
     }
      ++ p1;
     }
    } // the last constraint:
    if (im) {
       matind[pos]=p1;
       matval[pos]=l[im-1]*b[im-1];

       ++ pos; ++ p1;
    }
   ++ col;
   }
   col=0;
   for (i=0;i<=m;++i) { // filling columns for w[ij]
     for (j=i+1;j<=m;++j) {
       obj[col] = xVRP[i][j];
       lb[col] = 0; ub[col]=1;
       ctype[col] = 'C';
       ++col; //  next col
     }
   }
   for (i=0;i<=m;++i) { // filling columns for y[i]
       obj[col] = 0;
       lb[col] = 0; ub[col]=(i>0)?1:0;
       ctype[col] = 'I';
       ++col; //  next col
   }
   pos=0;
   for (i=0;i<=m;++i) { // rhs
     for (j=i+1;j<=m;++j) {
       sense[pos] = 'G';
       rhs[pos] = 0;
       ++pos; //  next col
     }
   }
   sense[pos] = 'G';
   rhs[pos] = 0; // for now;

   int      solstat;
   double   objval;
   Vector<double>   x(m*(m+1)/2+m+1),slack(m*(m+1)/2+1);


   CPXLPptr      lp = NULL;
   int           status;
//   int           i, j;
   int           cur_numrows, cur_numcols;

   /* Initialize the CPLEX environment */
   CPXENVptr env //WrpCpx::env;
			= CPXopenCPLEX (&status);
	 

   /* Turn on output to the screen */

   status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
   if ( status ) {
      fprintf (stderr,
               "Failure to turn on screen indicator, error %d.\n", status);
      goto TERMINATE;
   }

   /* Create the problem. */

   lp = CPXcreateprob (env, &status, probname);

   /* A returned pointer of NULL may mean that not enough memory
      was available or there was some other problem.  In the case of
      failure, an error message will have been written to the error
      channel from inside CPLEX.  In this example, the setting of
      the parameter CPX_PARAM_SCRIND causes the error message to
      appear on stdout.  */

   if ( lp == NULL ) {
      fprintf (stderr, "Failed to create LP.\n");
      goto TERMINATE;
   }

   /* Now copy the problem data into the lp */

   status = CPXcopylp (env, lp, m*(m+1)/2+m+1, m*(m+1)/2+1, CPX_MIN, &obj[0], &rhs[0],
                       &sense[0], &matbeg[0], &matcnt[0], &matind[0], &matval[0],
                       &lb[0], &ub[0], NULL);

   if ( status ) {
      fprintf (stderr, "Failed to copy problem data.\n");
      goto TERMINATE;
   }

   /* Now copy the ctype array */

   status = CPXcopyctype (env, lp, &ctype[0]);
   if ( status ) {
      fprintf (stderr, "Failed to copy ctype\n");
      goto TERMINATE;
   }

   /* Finally, write a copy of the problem to a file. */

   col=0;
   char nm[100];
   for (i=0;i<=m;++i) { // filling columns for w[ij]
     for (j=i+1;j<=m;++j) {
       sprintf( nm, "w%d_%d", i, j);
       status = CPXchgname (env, lp, 'c', col, nm);
       ++col; //  next col
     }
   }
   for (i=0;i<=m;++i) { // filling columns for y[i]
       sprintf( nm, "y%d", i);
       status = CPXchgname (env, lp, 'c', col, nm);
       ++col; //  next col
   }

   cout << "Starting MIPs";
   cin >> i;

  int M;
  for (M=1;M<=GetHeurBnd()-1;++M) {

   int ind = m*(m+1)/2;
   double val = double(M*L+1);
   status = CPXchgrhs (env, lp, 1, &ind, &val);
   /* Optimize the problem and obtain solution. */

   status = CPXwriteprob (env, lp, "mipex1.lp", NULL);
   if ( status ) {
      fprintf (stderr, "Failed to write LP to disk.\n");
      goto TERMINATE;
   }

   status = CPXmipopt (env, lp);
   if ( status ) {
      fprintf (stderr, "Failed to optimize MIP.\n");
      goto TERMINATE;
   }

   solstat = CPXgetstat (env, lp);

   /* Write the output to the screen. */

   printf ("\nSolution status = %d\n", solstat);

   status = CPXgetmipobjval (env, lp, &objval);
   if ( status ) {
      fprintf (stderr,"No MIP objective value available.  Exiting...\n");
      goto TERMINATE;
   }

   printf ("Solution value  = %f\n\n", objval);

   if (objval < M+1 - 1e-3)
     break;

  }
  if (M >= GetHeurBnd()) {
    log_ln("No Rnd Capacity violated!!! Break the program");
    cin >> M;
    return;
  }
   /* The size of the problem should be obtained by asking CPLEX what
      the actual size is, rather than using what was passed to CPXcopylp.
      cur_numrows and cur_numcols store the current number of rows and
      columns, respectively.  */

   cur_numrows = CPXgetnumrows (env, lp);
   cur_numcols = CPXgetnumcols (env, lp);

   status = CPXgetmipx (env, lp, &x[0], 0, cur_numcols-1);
   if ( status ) {
      fprintf (stderr, "Failed to get optimal integer x.\n");
      goto TERMINATE;
   }

   status = CPXgetmipslack (env, lp, &slack[0], 0, cur_numrows-1);
   if ( status ) {
      fprintf (stderr, "Failed to get optimal slack values.\n");
      goto TERMINATE;
   }

   for (i = 0; i < cur_numrows; i++) {
//      printf (" Sl[%d]=%10f", i, slack[i]);
   }

   col=0;
   for (i=0;i<=m;++i) { // filling columns for w[ij]
     for (j=i+1;j<=m;++j) {
//      printf (" w%d_%d=%10f", i,j , x[col]);
       ++col; //  next col
     }
   }
   printf("\n");
   for (i=0;i<=m;++i) { // filling columns for y[i]
      printf (" y%d=%10f", i , x[col]);
       ++col; //  next col
   }



TERMINATE:

   /* Free up the problem as allocated by CPXcreateprob, if necessary */

   if ( lp != NULL ) {
      status = CPXfreeprob (env, &lp);
      if ( status ) {
         fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
      }
   }

   status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
   if ( status ) {
      fprintf (stderr,
               "Failure to turn on screen indicator, error %d.\n", status);
      goto TERMINATE;
   }
   /* Free up the CPLEX environment, if necessary */
	 if ( env != NULL )
		 status = CPXcloseCPLEX (&env);
		 
   /* Free up the problem data arrays, if necessary. */


//}  /* END main */

/// -MIP ======================

    //dbg_outn (1.3,
    log__( "Separating HP: rhs = "
      << //GetCurrentLPBnd()
      M+1 << ". Enter the RHS to add actually: ");
    double rh;
    cin >> rh;
    Vector<int> SS (m+1);
    log__(" Adding with rhs " << rh);
    col= m*(m+1)/2+1;
    for (i=1;i<=m;++i) { // filling columns for y[i]
       SS[i] = (x[col] > 1e-6);
       printf(" %d",SS[i]);
       ++col; //  next col
   }
   printf("\n");

    cuts1.push_back(new Capacity(SS, rh));
      //GetCurrentLPBnd())); // Local!!
}


 // all involved cuts in! no branching columns
double Problem::CalcRedCost
  (Column *c, CutList *cuts, const d_vec &d) {
// assert that this col not in the forbidden set?
  double sum=c->GetObj();
  for_each_in(c->id,iid,Column::iterator)
    sum -= d[iid->i]*iid->d;
  CutList::iterator ic;
  for_each_in(*cuts,ic,) (*ic)->ClearNonRec();
  int i=Dim();
  for_each_in(*cuts,ic,)
    sum -= d[i++] * (*ic)->Calc__(c);
  return sum;
}


void Problem::MakeColumn(Column * c,Pattern * p) {
    c->clear();
    for_each_in(p->ix,iix,Pattern::iterator)
      c->PushID(iix->i,iix->x);
//    c.id.insert(c.id.end(),id.begin(),id.end());
    c->Sort(); // !!!
    c->ofc = p->ofc;
    c->GetAddiInfo() = p->GetAddiInfo();
}


void Problem::MakePattern(Pattern *p, Column *c) {
    p->clear();
    for_each_in(c->id,iid,Column::iterator)
      if (iid->i < Dim())
        p->PushIX(iid->i, (int)iid->d);
    p->SetObj(c->GetObj());
    p->GetAddiInfo() = c->GetAddiInfo();
} ////////////


void Problem::InitLP()
{ BBCuts::nStepsMin = 0; llrv = -1e100; /*local*/ }


void Problem::UpdateLagrBnd(const d_vec &d) {
  CalcLagrRel(d);
  if (lrv1 > llrv) {
    llrv = lrv1;
//    lrb = lrb1;
  }  // Here current as "instant", "last",
  // though we may need it like LPBnd for BPC.
}


void Problem::CalcLagrRel(const d_vec &d) {
  if (GetCurrentLPValue()<1e50)
    lrv1 = llrv =
      GetCurrentLPValue() * GetBestRedCost()
      + VectorProduct(
        d.begin(), d.end(), ba.begin()); // NEED RHS in b0
  else
    lrv1 = llrv = -1e100;
//  lrv1 = GetCurrentLPValue() / (1-GetBestRedCost());
//  if (vFarl > lrv1) lrv1 = vFarl;
}

void Problem::PrintColumn(ostream& os,Column* c) {
  int i;
//  if (c->GetCutSlackCut())
    //os << "Cut slack: " << c->GetCutSlackCoef()<< ' ';
/* sparse output:
      for (i=0; i<c->id.size(); ++i)
        os << c->id[i].i  <<':'<<c->id[i].d<<' ';
      os << "obj " << c->GetObj();
*/
// dense output:
   Vector<int> b(Dim());
      for (i=0; i<c->id.size(); ++i)
        b[c->id[i].i] = c->id[i].d;
      for (i=0; i<b.size(); ++i)
        os <<' ' << setw(2) << b[i];
      os << "  obj=" << c->GetObj();
}


void CSP1::Init() { // the object itself
  int i;
//  Vector<int> l(m);
  l.resize(m);

  for (i=0;i<m;++i) l[i] = (int)pc[i].l; // these 'int' are used in DP: TODO

  bb.Reallocate(m);
  BBCuts::Reallocate(m);

  //veps = deps;
  zi = +INFINITY__;
  lpb=llrv=lrv1= -INFINITY__; // lpb1 in FindPrimal?


  FillSigns();
  // + filling b:
  b.resize(m); // reserved FMax(100,m*2) or so (reall)
  ba.resize(m);
  b_cg.resize(m);
  s_libi = 0;
  for (i=0;i<m;++i) {
    b_cg[i] = ba[i] = b[i] = pc[i].b;
    s_libi += (pc[i].l * pc[i].b);
  }
  fill(iwstat, iwstat + sizeof(iwstat)/sizeof(int), 0);

  /// Item cost perturbation:
  pc_cost.clear(); pc_cost.resize(m); // fill 0
  Random rnd(cost_perturb_rndinit);
  for (i=0;i<m;++i)
    pc_cost[i] = float(-rnd * cost_perturb_max);

// Algorithms:
  bb.m = m; bb.L = L;
  BBCuts::m = m;
  BBCuts::mc = 0;
  BBCuts::d__.resize(m);
//  BBCuts::mc = mc;
  BBCuts::L = L; // effective ?
  for (i=0;i<m;++i) {
    bb.pieces__[i].l = pc[i].l;
    bb.pieces__[i].b = pc[i].b;
    BBCuts::pieces__[i].l = pc[i].l;
    BBCuts::pieces__[i].b = pc[i].b;
  }
  bb.zMin = 1;
  bb.zLowerInitial = -INFINITY__;
  bb.eps = GetBBEps();


  BBCuts::m = m;
//  BBCuts::mc = mc;
  BBCuts::L = L; // effective ?

  // FORBIDDEN:
  BBCuts::forbidden__ = bb.forbidden__ = forbidden;
}

void CSP1::Done() {
  spprc_dp::Clear();
  spprc_dp1::Clear();
  bkp_dp::Clear();
  if (accumulate(iwstat, iwstat + sizeof(iwstat)/sizeof(int), 0))
    log_n(1, "\nFreq. of CVRPi: " <<iwstat[0]<<' '<<iwstat[1]<<' '<<iwstat[2]<<' '
    <<iwstat[3]<<' '<<iwstat[4]<<' '<<iwstat[5]);
}


void CSP1::FillSigns() {
  validsign.resize(m);
  if (fModelEquality) // not yet, review FSS
    fill_n(&validsign[0],m,0); // Ax = b, optimal simplex multipliers arbitrary
  else
    fill_n(&validsign[0],m,1); // Ax >= b, optimal simplex multipliers non-neg
}




// returns 2 if I/O err, 1 if bad format, 0 else
int CSP1::Read(istream &ifs, long long L_)
{
// Old format for multiple stock:
  bool multi = (L_ < 0);  int M;
  L0 = (size)L_;
  size m_or_L; // for big numbers
  ifs >> m_or_L; //m0;
  if (multi) ifs >> M;
  if (!ifs) return 2;
  if (!multi && (m_or_L>L0) && !nLm) Swap(m_or_L,L0);
  if (!multi && 1==nLm) Swap(m_or_L,L0);
  m=m0=(int)m_or_L;
  if ((m0<2) or ((L0<3) and !multi))
  { PRINT_ERROR("Bad data."); return 1; }
  pc0.resize(m0);
  int i;
	int nSubst=0;
  for (i=0;i<m0;i++) {
    ifs >> pc0[i].l;             // later: >> pc0[i].b;
		if (!ifs) {
			PRINT_ERROR(infile<<", instance "<<inst
			<<": item "<<i<<", m="<<m0<<": bad format");
			return 1;
		}
		if ((pc0[i].l<=0)) { // b < 0 !!
			PRINT_ERROR(infile<<", instance "<<inst
			<<": l["<<i<<"] = "<<pc0[i].l<<", m="<<m0);
			return 1;
		}
		if ((pc0[i].l > L0) && !multi) {
			PRINT_ERROR(infile<<", instance "<<inst
			<<": too large piece "<<i<<": size " << pc0[i].l);
			return 1;
		}
		string str;
			getline(ifs, str);
		istringstream iss(str);
		iss >> pc0[i].b;
		if (!iss) {
			pc0[i].b = 1;
			if (1 == ++nSubst)
				cout << "\nSome demands missing/bad format, assumed 1" << endl;
		}
    if ((pc0[i].b<0)) { // b < 0 !!
      PRINT_ERROR(infile<<", instance "<<inst
        <<": b["<<i<<"]<0, m="<<m0);
      return 1;
    }
  }
  if (multi) {
    int dm;
    ifs >> L0 >> dm;
    for (i=1;i<M;++i) ifs >> dm >> dm; // skipping other lengths
  }
  if ((!ifs) && (!ifs.eof())) return 2;
  InitProblem();
  if (m<2) return 1; // But valid for M>1
  return 0;
}


void CSP1::InitProblem() // data, after input
{
  int i;
  // Sorting & compressing pieces' list (optionally)
  // + effective rod length ?
  pc = pc0;
  L = L0; // maybe calc. effective
  //for (int i=0;i<m0;++i)
  //  pc[i].i_ = i;   // original numbers
  if (true/*fSortPieces*/) {
    for (i=0;i<m0;++i)
      pc[i].w = (float)pc[i].l;    // weights for sorting
    sort(pc.begin(),pc.end(),greater<Piece>());
    if (fMergePieces
      && !fSP2Relaxation) { // not there!


      int i,j;
      for (j=0,i=1; i<m0; ++i) {
        if ((pc[i].l != pc[j].l)
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
  if (fEffectiveL) {
    if (L > INT_MAX) {
      cout << "Bin capacity behind the scope of 'int', cannot compute effective bin size.\n"
        "Use an equivalent instance with smaller sizes! Proceeding with L as is." << endl;
    }
    else
    {
      Vector<int> sz(m), bi(m), rp1;
      for (i=0;i<m;++i) sz[i] = (int)pc[i].l; // IMPORTANT: w's
      for (i=0;i<m;++i) bi[i] = (int)pc[i].b;
      log_n_(1.3, "Constructing raster points (option fEffectiveL)...");
      ConstructRP(sz,bi,(double)L,rp1);
      size L1 = rp1[FindRPUnder((int)L,rp1)];
      log_n(1.3, " Old L="<<L<<", effective L="<<L1<<", N raster points="<<rp1.size());
      if (L1 != L)
        L=L1;
    }
  }

}


// ALL VARS in the beginning like PASCAL


void CSP1::FFSBasis
  (const PieceContainer &pc__,ColumnList &bas)
{
  if (nStartBasis) {
    int i,im;
    if (2==nStartBasis) { // diag. matrix
      for (i=0;i<m;++i) {
        Column col;
        col.PushID(i, 1 //IMin(L/pc__[i].l, pc__[i].b)
          );
        col.SetObj(1);
        bas.Add(col);
      }
      return;
    }
    if (3==nStartBasis) // infeasible slacks only
      return;

    // SVC:
    if (4==nStartBasis or 5==nStartBasis) {
      lpd.resize(m);
      br.resize(m);
      for (i=0;i<m;++i) {
        lpd[i] = pow(pc__[i].l,SVC::pow_l);
        br[i] = fSP2Relaxation? 1 : pc__[i].b;
      }
      if (OUTP_LEV__ >=2) log__(" Ini basis: SVC... ");
      SVC svc(this, br, lpd, -1e100, 1e100);
//      int iterMax__ = SVC::iterMax;
//      if (fFastRounding) SVC::iterMax = IMin(SVC::iterMax,3);
      double SVCres=svc.Run();
//      SVC::iterMax = iterMax__;
      if (OUTP_LEV__ >=2) log__(" = " << SVCres);
      // SAVE:
      patBest.insert
        (patBest.end(),svc.patBest.begin(),svc.patBest.end());
      xiBest.insert
        (xiBest.end(),svc.xBest.begin(),svc.xBest.end());
      // BASIS:
      for (i=0;i<patBest.size();++i) {
        Column * pc = &(bas.Add());
        MakeColumn(pc,&patBest[i]);
      }
      if (4==nStartBasis) return; // otherwise greedy
    }

    // GREEDY:
    for (i=0;i<m;++i) {
      bb.pieces__[i].d = (double)pc__[i].l;
      bb.pieces__[i].b = fSP2Relaxation ? 1 : pc__[i].b;
    }
    for (im=0;im<m;++im)
    {
      for (i=0;i<im;++i) bb.pieces__[i].b = 0;
      bb.pieces__[im].b -= 1;
      bb.L = L - pc__[im].l;
      bb.zMin =
        bb.zLowerInitial = -INFINITY__;
      bb.eps = GetBBEps();

      ColSet cs; Column col;
      col.PushID(im,1); // the beginning piece is already there
//  col.PushID(m,pc[im].w); // the width constraint -- here not.
      bb.colsetRes = &cs;
      bb.colBest = &col; // not forget clBest

      bb.Run();
//  if (not (bb.found)) continue;

      col.SetObj(1);
      if (not col.id.empty())
        bas.Add(col);
    }
  // RESTORE:
    for (i=0;i<m;++i)
      bb.pieces__[i].b = pc__[i].b;
    bb.L = L;
    return;
  }
  // maybe sort pc with randomized weights
  // then use ->i0 to restore indices
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
    a1.SetObj(1);
    int i;
    for (i=j; i<m; ++i) {
      int a=int(L1/pc[i].l);
      if (fSP2Relaxation)
        if (a>1)
	  a=1;
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
    CSP1::MakeColumn(&bas.Add(),&a1);
  } // next j
//  UpdateHeurBnd(sx); // -- not in this form of FFD.
} //____________________________________________________


void CSP1::GenCol(ColSet &cs,const d_vec &d) {
  assert(cuts);
  assert(forbidden);
  if (cuts->size())
    if (cuts->size()==1
      and cuts->front()==GetLevelCut())
        GenColWithCuts(cs,d); // Could be a spec proc
    else
      if (cuts->front()->Type() == -1 // BrVRP
        or cuts->front()->Type() == 345) // CutVRP
        GenColWithHP_VRP(cs,d); // hyperplanes
      else
        if (cuts->front()->Type() == -2)
	  GenColWithHP_AFF(cs,d);
	else
          GenColWithCuts(cs,d);
  else
    GenColPure(cs,d);
}


void CSP1::GenColForHeur
  (Column *col,const d_vec &d,Vector<double> &b) {
  int i, j;

if (0 == fCompareBySPPRC_DP) {
  bb.fHeur = true;
  for (i=0;i<m;++i) {
    bb.pieces__[i].d = d[i];
    bb.pieces__[i].b = fSP2Relaxation ? IMin(1,int(b[i]))
      : int(b[i]);
  }
  bb.zMin = -INFINITY__;
  bb.zLowerInitial = -INFINITY__;
  bb.eps = GetBBEps();


  ColSet cs; // dummy
  bb.colsetRes = &cs;
  bb.colBest = col;


  bb.Run();
// Addi info to all cols:
  col->SetObj(1);
}
else {
    static Vector<int> bbb;
    bbb.resize(m);
    static Vector<int> xxx;

    assertm(L <=INT_MAX, // TODO
            "Cannot use sizes above 'int' in DP at the moment, choose B&B-pricing");
    if (bkp_dp::d.empty()) // l remains filled !!!
      bkp_dp::InitBasicData(m,(int)L,&(l.front()));

    log_n_(3,"dr");
    if (fSP2Relaxation)
      for (i=0;i<m;++i)
        bbb[i] = IMin(1,int(b[i]));
    else
      for (i=0;i<m;++i)
        bbb[i] = int(b[i]);

    /// The duals: (indexation: 1..m)
    for (j=1;j<=m;++j)
      bkp_dp::d[j] = d[j-1];
    /// No arcs {ij}

    bkp_dp::Run(&(bbb.front()));

    bkp_dp::RestoreSolution((int)L, xxx); // to check

    col->clear();
    col->SetObj(1);

  for (i=1;i<=m;++i)
    if (xxx[i] > 0)
      col->PushID(i - 1, xxx[i]);

	}
	static int cnt=0;
	if (0 == ++cnt % 100)
	log_n_(1.2,'|');
}

// To be called in each B&B node
/*void CSP1::CopyForbCols(list<Column*> & fc) {
  list<Column*>::iterator ifc;
  bb.forbidden__.clear();
  BBCuts::forbidden__.clear();
  list<Vector<IX> > * pfc=
//    (cuts->size()) ?   &(BBCuts::forbidden__) :
      &(bb.forbidden__);
  Vector<IX> vix;
  for_each_in(fc,ifc,) {
    vix.clear();
//    vix.reserve((*ifc)->id.size());
    for_each_in((*ifc)->id,iid,Column::iterator) {
      vix.push_back(IX(iid->i,iid->d));
      assert(iid->i < Dim());
    }
    int i;
    for (i=1;i<(*ifc)->id.size();++i)
      assert((*ifc)->id[i].i != (*ifc)->id[i-1].i);
    pfc->push_back(vix);
  }
  BBCuts::forbidden__ = *pfc;
}*/


// To be called only by LP
void CSP1::GenColPure
  (ColSet &cs,const d_vec &d) {
  int i, j;

	log_n_(3, " MAX DUAL VAL = " << *max_element(d.begin(), d.end()) << ", deps = " << GetDEps() << ' ');
	
	static int cnt=0;   // for printouts
	++cnt;
	
if (0 == fCompareBySPPRC_DP or 2 == fCompareBySPPRC_DP
 or fLocalUB) { // because then DP imposs.
  bb.fHeur = false;
  for (i=0;i<m;++i) {
    bb.pieces__[i].d = d[i];
    bb.pieces__[i].b = fSP2Relaxation ? IMin(1,int(b_cg[i]))
      : int(b_cg[i]);
  }
  bb.zMin = 1;
  bb.zLowerInitial = -INFINITY__;
	bb.eps = GetBBEps();
	bb.deps = GetDEps();
	bb.colsetRes = &cs;
  bb.colBest = &clBest;
  clBest.clear();
// Cleaning info from the prev generation: (+assert)
  redCostBest = - INFINITY__;
  cs.cs.clear();

	if (0 == cnt % 100)
		log_n_(1.2,'b');
  bb.Run();
//   log_n_(1.2,'|');
  if (not (bb.found)) return; // what if eterm?

  log_n_(3, " BB RC = " << (1 - bb.z) << ' ');

// Addi info to all cols:
  for_each_in(cs.cs,it,ColSet::iterator)
    ((Column*)(&*it))->SetObj(1); // GNU iterators are const
  if (not bb.fETerm)
    redCostBest = 1 - bb.z;
  UpdateLagrBnd(d);  //  !!!!!!!!!!
}

  /// Control: solving by spprc_dp
  if ((not fLocalUB and fCompareBySPPRC_DP)
// 		and forbidden->empty()             
	) {                       // DP: not with forbidden columns.. But need to check it differently TODO
    static Vector<int> bbb;
    bbb.resize(m);
    static Vector<int> xxx;

    assertm(L <= INT_MAX, // TODO
            "Cannot use sizes above 'int' in DP at the moment, choose B&B-pricing");
		bkp_dp::deps = GetDEps();
    if (bkp_dp::d.empty()) // l remains filled !!!
      bkp_dp::InitBasicData(m,(int)L,&(l.front()));

		if (0 == cnt % 100)
			log_n_(1.2,'d');
    if (fSP2Relaxation)
      for (i=0;i<m;++i)
        bbb[i] = IMin(1,int(b_cg[i]));
    else
      for (i=0;i<m;++i)
        bbb[i] = int(b_cg[i]);

    /// The duals: (indexation: 1..m)
    for (j=1;j<=m;++j)
      bkp_dp::d[j] = d[j-1];
    /// No arcs {ij}

    bkp_dp::Run(&(bbb.front()));

		bkp_dp::RestoreSolution((int)L, xxx); // to check

    if (2 == fCompareBySPPRC_DP) // BB performed also
			assertm(fabs(bb.z - bkp_dp::GetObjValue((int)L))<1e-6,
      "B&B col.gen.="<<bb.z<<", BKP_DP col.gen.="
			<<  bkp_dp::GetObjValue((int)L));

    clBest.clear();
    clBest.SetObj(1);

  for (i=1;i<=m;++i)
    if (xxx[i] > 0)
      clBest.PushID(i - 1, xxx[i]);
  // clBest is sorted
  cs.cs.insert(clBest);

	redCostBest = 1 - bkp_dp::GetObjValue((int)L);

//     log_n_(1.2,'|');
  }
}


void CSP1::GenColWithCuts(ColSet &cs,const d_vec &d)
{
  BBCuts::m = m;
  BBCuts::L = L;
  BBCuts::foundInit=false;
  int i;
  BBCuts::pieces__.resize(m);
  for (i=0;i<m;++i) {
    BBCuts::pieces__[i].d = d[i];
    BBCuts::pieces__[i].b = int(b_cg[i]);
  }
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
  // setting up the list of involved cuts:
  BBCuts::cuts__.AddInvolved(cuts);


  BBCuts::d__ = d;


  BBCuts::zMin = 1; // + deps ?
  BBCuts::zLowerInitial = 0;
  BBCuts::eps = GetBBEps();
  BBCuts::nStepsMin = FMax( // increasing with ... ??
    BBCuts::nStepsMin, nStepsMin0)
    * nStepsMinInc; // ... each col gen!


  BBCuts::colsetRes = &cs;
  BBCuts::colBest = &clBest;
  clBest.clear(); // Make it empty to init CG const terms
  clBest.SetObj(1); // Necessary for the level cut
  // calculation during generation
// Cleaning info from the prev generation: (+assert)
  redCostBest =  - 1e100;
  cs.cs.clear();


  log_n_(1.2,'B');
  BBCuts::Init();
  BBCuts::Run();
  log_n_(1.2,'|');
  if (not (BBCuts::found)) return;


// Addi info to all cols:
  for_each_in(cs.cs,it,ColSet::iterator)
    ((Column*)(&*it))->SetObj(1); // GNU iterators are const
  if (!BBCuts::fETerm)
    redCostBest = 1 - BBCuts::z;
  UpdateLagrBnd(d);  //  !!!!!!!!!!
} //____________________________________________________

void CSP1::GenColWithHP_VRP(ColSet &cs,const d_vec &d)
{
  int i,j,k;
  static Vector<int> bbb;
  bbb.resize(m);
  static Vector<int> xxx, xxx1;

  if (fSP2Relaxation) // GetCG_RHS() !!!!!
    for (i=0;i<m;++i)
      bbb[i] = IMin(1,int(b_cg[i]));
  else
    for (i=0;i<m;++i)
      bbb[i] = int(b_cg[i]);

  clBest.clear(); // Make it empty to init CG const terms
  clBest.SetObj(1); // Necessary for the level cut
  // calculation during generation
// Cleaning info from the prev generation: (+assert)
  redCostBest =  - 1e100;
  cs.cs.clear();

/*
  if (spprc_dp::d.empty())  // workes only if Clear()ed
    spprc_dp::InitBasicData(m,L,&(l.front()));

  /// The duals: (in spprc_dp, INDEXATiON: 1..m)
  for (i=0;i<=m;++i)
  for (j=IMax(1,i);j<=m;++j)
    spprc_dp::d[i][j] = d[j-1];
  /// Adding arcs {ij}:
  for (k=0;k<cuts->size();++k)
  if ((*cuts)[k]->Type() == -1) {
    BrVRP* php = dynamic_cast<BrVRP*>((*cuts)[k]);
    spprc_dp::d[php->i][php->j] += d[Dim() + k];
  }
  else
  if ((*cuts)[k]->Type() == 345) { // Capacity
    Capacity* php = dynamic_cast<Capacity*>((*cuts)[k]);
    for (i=0;i<=m;++i)
      if (not php->S[i])
        for (j=IMax(1,i);j<=m;++j)
          if (php->S[j])
            spprc_dp::d[i][j] += d[Dim() + k];
  }


  log_n_(1.2, 'd');
  spprc_dp::Run(&(bbb.front()));

  spprc_dp::RestoreSolution(L, xxx);
  log_n_(1.2,'|');

  for (i=1;i<=m;++i)
    if (xxx[i] > 0)
      clBest.PushID(i - 1, xxx[i]);
  // clBest is sorted
  cs.cs.insert(clBest);

  redCostBest = 1 - spprc_dp::GetObjValue(L);
  UpdateLagrBnd(d);  //  !!!!!!!!!!
*/

/////////////////// =======================================

  int nModified = 0;

  assertm(L <=INT_MAX, // TODO
       "Cannot use sizes above 'int' in SPP-DP at the moment, use an equivalent instance");
  if (spprc_dp1::d.empty())  // l remains filled !!!
		spprc_dp1::InitBasicData(m,(int)L,&(l.front()));

  /// The duals: (indexation: 1..m)
  for (i=1;i<=m;++i) {
    spprc_dp1::d[i] = d[i-1];
    spprc_dp1::dI[i].clear();
  }
  /// Adding arcs {ij}: (in the cuts, they are ij; in dp1, ji!!!)
  for (k=0;k<cuts->size();++k)
  if ((*cuts)[k]->Type() == -1) {
    BrVRP* php = dynamic_cast<BrVRP*>((*cuts)[k]);
    if (fabs(d[Dim() + k]) > GetDEps()) {
      spprc_dp1::dI[php->j][php->i] += d[Dim() + k];
      log_n_(2.5," d["<<php->i<<','<<php->j<<"]+="<< d[Dim() + k]);
      ++ nModified;
    }
  }
  else // CAPACITY d's NOT READY here !!!!!!!!! ++++++++++++++++
  if ((*cuts)[k]->Type() == 345) { // Capacity
    Capacity* php = dynamic_cast<Capacity*>((*cuts)[k]);
    for (i=0;i<=m;++i)
      if (not php->S[i])
        for (j=i+1;j<=m;++j)
          if (php->S[j]) {
            spprc_dp1::dI[j][i] += d[Dim() + k];
            log_n_(2.5," dc["<<i<<','<<j<<"]+="<< d[Dim() + k]);
            ++ nModified;
          }
  }


//  if (nModified)
    log_n_(3, " N modified CVRP cuts="<<nModified<<'/'<<cuts->size());
  spprc_dp1::Run(&(bbb.front()));

	spprc_dp1::RestoreSolution((int)L, xxx1);
  log_n_(1.2,'|');
//////////////////////// ------------------------------

/*
//  if (fVRPFullGraph && fVRPSmallGraph)
  if (fabs(spprc_dp::GetObjValue(L) - spprc_dp1::GetObjValue(L))
    > GetDEps()) {
    log_ln("\nspprc_dp::GetObjValue(L)="<<spprc_dp::GetObjValue(L)<<
      "\nspprc_dp1::GetObjValue(L)="<<spprc_dp1::GetObjValue(L));
    log_ln("tuples (i,d,l,b) for i=0..m-1:");
    int i;
    for (i=0;i<m;++i)
      log__(" ("<<i<<','<<d[i]<<','<<l[i]<<','<<b[i]<<')');
    assert(0);
  }
*/
  for (i=1;i<=m;++i)
    if (xxx1[i] > 0)
      clBest.PushID(i - 1, xxx1[i]);
  // clBest is sorted
  cs.cs.insert(clBest);

	redCostBest = 1 - spprc_dp1::GetObjValue((int)L);
  UpdateLagrBnd(d);  //  !!!!!!!!!!
} //____________________________________________________

void CSP1::GenColWithHP_AFF(ColSet &cs,const d_vec &d)
{
  int i,k, nModified=0;
  static Vector<int> bbb;
  bbb.resize(m);
  static Vector<int> xxx1;

  if (fSP2Relaxation) // GetCG_RHS() !!!!!
    for (i=0;i<m;++i)
      bbb[i] = IMin(1,int(b_cg[i]));
  else
    for (i=0;i<m;++i)
      bbb[i] = int(b_cg[i]);

  clBest.clear(); // Make it empty to init CG const terms
  clBest.SetObj(1); // Necessary for the level cut
  // calculation during generation
// Cleaning info from the prev generation: (+assert)
  redCostBest =  - 1e100;
  cs.cs.clear();

  assertm(L <=INT_MAX, // TODO
      "Cannot use sizes above 'int' in DP at the moment, use an equivalent instance");
  if (bkp_dp::d.empty()) // l remains filled !!!
		bkp_dp::InitBasicData(m,(int)L,&(l.front()));

  /// The duals: (indexation: 1..m)
  for (i=1;i<=m;++i) {
    bkp_dp::d[i] = d[i-1];
    bkp_dp::pv[i].clear();
  }
  /// Adding arcs {ij}: (in the cuts, they are ij; in dp1, ji!!!)
  for (k=0;k<cuts->size();++k)
  if ((*cuts)[k]->Type() == -2) {
    BrAFF* php = dynamic_cast<BrAFF*>((*cuts)[k]);
    if (fabs(d[Dim() + k]) > GetDEps()) {
      bkp_dp::pv[php->i].push_back // duplicates ???
        (make_pair(php->p, d[Dim() + k])); // only the increase!
      log_n_(2.5," d["<<php->i<<','<<php->p<<"]+="<< d[Dim() + k]);
      ++ nModified;
    }
  }

  log_n_(3, " N modified AFF cuts="<<nModified<<'/'<<cuts->size());
  bkp_dp::Run(&(bbb.front()));

	bkp_dp::RestoreSolution((int)L, xxx1);
  log_n_(1.2,'|');

  for (i=1;i<=m;++i)
    if (xxx1[i] > 0)
      clBest.PushID(i - 1, xxx1[i]);
  // clBest is sorted
  cs.cs.insert(clBest);

  /// Should be done centrally:
  assert(fabs(CalcRedCost(&clBest, cuts, d)
	- (1-bkp_dp::GetObjValue((int)L))) < GetRCEps());

	redCostBest = 1 - bkp_dp::GetObjValue((int)L);
  UpdateLagrBnd(d);  //  !!!!!!!!!!
}

void CSP1::PrintColumn(ostream& os,Column* c) {
  int i;
//  if (c->GetCutSlackCut())
    //os << "Cut slack: " << c->GetCutSlackCoef()<< ' ';
/*      for (i=0; i<c->id.size(); ++i)
        os << c->id[i].i
          <<':'<<c->id[i].d<<' ';
      os << "obj " << c->GetObj(); */
//  sparse output:
      for (i=0; i<c->id.size(); ++i)
        os << c->id[i].i  <<':'<<c->id[i].d<<' ';
      os << "obj " << c->GetObj();

// dense output:
//    Vector<int> b(Dim());
//       for (i=0; i<c->id.size(); ++i)
//         b[c->id[i].i] = c->id[i].d;
//       for (i=0; i<b.size(); ++i)
//         os <<' ' << setw(2) << b[i];
//       os << "  obj=" << c->GetObj();
}


LPCut * CSP1::GetLevelCut() {
  levelCut.lhs = GetLocalLPBnd();
  return &levelCut;
}


////////////////////////////////////////////////////////
/////////// Integer Rounding ///////////////////////////
////////////////////////////////////////////////////////


bool CSP1::ConstructIntSol() {
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


//double CSP1::GetXEps() { return xeps; }


void CSP1::SaveRoundedPart() {
  patBest.resize(pat.size()); patBest = pat;
  xiBest.resize(xi.size()); xiBest = xi;
}


  class J0Cmp { Vector<double> &xf; public:
    J0Cmp(Vector<double> &x_)  : xf(x_) { }
    bool operator () (int j1, int j2) const
    { return xf[j1] > xf[j2]; }
  };
bool CSP1::Round() {
  int i;
  xi.resize(lpx.size());
  xf.resize(lpx.size());    ///////// Rnd down:
	Vector<int> j0((int)xi.size()); // Indices; won't be used outs
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
  if (OUTP_LEV__ >=3) log__(" Rnd\\" << zr);
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
  if (OUTP_LEV__ >=3) log__("/" << zr<<flush);
//  } while (somefit)
  return true;
} //////////////////////////////// CSP1::Round


void CSP1::GetRHS(Vector<double> &br) {
  br.resize(m); int i;
  for (i=0;i<m;++i) br[i] = pc[i].b;
} //////////////////////////////// CSP1:: GetRHS

  CSP1::SVC::SVC(CSP1* p_, Vector<double> & b_,
      Vector<double> & d_, const double lb_,
      const double ub_)
    : Alg(p_), pr(p_), m(pr->m), b0(b_), d0(d_),
    lb(lb_), ub(ub_), kkk(0) { d.resize(m); bc.resize(m); }
  double CSP1::SVC::Run() {
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
  void CSP1::SVC::Execute() {   result = 1e+100;
    k=0;
    ConvertValues();
    do {
      try {
        zz = 0; bc = b0; pat.clear(); x.clear();
        while (GenPat()) {
          CorrectValues();
        }
        ControlSolution();
        if (pr->ExtractSolutionPart(pat, x))
          return;
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
//  d_vec dddd;
  void CSP1::SVC::ConvertValues() { // bec. of scale
    int i; // ATTENTION: when all duals of remaining
          //  pieces very small, our bb procedure doesnt
    for (i=0; i<m; ++i) {
      d[i] = d0[i] * double(pr->L);
      if (d[i] < 1) d[i] =1;
      d[i] = pow(d[i], pow_l);
    }
//    dddd = d;
  }
//  d_vec ddd;
  bool CSP1::SVC::GenPat() {
    double h = 0;  int i; double bmax=0;
    for (i=0; i<m; ++i) {
      h += double(pr->pc[i].l) * bc[i];
      if (bc[i] > bmax) bmax = bc[i];
    }
    if ( 0==h )
      return false; // this is simple but was  long unused
    if ( (!fSP2Relaxation || bmax <= 1) && h <= pr->L ) {
      FillLastPattern();
      return false;
    }
    pat.push_back(Pattern()); Column col;
    pr->GenColForHeur(&col,d,bc);
    if (col.id.empty()) {
      cout << "\nSVC could not generate patterns. Values: ";
      PrintVec(cout, d);
      cout << endl;
      cout << "rhs: ";
      PrintVec(cout, bc);
      cout << endl;
      assert(0);
    }
    pr->MakePattern(&pat.back(), &col);
    xMin = INT_MAX; // pattern intensity:
    waste = pr->L;  // original sorting ???
    Pattern::iterator iix;
    for_each_in(pat.back().ix,iix,) {
      assert(iix->x);
      assert(bc[iix->i] > 0);
      if (bc[iix->i] < double(xMin)*iix->x)
        xMin = int(double(bc[iix->i]) / iix->x); // floor
      waste -= pr->pc[iix->i].l * iix->x;
    }
    if (patternUseRatio > 0 and patternUseRatio <= 1)
      xMin = (int)ceil (double(xMin) * patternUseRatio);
    for_each_in(pat.back().ix,iix,)
       bc[iix->i] -= double(xMin)*iix->x;   // DECR. DMND
    zz += pat.back().GetObj() * xMin;
    x.push_back(xMin);
    return true;
  }
  void CSP1::SVC::CorrectValues() {
    SW = 1.1 + fabs( double( ++ kkk % 40 - 15)) / 10;
    assert ( waste < pr->L );
    double LLh = double(pr->L) / double(pr->L-waste);
    for_each_in(pat.back().ix,iix,Pattern::iterator) {
      d[iix->i] = ( LLh * pow(pr->pc[iix->i].l, pow_l) * iix->x
        + d[iix->i] * SW * (b0[iix->i] +bc[iix->i]))
        / (SW * (b0[iix->i] +bc[iix->i]) + iix->x);
//      assert((SW * (b0[iix->i] +bc[iix->i]) + iix->x) > 0);
    //  assert (d[iix->i]*0 == 0);
    }
//    ddd = d;
  }
  void CSP1::SVC::ControlSolution() { // the last one
    assert(pat.size() == x.size());
    Vector<double> s_aij(m); // fillled 0 ?
    double zzz = 0; int i;
    fill_n(&s_aij[0],m,0);
    for (i=0;i<pat.size(); ++i) {
      for_each_in(pat[i].ix,iix,Pattern::iterator) {
        s_aij[iix->i] += iix->x * x[i];
        assert(iix->x <= b0[iix->i]);
      }
      zzz += pat[i].GetObj() * x[i];
    }
    assertm(zzz == zz, "zzz="<<zzz<<", zz="<<zz);
    for (i=0;i<m;++i)
      assertm(s_aij[i] >=b0[i],"s_aij[i]="<<s_aij[i]<<", b0[i]"<<b0[i]);
}
  void CSP1::SVC::SaveSolution()
  { patBest = pat; xBest = x; }
  void CSP1::SVC::FillLastPattern() {
    pat.push_back(Pattern());
    int i;
    for (i=0;i<m;++i) if (bc[i] != 0)
      pat.back().PushIX(i,(int)bc[i]);
    pat.back().SetObj(1);
    zz += pat.back().GetObj();
    x.push_back(1);
  }


int CSP1::SVC::iterMax=20;


double CSP1::SVC::patternUseRatio=0.5;
double CSP1::SVC::pow_l=1.3;
double CSP1::SVC::outputLevel=0;


opt::OptContainer CSP1::SVC::Options() {
  opt::OptContainer oc;
  oc
    << opt::MakeOpt(&iterMax, 20,
      "iterMax", 0)
    << opt::MakeOpt(&patternUseRatio, 1,
      "patternUseRatio",
      "Usage intensity of a new pattern rel. to max")
    << opt::MakeOpt(&pow_l, 1.02,
      "pow_l", "value [i] ~ pow(l[i], pow_l);")
    << opt::MakeOpt(&outputLevel, DEF_OUTP_LEVEL,
      "outputLevel", "0-5");
  return oc;
} //____________________________________________________
opt::OptSection CSP1::SVC::opt
  ("CSP1_SVC", "SVC heur for the 1D CSP",
  SVC::Options(), opt::SolverCfg(), 5500);
//} ////////////////////// namespace { }


bool CSP1::CompleteIntSol() {
  int i;
  double LL=0;
///// Check if all fit in a single rod:
  for (i=0;i<m;++i) LL += (double(pc[i].l) * br[i]);
  if (0 == LL) {
    UpdateHeurBnd(zr);
    SaveRoundedPart();
    return true;
  }
  if (OUTP_LEV__ >=4)
    log__(" Adding to rounded value: " << zr);
  SVC svc(this, br, lpd, GetLPBnd()-zr, GetHeurBnd()-zr);
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
  }
//  if (Optimum())
  return true;
} //////////////////////////////// CSP1::CompleteIntSol


void CSP1::AddResidualSolution(CSP1::SVC & svc) {
  patBest.insert
    (patBest.end(),svc.patBest.begin(),svc.patBest.end());
  xiBest.insert
    (xiBest.end(),svc.xBest.begin(),svc.xBest.end());
}


// CSP 1D only:
void CSP1::ControlSolution() { // the best one
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
    if (fModelEquality)
      for (i=0;i<m;++i) assert(s_aij[i] == pc[i].b);
}


bool CSP1::VaryResProblem() {
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
  -- xi[i] ;  zr -= pat[i].GetObj();
  for_each_in(pat[i].ix, iix, Pattern::iterator)
    br[iix->i] += iix->x;
  return true;
} //////////////////////////////// CSP1::VaryResProblem


void CSP1::PrintProblem(ostream&os) {
  os
    //<< "CSP1 instance name: "
    << prName//<<", No. "<<inst<<
    //" from file "<<infile
    <<'\n'<<
    m<<'\n'<<L<<'\n';
  for (int i=0;i<m;++i)
    os << setw(6) << pc[i].l<<' '<<setw(6)<<pc[i].b<<'\n';
  os << endl;
}

void CSP1::PrintLog(ostream &os) {
  if (!fTestMSVC) return;
  d_vec lpd(m);
  int i;
  for (i=0;i<m;++i)
    if (double(pc[i].l) * 100 >L) // to prevent too long bb
      lpd[i] = double(pc[i].l)/double(L);
  Timer timer;
  timer.start();
  MSVC msvc(this, b, lpd, -1e100, // no bounds
    1e100);
  //double MSVCres=
    msvc.Run();
  timer.stop();
  MSVC::SolContainer::iterator its, itsOS, itsNN, itsEV;
  itsOS = itsNN = msvc.sols.end();
  itsEV = msvc.sols.begin(); // equivocal
  assert(not msvc.sols.empty());
  // SELECTING THE 3 EXTREME SOLUTIONS:
  for_each_in(msvc.sols, its,) {
    if (msvc.nOpenMaxMin == its->no) {// minOS with best NN
      if (msvc.sols.end() == itsOS)
        itsOS = its;
      else {
        if (its->nz < itsOS->nz)
          itsOS = its;
        }
    }
    if (msvc.zzMin == its->nz) { // minNN with best OS
      if (msvc.sols.end() == itsNN)
        itsNN = its;
      else {
        if (its->no < itsNN->no)
          itsNN = its;
      }
    }
    if (its->nz * msvc.nOpenMaxMin + its->no * msvc.zzMin
      < itsEV->nz * msvc.nOpenMaxMin + itsEV->no * msvc.zzMin)
      itsEV = its;
  }
  // + STATISTICS: (AVE ND...)
  sta[0] ++; // N tests
  sta[1] += timer.userTime();
  if ((itsNN->nz - zi)/zi > sta[2])
    sta[2] = (itsNN->nz - zi)/zi;
  sta[3] += (itsNN->nz - zi)/zi;
  int ii = 3;
  sta[++ii] += itsNN->nz;
  sta[++ii] += itsNN->no;
  sta[++ii] += itsNN->nd;
  sta[++ii] += itsNN->iter;
  sta[++ii] += itsOS->nz;
  sta[++ii] += itsOS->no;
  sta[++ii] += itsOS->nd;
  sta[++ii] += itsOS->iter;
  sta[++ii] += itsEV->nz;
  sta[++ii] += itsEV->no;
  sta[++ii] += itsEV->nd;
  sta[++ii] += itsEV->iter;
  os.precision(7);
  os
    // Quality of the BCP sol:
    << " Sol0: (" << msvc.sol0.no <<
    ' ' << msvc.sol0.nd << ' ' << msvc.sol0.nz << ')'
    ;
  for_each_in(msvc.sols, its,)
    os
    << " ("<< its->no <<
    ' '<< its->nd << ' ' <<its->nz << ')'
    ;
}

  /// STATISTICS:
void CSP1::InitStat() {
  sta.clear();
  sta.resize(50);
//  cout << "CSP1: Init Stat." << endl;
}
void CSP1::PrintStat(ostream& os) {
  if (!fTestMSVC or not sta[0]) return;
  int i;
  os << "MSVC stat. N="<<sta[0]<<", t,gapmax,gapave,(nz no nd iter) for minNN,OS,EV\n";
  os << sta[1]/sta[0] << ' ';
  os << sta[2] << ' ';
  for (i=3;i<16;++i)
    os << sta[i]/sta[0] << ' ';
//  os << '\n';
}

Vector<double> CSP1::sta(50);

bool CSP1::fModelEquality=true;
bool CSP1::fSortPieces=true;
bool CSP1::fMergePieces=true;
bool CSP1::fEffectiveL=true;
int CSP1::nStartBasis=1;
double CSP1::nStepsMin0=65536;
double CSP1::nStepsMinInc=1.05;
double CSP1::deps=1e-6, CSP1::bb_eps=1e-6, CSP1::xeps=1e-6, CSP1::rceps=1e-6, CSP1::veps=1e-6;
int CSP1::RPEMax = 10;
bool CSP1::fTestMSVC = 0;
bool CSP1::fSP2Relaxation = 0;
int CSP1::nLm;
double CSP1::outputLevel=DEF_OUTP_LEVEL;
bool CSP1::fCompareBySPPRC_DP=false;
bool CSP1::fInteractiveCapacity=0;
float CSP1::CVRPfrac=0.5;
int CSP1::iwstat[6];

float CSP1::cost_perturb_rndinit, CSP1::cost_perturb_max;
int CSP1::nHP;

opt::OptContainer CSP1::Options() {
  opt::OptContainer oc;
  oc
    << opt::MakeOpt(&fModelEquality, false,
      "fModelEquality", "Constraints Ax = b (take diag. start basis), otherwise Ax >= b (best)")
    << opt::MakeOpt(&fSortPieces, true,
      "fSortPieces", 0)
    << opt::MakeOpt(&fMergePieces, true,
      "fMergePieces",
      "Bool: Merge equal piece types. Set =0 if fSP2Relaxation")
    << opt::MakeOpt(&fEffectiveL, true,
      "fEffectiveL",
      "Bool: Set L to the largest integer combination of piece sizes not larger than initially. Care: for large L takes much mem.")
    << opt::MakeOpt(&nStartBasis, 4,
      "nStartBasis",
      "0:FFD, 1: Greedy (very ok but not with many small pieces), 2: IDENTITY mattrix. "
      "3: empty (best?"
      // when no artif cols: "USE IT ALSO when local upper bounds in BCP __OR__ Ax = b formulation."
      " But sometimes too long col gen when all items small); 4: SVC; 5: SVC+Greedy")
    << opt::MakeOpt(&nStepsMin0, 8192,
      "nStepsMin0",
      "Initial min N steps in B&B with cuts")
    << opt::MakeOpt(&nStepsMinInc, 1.01,
      "nStepsMinInc", "Incr ratio (each generation)")
    << opt::MakeOpt(&deps, 1e-9,
      "deps", "eps for dual multipliers (their difference from 0 and comparison)")
    << opt::MakeOpt(&rceps, 1e-8,
      "rceps", "eps for reduced costs")
    << opt::MakeOpt(&veps, 1e-5,
      "veps", "integrality eps for obj value. Needs to be even more sometimes! But then might lead to more nodes (just cautios opt test) See the README")
    << opt::MakeOpt(&xeps, 1e-8,
      "xeps", "eps for the variables")
    << opt::MakeOpt(&bb_eps, 1e-6,
      "bb_eps", "eps for b&b col gen (optimality test), setting it small (down to 1e-12) can help with small items)")
    << opt::MakeOpt(&RPEMax, 10,
      "RPEMax", "Max residual problem extensions")
    << opt::MakeOpt(&fTestMSVC, 0,
      "fTestMSVC", "Run spread/open stacks minimization SVC after main problem")
    << opt::MakeOpt(&fSP2Relaxation, 0,
      "fSP2Relaxation", "Relax. of 2D Strip Packing by setting b[i]=1 in col. gen.")
    << opt::MakeOpt(&nLm, 0,
      "nLm", "Choice order between L & m in the input. =0: choose m = min(both), "
        "=1: first m, =2: first L")
    << opt::MakeOpt(&fCompareBySPPRC_DP, false,
      "fRootPricing",
      "Pricing and SVC pattern generation in the root: 0: B&B, 1: Bounded DP Solver, if no Integer Bounding (fLocalUB) in the main algorithm "
      "-- can be faster sometimes; 2: both")
    << opt::MakeOpt(&fInteractiveCapacity, 0,
      "fInteractiveCapacity", "Interactive Capacity cuts")
    << opt::MakeOpt(&CVRPfrac, 0.5,
      "CVRPfrac", "between 0 and 1: attracting fractional part for branching choice when BrOnCVRP")
    << opt::MakeOpt(&nHP, 1,
      "nHP", "which hyperplanes: 0 - CVRP (correct for BPP only), 1 - AFF (see README/setting_params.txt)")
    << opt::MakeOpt(&cost_perturb_max, 1e-6,
      "cost_perturb_max", "product costs are -rnd(1)*value")
    << opt::MakeOpt(&cost_perturb_rndinit, 0.55,
      "cost_perturb_rndinit", "initializer for rnd")
    << opt::MakeOpt(&outputLevel, 5,
      "outputLevel", "0-5");


  return oc;
} //____________________________________________________
opt::OptSection CSP1::opt
  ("CSP1", "The 1D Cutting Stock Problem",
  CSP1::Options(), opt::SolverCfg(), 500);

  void CSP1::MSVC::Execute() {   result = 1e+100;
    k=0;
    //ConvertValues();
    int i;
    d.resize(m);
    for (i=0;i<m;++i)
      d[i] = pow(pr->pc[i].l,pow_l);
    bb.m = m;
    bb.Reallocate(m);
    bb.fOpen.resize(m);
    dd = d;
    mosR = mosR__;
    msR = msR__;
/*    nOpenMaxMin = 1e100;
    nDiffPatMin = 1e100;
    spreadMaxMin = 1e100;
    zzMin = 1e100;*/ // INIT only in constructor
    nOpenMaxTarget = m;
    spreadMaxTarget = pr->GetHeurBnd();
    nOpenMaxNext = (nOpenMaxInitial) ? nOpenMaxInitial: INT_MAX;
    do {
      try {
//        if (not k%5)
  //        pow_l = (1.3 - 1.01) * rnd + 1.01; // BAD
        StartNewPlan();
        zz = 0; bc = b0; pat.clear(); x.clear();
        while (GenPat()) {
          CheckPattern();
          CorrectValues();
        }
        ControlSolution();
//        if (pr->ExtractSolutionPart(pat, x))
    //      return;
        CheckSolution(); // pareto-opt
      } catch (...) {
        if (OUTP_LEV__ >=3) PRINT_ERROR("MSVC: Some error.");
      } // Let it be
    } while (++k < iterMax );
  }
  int CSP1::MSVC::GenPat() {
    h = 0;  int i;
    for (i=0; i<m; ++i) h += double(pr->pc[i].l) * bc[i];
    if ( 0==h )
      return false; // this is simple but was  long unused
    if (h <= pr->L) {
      FillLastPattern();
      return false;
    }
    pat.push_back(Pattern()); Column col;
/*    if (!fRestrictNOpen)
      pr->GenColForHeur(&col,d,bc);
    else */ // Always
      GenRestrCol(&col,d,bc);
    pr->MakePattern(&pat.back(), &col);
    xMin = INT_MAX; // pattern intensity:
    waste = pr->L;  // original sorting ???
    Pattern::iterator iix;
    for_each_in(pat.back().ix,iix,) {
      assert(iix->x);
      assert(bc[iix->i] > 0);
      if (bc[iix->i] < double(xMin)*iix->x)
        xMin = int(double(bc[iix->i]) / iix->x); // floor
      waste -= pr->pc[iix->i].l * iix->x;
    }
    if (patternUseRatio > 0 and patternUseRatio <= 1)
      xMin = (int)ceil (double(xMin) * patternUseRatio);
    for_each_in(pat.back().ix,iix,)
       bc[iix->i] -= double(xMin)*iix->x;   // DECR. DMND
    zz += pat.back().GetObj() * xMin;
    x.push_back(xMin);
    return true;
  }
  void CSP1::MSVC::GenRestrCol
  (Column *col,const d_vec &d,Vector<double> &b) {
  int i;

  for (i=0;i<m;++i)
    bb.fOpen[i] = (bc[i] > 0 && bc[i] < b0[i]);
  if (fRestrictNOpen) {
    bb.nOpenMax = nOpenMaxInitial;
    bb.fDominance = false; // ???
  } else {
    bb.nOpenMax = INT_MAX;
    bb.fDominance = 1; // ???
  }

  //bb.nOpenMax = nOpenMaxNext; // when the iterative decr.

  bb.fHeur = true;
  bb.L = pr->L;
  for (i=0;i<m;++i) {
    bb.pieces__[i].d = d[i];
    bb.pieces__[i].b = int(b[i]);
    bb.pieces__[i].l = int(pr->pc[i].l);
  }
  bb.zMin = -INFINITY__;
  bb.zLowerInitial = -INFINITY__;
  bb.eps = pr->GetBBEps()*double(pr->L); // HEREEEEEEEEEEEEEEEEEEEEE !!!!!!! ???????????


  ColSet cs; // dummy
  bb.colsetRes = &cs;
  bb.colBest = col;

  bb.Run();
  assert(not col->id.empty());
// Addi info to all cols:
  col->SetObj(1);
}
  void CSP1::MSVC::StartNewPlan() {
        zz = 0; bc = b0; pat.clear(); x.clear();
        nPatLast = 0;
        nOpenMax = 0;
        spreadMax = 0;
        spread.clear(); spread.resize(m,0);
        d = dd; // take original values for begin
        mosR *= mosRR;
        msR *= msRR;
        if (rndCMCM < 1) rndCMCM = 1;
        rndCM = (rndCMCM-1) * rnd +1;

 //       if (k >= kNOpenMaxNext)
     //     nOpenMaxNext = nOpenMaxMin - dNOpen;
        log_n_(3,"MSVC: iter="<<k);
  }
  void CSP1::MSVC::CheckPattern() {
    int i;
    Pattern::iterator iix;
    nStillOpen = 0;
    nOpen = 0;
		nPatLast = (int)pat.size();
    // Counting what is related to spread/N open stacks
    // and fixing maximums:
    for_each_in(pat.back().ix,iix,) {
      if (0==bc[iix->i]) { // just finished
        nOpen ++;
        spread[iix->i] += x.back()/double(pr->pc[iix->i].l)/b0[iix->i];
      }
    }
    for (i=0;i<m;++i) {
      if (b0[i] > bc[i] && 0 < bc[i]) // started but not finished
        spread[i] += x.back()/double(pr->pc[i].l)/b0[i];
      if (spread[i] > spreadMax)
        spreadMax = spread[i];
      nStillOpen += ((0 < bc[i]) && (b0[i]>bc[i]));
    }
    nOpen += nStillOpen;
    if (nOpen > nOpenMax)
      nOpenMax = nOpen;
    if (fRestrictNOpen)
      assertm(nOpen <= nOpenMaxInitial || h <= pr->L,
      nOpen<<"(nOpen) > (NOpenMax)"<<nOpenMaxInitial);
    // oder ..Next
  }
  void CSP1::MSVC::CorrectValues() {
    int i;
    Pattern::iterator iix;
    // correcting 'old' values in dd:
    SW = ((2-2/rndCM)*rnd + 1/rndCM);
      //1;//.1 + fabs( double( ++ kkk % 40 - 15)) / 10;
    assert ( waste < pr->L );
    double LLh = double(pr->L) / double(pr->L-waste);
    double weight1, weight2;
    for_each_in(pat.back().ix,iix,) {
      switch (weighScheme) {
      case 1:
        weight1 = iix->x;
        weight2 = SW * bc[iix->i];
        break;
      case 2:
        weight1 = iix->x;
        weight2 = SW * b0[iix->i];
        break;
      case 3:
        weight1 = iix->x;
        weight2 = SW * (b0[iix->i] - iix->x);
        break;
      case 4:
        weight1 = rnd;
        weight2 = rnd;
        break;
      default:
        weight1 = iix->x;
        weight2 = SW * (b0[iix->i] +bc[iix->i]);
      }
      dd[iix->i] = ( LLh * pow(pr->pc[iix->i].l, pow_l) * weight1
        + dd[iix->i] * weight2)
        / (weight1 + weight2);
    }
    // assigning actual values in d:
    double mOS = mosR * pow(FMax(nStillOpen/nOpenMaxTarget,1),mosP);
    assert(dd.size() >= m);
    for (i=0;i<m;++i) {
      if (bc[i] > 0 && bc[i] < b0[i]) {
        double mSi = msR * pow(FMax(spread[i]/spreadMaxTarget,1),msP);
        d[i] = dd[i] * (nObjective * mSi + (1-nObjective) * mOS);
      }
      else d[i] = dd[i];
    // + RANDOM for SW, see above!
      //d[i] = d[i] * ((rndCM-1/rndCM)*rnd + 1/rndCM);
    }
  }
  void CSP1::MSVC::CheckSolution() {
    bool fPrintSol=0;
    bool fBetterSol = sols.empty();
    SolContainer::iterator its, its1;
    if (pat.size() != nPatLast) {
      fill(bc.begin(),bc.end(),0); // all finished
      CheckPattern();
    }
    int nDiffPat = (int)set<Pattern>(pat.begin(),pat.end()).size();
    // whether it is pareto-better:
    fBetterSol = 1;
    for_each_in(sols, its, ) {
      if (its->no <= nOpenMax
        and its->nz <= zz and its->nd <= nDiffPat)
        fBetterSol = 0;
    }
    if (fBetterSol) { // deleting worse:
    its = sols.begin();
    while (its != sols.end()) {
      its1= its;
      ++ its1;
      if (its->no >= nOpenMax
        and its->nz >= zz and its->nd >= nDiffPat) {
        sols.erase(its);
      }
      its = its1;
    }
    }
    if (nOpenMax < nOpenMaxMin)
      nOpenMaxMin = nOpenMax;
    if (zz < zzMin)
      zzMin = zz;
    if (nDiffPat < nDiffPatMin)
      nDiffPatMin = nDiffPat;
    if (fBetterSol) {
      fPrintSol = 1;
      Solution s;
      s.nz = zz; s.no = nOpenMax; s.nd = nDiffPat;
      s.iter = k;
      sols.insert(s);
      if (!k) sol0 = s;
      log_n(2,"Iter "<<k<<". Sol. found: ("<<nOpenMax<<' '
        <<nDiffPat<<'['<<pat.size()<<"] "<<zz<<')');
    }
    if (spreadMax < spreadMaxMin)
      spreadMaxMin=spreadMax;
    nOpenMaxTarget = nOpenMaxMin * mosReduc;
    spreadMaxTarget = spreadMaxMin * msReduc;
      if (fPrintSol && OUTP_LEV__ >= 4) {
        int i;
        for (i=0;i<pat.size();++i) {
          log__(x[i]<<':');
          for_each_in(pat[i].ix,iix,Pattern::iterator)
            log__(' '<<iix->i<<'/'<<iix->x);
          log__('\n');
        }
      }
  }

int CSP1::MSVC::iterMax=20;

bool CSP1::MSVC::fRestrictNOpen=0;
int CSP1::MSVC::restrictDelay=5;
int CSP1::MSVC::dNOpen=1;
int CSP1::MSVC::nOpenMaxInitial=10;
double CSP1::MSVC::patternUseRatio=0.5;
double CSP1::MSVC::outputLevel=4;

double CSP1::MSVC::pow_l;
int CSP1::MSVC::weighScheme=0;

   double CSP1::MSVC::nObjective; // 0=open stacks min,
     // 1=spread min, between?
   double CSP1::MSVC::mosReduc;
   double CSP1::MSVC::mosP;
   double CSP1::MSVC::mosR__;
   double CSP1::MSVC::mosRR;
   double CSP1::MSVC::msReduc;
   double CSP1::MSVC::msP;
   double CSP1::MSVC::msR__;
   double CSP1::MSVC::msRR;

// randomization for each plan?
   double CSP1::MSVC::rndCMCM;
   // rndCM in [1,1+rnd*(rndCMCM-1)]



opt::OptContainer CSP1::MSVC::Options() {
  opt::OptContainer oc;
  oc
    << opt::MakeOpt(&iterMax, 200,
      "iterMax", 0)
    << opt::MakeOpt(&fRestrictNOpen, 0,
      "fRestrictNOpen",
      "whether to restrict N open stacks at each new pattern generation."
      " If yes, then remove value modifiers: set m(o)sR(R)=1 and m(o)sP=0")
/*    << opt::MakeOpt(&restrictDelay, 5,
      "restrictDelay",
      "NOT DONE YET. after how many sols. to set NOpenMax <= NOpenBest - dNOpen"
      " to let the method find better material consumptions")
    << opt::MakeOpt(&dNOpen, 1,
      "dNOpen", "s. above")*/
    << opt::MakeOpt(&nOpenMaxInitial, 10,
      "NOpenMaxInitial", "Max. NOpen. NOTE: if a solution appears with more,"
      " then this is in the last pattern with intensity =1")
    << opt::MakeOpt(&patternUseRatio, 1,
      "patternUseRatio",
      "Usage intensity of a new pattern rel. to max."
      " Setting=1 ensures that at least 1 stack closes after a pattern")
    << opt::MakeOpt(&pow_l, 1.02,
      "pow_l", "value [i] ~ pow(l[i], pow_l);")
    << opt::MakeOpt(&weighScheme, 0,
      "weighScheme", "weigh new values: 0=classic, ..., 4= random")
    << opt::MakeOpt(&nObjective, 0,
      "nObjective", "0=open stacks min, 1=spread min, between: "
      " d[i] = dd[i] * (nObjective * mSi + (1-nObjective) * mOS); (for all open types)")
    << opt::MakeOpt(&mosReduc, 0.7,
      "mosReduc", "nOpenMaxTarget = nOpenMaxMin * mosReduc")
    << opt::MakeOpt(&mosP, 0,
      "mosP", "mOS = mosR * pow(FMax(nStillOpen/nOpenMaxTarget,1),mosP);")
    << opt::MakeOpt(&mosR__, 1.02, "mosR__", "Initial mosR")
    << opt::MakeOpt(&mosRR, 1.001, "mosRR", "increase fo mosR in each new sol.")

    << opt::MakeOpt(&msReduc, 0.7,
      "msReduc", "spreadMaxTarget = spreadMaxMin * msReduc")
    << opt::MakeOpt(&msP, 0,
      "msP", "mSi = msR * pow(FMax(spread[i]/spreadMaxTarget,1),msP);")
    << opt::MakeOpt(&msR__, 1.02, "msR__", "Initial msR")
    << opt::MakeOpt(&msRR, 1.001, "msRR", "increase fo msR in each new sol.")
    << opt::MakeOpt(&rndCMCM, 1.5,
      "rndCMCM", " rndCM = max(rndCMCM-1,0) * rnd +1 for each solution; "
      " SW = ((2-2/rndCM)*rnd + 1/rndCM);")
    << opt::MakeOpt(&outputLevel, DEF_OUTP_LEVEL,
      "outputLevel", "0-5");
  return oc;
} //____________________________________________________
opt::OptSection CSP1::MSVC::opt
  ("CSP1_MSVC", "SVC heur for 1D MOSP",
  MSVC::Options(), opt::SolverCfg(), 5500);
//} ////////////////////// namespace { }


void CSP1::SSVC::Execute() { 
//  nStartUsed = 0;
  MSVC::Execute();
}
void CSP1::SSVC::StartNewPlan() {
  MSVC::StartNewPlan();
  fUsed.clear(); fUsed.resize(x0.size());
//  if (++ nStartUsed >= pat0.size())
  //  nStartUsed = 0; // each time a new
}

  int CSP1::SSVC::GenPat() {
    int i;
    double zkp;
    Pattern::iterator iix;
    iNext = -1; zkpBest = -1e100;
    if (pat.empty())
      iNext = int(rndStart*(double)pat0.size());
    else
    for (i=0;i<pat0.size();++i)  // "Generating": choose
      if (not fUsed[i]) {
        zkp = 0;
        for_each_in(pat0[i].ix, iix, )
          zkp += d[iix->i] * iix->x;
        if (zkp > zkpBest) {
          iNext = i; zkpBest = zkp;
        }
      }
    if (-1 == iNext) return 0;
    fUsed[iNext] = 1;
    pat.push_back(pat0[iNext]);
    int xMin = (int)x0[iNext]; // pattern intensity:
    waste = pr->L;  // original sorting ???
    for_each_in(pat.back().ix,iix,) {
      assert(iix->x);
//      assert(bc[iix->i] > 0);
      waste -= pr->pc[iix->i].l * iix->x;
    }
    for_each_in(pat.back().ix,iix,)
       bc[iix->i] -= double(xMin)*iix->x;   // DECR. DMND
    zz += pat.back().GetObj() * xMin;
    x.push_back(xMin);
    return true;
  }

  void CSP1::SSVC::CheckSolution() {
    bool fPrintSol=0;
    bool fBetterSol = sols.empty();
    SolContainer::iterator its, its1;
    if (pat.size() != nPatLast) {
      fill(bc.begin(),bc.end(),0); // all finished
      CheckPattern();
    }
    int nDiffPat = (int)pat0.size();
      //set<Pattern>(pat.begin(),pat.end()).size();
    // whether it is pareto-better:
    fBetterSol = 1;
    for_each_in(sols, its, ) {
      if (its->no <= nOpenMax
        and its->nz <= zz and its->nd <= nDiffPat)
        fBetterSol = 0;
    }
    if (fBetterSol) { // deleting worse:
    its = sols.begin();
    while (its != sols.end()) {
      its1= its;
      ++ its1;
      if (its->no >= nOpenMax
        and its->nz >= zz and its->nd >= nDiffPat) {
        sols.erase(its);
      }
      its = its1;
    }
    }
    if (nOpenMax < nOpenMaxMin)
      nOpenMaxMin = nOpenMax;
    if (zz < zzMin)
      zzMin = zz;
    if (nDiffPat < nDiffPatMin)
      nDiffPatMin = nDiffPat;
    if (fBetterSol) {
      fPrintSol = 1;
      Solution s;
      s.nz = zz; s.no = nOpenMax; s.nd = nDiffPat;
      s.iter = k;
      sols.insert(s);
      if (!k) sol0 = s;
      log_n(2,"Iter "<<k<<". Sol. found: ("<<nOpenMax<<' '
        <<nDiffPat<<'['<<pat.size()<<"] "<<zz<<')');
    }
    if (spreadMax < spreadMaxMin)
      spreadMaxMin=spreadMax;
    nOpenMaxTarget = nOpenMaxMin * mosReduc;
    spreadMaxTarget = spreadMaxMin * msReduc;
      if (fPrintSol && OUTP_LEV__ >= 4) {
        int i;
        for (i=0;i<pat.size();++i) {
          log__(x[i]<<':');
          for_each_in(pat[i].ix,iix,Pattern::iterator)
            log__(' '<<iix->i<<'/'<<iix->x);
          log__('\n');
        }
      }
  }


}                                                // SS_END_NAMESPACE__
