// wrpcpx.cpp: CPLEX wrapper using Callable Library

#include <unistd.h>

#include "stdafx.h"
#include "wrpcpx.h"
#include "lasthdr.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

// only CPLEX implemented
LP * LP::CreateLP(int )
{ return new WrpCpx(); }

// CPXENVptr WrpCpx::env = NULL;
// int nEnv=0;

WrpCpx::WrpCpx()
: //env(NULL),
 lp(NULL) { 
 	 env = NULL; 
}

void WrpCpx::Open()
{
Try:
   if (NULL == env)
		env = CPXopenCPLEX (&status);

// 	 ++ nEnv;
   /* If an error occurs, the status value indicates the reason for
      failure.  A call to CPXgeterrorstring will produce the text of
      the error message.  Note that CPXopenCPLEX produces no output,
      so the only way to see the cause of the error is to use
      CPXgeterrorstring.  For other CPLEX routines, the errors will
      be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON.  */

   if ( env == NULL ) {
      char  errmsg[1024];
      CPXgeterrorstring (env, status, errmsg);
      fprintf (stderr, "Could not open CPLEX environment.\n");
      fprintf (stderr, "%s", errmsg);
      fprintf (stderr, "  Waiting 50 seconds...\n");
      sleep(50);
      goto Try;
      //fprintf (stderr, "%s", errmsg);
      assertm(false,errmsg);
   }

   /* Turn on output to the screen */
/*
   status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
   if ( status ) {
      fprintf (stderr, 
        "Failure to turn on screen indicator, error %d.\n", status);
      assert(false);
   }
   // OR CPXsetlogfile()
*/
   /* Turn on data checking */
	 
	 status
	 = CPXsetintparam (env, CPX_PARAM_DATACHECK, CPX_ON);
	 if ( status ) {
		 fprintf (stderr, 
							"Failure to turn on data checking, error %d.\n", status);
		 assert(false);
	 }
	 
	 
   /* Create the problem. */

   lp = CPXcreateprob (env, &status, "__1D-CSP");

   /* A returned pointer of NULL may mean that not enough memory
      was available or there was some other problem.  In the case of 
      failure, an error message will have been written to the error 
      channel from inside CPLEX.  In this example, the setting of
      the parameter CPX_PARAM_SCRIND causes the error message to
      appear on stdout.  */

   if ( lp == NULL ) {
      fprintf (stderr, "Failed to create LP.\n");
      assertm(false,"Failed to create LP.\n");
   }

// OBJECTIVE::
   CPXchgobjsen (env, lp, CPX_MIN);
/*   status = CPXsetintparam (env, CPX_PARAM_SCAIND, -1);
   if (!status)
     fprintf(stderr, "CPX: Failed to cancel scaling.\n");*/
}

WrpCpx::~WrpCpx() { Close(); }

void WrpCpx::Close()
{
   /* Free up the problem as allocated by CPXcreateprob, if necessary */

   if ( lp != NULL ) {
      status = CPXfreeprob (env, &lp);
      if ( status ) {
         fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
      }
   }

   /* Free up the CPLEX environment, if necessary */

// 	 -- nEnv;
   if ( env != NULL
// 		 && !nEnv 
	) {
      status = CPXcloseCPLEX (&env);

      /* Note that CPXcloseCPLEX produces no output,
         so the only way to see the cause of the error is to use
         CPXgeterrorstring.
         For other CPLEX routines, the errors will
         be seen if the
         CPX_PARAM_SCRIND indicator is set to CPX_ON. */

      if ( status > 0 ) {
         char  errmsg[1024];
         fprintf (stderr, "Could not close CPLEX environment.\n");
         CPXgeterrorstring (env, status, errmsg);
         fprintf (stderr, "%s", errmsg);
      }
   }


}

int WrpCpx::Dim()
// { return m0 + cuts.size(); }
 { return CPXgetnumrows(env,lp); }
int WrpCpx::NCols()
{ return CPXgetnumcols(env,lp); }

void WrpCpx::SolvePrimal() {
// Again:
	/// These tolerances helped to obtain the exact LP value for instances 402_*:
	/// Having them here because for MIP it's better to have defaults:
	status
	= CPXsetdblparam (env, CPX_PARAM_EPOPT, epsOpt);
	if ( status ) {
		fprintf (stderr, 
						 "Failure to change CPLEX LP optimality tolerance to %f, error %d.\n", epsOpt, status);
		assert(false);
	}
	
	status
	= CPXsetdblparam (env, CPX_PARAM_EPRHS, epsFeas);
	if ( status ) {
		fprintf (stderr, 
						 "Failure to change CPLEX LP feas-ty tolerance to %f, error %d.\n", epsFeas, status);
		assert(false);
	}
	status = CPXprimopt(env,lp);
  if (status) {
       CPXwriteprob (env, lp, "PrimFailed.lp", NULL);
       char  errmsg[1024];
       fprintf (stderr, "Could not prim. optimize.\n");
       CPXgeterrorstring (env, status, errmsg);
       fprintf (stderr, "%s", errmsg);
  }
/*  assertm(!status,
    "Failed to prim. optimize. See PrimFailed.lp");*/
  GetStatus(LP::status);
// The following not good for a BCP where a node can inf
/*  int lpstat = CPXgetstat (env, lp);
  if (CPX_INFEASIBLE == lpstat
    || CPX_OPTIMAL_INFEAS == lpstat) {*/
  /* But unscaled infeasibilities: */
  int lpstat = CPXgetstat (env, lp);
  if (CPX_STAT_OPTIMAL_INFEAS == lpstat
	  or status) {
    log_ln("Solution infeasible...");
    // Trying to change scaling:
/*    if (0==iScale) iScale = 1; // more aggressive
    else if (1==iScale) iScale = -1;
    else goto ReallyInfeas;*/
    log_ln("Re-optimizing LP from scratch...");
    CPXsetintparam (env, CPXPARAM_Advance, 0);
//     log_ln("Turning off preprocessing...");
//     status = CPXsetintparam (env, CPX_PARAM_PREIND, 0);
//     assertm(!status,"Could not cancel preproc");
    status = CPXlpopt(env, lp);              // dual?
    if (status) CPXwriteprob (env, lp, "LP_resolve_failed.lp", NULL);
    assertm(!status,
      "Failed to resolve LP. See LP_resolve_failed.lp");
//     status = CPXsetintparam (env, CPX_PARAM_PREIND, 1);
//     assertm(!status,"Could not set preproc");
    GetStatus(LP::status);
	int lpstat = CPXgetstat (env, lp);
    CPXsetintparam (env, CPXPARAM_Advance, 1);       // returning to use advanced basis info
	if (CPX_STAT_OPTIMAL_INFEAS == lpstat
		or status) {
	log_ln("Solution infeasible...");
	// Trying to change scaling:
	/*    if (0==iScale) iScale = 1; // more aggressive
	else if (1==iScale) iScale = -1;
	else goto ReallyInfeas;*/
	log_ln("Turning on aggressive scaling and calling dual simplex method...");
	CPXsetintparam (env, CPX_PARAM_SCAIND, 1);
	//     log_ln("Turning off preprocessing...");
	//     status = CPXsetintparam (env, CPX_PARAM_PREIND, 0);
	//     assertm(!status,"Could not cancel preproc");
	status = CPXdualopt(env, lp);              // dual?
	if (status) CPXwriteprob (env, lp, "DualFailed_AggrScaling.lp", NULL);
	assertm(!status,
		"Failed to re-optimize. See DualFailed_AggrScaling.lp");
	//     status = CPXsetintparam (env, CPX_PARAM_PREIND, 1);
	//     assertm(!status,"Could not set preproc");
	GetStatus(LP::status);
	int lpstat = CPXgetstat (env, lp);
	if (CPX_STAT_OPTIMAL_INFEAS == lpstat
			or status) {
		log_ln("Solution infeasible again...");
		// Trying to change scaling:
		/*    if (0==iScale) iScale = 1; // more aggressive
		else if (1==iScale) iScale = -1;
		else goto ReallyInfeas;*/
		log_ln("Turning off scaling and calling dual simplex method...");
		CPXsetintparam (env, CPX_PARAM_SCAIND, -1);
		//     log_ln("Turning off preprocessing...");
		//     status = CPXsetintparam (env, CPX_PARAM_PREIND, 0);
		//     assertm(!status,"Could not cancel preproc");
		status = CPXdualopt(env, lp);              // dual?
		if (status) CPXwriteprob (env, lp, "DualFailed_NoScaling.lp", NULL);
		assertm(!status,
			"Failed to re-optimize. See DualFailed_NoScaling.lp");
		//     status = CPXsetintparam (env, CPX_PARAM_PREIND, 1);
		//     assertm(!status,"Could not set preproc");
		GetStatus(LP::status);
	}
    }
  }

// ReallyInfeas:;
// NOW SAVING VALUES TO ENABLE CONVERSION TO MIP?
    // Not, solve() after MIP.
}

void WrpCpx::SolveDual() {
	/// These tolerances helped to obtain the exact LP value for instances 402_*:
	status
	= CPXsetdblparam (env, CPX_PARAM_EPOPT, epsOpt);
	if ( status ) {
		fprintf (stderr, 
						 "Failure to change CPLEX LP optimality tolerance to %f, error %d.\n", epsOpt, status);
		assert(false);
	}
	
	status
	= CPXsetdblparam (env, CPX_PARAM_EPRHS, epsFeas);
	if ( status ) {
		fprintf (stderr, 
						 "Failure to change CPLEX LP feas-ty tolerance to %f, error %d.\n", epsFeas, status);
		assert(false);
	}
	//  try {
//    cplex.exportModel("rmp.lp");
  status = CPXdualopt(env,lp);
/*  if (status) {
    CPXwriteprob (env, lp, "DualFailed.lp", NULL);
         char  errmsg[1024];
         fprintf (stderr, "Could not prim. optimize.\n");
         CPXgeterrorstring (env, status, errmsg);
         fprintf (stderr, "%s", errmsg);
  }
  assertm(!status,
    "Failed to dually optimize. See DualFailed.lp");*/
  GetStatus(LP::status);
}

void WrpCpx::AddCol(double obj,int nnz,
  int *cmatind,double *cmatval,double lb,double ub) {

  int cmatbeg = 0;
/*  status = CPXcheckaddcols(env,lp,1,nnz,&obj,&cmatbeg,
      cmatind,cmatval,&lb,&ub,NULL);
  assertm(!status,"Could not check addcols"); */
  status = CPXaddcols(env,lp,1,nnz,&obj,&cmatbeg,
      cmatind,cmatval,&lb,&ub,NULL);
  if (status) {
    cout << " nnz="<<nnz<<" obj="<<obj<<" lb="<<lb<<" ub="<<ub<<" NCols=" << NCols() << " Dim="<<Dim()<<endl;
    for (int i=0;i<nnz;++i)
      cout << '(' << cmatind[i] << ':' << cmatval[i] << ')';
  }
  assertm(!status,"Could not addcols: status="<<status);
}

// Call this only after primal
void WrpCpx::GetStatus(StatusValues &s) {
   int lpstat = CPXgetstat (env, lp);
   switch (lpstat)
   {
//   case CPX_OPTIMAL_INFEAS :
   case CPX_STAT_OPTIMAL : s=LP::opt; break;
   //case IloAlgorithm::Unknown: s=LP::none; break;
   case CPX_STAT_INFEASIBLE : s=LP::infeas; break;
   case CPX_STAT_UNBOUNDED : s=LP::unbounded; break;
   default:
//      mylog << "An error occurred during the solution process";
     char buf[1024];
     CPXgetstatstring(env,lpstat,buf);
     log_ln("\n Unexpected CPX ERROR: "<<buf);
     s=LP::error; break;
   }
/*   if (lpstat != CPX_OPTIMAL) {
     char buf[1024];
     CPXgetstatstring(env,lpstat,buf);
     log_ln("\n See LPNotOpt.lp. CPX ERROR: "<<buf);
     status = CPXwriteprob (env, lp, "LPNotOpt.lp", NULL);
   }*/
}

double WrpCpx::GetValue()
{
  double v;
  status = CPXgetobjval (env, lp, &v);
//  assert2m(!status,"CPX: No solution exists.");
      if ( status > 0 ) {
         char  errmsg[1024];
         fprintf (stderr, "Could not get LP value.\n");
         CPXgeterrorstring (env, status, errmsg);
         fprintf (stderr, "%s", errmsg);
      }
  if (status)
    v = 1e100; // for min. pr.
  return v;
}

void WrpCpx::GetMultipliers(Vector<double> &d)
{
//  Solve();
  if (d.size() < Dim()) d.resize(Dim());
  status = CPXgetpi (env, lp, &d[0], 0, Dim()-1); 
  if (status) {
       CPXwriteprob (env, lp, "GetMultsFailed.lp", NULL);
       char  errmsg[1024];
       fprintf (stderr, "Could not get simplex mults.\n");
       CPXgeterrorstring (env, status, errmsg);
       fprintf (stderr, "%s", errmsg);
  }
  assertm(!status,"CPX: getting multipliers.");
}

void WrpCpx::GetRedCosts(d_vec &d)
{
//  Solve();
  d.resize(NCols());
  status = CPXgetdj (env, lp, &d[0], 0, NCols()-1);
  if (status) {
       CPXwriteprob (env, lp, "GetRedCostsFailed.lp", NULL);
       char  errmsg[1024];
       fprintf (stderr, "Could not get red costs.\n");
       CPXgeterrorstring (env, status, errmsg);
       fprintf (stderr, "%s", errmsg);
  }
  assertm(!status,"CPX: getting red costs.");
}

void WrpCpx::GetSolution
  (Vector<LP::ColStatus> *cs,Vector<double> *x)
{
  if (cs) {
    cs->resize(NCols());
    assert(sizeof(int) == sizeof(LP::ColStatus));
    status = CPXgetbase (env, lp, (int*)(&(*cs)[0]), NULL);
    assertm(!status,"CPX: could not get basis");
    int i;
    for (i=NCols()-1;i>=0;--i)
    switch(*(int*)&((*cs)[i])) {
    case CPX_AT_LOWER: (*cs)[i] = LP::atLower; break;
    case CPX_BASIC: (*cs)[i] = LP::basic; break;
    case CPX_AT_UPPER: (*cs)[i] = LP::atUpper; break;
    case CPX_FREE_SUPER: (*cs)[i] = LP::free_super; break;
    default: assert(false);
    }
  }
  if (x) {
    x->resize(NCols());
    status = CPXgetx (env, lp, &(*x)[0], 0, NCols()-1);
    assertm(!status,"CPX: could not get solution");
  }
}

void WrpCpx::InitConstr
  (Vector<double> &b, const char* pSense)
{
	status = CPXnewrows (env, lp, (int)b.size(),
           &b[0], pSense, NULL, NULL);
  assertm(!status,"CPX: Init constraints: status="<<status);
}

// x.size() must be NCols() 'cause x[i] => var non-int
double WrpCpx::SolveMIP
  (Vector<double> &x, double tm, double UB) {
  int i;
  double res;
  char t = CPX_INTEGER;

  try {

  status = CPXchgprobtype (env, lp, CPXPROB_MILP);
  assertm(!status,"Could not convert problem to MILP");
  for (i=0;i<NCols();++i) {
    if (x[i])
      status = CPXchgctype (env, lp, 1, &i, &t);
    // save var bounds
/*  char bndC='L';
  status = CPXchgbds(env,lp,1,&ic,&bndC,&lb);
  assert(!status);
  bndC = 'U';
  status = CPXchgbds(env,lp,1,&ic,&bndC,&ub);
  assert(!status);*/
  }
//    status = CPXsetintparam
//      (env, CPX_PARAM_STARTALG, CPX_ALG_PRIMAL);
//    assertm(!status,"MIP: Setting start algorithm");

	 status = CPXsetintparam
	 (env, CPX_PARAM_THREADS, nMIPThreads);
	 assertm(!status,"MIP: Setting nThreads: status="<<status);
	 
   status = CPXsetdblparam(env, CPX_PARAM_TILIM, tm);
//    status = CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 2);
	 
	 if (UB<1e75) {
		 status = CPXsetdblparam(env, CPX_PARAM_CUTUP, UB);    // +1e-5: CPLEX has own tolerances
		 assertm(!status,"MIP: Setting CUTUP to " << UB);
	 }

	 /* Turn on output to the screen */
		   status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
		   if ( status ) {
		      fprintf (stderr, 
		        "Failure to turn on screen indicator, error %d.\n", status);
		      assert(false);
		   }
	// OR CPXsetlogfile()
	
	 
	 
	 status = CPXwriteprob (env, lp, "rmpMIP.lp", NULL);
	 if (status) {
		 cerr << "\nCould not export rmpMIP.lp, status=" << status << endl;
	 }
	 
   status = CPXmipopt (env, lp);
   assertm(!status,"MIP failed, see rmpMIP.lp");

   status = CPXsetdblparam(env, CPX_PARAM_TILIM, 1e+75);
   int   lpstat = CPXgetstat (env, lp);

   if ((lpstat!=CPXMIP_OPTIMAL)
     && (lpstat!=CPXMIP_TIME_LIM_FEAS))
   { res = 1e100; goto Clean; }

  status = CPXgetmipx (env, lp, &x[0], 0, NCols()-1);

  status = CPXgetmipobjval (env, lp, &res);
  if (status) res = 1e100;
  } catch (const exception& e) {
    log_ln(e.what());
  }
Clean:
  status = CPXchgprobtype (env, lp, CPXPROB_LP);                    // NEED THIS?
  assertm(!status,"Could not convert problem to LP");
	status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
// 	status
//    = CPXprimopt (env, lp); // to reinstall lp solution
//   assertm(!status,"Could not resolve reinstalled LP");
	status = CPXsetdblparam(env, CPX_PARAM_CUTUP, 1e75);
	assertm(!status,"MIP: Resetting CUTUP to " << 1e75);
	status = CPXsetintparam (env, CPX_PARAM_THREADS, 0);
	return res;
}

void WrpCpx::GetBasisInverse
  (Vector<Vector<double> > &bi,Vector<double> &xx) {

  int i; //,j;

  cstat.resize(NCols());
  status = CPXgetbase (env, lp, &cstat[0], NULL);
  assertm(!status,"CPX: could not get basis");

  Vector<int> ibas;
  ibas.reserve(Dim());
  xx.clear(); xx.reserve(Dim());
  theX.resize(NCols());
  lb.resize(NCols());
  ub.resize(NCols());
  status = CPXgetx (env, lp, &theX[0], 0, NCols()-1);

  d_vec zeros(Dim()), oldobj(Dim());
  for (i=0;i<NCols();++i)
    if (CPX_BASIC == cstat[i]) {
//      bascol.push_back(GetCol(i));
      CPXgetobj(env, lp, &oldobj[ibas.size()], i, i);
      xx.push_back(theX[i]);
      ibas.push_back(i);
    }
    else {
      GetVarBnds(i, lb[i], ub[i]);
      if (CPX_AT_LOWER == cstat[i])
        ChangeVarBnds(i, lb[i], lb[i]);
      else
        ChangeVarBnds(i, ub[i], ub[i]);
    }

  status = CPXchgobj (env, lp,
    (int)ibas.size(), &ibas[0], &zeros[0]);
  assertm(!status,"Change bas obj coefs to 0");

  /// Extracting lines of the inverse:
  bi.resize(ibas.size());
  for (i=0;i<ibas.size();++i) {
    status = CPXchgcoef (env, lp, -1, ibas[i], 1);
    if (i>0)
      status = CPXchgcoef (env, lp, -1, ibas[i-1], 0);
    assert(!status);
    status = CPXlpopt(env, lp);
    if (status) {
      char  errmsg[1024];
      CPXgeterrorstring (env, status, errmsg);
      fprintf (stderr, "Could not solve the help problem"
      " to obtain basis inverse.\n"
//        "Writing the basis model to basinv.lp"
        " and the main model to basinvA.lp.\n");
//      CPXwriteprob (env, lp1, "basinv.lp", NULL);
      CPXwriteprob (env, lp, "basinvA.lp", NULL);
      __asErtm(!status,errmsg);
    log_ln("Turning off preprocessing...");
    status = CPXsetintparam (env, CPX_PARAM_PREIND, 0);
    assertm(!status,"Could not cancel preproc");
    status = CPXlpopt(env,lp);
    if (status) CPXwriteprob (env, lp, "HelpFailed.lp", NULL);
    assertm(!status,
      "Failed to prim. optimize. See HelpFailed.lp");
    status = CPXsetintparam (env, CPX_PARAM_PREIND, 1);
    assertm(!status,"Could not set preproc");
    }
//    if (cplexB.getNiterations())
//      cout << "N iter " << cplexB.getNiterations() << endl;
    bi[i].resize(Dim());
    status = CPXgetpi (env, lp, &bi[i][0], 0, Dim()-1); 
    assertm(!status,
      "CPX: getting multipliers for the basis inverse.");
  }


  CPXchgobj (env, lp,
    (int)ibas.size(), &ibas[0], &oldobj[0]);
  for (i=0;i<NCols();++i)
    if (CPX_BASIC != cstat[i])
      ChangeVarBnds(i, lb[i], ub[i]);
}

// NO FACILITY TO ADD MANY AT A TIME -- wantn't to.
// THE SLACK FOR THE CUT ADDED SEPARATELY
void WrpCpx::AddRow // equality
(int nnz,int *rmatind, double *rmatval,double lub) {
//  int i;
  int rmatbeg = 0;
/*  status = CPXcheckaddrows
    (env, lp, 0, 1, nnz, &lub,
    NULL, &rmatbeg, rmatind, rmatval,
    NULL, NULL);
  assertm(!status,"Failed to check add cut."); */
  status = CPXaddrows
    (env, lp, 0, 1, nnz, &lub,
    NULL, &rmatbeg, rmatind, rmatval,
    NULL, NULL);
  assertm(!status,"Failed to add cut.");
}

// DELS ALSO THE SLACK
void WrpCpx::DelRow(int ir) {
  status = CPXdelrows (env, lp, ir, ir);
  assertm(!status,"CPX: Del row");
}

void WrpCpx::DelCol(int ic) {
  status = CPXdelcols (env, lp, ic, ic);
  assertm(!status,"CPX: Del col");
}

void WrpCpx::ChangeRHS(int ir,double lub) {
  status = CPXchgcoef (env, lp, ir, -1, lub);
  assert(!status);
}

void WrpCpx::ChangeVarBnds(int ic,double lb,double ub) {
  char bndC='L';
  status = CPXchgbds(env,lp,1,&ic,&bndC,&lb);
  assert(!status);
  bndC = 'U';
  status = CPXchgbds(env,lp,1,&ic,&bndC,&ub);
  assert(!status);
}

void WrpCpx::GetVarBnds(int ic,double &lb,double &ub) {
 // char bndC='L';
  status = CPXgetlb(env,lp,&lb,ic,ic);
  assert(!status);
//  bndC = 'U';
  status = CPXgetub(env,lp,&ub,ic,ic);
  assert(!status);
}

void WrpCpx::LoadBasis(Vector<LP::ColStatus> &cstat) {
    assert(cstat.size() >= NCols());
    Vector<int> cs((int)cstat.size());
    Vector<int> rs(Dim());
    fill_n(rs.begin(),rs.size(),CPX_AT_LOWER);
    int i;
    for (i=NCols()-1;i>=0;--i)
    switch(cstat[i]) {
    case LP::atLower: cs[i] = CPX_AT_LOWER; break;
    case LP::basic: cs[i] = CPX_BASIC; break;
    case LP::atUpper: cs[i] = CPX_AT_UPPER; break;
    case LP::free_super: cs[i] = CPX_FREE_SUPER; break;
    default: assert(false);
    }
    status = CPXcopybase (env, lp, (&cs[0]), &rs[0]);
      if ( status > 0 ) {
         char  errmsg[1024];
         fprintf (stderr, "Could not set basis.\n");
         CPXgeterrorstring (env, status, errmsg);
         fprintf (stderr, "%s", errmsg);
      }
    assertm(!status,"CPX: could not set basis");
}

void WrpCpx::GetObjCoef(int ic,double & c) {
  c=GetLPCoef(-1, ic);
}
void WrpCpx::SetObjCoef(int ic,double c) {
  status = CPXchgcoef (env, lp, -1, ic, c);
  assertm(!status, "CPX: chg obj coef");
}

void WrpCpx::GetFullCol
(int j,Vector<double> &col,double &obj) {
  if (j<0 || j>NCols()) {
    status = CPXgetrhs (env, lp, &col[0], 0, Dim()-1);
    obj = 0; return;
  }

  col.clear(); col.resize(Dim()); // zeros
  int nzcnt, surplus;
  int cmatbeg;
  Vector<int> cmatind(Dim());
  Vector<double> cmatval(Dim());
  status = CPXgetcols
    (env, lp, &nzcnt, &cmatbeg, &cmatind[0],
    &cmatval[0], Dim(), &surplus, j, j);
  assert(surplus >= 0);
  for (int i=0;i< nzcnt; ++i)
    col[cmatind[i]] = cmatval[i];

}

double WrpCpx::GetLPCoef(int i, int j) {
  double t;
  status = CPXgetcoef(env,lp,i,j,&t);
  assert(!status);
  return t;
}

void WrpCpx::WriteModel(char *fln) {
     status = CPXwriteprob (env, lp, fln, NULL);
}



opt::OptContainer WrpCpx::Options() {
	opt::OptContainer oc;
	oc
	<< opt::MakeOpt(&epsOpt, 1e-8,
									"epsOpt", "Optimality tolerance in simplex method (for reduced costs). Default 1e-6 in CPLEX but 1e-8"
									" helps to get exact LP values on big/small-item instances like 402_*")
	<< opt::MakeOpt(&epsFeas, 1e-7,
									"epsFeas", "Feasibility tolerance in simplex method (for constraints). Default was 1e-6 in CPLEX")
	<< opt::MakeOpt(&nMIPThreads, 4,
									"nMIPThreads", "Max N threads when solving MIP (0: auto)")
	;
	return oc;
} //____________________________________________________

opt::OptSection WrpCpx::opt
("CPLEX", "CPLEX Optimizer",
 WrpCpx::Options(), opt::SolverCfg(), 6000);

double WrpCpx::epsOpt = 1e-8;
double WrpCpx::epsFeas = 1e-7;
int WrpCpx::nMIPThreads = 1;

}                                                // SS_END_NAMESPACE__
