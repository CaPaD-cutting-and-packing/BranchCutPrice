// probl_pmp1.cpp: Implementierung der Klasse PMP1.
//
//////////////////////////////////////////////////////////////////////


#include "stdafx.h"
#include "probl_csp1.h"
#include "probl_pmp1.h"
#include "solver.h"
#include "lasthdr.h"


namespace bgn {                     // SS_BEGIN_NAMESPACE__

PMP1::PMP1(Env *_penv, const char *_in,int _i,const char *_n)
   :CSP1(_penv, _in, _i, _n), csp(_penv, _in, _i, _n),
   // m=0 before SSVC
   ssvc(this, b, lpd, -1e100, // no bounds
    1e100)
{ }
//    :CSP1(_penv,_in,_i,_n) { }

void PMP1::Init() {
  CSP1::Init();
  // Create another object equal to this:
  csp = *this;
  csp.zi = cspIP;
  csp.lpb = cspLP; // is this enough for rounding?

  b.resize(m+1); // reserved FMax(100,m*2) or so (reall)
  ba.resize(m+1);
  b_cg.resize(m+1);

//  assertm(deltaK > 1e50, "deltaK > 1e50 int the first version.");
  if (Solver::problemType != 6
    or Solver::minZ >= 1e50) // if not solving the whole
    deltaK = ceil(cspIP * deltaKpercent/100);
  else {
    deltaK = Solver::minZ - cspIP;
    log_n(1,"PMP WHOLE PROBLEM: Setting deltaK = "<<deltaK);
  }

  b_cg[m] = ba[m] = b[m] = cspIP + deltaK;
  // TO REMOVE NUMERICAL ... :
//  pF = double (int (pF/GetVEps()) ) * GetVEps();
//  pV = double (int (pV/GetVEps()) ) * GetVEps();
  if (0 != pF*pV)
    gcdPP = __gcd(pF, pV);
  else gcdPP = FMax(pF, pV);  ////////////////
  if (0==pF and 0==pV)
  {pF=1; pV=1;}
  dbg_outn(4, " GCD(pF, pV)="<<gcdPP);
}

void PMP1::FillSigns() {
  CSP1::FillSigns(); // just like there
  validsign.resize(Dim());
//  fill_n(&validsign[0],m,1); // Ax = b, optimal simplex multipliers arbitrary
  validsign[m] = -1; // 1x <= K
//  dbg_outn(1," WARNING: Ax>=b");
}

void PMP1::FFSBasis
(const PieceContainer &pc__,ColumnList &bas) {
  int i;
  ColumnList::iterator ic;
  if (//deltaK<1e50// and not patBest.empty()
    nStartBasis
    ) {
/*    Column col(Dim()); // infeas col: seems not so super
    for (i=0;i<m;++i) col.PushID(i,int(b[i])); // or 0 last entry?
    col.PushID(m,-1);
    col.SetObj(100000000.0 * cspIP); // ???
    col.d = 1; // slack, not to branch on
    col.fInfeas = 1; //infeas
    bas.Add(col); */
  assertm(not patBest.empty(),
    "Need the best CSP solution for start.");
    for (i=0; i<xiBest.size();++i) {
      Column col;
      MakeColumn(&col, &patBest[i]);
      bas.Add(col);
    }
  for_each_in(bas.cl,ic,) {
    ic->Sort();
//    if (m!=ic->id.back().i)
    //  ic->PushID(m,1); // for K-constr
//    SetObj(&*ic);
    int lub = INT_MAX;
    int i;
    for (i=0;i<ic->id.size();++i)
      if (double(lub)*ic->id[i].d > b_cg[ic->id[i].i])
        lub = int(b_cg[ic->id[i].i] / ic->id[i].d);
    assertm(lub > 0, " ??? upper bound for column from CSP = 0!!! unproper ");
    ic->SetObj(pV + pF / lub);
  }
  }
  return;
//  assert(nStartBasis
  // STILL ARBITRARY (deltaK = inf)
  // but adding (m,1):
/*  int nsb = nStartBasis;
  nStartBasis = 2; // !!!!!!!!!!!!!! infeas otherwise
  CSP1::FFSBasis(pc__,bas);
  nStartBasis = nsb;
  for_each_in(bas.cl,ic,) {
    ic->Sort();
    if (m!=ic->id.back().i)
      ic->PushID(m,1); // for K-constr
//    SetObj(&*ic);
    int lub = INT_MAX;
    int i;
    for (i=0;i<ic->id.size();++i)
      if (double(lub)*ic->id[i].d > b_cg[ic->id[i].i])
        lub = int(b_cg[ic->id[i].i] / ic->id[i].d);
    ic->SetObj(pV + pF / lub);
  }*/
}


void PMP1::GenColPure
  (ColSet &cs,const d_vec &d) {
  int i;
  bb.fHeur = false;

  for (i=0;i<m;++i) {
    bb.pieces__[i].d = d[i];
    bb.pieces__[i].b = int(b_cg[i]);
  }
  bb.zLowerInitial = -INFINITY__;
  redCostBest = -1e100;
  bb.eps = GetBBEps();

  int q0 = 1, q0maxa;
  int q0lim = (q0max < INT_MAX) ? (int)q0max : INT_MAX-1;
  if (double(q0lim) > round(b_cg[m]))
    q0lim = int(round(b_cg[m])); // local bound
  double rc = 0;
  dbg_outn_(3,'c');

  for (;q0 <= (int)q0lim; ++q0) { // q0max changed in node init
    dbg_outn_(4,'q');
    for (i=0;i<m;++i)
      bb.pieces__[i].b = int(b_cg[i]/q0);
    bb.zMin = -rc - d[m]+pF/q0max + pV;

    Column col;
    col.PushID(m,1); // the K constr. If K=inf better weg?
    ColSet colset;
    bb.colBest = &col;
    bb.colsetRes = &colset;

    bb.Run();
    if (not (bb.found)) break; // cannot be better

    q0maxa = (int)q0max; // CALC THE ACTUAL q0(a) to calc the rc(a)
    for (i=0;i<col.id.size();++i)
      if (double(q0maxa*col.id[i].d) > b_cg[col.id[i].i])
        q0maxa = int(b_cg[col.id[i].i] / col.id[i].d);
    q0 = q0maxa;
    if ( pV + pF/q0 - bb.z - d[m] >= rc)
      // no actual improvement
      continue; // and what in LP opt? How long?
    // MAY BE not continue but check whether
    // a larger q0next is possible at next iteration?

    rc = pV + pF/q0 - bb.z - d[m];
    col.SetObj(pV + pF/q0); // integer bounding is necessary
    clBest.clear();
    clBest = col;
    cs.cs.clear();
    cs.Add(col);
  }
  if (cs.cs.empty()) return;

  redCostBest = rc;
//  UpdateLagrBnd(d);  //  !!!!!!!!!!
}

bool PMP1::ConstructIntSol() {
/*
Just change CSP1::xUse
and recalc the original objective.
Later: save solutions with different \Delta K
*/
//  double txUse = CSP1::SVC::patternUseRatio;
  int i;
//  CSP1::SVC::patternUseRatio = 1; // ???
//  csp.fFastRounding = fFastRounding;
//  csp.pat = pat;
  for (i=0; i<pat.size(); ++i)
    pat[i].SetObj(1); // !!!
//  csp.lpx = lpx;
//  csp.lpd = lpd;
  // BUT WE SHOULD ALSO EXTRACT DIFFERENT SOLUTIONS
  // NOT ONLY ACCORDING TO CSP BOUND. THUS:
//  csp.zi = 1e100;
  /*csp.*/CSP1::ConstructIntSol();
//  CSP1::SVC::patternUseRatio = txUse;
//  ExtractBestSolution(&csp);
  return Optimum();
}

// cspIP must be set before
void PMP1::ExtractBestSolution(Problem* pr) {
  ExtractSolution(pr->patBest, pr->xiBest);
}

void PMP1::ExtractSolution(Vector<Pattern>&pa_,d_vec& xi_) {
  set<Pattern> pats;
  set<Pattern>::iterator it;
  int i; double sx=0;
  for (i=0;i<xi_.size();++i)
  if (xi_[i] > 1e-6) {  // no zero multiplies
    sort(pa_[i].ix.begin(), pa_[i].ix.end());
    pa_[i].x = round(xi_[i]); // round: saves from rnd errs?
    sx += round(xi_[i]);
    it = pats.find(pa_[i]);
    if (pats.end() == it)
      pats.insert(pa_[i]);
    else
      it->x += round(xi_[i]);
  }
  assertm(sx >= cspLP - GetXEps(), "PMP rounding: sx="
    <<sx<<" < cspLP="<<cspLP);
//  assert(sx <= cspIP); // no
//  assert(fabs(sx - pr->zi) < 1e-6);
  dbg_outn(4, "Solution found with "<<pats.size()<<" diff patterns and "
    <<sx<<" stock sheets used.");
  double znew = pF * (double)pats.size() + pV * sx;
  dbg_outn_(4,"ObjVal = "<<znew);
  // Sequencing to get MOSP:
  if (znew < zi-1e-6 // better sol
    and sx <= cspIP + deltaK + GetXEps()) { // and feasible
    dbg_outn(2, "Better solution found with "<<pats.size()<<" diff patterns and "
      <<sx<<" stock sheets used.");
    dbg_outn_(2,"ObjVal = "<<znew);
    UpdateHeurBnd(znew);
    patBest.clear();
    xiBest.clear();
    for_each_in(pats,it,) {
      patBest.push_back(*it);
      xiBest.push_back(it->x);
    }
    //ControlSolution();
    Vector<double> s_aij(m); // fillled 0 ?
///    double zzz = 0;
    int i;
    fill_n(&s_aij[0],m,0);
    for (i=0;i<patBest.size(); ++i) {
      for_each_in(patBest[i].ix,iix,Pattern::iterator)
        s_aij[iix->i] += iix->x * xiBest[i];
//      zzz += patBest[i].GetObj() * xiBest[i];
    }
//    assert(zzz <=);
    for (i=0;i<m;++i)
      assertm(s_aij[i] >= pc[i].b,
      "s_aij[i] - pc[i].b = "<< s_aij[i] - pc[i].b
      <<", pc[i].b = "<<pc[i].b);
    if (fModelEquality)
      for (i=0;i<m;++i)
      assertm(s_aij[i] == pc[i].b,
      "s_aij[i] - pc[i].b = "<< s_aij[i] - pc[i].b
      <<", pc[i].b = "<<pc[i].b);
    if (fSequence && b.size()) goto SSVC; // Check all best
  } // CONTROL ?
  if (fSequence && rndSeq <= seqFreq && b.size()
      and sx <= cspIP + deltaK + GetXEps()  // no, only the same material
      // znew < zi - 1e-6 + 3 // not significantly worse sol. ???
//    and sx <= cspIP + deltaK + GetXEps() // and feasible
 ) {
SSVC:
    if (ssvc.b0.empty())
      (ssvc.b0) = b; // SSVC
    ssvc.d0.resize(m); // Init by li's ? In MSVC.
    ssvc.m = m;

    ssvc.pat0.assign(pats.begin(), pats.end());
    ssvc.x0.resize(pats.size());
    // ssvc.lpd = 
    for (i=0;i<pats.size();++i)
      ssvc.x0[i] = ssvc.pat0[i].x;
    ssvc.Run();
  }
}
  // returns 1 if opt
  // called from SVC
bool PMP1::ExtractSolutionPart
(Vector<Pattern>& pa_, d_vec& xi_) {
  Vector<Pattern> pa1 = pat; // rounded;
  d_vec xi1 = xi;
  pa1.insert(pa1.end(), pa_.begin(), pa_.end());
  xi1.insert(xi1.end(), xi_.begin(), xi_.end());
  ExtractSolution(pa1, xi1);
  return Optimum();
}


bool PMP1::CompleteIntSol() {
  int i;
  double LL=0;
///// Check if all fit in a single rod:
  for (i=0;i<m;++i) LL += (double(pc[i].l) * br[i]);
  if (0 == LL) {
    ExtractSolution(pat, xi);
    return true;
  }
  if (OUTP_LEV__ >=4)
    log__(" PMP1: Adding to CSP1 rounded value " << zr);
  SVC svc(this, br, lpd, -1e100, cspIP +deltaK +1 - zr);
  int iterMax__ = SVC::iterMax;
//  if (fFastRounding)
    SVC::iterMax = IMin(SVC::iterMax,4);
  // double SVCres=
    svc.Run();
  SVC::iterMax = iterMax__;
//  if (Optimum())
  return true;
} //////////////////////////////// CSP1::CompleteIntSol



void PMP1::MakeColumn(Column * c,Pattern * p) {
    c->clear();
    for_each_in(p->ix,iix,Pattern::iterator)
      c->PushID(iix->i,iix->x);
    c->PushID(m,1); //deltaK
    c->Sort();
//    c.id.insert(c.id.end(),id.begin(),id.end());
    c->ofc = p->ofc;
    c->GetAddiInfo() = p->GetAddiInfo();
}


void PMP1::MakePattern(Pattern *p, Column *c) {
    p->clear();
    for_each_in(c->id,iid,Column::iterator)
      if (iid->i < m)
        p->PushIX(iid->i, (int)iid->d);
    p->SetObj(c->GetObj());
    p->GetAddiInfo() = c->GetAddiInfo();
} ////////////

void PMP1::PrintColumn(ostream& os,Column* c) {
  int i;
//  if (c->GetCutSlackCut())
    //os << "Cut slack: " << c->GetCutSlackCoef()<< ' ';
      for (i=0; i<c->id.size(); ++i)
        os << c->id[i].i  <<':'<<c->id[i].d<<' ';
      os << "obj " << c->GetObj();
}

void PMP1::PrintLog(ostream &os) {
  if (!fSequence) return;
  MSVC::SolContainer::iterator its;
  MSVC::SolContainer::iterator itsOS, itsND, itsEV;
  itsOS = itsND = ssvc.sols.end();
  itsEV = ssvc.sols.begin(); // equivocal
  assert(not ssvc.sols.empty());
  // SELECTING THE 3 EXTREME SOLUTIONS:
  for_each_in(ssvc.sols, its,) {
    if (ssvc.nOpenMaxMin == its->no) { // minOS with best ND
      if (ssvc.sols.end() == itsOS)
        itsOS = its;
      else {
        if (its->nd < itsOS->nd)
          itsOS = its;
      }
    }
    if (ssvc.nDiffPatMin == its->nd) { // minND with best OS
      if (ssvc.sols.end() == itsND)
        itsND = its;
      else {
        if (its->no < itsND->no)
          itsND = its;
      }
    }
    if (its->nd * ssvc.nOpenMaxMin + its->no * ssvc.nDiffPatMin
      < itsEV->nd * ssvc.nOpenMaxMin + itsEV->no * ssvc.nDiffPatMin)
      itsEV = its;
  }
  // + STATISTICS: (AVE ND...)
  sta[0] ++; // N tests
//  sta[1] += timer.userTime();
  if ((itsND->nz - zi)/zi > sta[2])
    sta[2] = (itsND->nz - zi)/zi;
  sta[3] += (itsND->nz - zi)/zi;
  int ii = 3;
  sta[++ii] += itsND->nz;
  sta[++ii] += itsND->no;
  sta[++ii] += itsND->nd;
  sta[++ii] += itsND->iter;
  sta[++ii] += itsOS->nz;
  sta[++ii] += itsOS->no;
  sta[++ii] += itsOS->nd;
  sta[++ii] += itsOS->iter;
  sta[++ii] += itsEV->nz;
  sta[++ii] += itsEV->no;
  sta[++ii] += itsEV->nd;
  sta[++ii] += itsEV->iter;
  os.precision(7);
/*  os 
    // Quality of the BCP sol:
    << " Sol0: (" << msvc.sol0.no <<
    ' ' << msvc.sol0.nd << ' ' << msvc.sol0.nz << ')'
    ;*/
  for_each_in(ssvc.sols, its,)
    os
    << " ("<< its->no <<
    ' '<< its->nd << ' ' <<its->nz << ')'
    ;
}

void PMP1::PrintStat(ostream& os) {
  if (!fSequence or not sta[0]) return;
  int i;
  os << "SSVC stat. N="<<sta[0]<<", t,gapmax,gapave,(nz no nd iter) for minND,OS,EV\n";
  os << sta[1]/sta[0] << ' ';
  os << sta[2] << ' ';
  for (i=3;i<16;++i)
    os << sta[i]/sta[0] << ' ';
//  os << '\n';
}



double PMP1::deltaKpercent=1; // upper bound on addi stock.
   double PMP1::pF; // fixed (setup) cost
   double PMP1::pV; // variable cost (per stock sheet)

   int PMP1::nStartBasis;

   double PMP1::CSPOutpLev;
   double PMP1::CSPTiLim;
   double PMP1::BPPTiLim;

   bool PMP1::fSequence=0;
   double PMP1::seqFreq=0.1;

   double PMP1::outputLevel;

opt::OptContainer PMP1::Options() {
  opt::OptContainer oc;
  oc
    << opt::MakeOpt(&deltaKpercent, 1,
      "deltaKpercent", "Upper bound on addi stock = rndUP(dK%*cspIP).")
    << opt::MakeOpt(&pF, 1,
      "pF", "fixed (setup) cost")
    << opt::MakeOpt(&pV, 0,
      "pV", "variable cost (per stock sheet). NOTE: to avoid numerics, use integers")
    << opt::MakeOpt(&nStartBasis, 1,
      "nStartBasis", "0: empty, otherwise from CSP optimum")
    << opt::MakeOpt(&CSPOutpLev, 1,
      "CSPOutpLev", "outplevel when solving CSP/BPP")
    << opt::MakeOpt(&CSPTiLim, 30,
      "CSPTiLim", "")
    << opt::MakeOpt(&BPPTiLim, 30,
      "BPPTiLim", "")
    << opt::MakeOpt(&fSequence, 0,
      "fSequence", "Whether to sequence soltuons for min open stacks")
    << opt::MakeOpt(&seqFreq, 0.1,
      "seqFreq", "probability how often to sequence a solution")
/*    << opt::MakeOpt(&nStartBasis, 1,
      "nStartBasis",
      "0:FFD, 1: Greedy, 2: diag. mattrix with large obj coefs")
    << opt::MakeOpt(&nStepsMin0, 8192,
      "nStepsMin0",
      "Initial min N steps in B&B with cuts")
    << opt::MakeOpt(&nStepsMinInc, 1.01,
      "nStepsMinInc", "Incr ratio (each generation)")
    << opt::MakeOpt(&deps, 1e-8,
      "deps", "eps for dual multipliers")
    << opt::MakeOpt(&bb_eps, 1e-8,
      "bb_eps", "eps for b&b col gen")
    << opt::MakeOpt(&RPEMax, 10,
      "RPEMax", "FMax residual problem extensions")*/
    << opt::MakeOpt(&outputLevel, DEF_OUTP_LEVEL,
      "outputLevel", "0-5");

  return oc;
} //____________________________________________________
opt::OptSection PMP1::opt
  ("PMP1", "The 1D Pattern Minimization",
  PMP1::Options(), opt::SolverCfg(), 500);



}                                                // SS_END_NAMESPACE__

