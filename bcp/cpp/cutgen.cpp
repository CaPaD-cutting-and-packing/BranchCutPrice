// FILE: cutgen.cpp
// SuperAdditive Cut Generator
// AUTHOR: Gleb <belov@web.de>

#include "stdafx.h"
#include "cutgen.h"
#include "lasthdr.h"

// Be able to print everything
// SE = slack elimination

namespace bgn {                     // SS_BEGIN_NAMESPACE__

// With col gen, cols with bounds are forbidden (also lb)
// ; they are in LP already => not consider cut calc for
// ATTENTION: clear fVisited of all involved before recurs

#define CGeps (-1e-6) // care! numer. diffic. at cg if sma
#define MIeps 0.0 // 1e-6
// VERY IMPORTNAT TO INCREASE RHS raw sum by 3e-6
// because we choose 1e-6<=alfa<=0.999999
// or not ? don't write (a + b? c:d) which is equiv to ((a+b)?..)

#undef dbgpr
#define dbgpr(n,e) //cout<<e;

// : here, c is the rhs
double SACutSE::CalcRHS__(d_vec &b) {
  if (fNodeVisited) {
//    cout << "rhs for cut "<<this<<": pre-claculated "
//  <<value<<endl;
    return rhs;
  }
  fNodeVisited = true;
  sua = suaSE = 0;
  assert(b.size() >= u.size());
  assert(u.size());
  int i;
  for (i=0;i<u.size();++i) {
    sua += u[i] * b[i];
//    suaSE += uSE[iid->i] * iid->d;
  }
  value = sua; value2 = suaSE;
  for_each_in(dep, id, DepContainer::iterator) {
    value += id->u * id->c->CalcRHS__(b);
//    value2 += id->uSE * id->c->CalcLRHS__(c);
  }
// variable bounds:
  BndContainer::iterator ib;
  for_each_in(bnds,ib,) {
    assert(rawCoef.size() > ib->j);
    value -= rawCoef[ib->j] * ib->bnd;
  }
  value = F_Alfa(value - (cutType?(CGeps*2):(MIeps*2)));
  for_each_in(bnds,ib,)
    value += GetCoef(ib->j) * ib->bnd;
  rhs = value; /// +1e-6 ?
//    cout << "rhs for cut "<<this<<": claculated "<<value<<endl;
  return rhs;
}
double SACutSE::F_Alfa(double v) {
  if (cutType) {
//    if (fabs(frac(v)-CGeps-1) < 1e-6) cout << '!';
    return floor(v-CGeps); // Gomory fractional but not for
  }
    // cp22 level cut calc.
  v -= MIeps;
  double abru=floor(v);
  double frac=v-abru;
  if (frac > alfa)
    return abru + (frac-alfa) / (1.0-alfa);
  return abru;
} //____________________________________________________
double SACutSE::F_Bar(double v) {
  if (!cutType) {
//    return F_Alfa(v); // ...
    v -= MIeps;
    if (v<0) return v/(1.0-alfa);
    return 0;
  } else // CG cuts
    return F_Alfa(v);
}
// for cut checking:
double SACutSE::F_Alfa(double alfa,double v,int cutType) {
  if (cutType) return floor(v-CGeps); // Gomory fractional but not for
    // cp22 level cut calc.
  v -= MIeps;
  double abru=floor(v);
  double frac=v-abru;
  if (frac > alfa)
    return abru + (frac-alfa) / (1.0-alfa);
  return abru;
} //____________________________________________________
double SACutSE::F_Bar(double alfa,double v,int cutType) {
  if (!cutType) {
//    return F_Alfa(alfa, v,cutType); // ...
    v -= MIeps;
    if (v<0) return v/(1.0-alfa);
    return  0;
  } else // CG cuts
    return F_Alfa(alfa,v, cutType);
}
void SACutSE::CalcIntermSums(Column * col) {
  Column::iterator iid; sua = suaSE = 0;
  for_each_in (col->id,iid,) {
    sua += u[iid->i] * iid->d;
// SUPPOSING original constr. have integer slacks ?:
//    suaSE += uSE[iid->i] * iid->d;
  }
}
// NO CUT SLACKS POSSIBLE:
double SACutSE::CalcUsingIntermSums() {
// ex. during branching
  if (fNodeVisited) return value;
  fNodeVisited = true;
// KEEPING sua for branching
  value = sua; value2 = suaSE;
  dbgpr(5," Cut "<<this<<": sua="<<sua);
  for_each_in(dep, id, DepContainer::iterator) {
    value += id->u * id->c->CalcUsingIntermSums();
//    value2 += id->uSE * id->c->CalcUsingIntermSums();
  }
  dbgpr(5," Cut "<<this<<": sua="<<value);
  value = F_Alfa(value) + value2;
  dbgpr(5," value="<<value);
  return value;
} //____________________________________________________
// only col gen, no cols with upper bounds!
double SACutSE::Calc__(Column *c) {
  if (fNodeVisited) return value;
  fNodeVisited = true;
  assertm(!c->id.empty(),"no cut slacks please");
  CalcIntermSums(c);
// KEEPING sua for branching
  value = sua; value2 = suaSE;
  for_each_in(dep, id, DepContainer::iterator) {
    value += id->u * id->c->Calc__(c);
//    value2 += id->uSE * id->c->Calculate__(c);
  }
  value = F_Alfa(value) + value2;
  return value;
}
double SACutSE::GetCutSlackCoef__(LPCut *pcut) {
//  if (fNodeVisited) return value; // HERE NO USE ONLY .
//  fNodeVisited = true;
  dbgpr(5," cut slack "<<pcut);
  if (pcut == this)
    return (value = GetSlackCoef());
  const Dependence
    * pd = invCuts.Find(Dependence(pcut)); // cmp pcut
  if (NULL == pd) {
    dbgpr(5," val=0");
    return (value=0); // not involved
  }
  assertm(pd->u > -1e50, "Cut slack coef must be precalculated!");
  return (value = pd->u);
}

// coef in this cut for the slack of cut *pcut
double SACutSE::CalcCutSlackCoef__(LPCut *pcut) {
  if (fNodeVisited) return value;
  fNodeVisited = true;
  dbgpr(5," cut slack "<<pcut);
  if (pcut == this)
    return (value = GetSlackCoef());
  const Dependence
    * pd = invCuts.Find(Dependence(pcut)); // cmp pcut
  if (NULL == pd) {
    dbgpr(5," val=0");
    return (value=0); // not involved
  }
  if (pd->u > -1e50) {
    dbgpr(5," val precalc="<<pd->u);
    return (value=pd->u);
  }
  value = value2 = 0;
  for_each_in(dep, id, DepContainer::iterator) {
    double sc = id->c->CalcCutSlackCoef__(pcut);
    value += id->u * sc;
//    cout << " SL_COEF["<<id->c<<"]="<<sc<<'*'<<id->u;
//    value2 += id->uSE * id->c->Calculate__(c);
  }
  if (pcut->IntegerSlack())
    value = F_Alfa(value) + value2;
  else
    value = F_Bar(value) + value2;
//  cout << "\n => CUT ["<<this<<"]: SLCOEF["<<pcut<<"]="<<value
//    <<", pd->u was "<<pd->u<<endl;
  return (pd->u = value);
}
void SACutSE::ClearNonRec() {
  CheckID();
  fNodeVisited = false;
  value = value2 = 0; // value2 needed for ApprError
} //____________________________________________________
double SACutSE::CalcApprCoef(int k) {
  if (fNodeVisited) return value;
  fNodeVisited = true;
  for_each_in(dep, id, DepContainer::iterator)
    value += (id->u + id->uSE)
      * id->c->CalcApprCoef(k);
  value += u[k] + uSE[k];
  return value;
} //____________________________________________________
double SACutSE::CalcApprErrorL() { // ind. proof
  if (value) return value;
  for_each_in(dep, id, DepContainer::iterator) {
    if (id->u >= 0)
      value += (id->u + id->uSE)

        * id->c->CalcApprErrorL();
    else
      value += (id->u + id->uSE)
        * id->c->CalcApprErrorU();
  }
  value -= alfa; //  - 1e-15; // minus ???????? !!!!!!! cause
  return value;
} //____________________________________________________
double SACutSE::CalcApprErrorU() {
  if (value2) return value2; // visited flag
  value2 = 0; // 1e-15;
  for_each_in(dep, id, DepContainer::iterator) {
    if (id->u >= 0)
      value2 += (id->u + id->uSE)
        * id->c->CalcApprErrorU();
    else
      value2 += (id->u + id->uSE)
        * id->c->CalcApprErrorL();
  }
  // No alfa SUBTRACTING
  return value2;
} //____________________________________________________
void SACutSE::ClearSums() { // not used
  assert(false);
/*  if (fNodeVisited) return;
  fNodeVisited = true;
  sua = suaSE = 0;
  for_each_in(dep, id, DepContainer::iterator)
    id->c->ClearSums();*/
}
void SACutSE::AddToSums(int k, int x) {
  if (fNodeVisited)
    return;
  fNodeVisited = true;
  sua += u[k] * x;
//  suaSE += uSE[k] * x;
  for_each_in(dep, id, DepContainer::iterator)
    id->c->AddToSums(k,x);
} //____________________________________________________
void SACutSE::Sub1FromSums(int k) {
  if (fNodeVisited)
    return;
  fNodeVisited = true;
  sua -= u[k];
//  suaSE -= uSE[k];
  for_each_in(dep, id, DepContainer::iterator)
    id->c->Sub1FromSums(k);
} //____________________________________________________


void SACutSE::CalcConstTerms(Column* c) {
  if (fNodeVisited)
    return;
  fNodeVisited = true;
  sua0 = suaSE0 = 0; // - CGEps*2.0; // !!! would it help??
// ----- Fixing all elements set in *c:
  for_each_in(c->id,iid,Column::iterator) {
    sua0 += u[iid->i] * iid->d;
//    suaSE0 += uSE[iid->i] * iid->d;
  }
// ----- Fixing all constant cuts:
// THIS NOT.
/*  for_each_in(dep, id, DepContainer::iterator)
    if (id->c->is_const()) {
      sua0 += id->u * id->d;
      suaSE0 += uSE[iid->i] * id->d;
    }*/
  for_each_in(dep, id, DepContainer::iterator)
    id->c->CalcConstTerms(c);
} //____________________________________________________
void SACutSE::AssignConstTerms() {
  if (fNodeVisited) return;
  fNodeVisited = true;
  sua = sua0; suaSE = suaSE0;
  for_each_in(dep, id, DepContainer::iterator)
    id->c->AssignConstTerms(); // may be repeating !

} //____________________________________________________
void SACutSE::Print(ostream&os) {
//  ostream& os = cout;
  int i;
  DepContainer::iterator id;
  os<<"alfa"<<alfa<<" u's:\n";
  for (i=0;i<u.size();++i)
    os<<u[i]<<' ';
  os<<"Depend:\n";
  for_each_in(dep,id,)
    os<<'{'<<id->c<<':'<<id->u<<"} ";
  os<<"Involved:\n";
  for_each_in(invCuts,iic,Pool<Dependence>::iterator)
    os<<'{'<<iic->c<<':'<<iic->u<<"} ";
/* os<< "\nuSE's:\n";
  for (i=0;i<uSE.size();++i)
    os<<uSE[i]<<' ';
  os<<"Depend:\n";
  for_each_in(dep,id,)
    os<<'{'<<id->c<<':'<<id->uSE<<"} ";*/
  BndContainer::iterator ib;
  os<<"Var bnds:\n";
  for_each_in(bnds,ib,)
    os<<'{'<<ib->j<<':'<<ib->upper<<':'<<ib->bnd<<"} ";
  os<<"Coefs:\n";
  for (i=0;i<rawCoef.size();++i)
    os <<i<<':'<<GetCoef(i)<<' ';
  os << "\nrhs:"<<rhs<<'\n';
}

void SACutSE::Number(int &nn) {
  if (-1 != no) return; // incr. order, involved first
  for_each_in(dep,id,DepContainer::iterator)
    id->c->Number(nn);
  no = nn++;
}

void SACutSE::ProduceListOfInvolved
(Vector<LPCut*>& cuts) const {
  if (fNodeVisited) return;
  fNodeVisited = true;
  for_each_in(dep,id,DepContainer::const_iterator)
    id->c->ProduceListOfInvolved(cuts);
  cuts.push_back((LPCut*)this); // incr. order, involved first
}

bool SACutSE::operator<(const LPCut &cut) const
{
  if (Type()!=cut.Type())
    return Type()<cut.Type();
  // otherwise lexicogr. compare:
  int i;
  /// CHANGING COMPARISON TO FLOATS BECAUSE OF MESS WITH eps:
  /// double du;
  float du;
  const SACutSE * pcut =
    dynamic_cast<const SACutSE *>(&cut);
  SACutSE::DepContainer::iterator id1,id2;
  for (i=0;i<u.size(); ++ i) {
    du = float(u[i]) - float(pcut->u[i]);
    if (fabs(du) > 1e-8)
      return du < -1e-8;
  }
  return lexicographical_compare(
    dep.begin(),dep.end(),
    pcut->dep.begin(),pcut->dep.end());
} // compare bnds also ?

// no SE!

int SACutSE::GetNCoefs() { return (int)rawCoef.size(); }
/// calling only for iCol <= size
double SACutSE::Calc__(Column *c, int iCol) {
  assertm(!fSE,"Not programmed for slack elim!");
  assert(iCol>=0);
  if (GetNCoefs() > iCol)
    return GetCoef(iCol);
  CalcRawCoef(c,iCol);
  return GetCoef(iCol);
}
double SACutSE::GetCoef(int iCol) {
  assert(GetNCoefs() > iCol);
  const VarBnd *bnd = bnds.Find(VarBnd(iCol,true));
  if (bnd) // i.e. upper
    return - F_Alfa(-rawCoef[iCol] //- (cutType?(CGeps*2.0):(MIeps*2.0))
      );
  return F_Alfa(rawCoef[iCol]);
}
// can also calc cut slacks.
// returns the raw value.
double SACutSE::CalcRawCoef(Column *c,int iCol) {
  assert(not c->id.empty()); // no cut slacks
  CalcIntermSums(c);
// KEEPING sua for branching
  value = sua;// value2 = suaSE;
  for_each_in(dep, id, DepContainer::iterator) {
    value += id->u * id->c->Calc__(c,iCol);
//    value2 += id->uSE * id->c->Calculate__(c);
  }
  assert(rawCoef.size()==iCol); // calc only once???
  rawCoef.push_back(value);
  return value;
}

double SACutSE::CalcSlackValue
(i_vec &iNZ, d_vec &xNZ, map<LPCut*,double> &slVal) {
  map<LPCut*,double>::iterator
    iter = slVal.find(this);
  if (slVal.end() != iter)
    return iter->second;
  double sv = GetRHS();
  int i;
//  log__("Slack for cut "<<this<<": rhs="<<sv);
  for (i=0;i<iNZ.size();++i) {
    sv -= GetCoef(iNZ[i]) * xNZ[i];
//    log__(" x["<<iNZ[i]<<"]="<<xNZ[i]<<'*'<<GetCoef(iNZ[i]));
  }
  Pool<Dependence>::iterator idep;
  for_each_in(invCuts,idep,) {
     double sc = GetCutSlackCoef__(idep->c);
     double csv = idep->c->CalcSlackValue(iNZ,xNZ,slVal);
//     log__(" sv["<<idep->c<<"]="<<csv<<'*'<<sc);
     sv -= sc * csv;
  }
//  log_ln(" .="<<sv);
  if (GetSlackCoef()<0) sv = -sv;
  slVal.insert(make_pair(this,sv));
  return sv;
}

////////////////////////////////////////////////////////
////////////// IMPLEMENTATION OF CutSet ////////////////
////////////////////////////////////////////////////////
double CutSet::CalcCoefsUsingIntermSums() {
  ClearValues();
  double sum=0;
  for_each_in(dep, id, DepContainer::iterator)
    sum += id->u * id->c->CalcUsingIntermSums();
  return sum;
} //____________________________________________________
double CutSet::CalcApprCoefs(int k) {
  ClearValues();
  double sum=0;
  for_each_in(dep, id, DepContainer::iterator)
    sum += id->u * id->c->CalcApprCoef(k);
  return sum;
} //____________________________________________________
double CutSet::CalcApprError() {
  ClearValues();
  assert(invCuts.size()>=dep.size());
  double sum=0;
  for_each_in(dep, id, DepContainer::iterator)
    sum += id->u > 0 ?
      id->u * id->c->CalcApprErrorU() :
      id->u * id->c->CalcApprErrorL();
  return sum;
} //____________________________________________________
void CutSet::ClearValues() {
  for_each_in(invCuts, ic, CutList::iterator)
    (*ic)->ClearNonRec();
} //____________________________________________________
void CutSet::PrepareRecursion() {
//  for_each_in(dep, id, DepContainer::iterator)
//    id->c->Clear(); // ... not Prepare Rec.
    //ClearValues();
  for_each_in(invCuts, ic, CutList::iterator)
    (*ic)->ClearNonRec();
} //____________________________________________________
void CutSet::AddToSums(int k, int x) {
  PrepareRecursion();
  for_each_in(dep, id, DepContainer::iterator)
    id->c->AddToSums(k,x);
} //____________________________________________________
void CutSet::Sub1FromSums(int k) {
  PrepareRecursion();
  for_each_in(dep, id, DepContainer::iterator)
    id->c->Sub1FromSums(k);
} //____________________________________________________
void CutSet::CalcConstTerms(Column *c) {
  PrepareRecursion();
  for_each_in(dep, id, DepContainer::iterator)
    id->c->CalcConstTerms(c);
} //____________________________________________________
void CutSet::AssignConstTerms() {
  PrepareRecursion();
  for_each_in(dep, id, DepContainer::iterator)
    id->c->AssignConstTerms(); // may be repeating !
} //____________________________________________________
void CutSet::AddInvolved(CutList * cuts) {
//  for_each_in(*cuts,ic,) { (*ic)->fNodeVisited = true; }
  int nn=0;
  CutList::iterator ic;
  for_each_in(*cuts,ic,) { (*ic)->no = -1; }
  for_each_in(dep,id,DepContainer::iterator)
    id->c->Number(nn);
  invCuts.clear();
  for_each_in(*cuts,ic,)
    if ((*ic)->no != -1)
      invCuts.push_back(*ic);
} //____________________________________________________


}                                                // SS_END_NAMESPACE__
