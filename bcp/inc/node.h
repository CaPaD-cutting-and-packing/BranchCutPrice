/*****************************************************
                          node.h  -  description
                             -------------------
    begin                : Don Jan 9 2003
    copyright            : (C) 2003 by Gleb Belov
    email                : belov@math.tu-dresden.de
 *****************************************************/

#ifndef __NODE_H__32
#define __NODE_H__32

//#include "lpcol.h"
#include "branch.h"
#include "cutgen.h"

namespace bgn {                     // SS_BEGIN_NAMESPACE__

typedef Pool<SACutSE> CutContainer; // why here? for Node

  // try to leave only active cuts in a node?
struct Node {
    Vector<int> iBas, iUpper; // other vars at lower
    Vector<LPCut*> slkBas; // which cut slacks in bas
    typedef Vector<VarBnd> BndCont;
    BndCont bnds; // the branches
//    LPCut* brHP;  a branching hyperplane - now in cutcont1
    Pool<VarBnd> rcBnds; // red costs bounds
    typedef Pool<SACutSE> CutContainer; // why here? for Node
    // - storing SA cuts here because they are not problem-dep
    CutContainer cutcont; // local SA cuts
    typedef list<//my_auto_ptr<
      LPCut*> CutCont1;
    typedef CutCont1::iterator CutIter1;
//    mutable // because of auto_ptr
      CutCont1 cutcont1;
    void CheckNew() const
    { assert(cutcont1.empty()); } // no pointers in added nd
    CutList cuts; // cuts present in formulation, incl. branching HP
    Node * parent;
    Node * sonLeft, * sonRight; // left, right son

    double no; // number
    int depth;
    int nDesc; // n ACTIVE sons
    double xValue; // value of the branching var
    int nLeft; // N left branches (xj<=..) from the root including oneself

    double llv, llb;
    bool fLPOpt; // if llb can be trusted
    enum State {open,current,processed,fathomed}
      state; // pruned if processed & no sons
    bool Active() const
      { return (open==state or current==state); }
    bool Open() const { return (open==state); }
    Node() : no(0), depth(0), nDesc(0),
      parent(NULL), sonLeft(NULL), sonRight(NULL), nLeft(0),
      state(open) { }
    ~Node() { // deleting non-SA cuts:
      for_each_in(cutcont1, it, CutIter1) delete *it;
    }
//    Node(const Node& n) { *this = n; }
    // 1 of the following 2 is called after bounding:
    void MarkProcessed() { state = processed; }
    // No fathomed state: just delete
//    void Fathom() { state = fathomed; }
    /// For the case if we store all nodes in a tree:
    bool operator<(const Node & n) const {
      if (depth < n.depth) return 1;
      if (depth > n.depth) return 0;
      return (no > n.no); // last nodes first
    }
};

/// Incapsulate the solution tree:
class SolutionTree {
public:
    typedef list<Node> NodeContainer; // could be a multiset
    typedef NodeContainer::iterator node_iterator;
    typedef node_iterator INode;
    typedef node_iterator PNode;
		
    bool SearchEmpty() const { return tree.empty(); }
    INode ExtractBestOpen() {
      assert(not SearchEmpty());
      INode in = *(tree.begin());
      tree.erase(in);
      return in;
    }
    INode GetRoot()
    { assert(0==nodes.begin()->no); return nodes.begin(); }
    NodeContainer::size_type GetSize()
      { return nodes.size(); }
    bool GetLastOpen(INode & in) {
      ISrch is = tree.end();
      if (not tree.empty()) {
        -- is; in = *is;
        return true;
      }
      return false;
    }
//    void PopFirst() { } // ????? needed ? or deleted by other means
     // currently BFS using DefaultNodeCmp::operator<
     // Need pure DFS??
//    INode GetEnd() { return nodes.end(); } no iteration in the list
    INode AddRoot(const Node &);
    INode AddOpen(const Node &);
//    void Del(const INode );  deletion by pruning only!
    void ReorderOpen(); // after changing the bounds
    void MarkFathomed(INode ); // problem: have only Node*
    void MarkFathomed(Node* pn) {
      MarkFathomed(Ptr2Iter<INode>(pn));
    }
    void DeleteFathomed() { flist.clear(); } // Del all fathomed

     // ordering at 1st acc. to rounded lpb!
  class DefaultNodeCmp {
  public:
    bool operator() (const node_iterator& m,
         const node_iterator& n) const {
      if (m->llb < n->llb) return 1;
      if (m->llb > n->llb) return 0;
      if (m->depth < n->depth) return 1;
      if (m->depth > n->depth) return 0;
      if (m->nLeft < n->nLeft) return 1;
      if (m->nLeft > n->nLeft) return 0;
 /// AND NOW comes the LP value, though it may be non-exact
 /// under tailing-off control:
      if (m->llv < n->llv) return 1;
      if (m->llv > n->llv) return 0;
 //     assert(m->no != n->no);
      return (m->no > n->no); // last node first.
       // this should provide uniqueness
    }
  };

    typedef set<node_iterator, DefaultNodeCmp> SearchTree;
    typedef SearchTree::iterator ISrch;
    SearchTree::size_type GetSearchSize()
      { return tree.size(); }
    ISrch GetBeginSearch() { return tree.begin(); } // in the s. tree
    void Inc(ISrch & is) {  ++ is; } // by reference _____
    ISrch GetEndSearch() { return tree.end(); }
    void DeleteSearch(const INode in)
    { tree.erase(in); }

    void clear() {
			nodes.clear();
			flist.clear();
			tree.clear();
		}
private:
  NodeContainer nodes,
    flist; // those marked fathomed
  SearchTree tree;
};

inline SolutionTree::INode
  SolutionTree::AddRoot(const Node & n) {
    assert(!n.no);
    nodes.push_back(n); // not to the tree
    n.CheckNew();
    return (--nodes.end());
}

inline SolutionTree::INode
  SolutionTree::AddOpen(const Node & n) {
    assert(n.no);
    nodes.push_back(n);
    n.CheckNew(); // assert pointers are 0
    INode in = nodes.end();
    -- in;
    pair<SearchTree::iterator, bool> res;
    res = tree.insert(in);
    assert(res.second);
    return in;
}
//inline void SolutionTree::Del(const SolutionTree::INode i) {
   // tree.erase(i);
//}
inline void SolutionTree::ReorderOpen() {
    SearchTree tr1(tree.begin(), tree.end());
    tree = tr1;
}
inline void SolutionTree::MarkFathomed(INode itF) {
  itF->state = Node::fathomed;
  flist.splice(flist.end(), nodes, itF);
//  tree.erase(itF); processed => extracted
}

}                                                // SS_END_NAMESPACE__

#endif // __NODE_H__32
