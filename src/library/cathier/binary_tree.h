#ifndef _BINARYTREE_H_
#define _BINARYTREE_H_

// includes from STL
#include <list>
#include <queue>
#include <set>

// includes from BOOST
#include <boost/array.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

// includes from TIL
#include "globalTraits.h"

namespace til
{

  //---------------------------------------------------------------------------
  
    //------------//
   //  NaryTree  //
  //------------//
  
  /// A class to represent an N-ary tree.
  /// An N-ary tree is a tree with N numbered branches at each node. Of course, 
  /// not all of these branches need to be there. But the fact that these branches
  /// are numbered (till N) has a meaning: a node having only one branch is only
  /// part of the story; whether this branch is number 0 or number 7 is also
  /// an information stored in the tree.
  /// From that point of view, a tree is not a graph -- or then, it is a graph
  /// with valued edges.
  /// The type of the value stored at each node is of type T. This could be
  /// anything -- including pointers on actual values.
  /// A tree is a container, but it is also much more than a container. Therefore
  /// the iterators can be used just like STL iterators, but have also other
  /// tree-specific functionalities.
  template < typename T, std::size_t N >
  class NaryTree
  {
  public: // exceptions -------------------------------------------------------
  
    class TreeRootAlreadyDefined : public std::exception {};
    
  public: // classes ----------------------------------------------------------
  
    class iterator;
    class const_iterator;
  
  private: // classes ---------------------------------------------------------
  
    template < typename X > class iterator_base;
    class Node;
    struct Destructor { void operator()(iterator i) const { delete &*i; } };
    
  public: // constructors & destructor ----------------------------------------
  
    /// Default constructor builds an emtpy tree.
    NaryTree() : m_root(), m_size() {}
    
    /// Destructor
    ~NaryTree();
  
  public: // set & get --------------------------------------------------------
  
    std::size_t size() const { return m_size; }
    
  public: // functions --------------------------------------------------------
  
    /// Return a const iterator on the root element of the tree.
    const_iterator
    root() const
    { return const_iterator(m_root); }
  
    /// Return a non_const iterator on the root element of the tree.
    iterator 
    root()
    { return iterator(m_root, 0, 0); }
    
    /// Add a child to a node.
    iterator 
    addChild
    (
      Node * parent
    , std::size_t childNumber
    , T value = T()
    );
  
    /// Add a child to a node.
    iterator
    addChild
    (
      iterator & parent
    , std::size_t childNumber
    , T value = T()
    )
    { return this->addChild(parent.node(), childNumber, value); }
       
  private: // data ------------------------------------------------------------
    
    Node * m_root;  
    std::size_t m_size;
  };
  

  //---------------------------------------------------------------------------

    //------------------//
   //  NaryTree::Node  //
  //------------------//
  
  /// Internal class used inside binary tree.
  /// This class is a private class, theoretically you should not have to
  /// know anything about it to use trees, the same way you probably haven't heard
  /// of _List_node even though you have been using a GNU implementation 
  /// of std::list.
  template < typename T, std::size_t N >
  class NaryTree<T,N>::Node
  {
  public: // typedefs ---------------------------------------------------------
  
    typedef T value_type; 
   
  public: // constructors -----------------------------------------------------
  
    /// Default constructor, set everything to zero.
    Node() :
        m_value()
      , m_branch()
    {
      m_branch.assign(0);
    }
  
    /// Constructor with a single value: collect value and set branches to zero.
    explicit Node(const T & value) :
        m_value(value)
      , m_branch()
    {
      m_branch.assign(0);
    }
    
  public: // set & get --------------------------------------------------------
  
    const T & operator()() const { return m_value; }
    T & operator()() { return m_value; }
  
    const T & get() const { return m_value; }
    T & get() { return m_value; }
  
    /// Returns a pointer on the n-th child.
    Node * child(std::size_t i) const
    {
      assert(i < til::size(m_branch));
      return m_branch[i];
    }
  
    /// Returns a pointer on the n-th child.
    Node * & child(std::size_t i)
    {
      assert(i < til::size(m_branch));
      return m_branch[i];
    }
    
    void operator=(const T & value) { m_value = value; }
  
  private: // data ------------------------------------------------------------
  
    // Index to the current value
    T m_value;
    // pointers to branches
    boost::array<Node*, N> m_branch;  
  };


  //---------------------------------------------------------------------------

    //---------------------------//
   //  NaryTree::iterator_base  //
  //---------------------------//
  
  /// A base for NaryTree iterators, that collects common code.
  /// Useless on its own.
  template < typename T, std::size_t N >
  // NB: It is understood that X is a Node
  template < typename X >
  class NaryTree<T,N>::iterator_base
  {
  public: // typedefs
    typedef iterator_base<X>                Self;
    typedef typename X::value_type          value_type;
    
  public: // constructors & destructor ----------------------------------------
  
    iterator_base() {}
    explicit iterator_base(X * pNode) : m_pNode(pNode) {}
  
  public: // functions --------------------------------------------------------
  
    /// Accessing values
    const value_type & operator*() const { return m_pNode->get(); }
    const value_type * operator->() const { return &(m_pNode->get()); }
  
    /// Accessing node
    operator X*() { return m_pNode; }
    X * node() const { return m_pNode; }
  
  public: // static functions -------------------------------------------------
  
    static std::size_t nchildren() { return N; }
  
  protected: // functions -----------------------------------------------------
  
    X * child(std::size_t i) const { return this->node()->child(i); }
  
  protected: // data ----------------------------------------------------------
  
    X * m_pNode;
  };
  
  //---------------------------------------------------------------------------
  
    //----------------------------//
   //  NaryTree::const_iterator  //
  //----------------------------//
  
  /// Const iterator on n-ary trees.
  /// On top of traditional, STL iterator features, the class offers the
  /// possibility to go to one of its children via the child() member function.
  template < typename T, std::size_t N >
  //class NaryTree<T,N>::const_iterator : public iterator_base<const Node>
  class NaryTree<T,N>::const_iterator : public NaryTree<T,N>::template iterator_base<const Node>
  {
  public: // typedefs ---------------------------------------------------------
  
    typedef const_iterator                  Self;
    typedef iterator_base<const Node>       Base;
    typedef typename Base::value_type       value_type;
  
  public: // constructors -----------------------------------------------------
  
    // NB: default constructor set const_iterator to zero
    const_iterator() : Base() {}
    explicit const_iterator(const Node * p) : Base(p) {}
  
  public: // functions --------------------------------------------------------
  
    /// Get an iterator on the i-th child.
    Self child(std::size_t i) const { return Self(this->Base::child(i)); }
  };
  
  //---------------------------------------------------------------------------
  
    //----------------------//
   //  NaryTree::iterator  //
  //----------------------//
  
  /// Non-const iterator on n-ary trees.
  /// On top of traditional, STL iterator features, the class offers the
  /// possibility to go to one of its children via the child() member function.
  /// It also keeps track of where it is coming from (its parent and its branch
  /// number), because this information is necessary for new allocations.
  // TODO: get rid of code duplication with const_iterator if possible.
  // TODO: the advertising of parent() and branch() is not cool -- this is more of
  // an internal stuff.
  template < typename T, std::size_t N >
  class NaryTree<T,N>::iterator : public NaryTree<T,N>::template iterator_base<Node>
  {
  public: // typedefs ---------------------------------------------------------
  
    typedef iterator_base<Node>             Base;
    typedef typename Base::value_type       value_type;
    
  public: // constructors -----------------------------------------------------
  
    // NB: default constructor set iterator to zero
    iterator() : Base(), m_parent() {}
    iterator(Node * p) : Base(p) {}
    iterator(Node * p, Node * parent, std::size_t branch) : Base(p), m_parent(parent), m_branch(branch) {}
    //iterator(const iterator &i) : TCollection::iterator(i) {}
  
  public: // functions --------------------------------------------------------
  
    /// Convertion into a Node pointer.
    //operator Node* () { return &**this; }
  
    iterator child(std::size_t i) const
    {
      // NB: no need to assert here, as it is done in Node
      return iterator(this->node()->child(i), this->node(), 0);
    }
  
    std::size_t branch() const { return m_branch; }
  
    Node * parent() const { return m_parent; }
  
    value_type & operator*() { return this->m_pNode->get(); }
    value_type * operator->() { return &(this->m_pNode->get()); }
  
  private: // data ------------------------------------------------------------
  
    // The private information of the non-const iterator contains information
    // necessary in case of a new allocation, i.e. its relation with its parent.
  
    // Pointer on the parent  
    Node * m_parent;
    // Branch number on which we lie.
    std::size_t m_branch;
  };
  
  template < typename T, std::size_t N >
  inline
  std::size_t
  size(const NaryTree<T,N> & ntree) { return ntree.size(); }
  

  //---------------------------------------------------------------------------
  
    //------------------------//
   //  functors and helpers  //
  //------------------------//


  
  namespace treeFunctors
  {
    // a dummy functor to test tree traversals.
    template < typename TTreeIterator >
    class Count
    {
    public: // constructors & destructor
      Count() : m_count() {}
      
    public: // functions
      static bool leftFirst() { return true; }
      static bool left() { return true; }
      static bool right() { return true; }
  
      void operator()(const TTreeIterator &) { ++m_count; }
      std::size_t operator()() { return m_count; }
  
    private:
      std::size_t m_count;
    };
  }

  
  //---------------------------------------------------------------------------
  
  /// Pre-order scan.
  /// Pre-order scanning means looking at a node, and then at its children in
  /// order. This has some similarity with a graph depth first search, except that
  /// children scanning order is important in a tree.
  template < typename TTreeIterator, typename TTreeFunctor >
  void pre_order_scan(TTreeIterator iTree, TTreeFunctor & f)
  {
    f(iTree);
    for (std::size_t i=0; i < iTree.nchildren(); ++i)
    {
      if (&*(iTree.child(i))) pre_order_scan(iTree.child(i), f);
    }
  }
  
  //---------------------------------------------------------------------------

  /// Post-order scan.
  /// Post-order scanning means looking at the children in order, and then at the
  /// node. An exemple of such visit order is the inverse-polish notation.
  template < typename TTreeIterator, typename TTreeFunctor >
  void post_order_scan(TTreeIterator iTree, const TTreeFunctor & f)
  {
    for (std::size_t i=0; i < iTree.nchildren(); ++i)
    {
      if (&*(iTree.child(i))) post_order_scan(iTree.child(i), f);
    }
    f(iTree);
  }
  
  //---------------------------------------------------------------------------

  /// In-order scan.
  /// In-order scanning means looking at the left children, then the node, then
  /// the right children. This works with binary trees only.
  template < typename TTreeIterator, typename TTreeFunctor >
  void
  in_order_scan(TTreeIterator iTree, const TTreeFunctor & f)
  {  
    while (iTree)
    {
      if (&*(iTree.child(0))) in_order_scan(iTree->child(0));
      f(iTree);
      iTree = iTree.child(1);
    }
  }
  
  //---------------------------------------------------------------------------

  namespace detail
  {
    template < typename TTreeIterator, typename TBFFunctor >
    void breadth_first(std::queue<TTreeIterator> q, const TBFFunctor & f)
    {
      std::queue<TTreeIterator> newQ;
      while (!q.empty())
      {
        f(q.front());
        for (std::size_t i = 0; i < q.front().nchildren(); ++i)
        {
          if (q.front().child(i)) newQ.push(q.front().child(i));
        }
        q.pop();
      }
      if (!newQ.empty())
      {
        // warn functor that we are about to go to the next level  
        f.nextLevel();
        // iterate
        breadth_first(newQ, f);
      }
    }
  } // namespace detail
  
  //---------------------------------------------------------------------------

  /// Apply a functor in a breadth-first traversal order, starting from the iterator
  /// given in input.
  /// Note that so far, the functor needs to be breadth_first specific.
  template < typename TTreeIterator, typename TBFFunctor >
  void breadth_first(TTreeIterator iTree, const TBFFunctor & f)
  {
    // Do nothing if nothing to be done.
    if (!&*iTree) return;
    // Here, we simply construct a queue containing a single element, and call
    // the appropriate version of breadth_first.
    std::queue<TTreeIterator> q;
    q.push(iTree);
    detail::breadth_first(q, f);
  }
  
  //---------------------------------------------------------------------------

  namespace detail
  {
    /// Breadth-first tree functor to print out tree.
    struct Print
    {
      /// Function called at 
      static void nextLevel() { std::cout << std::endl; }
  
      template < typename TIterator >
      void operator()(const TIterator & i) const
      {
        if (&*i)
        {
          std::cout << *i << " , ";
        }
        else 
        {
          std::cout << " [X] , ";
        }      
      }
    };
  } // namespace detail
  
  //---------------------------------------------------------------------------

  /// Print tree on stdout (very ugly and experimental).
  template < typename T, std::size_t N >
  void print(const NaryTree<T,N> & tree)
  {
    breadth_first(tree.root(), detail::Print());
  }

  //---------------------------------------------------------------------------

} // namespace til

// package include
#include "binary_tree.tpp"

#endif //_BINARYTREE_H_
