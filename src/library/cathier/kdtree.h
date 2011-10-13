#ifndef _KDTREE_H_
#define _KDTREE_H_

// include from STL
#include <ostream>

// include from BOOST
#include <boost/shared_ptr.hpp>
#include <boost/type_traits.hpp>

// includes from TIL
#include "til/numeric_array.h"

// includes from TIL
#include "binary_tree.h"
#include "index_collection.h"

namespace til
{
  //---------------------------------------------------------------------------

  // technical details to initialize result array to -1 only if index type
  // is integer
  template < typename _TIndex >
  inline _TIndex kdtreeDefaultValue() { return _TIndex(); }
  template < >
  inline std::size_t kdtreeDefaultValue<std::size_t>()
  { return std::size_t(-1); }
  
  //---------------------------------------------------------------------------

  namespace
  {
    /// A type used to build KDTree.
    /// It has no value for the user; the only reason it is not a private class
    /// of KDTree is that it can't because KDTree derives from NaryTree<KDNode>
    /// and it seems incompatible with KDNode being declared inside KDTree.
    template < typename T >
    class KDNode
    {
    public: // set & get ------------------------------------------------------
    
      int   dim() const { return m_dim; }
      int & dim()       { return m_dim; }
      
      T const & value() const { return m_value; }
      T       & value()       { return m_value; }
      
      const T & operator()() { return m_value; }
    
      operator T () const { return m_value; }
    
    private: // data ----------------------------------------------------------
    
      T m_value;
      // The dimension along which the split is done
      int m_dim;
    };
  }  
  
  //---------------------------------------------------------------------------
  
  template < typename T >
  inline std::ostream & operator<<(std::ostream & os, const KDNode<T> & n)
  {
    return os << n.value();
  }
  
  //---------------------------------------------------------------------------

    //----------//
   //  KDTree  //
  //----------//

  template < typename TIndex, typename TContainer >
  class KDTree
    : public NaryTree< KDNode< TIndex >, 2 >
    , public index_collection< const TContainer, TIndex >
  {
  public: // typedefs ---------------------------------------------------------
  
    typedef KDTree<TIndex, TContainer>  Self;
    // For safety reasons, because of multiple inheritance.
    typedef void                        Base;
    typedef TContainer Collection;
    
  public: // enums ------------------------------------------------------------
  
    enum { LEFT = 0, RIGHT };  
  
  public: // construtors ------------------------------------------------------
  
    KDTree() : 
        NaryTree<KDNode<TIndex>,2>()
      , index_collection<const TContainer, TIndex>()
    {}
    
    // templated to avoid compilation when used
    template < typename T_Container >
    KDTree(const T_Container & c)
      : NaryTree<KDNode<TIndex>,2>()
      , index_collection<const TContainer, TIndex>(c)
    {}
  };
    
  //---------------------------------------------------------------------------

    //-----------------//
   //  KDTreeFactory  //
  //-----------------//
    
  template < typename TKDTree >
  class KDTreeFactory
  {
  public: // typedefs ---------------------------------------------------------
    
    typedef KDTreeFactory<TKDTree>             Self;
    typedef typename TKDTree::index_type       index_type;
    typedef typename TKDTree::indexed_type     indexed_type;
    typedef typename TKDTree::Collection       Collection;
    // Of course the real deal is to have KDTreeIndex as a template parameter
    // of this class
  
  public: // constructors -----------------------------------------------------
  
    /// Initialize the factory with the collection of which we want to compute
    /// the KDTree.
    KDTreeFactory(const Collection & v) : m_v(v), m_compare(*this) {}
    
  public: // functions --------------------------------------------------------
  
    void build(TKDTree & res);

  private: // classes ---------------------------------------------------------
  
    /// Compares the n-th coordinate of two vectors    
    class Coordinate_comparison
    {
    public:
      Coordinate_comparison(const Self & parent) : m_parent(parent) {}
    
      /// When indices are pointers
      template < class TVector >
      inline bool operator()(const TVector *v1, const TVector *v2)
      { return (*v1)[m_parent.m_dim] < (*v2)[m_parent.m_dim]; }
  
      /// When indices are integers
      inline bool operator()(std::size_t i1, std::size_t i2)
      { return m_parent.m_v[i1][m_parent.m_dim] < m_parent.m_v[i2][m_parent.m_dim]; 
      }
    private:
      const Self & m_parent;
    };
  
  private: // functions -------------------------------------------------------
  
    /// An iteration for building KDTree.
    /// This function calls itself recursively to build the KDTree
    /// The main idea is to find the median of the data, and to recursively
    /// find the median for the elements below and above the median.
    /// More precisely, a call to this function will store in res[index] the
    /// median of the elements in the range [left, right] (yes, right is included
    /// in the range, otherwise it's getting more of a headache), sorted along
    /// dimension 'dim' -- and will continue sorting for lower and higher sets.
    void
    genKDIter
    (
     typename std::vector<index_type>::iterator left,   ///< left-most element of the range we want to sort
     typename std::vector<index_type>::iterator right,  ///< right-most element of the range we want to sort
     typename TKDTree::iterator iRes,                   ///< the kdtree being built
     int dim                                            ///< the dimension along which the sorting is done to find the median
    );
      
    /// If the index is of type T*, take the reference of the element
    inline
    void getIndex(const typename Collection::const_iterator & iV, indexed_type * & vIndex)
    { vIndex = const_cast<indexed_type*>(&*iV); }
  
    /// If the index is of type std::size_t, take the distance to the first
    /// element
    inline
    void getIndex(const typename Collection::const_iterator & iV, std::size_t & vIndex)
    { vIndex = std::distance(m_v.begin(), iV); }
        
    /// Return a vector containing the indices of the collection elements.
    /// E.g, if TIndex is an integer, the indices will be (0, 1, 2, 3...).
    /// If TIndex is T*, the indices will be (&v[0], &v[1], ...).
    void makeIndexVector(const Collection & v, std::vector<index_type> & vIndex);
    
  private: // variables -------------------------------------------------------
  
    /// A reference on the collection whose KDTree is to be computed
    const Collection & m_v;
    Coordinate_comparison m_compare;
    int m_dim;
    std::vector<index_type> m_vIndex;
    //KDTree<TIndex> m_kdt;
    TKDTree *m_pKdt;
  };
  
  //-------------------------------------------------------------------------------------
  
    //----------------//
   //  Find_closest  //
  //----------------//
  
  /// A class to find a KDTree point closest to an input point.
  // TODO: This should be renamed KDTFinder, or put inside a namespace or something.
  template < typename TPrecision, typename TKDTree >
  class Find_closest
  {
  public: // typedefs ---------------------------------------------------------
  
    typedef typename TKDTree::index_type index_type;
    typedef typename TKDTree::indexed_type indexed_type;
    
  public: // constructors -----------------------------------------------------
    
    /// Provide input kdtree.
    // TODO: remove reference. People can still put the reference in the template parameters if needed, I think.
    Find_closest(const TKDTree & kdtree) : m_pKdt(&kdtree), m_niter(0) {}
    
  public: // functions --------------------------------------------------------
  
    /// Get kdtree point closest to input point.
    index_type operator()(const indexed_type & vec);
    /// Return the number of iterations taken by the last search
    int niter() { return m_niter; }
    
  private: // functions -------------------------------------------------------
  
    void init();
    void lookMax(typename TKDTree::const_iterator iNode);
  
  private: // data ------------------------------------------------------------
  
    const TKDTree * m_pKdt;
    indexed_type m_vec;
    TPrecision m_min;
    index_type m_minIndex;
    int m_niter;
  };

  //---------------------------------------------------------------------------
  
  template < typename T, typename TIndex, typename TContainer >
  inline void makeKDTree(const std::vector<T> & v,
                         KDTree<TIndex, TContainer> & res)
  {
    KDTreeFactory<KDTree<TIndex, TContainer> > kf(v);
    kf.build(res);
  }
  
  //---------------------------------------------------------------------------

  /// Returns the closest point to vec in kdtree
  template < typename TPrecision, typename TKDTree >
  inline
  typename TKDTree::index_type
  find_closest(const TKDTree & kdtree, const typename TKDTree::indexed_type & vec)
  {
    Find_closest<TPrecision, TKDTree> fc(kdtree);
    return fc.doit(vec);
  }
    
} // namespace til

#include "kdtree.tpp"

#endif //_KDTREE_H_

