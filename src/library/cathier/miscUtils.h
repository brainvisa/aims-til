#ifndef _MISCUTILS_H_
#define _MISCUTILS_H_

// includes from STL
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdarg> // va_start
#include <numeric> // accumulate

// includes from BOOST
#include <boost/array.hpp>
#include <boost/type_traits.hpp>

// includes from AIMS
#include "aims/vector/vector.h"

// includes from TIL
#include "til/is_traits.h"
#include "til/TExpr.h"
#include "til/TExprOperators.h"

// includes from TIL
#include "globalTraits.h"

namespace til
{
  namespace policy
  {
    //---------------------------------------------------------------
    
    /// A policy label to indicate that radian range is ]-PI; PI]
    class Radian_CenteredRange {};
    
    //---------------------------------------------------------------

    /// A policy label to indicate that radian range is [0; PI[
    class Radian_PositiveRange {};

    //---------------------------------------------------------------

  } // namespace policy
  
  

  
  //---------------------------------------------------------------------------  
  /*  
  /// Compute the average of a list of points
  template < typename TVertexCollection >
  Point3dd
  getAverage(const TVertexCollection &vertices)
  {
    typename TVertexCollection::const_iterator iVertex;
    Point3dd res(0.0, 0.0, 0.0);
  
    // Add position of all vertices
    for (iVertex = vertices.begin(); iVertex != vertices.end(); ++iVertex)
    {
      res += *iVertex;
    }
    // Divide by the number of points
    res *= 1.0/size(vertices);
    return res;
  }
  */  
  
  //---------------------------------------------------------------------------

  /// Squared euclidean norm of vector (x,y,z)
  template < typename T >
  inline
  T
  norm2(T x, T y, T z)
  { return x*x+y*y+z*z; }
  
  //---------------------------------------------------------------------------

  /*
  template < typename T, typename TPoint3DFrom >
  inline
  typename boost::enable_if_c<is_3DPoint<TPoint3DFrom>::value>::type
  cast(const TPoint3DFrom & x, boost::array<T,3> & y)
  {
    y[0] = p[0];
    y[1] = p[1];
    y[2] = p[2];
  }
  */
  
  /*
  /// User-definable static-cast. By default, uses standard static_cast
  template < typename TTo, typename TFrom >
  typename boost::disable_if_c<
  (is_3DPoint<TFrom>::value && is_3DPoint<TTo>::value) ||
  (is_3DPoint<TFrom>::value && is_BoostArray_N<TTo,3>::value)
  , TTo>::type
  cast(const TFrom & v) { return static_cast<TTo>(v); }
  
  /// Standard convertion between 3D points
  template < typename TPoint3DTo, typename TPoint3DFrom >
  typename boost::enable_if_c<is_3DPoint<TPoint3DFrom>::value && is_3DPoint<TPoint3DTo>::value, TPoint3DTo && !is_BoostArray_N<TPoint3DTo>::value>::type
  cast(const TPoint3DFrom &p)
  {
    return TPoint3DTo(p[0], p[1], p[2]);
  }
  
  template < typename TBoostArray, typename TPoint3DFrom >
  typename boost::enable_if_c<is_3DPoint<TPoint3DFrom>::value && is_BoostArray_N<TBoostArray, 3>::value, TBoostArray>::type
  cast(const TPoint3DFrom & p)
  {
    TBoostArray tmp;
    tmp[0] = p[0];
    tmp[1] = p[1];
    tmp[2] = p[2];
    return tmp;
  }
  */
  
  /*
  namespace til
  {
    /// Mostly similar to std::copy, except that convert is used in the
    // copy process
    template < typename TIterator1, typename TIterator2 >
    void
    copy_convert
    (
     const TIterator1 & begin1,
     const TIterator1 & end1,
     const TIterator2 & begin2
    )
    {
      TIterator1 i1 = begin1;
      TIterator2 i2 = begin2;
      for (; i1 != end1; ++i1, ++i2)
      {
        convert(*i1, *i2);
      }
    }
  }
  */

    
  namespace detail
  {
    template < typename TContainer1, typename TContainer2 >
    typename boost::disable_if_c<
      is_container<typename TContainer1::value_type>::value &&
      is_container<typename TContainer2::value_type>::value
    >::type
    _allocate_sameSize
    (
     const TContainer1 & c1,    ///< The model
     TContainer2 & c2           ///< The container to be allocated
    )
    {
      if (size(c2) != size(c1))
      {
        c2.resize(size(c1));
      }
    }
  
    template < typename TContainer1, typename TContainer2 >
    typename boost::enable_if_c<
      is_container<typename TContainer1::value_type>::value &&
      is_container<typename TContainer2::value_type>::value
    >::type
    _allocate_sameSize
    (
     const TContainer1 & c1,    ///< The model
     TContainer2 & c2           ///< The container to be allocated
    )
    {
      if (size(c2) != size(c1))
      {
        c2.resize(size(c1));
      }
      typename TContainer1::const_iterator iC1 = c1.begin();
      typename TContainer2::iterator iC2 = c2.begin();
      for (; iC1 != c1.end(); ++iC1, ++iC2)
      {
        _allocate_sameSize<typename TContainer1::value_type, typename TContainer2::value_type>(*iC1, *iC2);
      }
    }
  }
  
  
  /// A simple function to allocate a container so that its size matches the size
  /// of a second container.
  /// The trick is that it is recursive, so that we can easily create a container
  /// of container (of container...) whose size matches the size of another
  /// container.
  /// NB: You have to be careful that this won't go farther than what you want.
  /// In particular if you are dealing with arrays of vectors of matrices, that
  /// might themselves be considered as containers.
  template < typename TContainer1, typename TContainer2 >
  void
  allocate_sameSize
  (
   const TContainer1 & c1,    ///< The model
   TContainer2 & c2           ///< The container to be allocated
  )
  { detail::_allocate_sameSize<TContainer1, TContainer2>(c1,c2); }
  
  /*
  template < typename TContainer1, typename TContainer2 >
  void
  allocate_sameSize
  (
   const TContainer1 & c1,    ///< The model
   TContainer2 & c2           ///< The container to be allocated
  )
  {
    if (size(c2) != size(c1))
    {
      c2.resize(size(c1));
    }
    if (is_container<typename TContainer1::value_type>::value &&
        is_container<typename TContainer2::value_type>::value)
    {
      typename TContainer1::const_iterator iC1 = c1.begin();
      typename TContainer2::iterator iC2 = c2.begin();
      for (; iC1 != c1.end(); ++iC1, ++iC2)
      {
        allocate_sameSize(*iC1, *iC2);
      }
    }
  }
  */
  
  
  namespace functor
  {
    
    /// Computes the squared Euclidean distance between two vectors.
    // TODO: I think boost::numeric::converter could be used there but maybe
    // without range checking -- this could really drag performance down.
    template < typename TPrecision, typename TVector3D  >
    inline 
    TPrecision
    diff2(const TVector3D &v1, const TVector3D &v2)
    { return static_cast<TPrecision>(norm2(v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2])); }
  
    /// Computes the Euclidean distance between two vectors.
    template < typename TPrecision, typename TVector3D  >
    inline 
    TPrecision
    diff(const TVector3D &v1, const TVector3D &v2)
    { return static_cast<TPrecision>(std::sqrt(diff2<double>(v1,v2))); }
      
    /// Computes the squared Euclidean distance between two vectors.
    struct Diff2
    {
      template < typename TVector3D, typename TPrecision >
      typename boost::enable_if_c<
        is_3DVector<TVector3D>::value
      >::type
      operator()(const TVector3D &v1, const TVector3D &v2, TPrecision &res) const
      { res = diff2<TPrecision>(v1,v2); }
    };
  
    /// Computes the Euclidean distance between two vectors.
    /// TODO: remove this
    struct Diff
    {
      template < typename TVector3D, typename TPrecision >
      void operator()(const TVector3D &v1, const TVector3D &v2, TPrecision &res) const
      { res = diff<TPrecision>(v1,v2); }
    };
  
    /// Normalize a vector with its Euclidean norm.
    template < typename TVector3D >
    struct Normalize
    {
      void operator()(const TVector3D &v1, const TVector3D &v2, TVector3D & res) const
      {
        sub(v1, v2, res);
        div(res, norm(res));
      }
    };
  }
 
  //---------------------------------------------------------------------------
  
  /// Computes the centroid of a list of vertices.
  template < typename TRes, typename TIterator >
  void mean(TIterator begin, TIterator end, TRes & res);
  
  //---------------------------------------------------------------------------

  /// Computes the centroid of a list of vertices.
  /// NB: this assumes that the default constructor of TVector inits to zero.
  /// NB: this is obsolete. Use mean instead.
  // TODO: change this: remove res = T()
  template < typename TRes, typename TCollection >
  void centroid(const TCollection & c, TRes & res)
    __attribute__((__deprecated__("use mean() instead")));

  template < typename TRes, typename TCollection >
  void centroid(const TCollection & c, TRes & res)
  {
    mean(c.begin(), c.end(), res);
  }
  
  //---------------------------------------------------------------------------
  
  template < typename TPoint, typename TPointCollection >
  void stdev(const TPointCollection & c, TPoint & res);
  
  //---------------------------------------------------------------------------
  
  /// Return true iff  both arguments have same sign.
  template < typename T1, typename T2 >
  inline bool 
  same_sign(const T1 & x, const T2 & y)
  {
    if (x >= 0) return (y >= 0);
    else        return (y <= 0);
  }
  
  //---------------------------------------------------------------------------

  /// Change the sign of x if y < 0
  template < typename T1, typename T2 >
  inline T1 sign(T1 x, T2 y)
  { return (y >= 0 ? x : -x); }
    
  //---------------------------------------------------------------------------

  template < typename T >
  char * sprintf(T format, ... )
  {
    static char *s = new char[1024];
    va_list vl;
    va_start(vl, format);
    vsprintf(s, format, vl);
    return s;
  }
  
  //---------------------------------------------------------------------------

  template < typename T >
  struct Lexicographical_compare
  {
    bool operator()(const T & x, const T & y) const
    { return std::lexicographical_compare(x.begin(), x.end(), y.begin(), y.end()); }
  };

  //---------------------------------------------------------------------------

  /// Invert an index array.
  /// This means, if the i-th element of the collection contains indices (j,k), then the j-th and k-th
  /// element of the result will contain indice i.
  /// This could correspond to inverting the direction of a directed graph.
  /// This can be useful also to reverse a look-up table, e.g. if a mesh has a collection of faces that
  /// indexes the vertex, then inverting this collection gives for each vertex the index of faces it belongs
  /// to.
  template < typename TIndexCollection >
  std::vector<std::list<std::size_t> >
  invertIndices(const std::vector<TIndexCollection> & c);
  
  //-----------------------------------------------------------------
  
  /// operator< on the address in memory of the pointed object
  template < typename TIterator >
  struct Lesser_PointeeAddress
  {
    bool operator()(TIterator i, TIterator j) const
    { return &*i < &*j; }
  };

  //-----------------------------------------------------------------
  
  /// operator< on the second member of a pair.
  template < typename T1, typename T2 >
  struct Lesser_Pair2
  {
    bool operator()(const std::pair<T1,T2> & p1, const std::pair<T1,T2> & p2) const
    { return p1.second < p2.second; }
  };
  
  //-----------------------------------------------------------------
  
  /// operator> on the second member of a pair.
  template < typename T1, typename T2 >
  struct Greater_Pair2
  {
    bool operator()(const std::pair<T1,T2> & p1, const std::pair<T1,T2> & p2) const
    { return p1.second > p2.second; }
  };
  
  //-----------------------------------------------------------------

  /// Returns first value of pair.
  template < typename T1, typename T2 >
  struct Return_pair1
  {
    inline T1 operator()(std::pair<T1, T2> const & p) const { return p.first; }
  };
    
  //-----------------------------------------------------------------

  /// Returns i0, where X_i0 is the greatesstructt of all Xi's.
  template < typename T >
  inline std::size_t max_index(T x0, T x1, T x2);

  //---------------------------------------------------------------------------

  /// Returns i0, where X_i0 is the smallest of all Xi's.
  template < typename T >
  inline std::size_t min_index(T x0, T x1, T x2);

  //---------------------------------------------------------------------------

  /// Check if x is NaN.
  template < typename T >
  inline bool is_nan(T x);
   
  //---------------------------------------------------------------------------

  template < typename TPrec >
  struct SquaredEuclideanDist
  {
    template < typename T, std::size_t D >
    TPrec operator() (const numeric_array<T,D> & a1, const numeric_array<T,D> & a2) const
    {
      return dist2(a1, a2, prec<TPrec>());
    }
  };
  
  //---------------------------------------------------------------------------

  /// Returns true if argument is of the form 2^m, m>=1.
  // TODO: add a test so that it works only for signed integers... or is >> working only on signed integers anyway?
  template < typename T >
  bool is_dyadic(T n);

  //---------------------------------------------------------------------------

  /// Returns the largest m so that n >= 2^m.
  template < typename T >
  inline int log2(T n);

  //---------------------------------------------------------------------------

  /// Returns 2^n
  template < typename T >
  inline T exp2(unsigned int n);
  
  //---------------------------------------------------------------------------

  /// Return the greater power of two inferior or equal to n.
  template < typename T >
  inline T lower_dyadic(T n);
  
  //-------------------------------------------------------------------------------------------
  
  /// A class which returns a value that is increased everytime it is returned.
  // TODO: shouldn't this guy be in functor or something? Or something like source or generator?
  template < typename T >
  class Incrementor
  {
  public: // constructors
    /// Initialize with value.
    Incrementor(T i0) : m_i(i0) {}
  public: // operator
    /// Return current value.
    /// Value is increased at each call.
    T operator()() { return m_i++; }
  private: // data
    T m_i;
  };
    
  //-------------------------------------------------------------------------------------------

    
} // namespace til

// package include
#include "misc_utils.tpp"


#endif //_MISCUTILS_H_
