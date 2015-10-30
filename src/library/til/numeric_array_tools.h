#ifndef TIL_NUMERIC_ARRAY_TOOLS_H
#define TIL_NUMERIC_ARRAY_TOOLS_H

/// \file Belongs to numeric array package.
/// Do not include directly, include til/numeric_array instead

namespace til
{
  
  //---------------------------------------------------------------------------
  
  // Technically, it's OK if newSize if simply <= D. But the problem is that
  // if we allow newSize < D, we have the surprising effect that doing a resize
  // won't affect a.size()...
  //TODO: I think this crap is completely obsolete, right? I hope so.
  template < typename T, std::size_t D >
  inline void resize(numeric_array<T,D> &, std::size_t
#ifndef NDEBUG
  newSize
#endif
  ) __attribute__((__deprecated__))
  {
    assert(newSize == D);
  }

  //---------------------------------------------------------------------------

  // TODO: keep this in a perf namespace and compare with current min version.
  /*
  template < typename T, std::size_t D >
  inline T min(const numeric_array<T,D> & x)
  { return min(x.whole_range()); }
  */

  template < typename TCollection >
  inline
  std::ostream & stream_collection(std::ostream & os, const TCollection & v)
  {
    os << "[ ";
    typename TCollection::const_iterator iV = v.begin();
    for (; iV != v.end(); ++iV)
    {
      os << *iV << ", ";
    }
    os << " ]";
    return os;
  }

  //---------------------------------------------------------------------------

  template < typename T, std::size_t D >
  inline std::ostream & operator<<(std::ostream & os, const numeric_array<T,D> & x)
  { return stream_collection(os,x); }

  /// Return the squared euclidean norm of a.
  // TODO: How to adapt this when storage is sparse?
  template < typename TPrec, typename T, std::size_t D >
  inline TPrec norm2(const numeric_array<T,D> & a, prec<TPrec>)
  {
    typedef TPrec prec_type;
    prec_type res(0);
    for (std::size_t i = 0; i < D; ++i)
    {
      res += square(static_cast<prec_type>(a[i]));
    }
    return res;
  }
  
  //---------------------------------------------------------------------------

  /// Return the squared euclidean norm of array.
  template < typename T, std::size_t D >
  inline T norm2(const numeric_array<T,D> & a)
  { return norm2(a, prec<T>()); }

  //---------------------------------------------------------------------------

  /// Return the euclidean norm of array.
  template < typename TPrec, typename T, std::size_t D >
  inline TPrec norm(const numeric_array<T,D> & a, prec<TPrec>)
  { return std::sqrt(norm2(a, prec<TPrec>())); }

  //---------------------------------------------------------------------------

  /// Return the squared euclidean norm of array.
  template < typename T, std::size_t D >
  inline T norm(const numeric_array<T,D> & a)
  { return norm(a, prec<T>()); }

  //---------------------------------------------------------------------------

  /// Return the squared Euclidean distance between two vectors, computed with
  /// a precision given as the first template parameter.
  template < typename TPrec, typename T1, typename T2, std::size_t D >
  inline 
  TPrec 
  dist2(const numeric_array<T1,D> & v1, const numeric_array<T2,D> & v2, prec<TPrec>)
  {
    typedef TPrec prec_type;
    prec_type res(0);
    // TODO: replace this by a template expression when have time  
    for (std::size_t i = 0; i < D; ++i)
    {
      // NB: I prefer to directly convert into TPrec in case
      // the content of vectors is unsigned.
      res += square(static_cast<prec_type>(v1[i]) - static_cast<prec_type>(v2[i]));
    }
    return res;
  }

  //---------------------------------------------------------------------------

  /// Return the squared Euclidean distance between two arrays.
  template < typename T, std::size_t D >
  inline 
  T 
  dist2(const numeric_array<T,D> & v1, const numeric_array<T,D> & v2)
  { return dist2(v1, v2, prec<T>()); }

  //---------------------------------------------------------------------------

  /// Return the Euclidean distance between two arrays.
  template < typename TPrec, typename T1, typename T2, std::size_t D >
  inline 
  TPrec 
  dist(const numeric_array<T1,D> & v1, const numeric_array<T2,D> & v2, prec<TPrec>)
  { return std::sqrt(dist2(v1, v2, prec<TPrec>())); }

  //---------------------------------------------------------------------------

  /// Return the Euclidean distance between two arrays.
  template < typename T, std::size_t D >
  inline T 
  dist(const numeric_array<T,D> & v1, const numeric_array<T,D> & v2)
  { return dist(v1, v2, prec<T>()); }
  
  //---------------------------------------------------------------------------

  /// Return the dot product of two vectors.
  template < typename TPrec, typename T1, typename T2, std::size_t D >
  inline TPrec 
  dot(const numeric_array<T1,D> & a1, const numeric_array<T2,D> & a2, prec<TPrec>)
  {
    typedef TPrec prec_type;
    prec_type res = 0;
    for (std::size_t i = 0; i < D; ++i) res += static_cast<prec_type>(a1[i]) * static_cast<prec_type>(a2[i]);
    return res;
  }

  //---------------------------------------------------------------------------

  /// Return the dot product of two vectors.
  template < typename T, std::size_t D >
  inline T dot(const numeric_array<T,D> & a1, const numeric_array<T,D> & a2)
  { return dot(a1, a2, prec<T>()); }
 
  //---------------------------------------------------------------------------

  namespace
  {
    template < typename T >
    inline T 
    cross_line(T a1, T a2, T b1, T b2) { return a1 * b2 - a2 * b1; }    
  }

  //---------------------------------------------------------------------------

  /// Return the cross product of two 3D vectors.
  template < typename TPrec, typename T1, typename T2 >
  inline numeric_array<TPrec,3> 
  cross(const numeric_array<T1,3> & vec1, const numeric_array<T2,3> & vec2, prec<TPrec>)
  {
    // TODO: each element is cast twice to TPrec -- this could be improved.
    return numeric_array<TPrec,3>
      (
      cross_line<TPrec>(vec1[1], vec1[2], vec2[1], vec2[2]),
      cross_line<TPrec>(vec1[2], vec1[0], vec2[2], vec2[0]),
      cross_line<TPrec>(vec1[0], vec1[1], vec2[0], vec2[1])
      );
  }

  //---------------------------------------------------------------------------

  /// Return the norm of the cross product of two 3D vectors.
  template < typename T1, typename T2 >
  inline 
  double
  cross_norm(const numeric_array<T1,3> & vec1, const numeric_array<T2,3> & vec2)
  {
    return std::sqrt(
      til::square(cross_line<double>(vec1[1], vec1[2], vec2[1], vec2[2])) +
      til::square(cross_line<double>(vec1[2], vec1[0], vec2[2], vec2[0])) +
      til::square(cross_line<double>(vec1[0], vec1[1], vec2[0], vec2[1]))
      );
  }

  //---------------------------------------------------------------------------

  /// Absolute value, element-wise.
  template < typename T, std::size_t D >
  inline numeric_array<T,D> 
  abs(const numeric_array<T,D> & a)
  {
    numeric_array<T,D> res;
    for (std::size_t i = 0; i < D; ++i) res[i] = std::abs(a[i]);
    return res;
  }

  //---------------------------------------------------------------------------

  /// Return the maximum value.
  template < typename T, std::size_t D >
  inline T max(const numeric_array<T,D> & a)
  {
    T res = a[0];
    for (std::size_t i = 1; i < D; ++i) max_helper(res, a[i]);
    return res;
  }

  //---------------------------------------------------------------------------

  /// Return the minimum value.
  template < typename T, std::size_t D >
  inline T min(const numeric_array<T,D> & a)
  {
    T res = a[0];
    for (std::size_t i = 1; i < D; ++i) min_helper(res, a[i]);
    return res;
  }

  //---------------------------------------------------------------------------

  /// Return the cross product of two 3D ararys.
  template < typename T >
  inline numeric_array<T,3> cross(const numeric_array<T,3> & vec1, const numeric_array<T,3> & vec2)
  { return cross(vec1, vec2, prec<T>()); }
 
  //---------------------------------------------------------------------------

  /// Return the infinity norm (i.e. max absolute value) of array.
  template < typename T, std::size_t D >
  inline T norm_inf(const numeric_array<T,D> & a)
  {
    T res = std::abs(a[0]);
    for (std::size_t i = 1; i < D; ++i) max_helper(res, std::abs(a[i]));
    return res;
  }
  
  //---------------------------------------------------------------------------

  /// Negate an array.
  template < typename T, std::size_t D >
  inline void neg(numeric_array<T,D> & vec)
  { for (std::size_t i = 0; i < D; ++i) vec[i] = -vec[i]; }
  
  //---------------------------------------------------------------------------

  /// Change precision.
  template < typename T, std::size_t D, typename TNewPrecision>
  struct change_precision< numeric_array<T,D> , TNewPrecision>
  {
    typedef numeric_array<typename change_precision<T, TNewPrecision>::type, D> type;
  private: // requirements
    typedef require_that<std::numeric_limits<TNewPrecision>::is_specialized> Require_is_numeric;
  };
  
  //---------------------------------------------------------------------------

  /// Check that numeric_array does not contain any NAN
  /// TODO: Actually, is_nan should be a functor, and we should use a check loop.
  template < typename T, std::size_t D >
  inline bool is_nan(const numeric_array<T,D> & a)
  {
    for (std::size_t i = 0; i < D; ++i) if (is_nan(a[i])) return true;
    return false;
  }
  
  //---------------------------------------------------------------------------

  template < typename T >
  inline T cross(const numeric_array<T,2> & x, const numeric_array<T,2> & y)
  { return x[0] * y[1] - x[1] * y[0]; }
 
 
  //---------------------------------------------------------------------------

  // This is just to speed up things. Not sure whether this is really necessary...
  // Actually, this does not have to know numeric_array explicitely, does it?
  template < int x, int y, int z, typename TIN, typename TOUT >
  inline void addTo
  (
   const numeric_array<TIN,3> & vIn,    ///< The input vector  
   numeric_array<TOUT,3> & vOut         ///< The output vector
  )
  {
    if (x) vOut[0] = vIn[0] + x;
    if (y) vOut[1] = vIn[1] + y;
    if (z) vOut[2] = vIn[2] + z;
  
  }
  
  //---------------------------------------------------------------------------

  /// Stores v.v^T in mat
  template < typename T >
  inline void
  tdot
  (
    numeric_array<T,3> const & v
  , SymMatrix3<T> & mat
  )
  {
    for (int j = 0; j < 3; ++j)
    {
      T vj = v[j];
      for (int i = j; i < 3; ++i)
      {
        mat(i,j) = v[i] * vj;
      }
    }
  }

  //---------------------------------------------------------------------------

  /// Stores v1.v2^T in mat
  template < typename T >
  inline void
  tdot
  (
    numeric_array<T,3> const & v1
  , numeric_array<T,3> const & v2
  , Matrix3<T> & mat
  )
  {
    for (int j = 0; j < 3; ++j)
    {
      T v1j = v1[j];
      for (int i = 0; i < 3; ++i)
      {
        mat(i,j) = v2[i] * v1j;
      }
    }
  }
  //---------------------------------------------------------------------------

} // namespace til


#endif


