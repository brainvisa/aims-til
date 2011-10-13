#ifndef POINT_TOOLS_H_
#define POINT_TOOLS_H_

/// \file Belongs to package Point
/// Do not include directly, include til/Point.h instead

// includes from TIL
#include "til/Accumulator.h"

namespace til
{
   
  template < typename T, std::size_t D, typename TStorage, typename TNewPrecision>
  struct change_precision< Point<T,D,TStorage> , TNewPrecision>
  {
    typedef Point<T,D, typename change_precision<TStorage, TNewPrecision>::type> type;
  private: // requirements
    typedef require_that<std::numeric_limits<TNewPrecision>::is_specialized> Require_is_numeric;
  };
    
  template < typename TAccumulationPoint, typename TPointCollection >
  TAccumulationPoint centroid(const TPointCollection & c)
  {
    MeanAccumulator<typename TPointCollection::value_type, TAccumulationPoint> m;
    m.accumulate(c);
    return m.get();
  }
  
  template < typename T, std::size_t D >
  T dist2(const Point<T,D> & p1, const Point<T,D> & p2)
  {
    return dist2(p1.data(), p2.data());
  }

  /// Return the squared Euclidean distance between two vectors, computed with
  /// a precision given as the first template parameter.
  template < typename TPrecision, typename T, std::size_t D >
  TPrecision
  dist2(const Point<T,D> & p1, const Point<T,D> & p2, prec<TPrecision>)
  {
    return dist2(p1.data(), p2.data(), prec<TPrecision>());
    /*
    TPrecision res(0);
    // TODO: replace this by a template expression when have time  
    typename Point<T,D>::const_iterator iV1 = v1.begin();
    typename Point<T,D>::const_iterator iV2 = v2.begin();
    for (; iV1 != v1.end(); ++iV1, ++iV2)
    {
      // NB: I prefer to directly convert into TPrecision in case
      // the content of vectors is unsigned.
      res += square(TPrecision(*iV1) - TPrecision(*iV2));
    }
    return res;
    */
  }

  /*
  template < typename TPrecision, typename T, std::size_t D >
  TPrecision
  dist2(const Point<T,D> & v1, const Point<T,D> & v2)
  {
    TPrecision res(0);
    // TODO: replace this by a template expression when have time  
    typename Point<T,D>::const_iterator iV1 = v1.begin();
    typename Point<T,D>::const_iterator iV2 = v2.begin();
    for (; iV1 != v1.end(); ++iV1, ++iV2)
    {
      // NB: I prefer to directly convert into TPrecision in case
      // the content of vectors is unsigned.
      res += square(TPrecision(*iV1) - TPrecision(*iV2));
    }
    return res;
  }
*/
  
  /// Return the squared Euclidean distance between two vectors.
  /// Simpy take T as the precision. NB: This can be dangerous say if
  /// T is unsigned.
  // TODO: maybe add this kind of checks in the code?
  /*
  template < typename T, std::size_t D >
  T dist2(const Point<T,D> & v1, const Point<T,D> & v2)
  {
    return dist2<T,T,D>(v1,v2);
  }
  */
  
  template < typename TPrecision, typename T, std::size_t D >
  TPrecision
  dist(const Point<T,D> & v1, const Point<T,D> & v2)
  {
    return TPrecision(std::sqrt(dist2(v1, v2, prec<double>())));
  }
  
      
} // namespace til

#endif /*POINT_TOOLS_H_*/
