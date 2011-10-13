#ifndef TIL_MISC_SORT_H_
#define TIL_MISC_SORT_H_

/// \file Small helper functions related to sorting.

// includes from STL
#include <utility>  // std::pair

namespace til
{
  //---------------------------------------------------------------------------
  
  /// returns vector (x, y) with x <= y
  template < typename TVector2D, typename T >
  inline TVector2D 
  sortedVector(T a, T b);

  //---------------------------------------------------------------------------  

  template < typename TVector3D, typename T >
  inline TVector3D 
  sorted_vector(T a, T b, T c);

  //---------------------------------------------------------------------------  

  /// NB: the argument is passed by value *on purpose*.
  template < typename TVector3D >
  inline TVector3D 
  sorted_vector(TVector3D v);
  
  //---------------------------------------------------------------------------

  /// Return (a, b) if a >= b, (b, a) else.  
  template < typename T >
  inline std::pair<T,T> 
  make_sorted_pair(T a, T b);

  //---------------------------------------------------------------------------
  
} // namespace til

// package include
#include "misc_sort.tpp"

#endif /*MISC_SORT_H_*/
