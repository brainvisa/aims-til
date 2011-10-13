#ifndef SPARSE_VECTOR_TOOLS_H_
#define SPARSE_VECTOR_TOOLS_H_

/// \file Belong to sparse_vector package.
/// Do not include directly, include "til/sparse_vector.h" instead.

namespace til
{
  /// Check whether there is an element of sparse array which is NaN.
  /// NB: This is to be used only in generic programming situations were we check for NaN on type T,
  /// not knowing whether T is a concrete type or a collection.
  /// If we know the type is a collection, simply use STL functions like find, which has the advantage of
  /// returning the position where the fault occurs.
//   template < typename T >
//   bool is_nan(const sparse_vector<T> & a)
//   {
//     typename sparse_vector<T>::Map::const_iterator i = a.getMap().begin();
//     typename sparse_vector<T>::Map::const_iterator end = a.getMap().end();
//     for (; i != end; ++i)
//     {
//       if (is_nan(i->second)) return true;
//     }
//     return false;
//   }
//   template < typename T >
  inline bool is_nan(const sparse_vector<double> & a)
  {
    for (sparse_vector<double>::Map::const_iterator i = a.getMap().begin(); i != a.getMap().end(); ++i)
    {
      if (is_nan(i->second)) return true;
    }
    return false;
  }
} // namespace til

#endif /*SPARSE_VECTOR_TOOLS_H_*/
