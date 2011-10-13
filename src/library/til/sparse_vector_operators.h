#ifndef TIL_SPARSE_VECTOR_OPERATORS_H_
#define TIL_SPARSE_VECTOR_OPERATORS_H_

/// \file Belongs to sparse_vector package.
/// Do not include directly, include "til/sparse_vector.h" instead

// includes from BOOST
#include <boost/call_traits.hpp>

namespace til
{

// NB: to be undef'd at the end of this file
#define TIL_DEFINE_SPARSE_OP(op)                                                                  \
  template < typename T >                                                                         \
  sparse_vector<T>                                                                                \
  operator op (const sparse_vector<T> & a, typename boost::call_traits<T>::param_type value)      \
  {                                                                                               \
    sparse_vector<T> res(a.size());                                                               \
    typename sparse_vector<T>::Map::const_iterator i = a.getMap().begin();                        \
    typename sparse_vector<T>::Map::const_iterator end = a.getMap().end();                        \
    for (; i != end; ++i)                                                                         \
    {                                                                                             \
      res[i->first] = i->second op value;                                                         \
    }                                                                                             \
    return res;                                                                                   \
  }                                                                                               \
                                                                                                  \
  template < typename T >                                                                         \
  sparse_vector<T>                                                                                \
  operator op (typename boost::call_traits<T>::param_type value, const sparse_vector<T> & a)      \
  {                                                                                               \
    sparse_vector<T> res(a.size());                                                               \
    typename sparse_vector<T>::Map::const_iterator i = a.getMap().begin();                        \
    typename sparse_vector<T>::Map::const_iterator end = a.getMap().end();                        \
    for (; i != end; ++i)                                                                         \
    {                                                                                             \
      res[i->first] = value op i->second;                                                         \
    }                                                                                             \
    return res;                                                                                   \
  }                                                                                               \
  
  TIL_DEFINE_SPARSE_OP( + );
  TIL_DEFINE_SPARSE_OP( - );
  TIL_DEFINE_SPARSE_OP( * );

#undef TIL_DEFINE_SPARSE_OP

} // namespace til


#endif /*SPARSE_VECTOR_OPERATORS_H_*/
