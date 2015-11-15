#ifndef DECLARATIONS_EXTERNAL_H_
#define DECLARATIONS_EXTERNAL_H_

#ifndef _GLIBCXX_BEGIN_NAMESPACE_CXX11
#define _GLIBCXX_BEGIN_NAMESPACE_CXX11
#define _GLIBCXX_END_NAMESPACE_CXX11
#endif

/// \file Forward declarations of external classes.

namespace boost
{
  template < typename T, std::size_t D >          class array;
  template < typename T >                         class shared_ptr;
}

namespace std
{
  // in gcc 5.2 list is defined in std::__cxx11 namespace
_GLIBCXX_BEGIN_NAMESPACE_CXX11
  template < typename T, typename TAlloc >        class list;
_GLIBCXX_END_NAMESPACE_CXX11

  template < typename TKey, typename TValue, typename TCompare, typename TAlloc >
                                                  class map;
  template < typename TKey, typename TValue, typename TCompare, typename TAlloc >
                                                  class multimap;
  template < typename TKey, typename TCompare, typename TAlloc >
                                                  class set;
  template < typename T, typename TAlloc >        class vector;
}

#endif /*DECLARATIONS_EXTERNAL_H_*/
