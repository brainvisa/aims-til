#ifndef DECLARATIONS_EXTERNAL_H_
#define DECLARATIONS_EXTERNAL_H_

/// \file Forward declarations of external classes.

namespace boost
{
  template < typename T, std::size_t D >          class array;
  template < typename T >                         class shared_ptr;
}

namespace std
{
  template < typename T, typename TAlloc >        class list;
  template < typename TKey, typename TValue, typename TCompare, typename TAlloc >
                                                  class std::map;
  template < typename TKey, typename TValue, typename TCompare, typename TAlloc >
                                                  class std::multimap;  
  template < typename TKey, typename TCompare, typename TAlloc >
                                                  class set;
  template < typename T, typename TAlloc >        class vector;
}

#endif /*DECLARATIONS_EXTERNAL_H_*/
