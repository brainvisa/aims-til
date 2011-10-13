#ifndef TIL_STD_WRAP_H_
#define TIL_STD_WRAP_H_

// includes from STL
#include <algorithm>

// includes from BOOST
#include "boost/call_traits.hpp"

// includes from TIL
#include "til/basic_range.h"
#include "til/ext_declarations.h"
#include "til/traits.h"

///\file Wrap of STL functions to take objects as input.
// I think it is important to pass objects rather than iterators to these functions.
// I mean, from what I understand, the reasons why STL is passing iterators around is for technical
// reasons (and I mean technical ease, here) -- one of them being compatibility with C-style
// arrays, others being I think const and reference type attribute headache.

// Okay, actually it's more suble: the iterator can contain part of the algorithm itself,
// e.g. reverse or insert or whatever. Such iterators actually extends the power of the
// algorithm to degrees the algorithm is not even aware of. That could be taken care of
// here by using extra template arguments reflecting the iterator types to use in the algo.

// The pb is that iterators generally know almost nothing about the objects they are pointing to.
// STL iterators are quite plain and thus manipulating iterators is okay. But when we try to do
// more complicated stuff, like compressed arrays, this hurts.
// I'm perfectly aware that precisely, compressed arrays and stuff should rather not be STL compliant,
// because if they were, it would be too easy to misuse them.
// Still: if we want to have a function that is generic and would accept both, we need these functions.

namespace til
{
  /// Size of a container.
  template < typename TContainer >
  inline
  //typename boost::enable_if<til::is_container<TContainer>, std::size_t>::type
  std::size_t
  size(const TContainer &c)
  {
    return c.size();
  }
  
  /// Wrapper of std::fill
  template < typename TContainer >
  inline void fill(TContainer & c, typename boost::call_traits<typename TContainer::value_type>::param_type value)
  {
    std::fill(c.begin(), c.end(), value);
  }

  // Commented out because I think that's already the case
  /*
  /// Use lexicographical order as default order on boost::array
  template<class T, std::size_t N>
  bool operator<(const boost::array<T,N> & x, const boost::array<T,N> & y)
  {
      return std::lexicographical_compare(x.begin(),x.end(),y.begin(),y.end());
  }
  */


  //------------------------------------------------------------------------------


  // range for std::vector
  template < typename T, typename TAlloc >
  struct range_of<std::vector<T,TAlloc> >
  { typedef basic_range<typename std::vector<T,TAlloc>::iterator> type; };
  template < typename T, typename TAlloc >
  struct const_range_of<std::vector<T,TAlloc> >
  { typedef basic_range<typename std::vector<T,TAlloc>::const_iterator> type; };

  template < typename T, typename TAlloc >
  typename range_of<std::vector<T,TAlloc> >::type
  whole_range(std::vector<T,TAlloc> & v)
  { return typename range_of<std::vector<T,TAlloc> >::type (v.begin(), v.end()); }
  template < typename T, typename TAlloc >
  typename const_range_of<std::vector<T,TAlloc> >::type
  whole_range(const std::vector<T,TAlloc> & v)
  { return typename const_range_of<std::vector<T,TAlloc> >::type (v.begin(), v.end()); }


  // range for std::list
  template < typename T, typename TAlloc >
  struct range_of<std::list<T,TAlloc> >
  { typedef basic_range<typename std::list<T,TAlloc>::iterator> type; };
  template < typename T, typename TAlloc >
  struct const_range_of<std::list<T,TAlloc> >
  { typedef basic_range<typename std::list<T,TAlloc>::const_iterator> type; };

  template < typename T, typename TAlloc >
  typename range_of<std::list<T,TAlloc> >::type
  whole_range(std::list<T,TAlloc> & v)
  { return typename range_of<std::list<T,TAlloc> >::type (v.begin(), v.end()); }
  template < typename T, typename TAlloc >
  typename const_range_of<std::list<T,TAlloc> >::type
  whole_range(const std::list<T,TAlloc> & v)
  { return typename const_range_of<std::list<T,TAlloc> >::type (v.begin(), v.end()); }


} // namespace til

/*
/// Size of a boost::array
template < typename T, std::size_t D >
inline
std::size_t
size(const boost::array<T,D> &) { return D; }
*/

#endif /*STD_WRAP_H_*/
