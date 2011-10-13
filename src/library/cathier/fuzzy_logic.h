#ifndef TIL_FUZZY_LOGIC_H
#define TIL_FUZZY_LOGIC_H

/// \file Basic fuzzy logic/floating point logic stuff.
/// OK, actually, it's not really fuzzy stuff since it is not continous
/// but tri-bool stuff...

// includes from STL
#include <limits>

// includes from BOOST
#include <boost/logic/tribool.hpp>


namespace til { namespace fuzzy
{

  //---------------------------------------------------------------------------

  /// Tests whether x is positive, with a +- delta uncertainty cushion.
  /// In other words: 
  /// if x > delta, x is definitely positive, and the result is true
  /// if x < -delta, x is definitely negative, and the result is false
  /// otherwise, the result is indeterminate.
  template < typename T >
  inline
  boost::logic::tribool
  is_positive(T x, T delta = 128*std::numeric_limits<T>::epsilon());

  //---------------------------------------------------------------------------

  /// Returns true iff both arguments have same sign, allowing some degree of
  /// imprecision.
  // NB: I think it is a good idea to keep the name different from the function above, rather than overloading 
  // using a disable_if on numeric_limits::is_integer. The mechanism is really different and we loose some properties
  // like transitivity by using this fuzzy version, so better keep the user aware of the change.
  template < typename T >
  inline 
  boost::logic::tribool
  same_sign(T x, T y, T delta = 128*std::numeric_limits<T>::epsilon());

  //---------------------------------------------------------------------------

}} // namespace til::fuzzy

// package include
#include "fuzzy_logic.tpp"

#endif

