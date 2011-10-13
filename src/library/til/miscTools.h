#ifndef TIL_MISCTOOLS_H
#define TIL_MISCTOOLS_H

// includes from STL
#include <cassert>
#include <iostream>
#include <limits>
#include <cstdio>		// sprintf

// includes from TIL library
#include "til/cat2type.h"
#include "til/traits.h"

namespace til
{

  //---------------------------------------------------------------------------

  /// Returns a number (a positive integer) with a different value than the input numbers.
  template < typename T >
  T findValueOtherThan(T v1, T v2)
  {
    T res;
    for (res = 0; res == v1 || res == v2; ++res) {}
    return res;
  }
  
  //---------------------------------------------------------------------------

  /// Returns a number (a positive integer) with a different value than the input numbers.
  template < typename T >
  T findValueOtherThan(T v1, T v2, T v3)
  {
    T res;
    for (res = 0; res == v1 || res == v2 || res == v3; ++res) {}
    return res;
  }
  
  //---------------------------------------------------------------------------

  /// Returns a number (a positive integer) with a different value than the input numbers.
  template < typename T >
  T findValueOtherThan(T v1, T v2, T v3, T v4)
  {
    T res;
    for (res = 0; res == v1 || res == v2 || res == v3 || res == v4; ++res) {}
    return res;
  }
  
  //---------------------------------------------------------------------------
  
  /// Returns a number (a positive integer) with a different value than the input numbers.
  template < typename T >
  T findValueOtherThan(T v1, T v2, T v3, T v4, T v5)
  {
    T res;
    for (res = 0; res == v1 || res == v2 || res == v3 || res == v4 || res == v5; ++res) {}
    return res;
  }

  //---------------------------------------------------------------------------


namespace
{
  template < typename TFrom, typename TTo >
  inline TTo castValue_impl(TFrom value, bool_type<true>)
  {
    // TODO: I'm afraid we have to use a cat2type to avoid a warning here... :(
    if (value < 0)
    //if (std::numeric_limits<TFrom>::is_signed && value < 0)
    {
      return TTo(value - TFrom(0.5));
    }
    else
    {
      return TTo(value + TFrom(0.5));
    }
  }

  template < typename TFrom, typename TTo >
  inline TTo castValue_impl(TFrom value, bool_type<false>)
  {
    return TTo(value);
  }
}

// Cast a numerical value from type TFrom to type TTo using the closest value
template < typename TFrom, typename TTo >
inline TTo castValue(TFrom value)
{
  // NB: I used that instead of a more readable cat2type, because it involved two tests, and cat2type would get more
  // complicated... Those two tests cannot be collapsed into one for cat2type because this test involves two template
  // parameters, while cat2type is templated over a single templated parameter...
  typedef bool_type<!std::numeric_limits<TFrom>::is_integer && std::numeric_limits<TTo>::is_integer> test_type;
  return castValue_impl<TFrom,TTo>(value, test_type());
/*
	if ( !std::numeric_limits<TFrom>::is_integer && 
		  std::numeric_limits<TTo  >::is_integer)
	{
    // TODO: I'm afraid we have to use a cat2type to avoid a warning here... :(
		if (value < 0)
    //if (std::numeric_limits<TFrom>::is_signed && value < 0)
		{
			return TTo(value - TFrom(0.5));
		}
		else
		{
			return TTo(value + TFrom(0.5));
		}
	}
	else
	{
		return TTo(value);
	}
  */
}

// NB : the assignment makes it possible to use ImageAxis as an array index
// following this precise convention.
// TODO: this image axis crap should go away.

enum ImageAxis { X_AXIS = 0, Y_AXIS = 1, Z_AXIS = 2 };

TIL_API ImageAxis operator++(ImageAxis &axis);


  //---------------------------------------------------------------------------

  /// Returns the min of a type.
  /// This is used as a work-around for std::numeric_limits::min(), which has a different meaning
  /// for integer and non-integer types.
  template < typename T >
  inline T type_min()
  {
  	return std::numeric_limits<T>::is_integer?
  		std::numeric_limits<T>::min() :
  	-std::numeric_limits<T>::max();
  }

  //---------------------------------------------------------------------------

  /// Shifts value to the right, i.e c = b, b = a.
  /// The value of c is lost; the value of a does not change. Note that a can be an rvalue.
  template < typename T >
  inline
  void shift_right(T a, T & b, T & c)
  {
    c = b;
    b = a;
  }
  
  //---------------------------------------------------------------------------

  /// Shifts value to the right, i.e d = c, c = b, b = a.
  /// The value of d is lost; the value of a does not change. Note that a can be an rvalue.
  template < typename T >
  inline void shift_right(T a, T & b, T & c, T & d)
  {
    d = c;
    c = b;
    b = a;
  }

  //---------------------------------------------------------------------------

  /// Compute the N-root of an integer.
  /// Right now, just a stupid, super-brute-force algorithm.
  /// NB: T can actually be a non integer type, this does not have to be a restriction, and at the same time
  /// the name of the class advertise quite clearly what it does, so I guess it's OK.
  template < typename T >
  class IntegerRoot
  {
  public: // exceptions
    class NoRoot : public std::exception {};
  public: // operators
    /// Return an N-root of integer n, that is, the integer i so that i^N = n.
    /// If several solutions are possible, returns the positive solution. This happens only if n > 0 and N%2 = 0, and
    /// in that case, the other solution is -n.
    /// If there exist no integer solution, throws an exception.
    T operator()(T n, std::size_t N);
  private: // functions
    T compute(T n, std::size_t N, label::Passed<is_signed>);
    T compute(T n, std::size_t N, label::Failed<is_signed>);
    T solve(T n, std::size_t N);
  };
  
  //---------------------------------------------------------------------------

  /// Return the n-th root of integer i.
  /// Throws an exception if unsuccessful
  template < typename T >
  T integer_root(T i, std::size_t n)
  { return IntegerRoot<T>()(i, n); }

  //---------------------------------------------------------------------------

  /// Sort a, b, c in decreasing order (a>=b>=c).
  template < typename T >
  inline void sort(T & a, T & b, T & c);
  
  //---------------------------------------------------------------------------

  // Like sprintf, but returns the string -- usefull to concatenate in other
  // functions.
  // TODO: check std::string for similar use... that would spare us a .o file
  TIL_API char * mySprintf(const char * format, ... );

  //---------------------------------------------------------------------------

  /// Little helper for something often used when looking for a maximum.
  template < typename T >
  inline void max_helper(T & maxx, T x)
  { if (maxx < x) maxx = x; }
    
  //---------------------------------------------------------------------------

  /// Little helper for something often used when looking for a minimum.
  template < typename T >
  inline void min_helper(T & minx, T x)
  { if (minx > x) minx = x; }

  //---------------------------------------------------------------------------

  template < typename TPIterator, typename TSIterator >
  inline std::size_t pos2offset(TPIterator pbegin, TPIterator pend, TSIterator sbegin)
  {
    std::size_t res = *pbegin;
    while (++pbegin != pend)
    {
      ++sbegin;
      res *= *sbegin;
      res += *pbegin;
    }
    return res;
  }

  //---------------------------------------------------------------------------

  template < std::size_t D >
  inline std::size_t pos2offset(const numeric_array<std::size_t,D> & pos, const numeric_array<std::size_t,D> & size)
  {
    std::size_t i = D-1;
    std::size_t res = pos[i];
    // NB: we avoid using a 'for' loop because i is unsigned.
    while (i != 0)
    {
      --i;
      res *= size[i];
      res += pos[i];
    }
    return res;
  }

  //---------------------------------------------------------------------------

} // namespace til

// package include
#include "misc_tools.tpp"

#endif

