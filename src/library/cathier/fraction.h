#ifndef FRACTION_H_
#define FRACTION_H_

/// \file The content of this file is entirely dedicated to... dividing one number by another.
/// I know, this sounds kind of dull. Doing a x/y shouldn't take that many lines of code...
/// The problem arise of course when y is zero. Then, shit begins, and usually the best way people find
/// to deal with it is just to ignore it. By using this predefined code, however, I ensure that I have
/// to deal with it, and then, all the crap is already programmed.

/// In short: when you are not sure x/y will never fail and have no easy way to deal with it, you will
/// need some lines of code to deal with it, and fortunately they are all there. If you can change your
/// formula so that x/y will never fail, then of course it's better and you can ignore this.


// includes from STL
#include <exception>

namespace til
{
  
  //---------------------------------------------------------------------------
  
  /// A class to divide one number by another. 
  /// The behavior to follow in case of a zero / zero is given by a template policy.
  template < typename T, typename ZeroByZeroPolicy >
  class Fraction
  {
  public: // constructors -----------------------------------------------------
    Fraction() : m_zeroByZeroPolicy() {}
    Fraction(ZeroByZeroPolicy p) : m_zeroByZeroPolicy(p) {}
  public: // functions --------------------------------------------------------
    /// Returns nom / denom. The 0/0 case is handled by the policy.  
    T operator()(T nom, T denom);
  private: // data ------------------------------------------------------------
    ZeroByZeroPolicy m_zeroByZeroPolicy;
  };
  
  //---------------------------------------------------------------------------

  /// Divides one number by another.
  /// In case of a zero / zero, do what the template policy says.
  template < typename ZeroByZeroPolicy, typename T >
  typename boost::enable_if_c<!std::numeric_limits<T>::is_integer, T>::type
  fraction(T nom, T denom);

  //---------------------------------------------------------------------------

} // namespace til

// package include
#include "fraction.tpp"
#include "fraction_policies.h"

#endif /*FRACTION_H_*/
