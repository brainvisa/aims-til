#ifndef TIL_FRACTION_POLICIES_H_
#define TIL_FRACTION_POLICIES_H_

/// \file Contains policies for the Fraction class.
/// Do not include this file, include fraction.h instead.

namespace til { namespace policy
{

  //---------------------------------------------------------------------------

  /// Throw an exception
  struct ZeroByZero_Throw
  {
    class DivisionByZero : std::exception {};
          
    template < typename T >
    T operator()(T&,T&)
    {
      throw DivisionByZero();
    }
  };
  
  //---------------------------------------------------------------------------

  /// Return zero
  struct ZeroByZero_Zero
  {
    template < typename T >
    T operator()(T&,T&)
    {
      return T(0);
    }
  };

  //---------------------------------------------------------------------------

}} // namespace til::policy


#endif /*FRACTION_POLICIES_H_*/
