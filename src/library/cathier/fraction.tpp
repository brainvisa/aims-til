

namespace til
{

  //---------------------------------------------------------------------------

  template < typename T, typename ZeroByZeroPolicy >
  T
  Fraction< T, ZeroByZeroPolicy >::
  operator()(T nom, T denom)
  {
    const T epsilon = 128*std::numeric_limits<T>::epsilon();
    if (denom < epsilon)
    {
      if (nom < epsilon)
      {
        return m_zeroByZeroPolicy(nom, denom);
      }
      return T(0);
    }
    return nom/denom;
  }

  //---------------------------------------------------------------------------

  template < typename ZeroByZeroPolicy, typename T >
  typename boost::enable_if_c<!std::numeric_limits<T>::is_integer, T>::type
  fraction(T nom, T denom)
  {
    Fraction<T, ZeroByZeroPolicy> f;
    return f(nom, denom);
  }

  //---------------------------------------------------------------------------
  
} // namespace til

