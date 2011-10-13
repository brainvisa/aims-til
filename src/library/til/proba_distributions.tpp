#ifndef PROBA_DISTRIBUTIONS_TPP_
#define PROBA_DISTRIBUTIONS_TPP_

// includes from STL
#include <cassert>
#include <cstdlib>  // std::rand
#include <limits>

namespace til
{
  
  //---------------------------------------------------------------------------
  
  template < typename T >
  UniformRandomDistribution<T>::UniformRandomDistribution(T min, T max)
    : m_min(min)
    , m_max(max)
  { assert(min <= max); }
  
  //---------------------------------------------------------------------------
 
  template < typename T >
  T UniformRandomDistribution<T>::operator()()
  {
    if (std::numeric_limits<T>::is_integer)
    {
      // Epsilon correspond to the safety margin we want on the side of the
      // continuous uniform distribution, to be sure that when rounding we
      // never get min-1 or max+1
      // NB: with epsilon == 0 , I once obtained max+1 !!
      // TODO: is this much better / efficient than two tests (on min and max)?
      const double eps = 32 * std::numeric_limits<double>::epsilon();
      T res = castValue<double, T>(std::rand() / double(RAND_MAX) * (m_max - m_min + 1 - 2*eps) + m_min - 0.5 + eps);
      assert(res >= m_min);
      assert(res <= m_max);
      return res;
    }
    else
    {
      // TODO: is there the same issue here as above?
      T res = castValue<double, T>(std::rand() / double(RAND_MAX) * (m_max - m_min) + m_min);
      assert(res >= m_min);
      assert(res <= m_max);
      return res;
    }
  }
  
  //---------------------------------------------------------------------------

  template < typename T >
  T rand(T min, T max)
  {
    UniformRandomDistribution<T> u(min, max);
    return u();
  }
  
  //---------------------------------------------------------------------------
  
} // namespace til

#endif /*PROBA_DISTRIBUTIONS_TPP_*/
