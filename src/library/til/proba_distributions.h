#ifndef TIL_PROBA_DISTRIBUTIONS_H_
#define TIL_PROBA_DISTRIBUTIONS_H_

// includes from STL

// includes from TIL
#include "til/miscTools.h"

namespace til
{

  //---------------------------------------------------------------------------

  // TODO: rajouter un namespace qqpart...
  
  // TODO: Je pense qu'au final il faudra bien avoir son propre generateur de nombre aleatoire,
  // parce qu'on veut pouvoir controler les seeds de maniere independante pour differentes instances
  // de generateurs.
  
  /// Uniform random distribution.
  template < typename T >
  class UniformRandomDistribution
  {
  public: // constructors
    /// Set min and max bounds.
    UniformRandomDistribution(T min, T max);
  public: // set & get
    T & min() { return m_min; }
    const T & min() const { return m_min; }
    T & max() { return m_max; }
    const T & max() const { return m_max; }
  public: // functions
    /// Draw a random number according to the current distribution.
    inline T operator()();
  private: // data
    T m_min;
    T m_max;
  };

  //---------------------------------------------------------------------------
 
  /// draw a random number according to the uniform distribution on [min, max].
  // TODO: rename into uniform_rand
  template < typename T >
  T rand(T min, T max);

  //---------------------------------------------------------------------------

} // namespace til

#include "til/proba_distributions.tpp"


#endif /*PROBA_DISTRIBUTIONS_H_*/
