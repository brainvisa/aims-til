#ifndef MATH_FUNCTIONS_H_
#define MATH_FUNCTIONS_H_


#include <functional>

// TODO: I guess all this crap should rather go into functor? I mean, why Exp in functor and Gaussian in Math?

namespace til { namespace math {

  //---------------------------------------------------------------------------

    //---------//
   //  huber  //
  //---------//

  /// Huber's penalty function.
  /// x^2 if |x| < K
  /// 2*K*|x| - K^2 if |x| > K
  template < typename T >
  inline
  T huber(T x, T K);

  //---------------------------------------------------------------------------

    //-----------//
   //  d_huber  //
  //-----------//

  /// Derivative of Huber's penalty function.
  /// 2*x if |x| < K
  /// +/-2*K if |x| > K
  template < typename T >
  inline
  T d_huber(T x, T K);

  //---------------------------------------------------------------------------

    //----------//
   //  shrink  //
  //----------//
  
  /// The shrink operator
  /// if |x| > K then x +/- K
  /// else 0
  template < typename T >
  inline
  T shrink(T x, T K);
  
  //---------------------------------------------------------------------------
  
    //------------//
   //  Gaussian  //
  //------------//

  /// Unnormalized Gaussian.
  /// The equation of this function is exp(-x^2 / ( 2 * sigma^2 ) ), where sigma is the standard deviation.
  /// It is unnormalized: there is no multiplicative coefficient and the value at 0 is always 1.
  template < typename TPrec >
  class Gaussian
    : public std::unary_function<TPrec, TPrec>
  {
  public: // typedefs
    
    typedef TPrec prec_type;

  public: // exceptions
    
    class StandardDeviationTooSmall : public std::exception {};
    
  public: // constructors
  
    /// Constructor with Gaussian standard deviation.
    explicit Gaussian(prec_type sigma);
    
  public: // set & get
  
    prec_type getSigma() const { return std::sqrt(m_2var / 2); }
    
    inline void setSigma(prec_type sigma);
    
  public: // operators
  
    /// Return the value of the Gaussian kernel at x.
    prec_type operator()(prec_type x) const;
    
    /// Return the value of the Gaussian kernel at x, the input being x^2.
    prec_type squared_input(prec_type x2) const;
    
  private: // data
    prec_type m_2var;
  };

  //---------------------------------------------------------------------------

    //---------------------------//
   //  IsotropicGaussianKernel  //
  //---------------------------//

  /// Unnormalized isotropic Gaussian kernel.
  // NB: as a reminder to myself: should be labeled 'kernel' a function of two (spatial) parameters.
  // So the 'gaussian kernel' really is exp(- ||x - y||^2)
  // as opposed to the standard scalar gaussian function exp (-x^2)
  // TODO: There should be a 'scalar matrix' class, with just a number. Then the GaussianKernel
  // would be templated over the matrix type. When this matrix type is a scalar matrix, we get
  // the fast isotropic gaussian kernel for free.
  template < typename TArray, typename TPrec >
  class IsotropicGaussianKernel
    : public std::binary_function<const TArray &, const TArray &, TPrec>
  {
  public: // typedefs
    //typedef typename TArray::value_type prec_type;
    typedef TPrec prec_type;
  public: // constructors
    explicit IsotropicGaussianKernel(prec_type sigma) : m_gaussian(sigma) {}
  public: // set & get
    void setSigma(prec_type sigma) { m_gaussian.setSigma(sigma); }
    prec_type getSigma() const { return m_gaussian.getSigma(); }
  public: // functions
    prec_type operator()(const TArray & x, const TArray & y) const
    {
      return m_gaussian.squared_input(dist2(x, y, prec<TPrec>()));
    }
  private: // data
    Gaussian<prec_type> m_gaussian;
  };
  
  //---------------------------------------------------------------------------
  
}} // namespace til::math

// package include
#include "math_functions.tpp"

#endif /*MATH_FUNCTIONS_H_*/
