#ifndef TIL_CUBIC_SPLINE_INTERPOLATION_H
#define TIL_CUBIC_SPLINE_INTERPOLATION_H

// includes from TIL library
#include "til/til_common.h"
#include "til/labels.h"

namespace til
{
  /// Cubic interpolation.
  template < typename T >
  class CubicSplineInterpolation
    : public Interpolator_label
  {
  public: // typedefs
  	typedef CubicSplineInterpolation<T> Self;
  	
  public: // static functions
  
  	/// Returns the cubic spline interpolation value of four numbers,
  	/// given their four spline coefficients. 
    /// The cubic spline coefficients
  	/// of the data to interpolate therefore must have been previoulsy
  	/// computed. To do that, one technique is to use recursive filtering,
  	/// see e.g. til::cubicSplineCoefficients for images.
  	/// The cubic spline is
  	/// 3/4*x^3-3/2*x^2+1,      0<x<1
  	/// -1/4*x^3+3/2*x^2-3*x+2, 1<x<2
  	// NB: with this method, I use 16 muls and 9 adds. By expanding the
  	// equation above and regrouping by coefficients of the polynomial in x,
  	// I get only to 15 muls and 11 adds -- not really better.
  	// Any better factorization?
    // TODO: use operator() instead
  	inline static T compute
  	(
    	T f1,	///< Cubic spline coefficient assumed to lie at position -1
    	T f2,	///< Cubic spline coefficient assumed to lie at position 0
    	T f3,	///< Cubic spline coefficient assumed to lie at position 1
    	T f4,	///< Cubic spline coefficient assumed to lie at position 2
    	T x		///< The interpolated value at position x, supposed to lie between 0 and 1 (no checking on that)
  	);
  };
} // namespace til


#include "til/cubic_spline_interpolation.tpp"


#endif

