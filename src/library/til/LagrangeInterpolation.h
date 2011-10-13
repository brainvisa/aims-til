#ifndef TIL_LAGRANGE_INTERPOLATION_H
#define TIL_LAGRANGE_INTERPOLATION_H

// includes from TIL library
#include "til/til_common.h"
#include "til/labels.h"

namespace til
{

  /// Interpolation using Lagrange polynomials of order 4.
  /// Computes the Lagrange interpolation of the sequence (f1, f2, f3, f4)
  /// (supposed to be equally spaced and at position (-1, 0, 1, 2)).
  /// 'x' is the position of the point where we want an interpolated value,
  /// and has to be in the range [0, 1].
  /// The Lagrange polynomial is the only third order polynomial fitting all four
  /// points.
  // NOTE: The templating has been done on the function members rather than on
  // the class itself, to make use of the automatic template deduction of functions.
  // It may be changed.
  template < typename T >
	class Lagrange4Interpolation
    : public Interpolator_label
	{
	public: // static members ---------------------------------------------------

		/// Returns the interpolated value of four numbers using Lagrange4 method.
		// NB: There is a priori no reason for f to be the same type as x.
		// However, forcing x and f to have the same type forces the
		// (necessary anyway) casting to be done once, before function
		// call.
		inline static T compute(
  		T f1,	///< Value assumed to lie at position -1
  		T f2,	///< Value assumed to lie at position 0
  		T f3,	///< Value assumed to lie at position 1
  		T f4,	///< Value assumed to lie at position 2
  		T x		///< The interpolated value at position x
		)
		{
			return (((3*(f2-f3)+f4-f1)/6.0*x+((f1+f3)/2.0-f2))*x+(6*f3-2*f1-3*f2-f4)/6.0)*x+f2;
		}
	};

}

#endif

