#ifndef TIL_KEYS_INTERPOLATION_H
#define TIL_KEYS_INTERPOLATION_H

// includes from TIL library
#include "til/til_common.h"
#include "til/labels.h"

namespace til
{
  
  //-----------------------------------------------------------------------------
  
    //---------------------//
   //  KeysInterpolation  //
  //---------------------//
  
	/// Interpolation using Keys polynomial.
	/// Computes the Keys interpolation of the sequence (f1, f2, f3, f4, f5, f6)
	/// (supposed to be equally spaced and at position (-2, -1, 0, 1, 2, 3)).
	/// 'x' is the position of the point where we want an interpolated value,
	/// and has to be in the range [0, 1].
	///
	/// The Keys spline is 4/3.x^3-7/3.x^2+1, -7/12.x^3 + 3.x^2-59/12.x+5/2,
	/// 1/12.x^3-2/3.x^2+7/4.x-3/2.
	/// It has been implemented using the factorization of (Meijering03) that
	/// needs the computation of only four cubic polynomial per point.
	//
	// NOTE: The templating has been done on the function members rather than on
	// the class itself, to make use of the automatic template deduction of functions.
	// It may be changed.
  template < typename T >
	class KeysInterpolation
    : public Interpolator_label
	{
	public: // typedefs ---------------------------------------------------------
		
    typedef KeysInterpolation Self;
	
	public: // static methods ---------------------------------------------------

		/// Returns the interpolated value of four numbers using
		/// the Keys method.
		// NB: There is a priori no reason for f to be the same type as x.
		// However, forcing x and f to have the same type forces the
		// (necessary anyway) casting to be done once, before function
		// call.
		inline static T compute(
  		T f1,	///< Value assumed to lie at position -2
  		T f2,	///< Value assumed to lie at position -1
  		T f3,	///< Value assumed to lie at position 0
  		T f4,	///< Value assumed to lie at position 1
  		T f5,	///< Value assumed to lie at position 2
  		T f6,	///< Value assumed to lie at position 3
  		T x		///< The interpolated value at position x
		);
	};

} // namespace til

#include "til/keys_interpolation.tpp"

#endif

