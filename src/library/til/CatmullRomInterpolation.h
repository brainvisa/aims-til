#ifndef TIL_CATMULL_ROM_INTERPOLATION_H
#define TIL_CATMULL_ROM_INTERPOLATION_H

// includes from TIL library
#include "til/til_common.h"
#include "til/labels.h"

namespace til
{

  //----------------------------------------------------------------------------------
  
    //---------------------------//
   //  CatmullRomInterpolation  //
  //---------------------------//

  /// Interpolation using the Catmull-Rom polynomial.
  /// Computes the Catmull-Rom interpolation of the sequence (f1, f2, f3, f4)
  /// (supposed to be equally spaced and at position (-1, 0, 1, 2)).
  /// 'x' is the position of the point where we want an interpolated value,
  /// and has to be in the range [0, 1].
  ///
  /// The Catmull-Rom spline is 3/2.x^3-5/2.x^2+1, -1/2.x^3 + 5/2.x^2-4.x+2.
  /// It it the only 3rd order interpolating spline of support of length 4 that
  /// is C^1 and has no 2nd order moment.
  /// It has been implemented using the factorization of (Meijering03) that
  /// needs the computation of only two cubic polynomials per point.
  // NOTE: The templating has been done on the function members rather than on
  // the class itself, to make use of the automatic template deduction of functions.
  // It may be changed.
  
  template < typename T >
  class CatmullRomInterpolation
    : public Interpolator_label
  {
  public: // typedefs
  	typedef CatmullRomInterpolation Self;
  
  public: // static method
  
  	/// Returns the interpolated value of four numbers using
  	/// the Catmull-Rom method.
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
  	);
  };
    
  //----------------------------------------------------------------------------------
  
} // namespace til

#include "til/catmull_rom_interpolation.tpp"

#endif

