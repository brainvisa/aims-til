#ifndef TIL_LINEAR_INTERPOLATION_H
#define TIL_LINEAR_INTERPOLATION_H
// TODO: this file should be rename "Linear interpolation"

// includes from TIL library
#include "til/til_common.h"
#include "til/labels.h"

namespace til
{

  /// Linear interpolation between two values
  template < typename T >
  class LinearInterpolation
    : public Interpolator_label
  {
  public: // typedefs
  
  	typedef LinearInterpolation Self;
  	
  public:
  
  	/// Returns the linear interpolated value of two numbers
  	// NB: There is a priori no reason for f to be the same type as x.
  	// However, forcing x and f to have the same type forces the
  	// (necessary anyway) casting to be done once, before function
  	// call.
  	inline static T compute
  	(
    	T f1,	///< Value assumed to lie at position 0
    	T f2,	///< Value assumed to lie at position 1
    	T x		///< The interpolated value at position x
  	)
  	{
  		//return f2*x+f1*(1-x);
  		// this saves one multiplication!
  		return f1 + x*(f2-f1);
  	}
  };

} // namespace til


#endif

