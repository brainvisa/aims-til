#ifndef TIL_MATH_H
#define TIL_MATH_H

// include from TIL library
#include "til/til_common.h"

namespace til {

	/// Namespace for misc. math definitions
	namespace math {


		/// sinus cardinal
		INLINE double sinc(double x)
		{
			const double EPSILON = 128 * std::numeric_limits<double>::epsilon();
			// If x is too small, return assymptotic value
			if (fabs(x) <= EPSILON) return 1.0;
			// formula
			return sin(x) / x;
		}


		/// first derivative of sinc
		INLINE double dsinc(double x)
		{
			const double EPSILON = 128 * std::numeric_limits<double>::epsilon();
			// If x is too small, return assymptotic value
			if (fabs(x) <= EPSILON) return 0.0;
			// formula
			return cos(x) / x - sin(x) / (x*x);
		}

		/// second derivative of sinc
		INLINE double d2sinc(double x)
		{
			const double EPSILON = 128 * std::numeric_limits<double>::epsilon();
			// if x is too small, return assymptotic value
			if (fabs(x) <= EPSILON) return -1.0 / 3.0;	
			// formula
			return 2*sin(x)/(x*x*x) - 2*cos(x)/(x*x) - sin(x)/x;
		}

	} // namespace math
} // namespace til

#endif

