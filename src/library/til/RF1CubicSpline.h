#ifndef TIL_RECURSIVE_FILTER_CUBIC_SPLINE_H
#define TIL_RECURSIVE_FILTER_CUBIC_SPLINE_H

// includes from STL
#include <cmath>	// sqrt

// includes from TIL library
#include "til/til_common.h"
#include "til/rf1.h"


namespace til {

	/// Replaces im by a pre-processed version of im for cubic spline
	/// interpolation.
	/// NB: out can be the same image as in.
	// TODO: so far there is no choice in extrapolation. However recursive
	// filters propose a small choice. We could enable this choice, but then the
	// size of im might change! Might be a problem here...
	// Maybe a wrapper on an image class with templating on the shift might just
	// help, WindowedImage<TImage, x0,y0,z0,x1,y1,z1>?
	// We may also enable to be able to translate from the standard image extrapolation
	// classes. In this case I see yet another trait coming...
	template < typename TImage >
	void cubicSplineCoefficients(const TImage &im, TImage &coeffs)
	{
		typedef typename TImage::value_type value_type;
		const value_type z1 = ::sqrt(3.0) - 2;

		// Check that im and coeffs are allocated and have the same size
		similarityCheck(im, coeffs);

		RecursiveFilter<value_type, +1, 1> filter1;
		RecursiveFilter<value_type, -1, 1> filter2;

		filter1.setFilter(1, -z1);
		filter2.setFilter(-z1, -z1);

	
		for (ImageAxis axis = X_AXIS; axis <= Z_AXIS; ++axis)
		{
			// Apply filter only if image is not flat in that direction
			if (im.dim()[axis] > 1)
			{
				filterAlongAxis(im, coeffs, axis, filter1);
				filterAlongAxis(im, coeffs, axis, filter2);
			}
		}
	}
} // namespace

#endif

