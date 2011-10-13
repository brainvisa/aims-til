#ifndef TIL_SCALING_MAP_TOOLS_H
#define TIL_SCALING_MAP_TOOLS_H

/// \file Belongs to ScalingMap package -- do not include directly, include til/ScalingMap.h instead

// includes from TIL library
#include "til/til_declarations.h"


namespace til {


	/// Computes the mapping that transforms im1 into im2, assuming that
	/// they cover exactly the same volume in space.
	/// This computation hence discards any information from world coordinates.
	template < class TImage1, class TImage2 >
	ScalingMap<double>
	getScalingBetween(const TImage1 & from, const TImage2 & to)
	{
		allocationCheck(from, to);
		// To prevent from dividing by zero below
		if (from.size() == 0)
		{
			throw std::invalid_argument("Input image has no elements");
		}
		numeric_array<double,3> mulFactor;
    convert(mulFactor, to.dim());
    mulFactor /= from.dim();
		return ScalingMap<double>(mulFactor, 0.5*mulFactor - 0.5);
	}


} // namespace

#endif

