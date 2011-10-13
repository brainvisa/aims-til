#ifndef TIL_CHAMFERDISTANCE_H
#define TIL_CHAMFERDISTANCE_H

#include <vector>

#include "til/til_common.h"

#include "til/ImageExtrapolator.h"

namespace til
{
	
	template < typename TImage >
	void chamferDistance_3
	(
 	 TImage &im,
	 const std::vector<typename TImage::value_type> &mask,
	 typename TImage::value_type foreground,
	 typename TImage::value_type background
	)
	{
		typedef typename TImage::value_type value_type;
		/// TODO: this probably has to be changed
		/// First, k_infinity could be set to infinity if value_type has a representation
		/// for it (use has_infinity)
		/// Second, big (ok, theoretical) problem if value_type is unsigned float: then
		/// infinity is max - 1... which is numerically certainly equal to 1 due to rounding.
		const value_type k_doNotProcess = std::numeric_limits<value_type>::is_signed ? -1 : std::numeric_limits<value_type>::max();
		const value_type k_infinity = std::numeric_limits<value_type>::is_signed ? std::numeric_limits<value_type>::max() : std::numeric_limits<value_type>::max()-1;

		// Initialize the image by putting initial distances
		{
			typename Iterator<TImage>::Linear iIm(im);
			for (; !iIm.isAtEnd(); ++iIm)
			{
				if (*iIm == foreground)
				{
					// Distance to object within object is zero
					*iIm = 0;
				}
				else if (*iIm == background)
				{
					// Distance to object is infinity
					*iIm = k_infinity;
				}
				else
				{
					// Special label for other objects present in the image
					*iIm = k_doNotProcess;
				}
			}
		}

		//
		{
			typename Iterator<TImage>::Volumetric iIm(im);
	
			value_type x = mask[0];
			value_type y = mask[1];
			value_type z = mask[2];

			//ConstantExtrapolator

			for (; !iIm.isAtEnd(); ++iIm)
			{
				//*iIm = MIN(*iIm, iIm.getValue(-1,0,0) * x);
			}
		}
	}


} // namespace


#endif

