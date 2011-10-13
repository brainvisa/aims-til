#ifndef TIL_CONVOLUTION_H
#define TIL_CONVOLUTION_H

// local includes

#include "til/til_common.h"

namespace til
{

	// NOTE: this is probably not the way to go.
	// Would rather define a Action class to do point-wise correlation
	// Thus I could be able to have convolution method independant of
	// say boundary conditions


	/// Correlation between an image and a mask.
	/// Similar to convolution, except that the coordinates of the mask
	/// are not inverted. For symmetric masks, convolution and correlation are
	/// strictly the same.
	template < typename Extrapolation, typename TImageIn, typename TImageMask, typename TImageOut >
	void correlation(
	const TImageIn &im,				///< the input image
	const TImageMask &mask,			///< the correlation mask
	TImageOut &out					///< the result image, already allocated
	)
	{
		typename TImageIn::ConstVolumetricIterator iIm(im);
		typename TImageMask::ConstVolumetricIterator iMask(mask);
		typename TImageOut::ConstVolumetricIterator iOut(out);

		for (; !iIm.isAtEnd(); ++iIm, ++iOut)
		{
		}
	}

	/// So-called 'Symmetric Neighbor' correlation.
	/// For each couple of symmetric points, only the one having an intensity
	/// closest to the intensity of the center point is used. Therefore
	/// only half of the mask points are used; one should then use masks
	/// slightly bigger than usual (theoretically, two times bigger in volume)
	template < typename Extrapolation, typename TImageIn, typename TImageMask, typename TImageOut >
	void correlationSN(
	const TImageIn &im,
	const TImageMask &mask,
	TImageOut &out)
	{
	}

} // namespace


#endif

