#ifndef TIL_SUBIMAGEEXTRACTOR_H
#define TIL_SUBIMAGEEXTRACTOR_H


// Standard library includes

#include <iostream>


// Local includes

#include "til/til_common.h"

#include "til/Range.h"
#include "til/imageTools.h"
#include "til/miscTools.h"


// Extract a sub-image out of an original image
// If the crop box is out of the image range, the image is extrapolated


// Namespace 

namespace til {


template <class Extrapolator, class TImageIn, class TImageOut>
void extractSubImage(const TImageIn & im, 
					 TImageOut & subIm,
					 const Range<int,3> & cropbox)
{
	typedef typename TImageIn::value_type TPixelIn;
	typedef typename TImageOut::value_type TPixelOut;

	// Check whether input image is allocated
	allocationCheck(im);

	// Allocate output image
	subIm.init(cropbox.dims(), im.vdim());
  /*
		cropbox.getSizeX(),
		cropbox.getSizeY(),
		cropbox.getSizeZ(),
		im.vdim()[0],
		im.vdim()[1],
		im.vdim()[2]);
    */
	
	// set origin of subimage
	//subIm.setOrigin(im.getOrigin() - cropbox.min_bounds());
		
	typename Iterator<TImageOut>::Linear iSubIm(subIm);

	// extract sub image
	numeric_array<int,3> p;
	for (p[2] = cropbox.min_bounds()[2]; p[2] <= cropbox.max_bounds()[2]; ++p[2])
	{
		for (p[1] = cropbox.min_bounds()[1]; p[1] <= cropbox.max_bounds()[1]; ++p[1])
		{
			for (p[0] = cropbox.min_bounds()[0]; p[0] <= cropbox.max_bounds()[0]; ++p[0])
			{
				*iSubIm = castValue<TPixelIn, TPixelOut>(Extrapolator::getValue(im, p));
				/*
				if (contains(im, i, j, k))
				{
					*iSubIm = castValue<TPixelIn, TPixelOut>(im.getValue(i,j,k));
				}
				else
				{
					*iSubIm = castValue<TPixelIn, TPixelOut>(Extrapolator::getExtrapolatedValue((const TImageIn*)im, i, j, k));
				}
				*/
				++iSubIm;
			}
		}
	}
}

} // namespace


#endif

