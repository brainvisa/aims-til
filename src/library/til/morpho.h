#ifndef TIL_MORPHO_H
#define TIL_MORPHO_H

// includes from STL library
#include <memory>
#include <vector>

// includes from TIL library
#include "til/til_common.h"
#include "til/miscTools.h"
#include "til/numeric_array.h"

namespace til
{

/*
/// This has to be rewritten (first, to template with neighborhood, and then,
/// because we probably don't need a function for this).
template < typename TImage >
std::auto_ptr<std::vector<numeric_array<int,3> > >
findBorderPoints6(const ConstPtr<TImage> &im,
				  typename TImage::value_type foreground)
{
	Ptr<PointList<int> > borderPoint = new PointList<int>;

	typename Iterator<TImage>::ConstVolumetric iIm(im);
	
	while (!iIm.isAtEnd())
	{
		if (*iIm == foreground)
		{
			if ((iIm.getX() > 0				&& iIm(-1,0,0) != foreground) ||
				(iIm.getX() < im.getX()-1	&& iIm(+1,0,0) != foreground) ||
				(iIm.getY() > 0				&& iIm(0,-1,0) != foreground) ||
				(iIm.getY() < im.getY()-1	&& iIm(0,+1,0) != foreground) ||
				(iIm.getZ() > 0				&& iIm(0,0,-1) != foreground) ||
				(iIm.getZ() < im.getZ()-1	&& iIm(0,0,+1) != foreground))
			{
				borderPoint->push_back(iIm.pos());
			}
		}

		++iIm;
	}
	
	return borderPoint;
}
*/


template < typename TImage >
void oneStepDilation(TImage &im,
                     typename TImage::value_type foreground,
                     typename TImage::value_type background)
{

	int xs = im.dim()[0];
	int ys = im.dim()[1];
	int zs = im.dim()[2];

	int xy = xs*ys;

	typename TImage::value_type color = findValueOtherThan(foreground, background);

	typename Iterator<TImage>::Volumetric iIm(im);

	for (; !iIm.isAtEnd(); ++iIm)
	{
		if (*iIm == background)
		{	
			if (
        // TODO: this has to change to, so that we don't have to write this sort of code
        (containsNeighbor<-1, 0, 0>(iIm) && (iIm.template getUnsafeValue<-1,0,0>() == foreground)) ||
				(containsNeighbor<+1, 0, 0>(iIm) && (iIm.template getUnsafeValue<+1,0,0>() == foreground)) ||
				(containsNeighbor< 0,-1, 0>(iIm) && (iIm.template getUnsafeValue<0,-1,0>() == foreground)) ||
				(containsNeighbor< 0,+1, 0>(iIm) && (iIm.template getUnsafeValue<0,+1,0>() == foreground)) ||
				(containsNeighbor< 0, 0,-1>(iIm) && (iIm.template getUnsafeValue<0,0,-1>() == foreground)) ||
				(containsNeighbor< 0, 0,+1>(iIm) && (iIm.template getUnsafeValue<0,0,+1>() == foreground)))
			{
				*iIm = color;	
			}
		}
	}

	//IF_THEN(im, im[n] == color, im[n] = foreground);
	{
		typename Iterator<TImage>::Linear iIm(im);
		for (; !iIm.isAtEnd(); ++iIm)
		{
			if (*iIm == color) *iIm = foreground;
		}
	}
}



template < typename TImage >
void diamondDilation(TImage &im, int nsteps,
					 typename TImage::value_type foreground,
					 typename TImage::value_type background)
{
	for (int i = 0; i < nsteps; ++i)
		oneStepDilation(im, foreground, background);
}



template < typename TImage >
void diamondErosion(TImage &im, int nsteps,
					typename TImage::value_type foreground,
					typename TImage::value_type background)
{
	for (int i = 0; i < nsteps; ++i)
		oneStepDilation(im, background, foreground);
}



template < typename TImage >
void diamondClosure(TImage &im,
					int nsteps,
					typename TImage::value_type foreground,
					typename TImage::value_type background)
{
	diamondDilation(im, nsteps, foreground, background);
	diamondErosion(im, nsteps, foreground, background);
}

} // namespace


#endif

