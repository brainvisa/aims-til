#ifndef TIL_DISTANCE_TRANSFORM_H
#define TIL_DISTANCE_TRANSFORM_H

// Local includes
#include "til/til_common.h"

namespace til
{

	// NB: both foreground and background are needed, since it might be that
	// other objects are present inside the image, inside which the distance
	// is not computed. However they are important, since the geodesics
	// goes around these objects.

	template < typename TImage >
	void
	cityBlock(TImage &seg,				///< The object image, to be overwritten
	typename TImage::value_type foreground,	///< Object label
	typename TImage::value_type background  ///< Background label
	)
	{
		{
		typename TImage::LinearIterator iSeg(seg);
		
		const TImage::value_type INFINITY		= std::numerical_limits<T>::max() - 1;
		const TImage::value_type OTHER_OBJECT	= std::numerical_limits<T>::max();

		// Initialize seg for distance map
		for (; !iSeg.isAtEnd(); ++iSeg)
		{
			// Object has a distance of zero
			if (*iSeg == foreground)		*iSeg = 0;
			// Background has an initial distance of as close as infinity we can get
			else if (*iSeg == background)	*iSeg = INFINITY;
			// Other objects have a special value
			else							*iSeg = OTHER_OBJECT;
		}
		

		{			
			// forward loop
			typename TImage::LinearIterator iSeg(seg);
			for(; !iSeg.isAtEnd(); ++iSeg)
			{
				if (*iSeg == OTHER_OBJECT) continue;
				*iSeg = min(*iSeg, min(min(min(iSeg.getValue<-1,0,0>(), iSeg.getValue<0,-1,0>()), iSeg.getValue<0,0,-1>) + T(1), INFINITY)));
			}
		}
		{			
			// backward loop
			typename TImage::LinearIterator iSeg(seg);
			for(; !iSeg.isAtEnd(); --iSeg)
			{
				if (*iSeg == OTHER_OBJECT) continue;
				*iSeg = min(*iSeg, min(min(min(iSeg.getValue<+1,0,0>(), iSeg.getValue<0,+1,0>()), iSeg.getValue<0,0,+1>) + T(1), INFINITY)));
			}
		}
	}


} // namespace

#endif

