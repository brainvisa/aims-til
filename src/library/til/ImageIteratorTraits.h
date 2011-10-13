#ifndef TIL_IMAGE_ITERATORS_COMMON_H
#define TIL_IMAGE_ITERATORS_COMMON_H

// til include
#include "til_common.h"



namespace til
{


	template < class T >
	struct IteratorOf<ImageC<T> >
	{
		typedef ConstImageCLinearIterator<T>		ConstLinear;
		typedef ConstImageCVolumetricIterator<T>	ConstVolumetric;
		typedef ImageCLinearIterator<T>				Linear;
		typedef ImageCVolumetricIterator<T>			Volumetric;
	};


}

#endif

