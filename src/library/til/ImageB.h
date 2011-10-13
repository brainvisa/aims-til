#ifndef TIL_IMAGEB_H
#define TIL_IMAGEB_H

namespace til
{
	
	/// A 3D image class that stores data in cubes of size NX, NY, NZ
	template < typename _TPixel, int NX, int NY, int NZ >
	class ImageB : public SmartObject, public ImageBase
	{
	};


} // namespace til

#endif

