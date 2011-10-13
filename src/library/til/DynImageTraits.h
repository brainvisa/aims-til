#ifndef TIL_IMAGE_OF_TYPE_H
#define TIL_IMAGE_OF_TYPE_H

// include from TIL library
#include "til/til_common.h"

namespace til
{
	/// Traits for dynamic images.
	/// Its role is to help to know how to construct a type 'TImage'
	/// of pixel type 'TPixel'.
	/// For template image types with a single template argument, 
	/// the default behavior is to return a TImage<TPixel>.
	/// If your image class behaves differently (say, its template argument
	/// is not directly the pixel type), simply specialize this traits class
	/// for you image class to suit your needs.
	template < template <typename> class TImage, typename TPixel >
	struct ImageOfType_OneTemplateArgument
	{
		typedef TImage<TPixel> Type;
	};

	template < template < typename > class TImage, typename TPixel >
	ImageOfType_OneTemplateArgument<TImage, TPixel>
	imageOfType()
	{
		return ImageOfType_OneTemplateArgument<TImage, TPixel>;
	}


} // namespace til

#endif

