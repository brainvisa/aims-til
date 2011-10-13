#ifndef TIL_IMAGE_COMMON_H
#define TIL_IMAGE_COMMON_H

// includes from TIL library
#include "til/til_common.h"
#include "til/ImageBase.h"


namespace til {
	
	/// Collects image information to create similar images.
	/// This structure contains parameters that help define similar images
	/// These parameters include dimension, voxel size, origin...
	/// This class is only design to communicate these data from
	/// one image to the other, in particular to instanciate a new
	/// image with identical parameters as an existing one
	/// (cf im1 = new TImage(param(im2)); )

	struct ImageParameter
	{
		/// Image dimensions
    numeric_array<int,3>    m_dim;
		
    /// Voxel size
    numeric_array<t_voxsize,3>     m_vDim;	

		/// Simple constructor
		ImageParameter(numeric_array<int,3> dim, numeric_array<t_voxsize,3> vDim)
			: m_dim(dim), m_vDim(vDim){};
	};

	/// Create an ImageParameter structure out of an Image.
	template < typename TImage >
	ImageParameter param(const TImage &im)
	{
		allocationCheck(im);
		return ImageParameter(im.dim(), im.vdim());
	}

} // namespace


#endif

