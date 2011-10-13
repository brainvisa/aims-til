#ifndef TIL_RECURSIVE_FILTERS_H
#define TIL_RECURSIVE_FILTERS_H

// includes from STL
#include <math.h>

// includes from TIL library
#include "til/til_common.h"
#include "til/convert.h"
#include "til/imageArith.h"
#include "til/miscTools.h"
#include "til/RF4Gaussian.h"
#include "til/rf4Tools.h"


// Namespace 
namespace til {


// Assumption: the sigma entered is the largest sigma possible
// i.e. it corresponds to the axis with the smallest voxel size.
// The motivation is that, e.g. for CT images of different
// voxel size in the Z dimension, at least the smoothing
// on X-Y planes remains the same.
/*
template < typename T >
T modifySigmaAccordingToVoxelSize
(T sigma,
 float vx, float vy, float vz,
 ImageAxis axis)
{
	
	float vMin = min(vx, min(vy, vz));
	
	if (axis==X_AXIS) return T(sigma * vMin/vx);
	if (axis==Y_AXIS) return T(sigma * vMin/vy);
	if (axis==Z_AXIS) return T(sigma * vMin/vz);

	throw std::invalid_argument("Invalid axis");
}
*/
template < typename T >
T modifySigmaAccordingToVoxelSize
(T sigma, const numeric_array<float,3> & vdim, ImageAxis axis)
{
	
	float vMin = min(vdim);
	
	if (axis==X_AXIS) return T(sigma * vMin/vdim[0]);
	if (axis==Y_AXIS) return T(sigma * vMin/vdim[1]);
	if (axis==Z_AXIS) return T(sigma * vMin/vdim[2]);

	throw std::invalid_argument("Invalid axis");
}

// Recursive Gaussian filtering of size sigma
// Warning: pointer im is changed

template < class TImage >
void recursiveGaussianSmoothing(TImage &im, typename TImage::value_type sigma)
{
	typedef typename TImage::value_type T;
	typedef RecursiveFilterSum < RecursiveFilter<double, +1, 4>, RecursiveFilter<double, -1, 4> > RecGaussian;

	T sigmaOriginal = sigma;
	
	for (ImageAxis axis = X_AXIS; axis <= Z_AXIS; ++axis)
	{
		// Apply filter only if image is not flat in that direction
		if (im.dim()[axis] > 1)
		{
			sigma = modifySigmaAccordingToVoxelSize(sigmaOriginal, im.vdim(), axis);
			RecGaussian rf4 = RF4Gaussian(sigma);
			filterAlongAxis(im, im, axis, rf4);
		}
	}
}




// Recursive Gaussian filtering of size sigma

template < typename TImage >
void recursiveGaussianDerivative(TImage &im, typename TImage::value_type sigma,
								 ImageAxis axis)
{	
	typedef typename TImage::value_type T;
	typedef RecursiveFilterSum < RecursiveFilter<double, +1, 4>, RecursiveFilter<double, -1, 4> > RecGaussian;

	T sigmaOriginal = sigma;	

	for (ImageAxis i=X_AXIS; i<=Z_AXIS; ++i)
	{	
		if (i == axis)
		{
			sigma = modifySigmaAccordingToVoxelSize(sigmaOriginal, im.vdim(), i);
			RecGaussian rf4 = RF4GaussianDerivative<T>(sigma, im.vdim()[i]);
			filterAlongAxis(im, im, i, rf4);
		}
		else 
		{
			sigma = modifySigmaAccordingToVoxelSize(sigmaOriginal, im.vdim(), i);
			RecGaussian rf4 = RF4Gaussian(sigma);
			filterAlongAxis(im, im, i, rf4);
		}
	}
}



template < typename TImage >
void recursiveGaussianSecondDerivative(TImage &im, typename TImage::value_type sigma, ImageAxis axis1, ImageAxis axis2)
{
	typedef typename TImage::value_type T;
	typedef RecursiveFilterSum < RecursiveFilter<double, +1, 4>, RecursiveFilter<double, -1, 4> > RecGaussian;

	T sigmaOriginal = sigma;
	
	if (axis1 == axis2)
	{
		for (ImageAxis i=X_AXIS; i<=Z_AXIS; ++i)
		{	
			if (i == axis1)
			{
				sigma = modifySigmaAccordingToVoxelSize(sigmaOriginal, im.vdim(), i);
				RecGaussian rf4 = RF4GaussianSecondDerivative(sigma);
				filterAlongAxis(im, im, i, rf4);
			}
			else 
			{
				sigma = modifySigmaAccordingToVoxelSize(sigmaOriginal, im.vdim(), i);
				RecGaussian rf4 = RF4Gaussian(sigma);
				filterAlongAxis(im, im, i, rf4);
			}
		}
	}
	else
	{
		for (ImageAxis i=X_AXIS; i<=Z_AXIS; ++i)
		{	
			if (i == axis1 || i == axis2)
			{
				sigma = modifySigmaAccordingToVoxelSize(sigmaOriginal, im.vdim(), i);
				RecGaussian rf4 = RF4GaussianDerivative(sigma);
				filterAlongAxis(im, im, i, rf4);
			}
			else 
			{
				sigma = modifySigmaAccordingToVoxelSize(sigmaOriginal, im.vdim(), i);
				RecGaussian rf4 = RF4Gaussian(sigma);
				filterAlongAxis(im, im, i, rf4);
			}
		}
	}
}



// Compute image gradient using recursive Gaussian derivative filtering

template <class TImage1, class TImage2>
void recursiveGaussianGradient(const TImage1 &im,
							   TImage2 &dx,
							   TImage2 &dy,
							   TImage2 &dz,
							   typename TImage2::value_type sigma)
{
	
	convert_im(im, dx);
	copy(dx, dy);
	copy(dx, dz);
	
	recursiveGaussianDerivative(dx, sigma, X_AXIS);
	recursiveGaussianDerivative(dy, sigma, Y_AXIS);
	recursiveGaussianDerivative(dz, sigma, Z_AXIS);

}

template < class TImage1, class TImage2 >
void recursiveGaussianHessian(const TImage1 &im,
							  TImage2 &dxx,
							  TImage2 &dxy,
							  TImage2 &dxz,
							  TImage2 &dyy,
							  TImage2 &dyz,
							  TImage2 &dzz,
							  typename TImage2::value_type sigma)
{
	
	convert_im(im, dxx);
	copy(dxx, dxy);
	copy(dxx, dxz);
	copy(dxx, dyy);
	copy(dxx, dyz);
	copy(dxx, dzz);
	
	recursiveGaussianSecondDerivative(dxx, sigma, X_AXIS, X_AXIS);
	recursiveGaussianSecondDerivative(dxy, sigma, X_AXIS, Y_AXIS);
	recursiveGaussianSecondDerivative(dxz, sigma, X_AXIS, Z_AXIS);
	
	recursiveGaussianSecondDerivative(dyy, sigma, Y_AXIS, Y_AXIS);
	recursiveGaussianSecondDerivative(dyz, sigma, Y_AXIS, Z_AXIS);
	
	recursiveGaussianSecondDerivative(dzz, sigma, Z_AXIS, Z_AXIS);
}


} // namespace

#endif

