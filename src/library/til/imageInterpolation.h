#ifndef TIL_IMAGEINTERPOLATION_H
#define TIL_IMAGEINTERPOLATION_H

// includes from BOOST library
//#include <boost/utility/result_of.hpp>

// include from TIL library
#include "til/til_common.h"
#include "til/til_declarations.h"
#include "til/imageTools.h"
#include "til/PointerTraits.h"
#include "til/recursiveFilters.h"


// Namespace 
namespace til {


template < typename Interpolation, typename Extrapolation, typename Mapping, typename TImageIn, typename TImageOut >
typename enable_if<is_Mapping<Mapping> >::type
resample
(
 const TImageIn &in,	///< The input image
 TImageOut &out,		///< The output image
 const Mapping &map				///< The mapping (from out to in) used for resampling
 )
{
	typedef typename TImageIn::value_type TIn;
	typedef typename TImageOut::value_type TOut;

	notSameImageCheck(in, out);

	// Add extrapolation capability to image
	//ConstPtr<ExtrapolableImage<Extrapolation, TImageIn> > eIm = extrapolable<Extrapolation>(in);

	typename Iterator<TImageOut>::Volumetric iOut(out);
	for (; !iOut.isAtEnd(); ++iOut)
	{
		//*iOut = castValue<TIn, TOut>(Interpolation::template compute<Extrapolation, TImageIn>(in, map(iOut.pos())));
		*iOut = castValue<typename Interpolation::ReturnType, TOut>(Interpolation::template compute<Extrapolation, TImageIn>(in, map(iOut.pos())));
	}

	// TODO: what about the origin??
}


/// Resample an image to match given dimension.
/// If the output image is not allocated, it is allocated to match given dimension.
/// If not, it has to have the same size as the one given in 'dim'. 'Dim' is given
/// as an explicit parameter, and not through 'out', to make it clear that even
/// if out is allocated, things might get modified in there (in particular, 
/// voxel sizes, that are deduced automatically).
// NB: even if dim == in.getDim(), we still resample the image. The reason is that
// in can be a coefficient image (e.g. cubic spline coefficients) that would
// need resampling anyway.
// NB2: YES, we need Extrapolation even to resample an image to another size! Actually
// this is necessary only for supersampling, but I figured out it would be better not to
// split these two cases. TODO: is it really so? Maybe it is better to have an subsampling
// function, that could for example also take care of smothing before sampling (this
// function does not do that).
template < typename Interpolation, typename Extrapolation, typename TImageIn, typename TImageOut >
void
resample
(
 const TImageIn	&in,		///< the input image
 TImageOut				&out,		///< the result of interpolation
 const numeric_array<int,3>			&dim		///< the new dimension
 )
{
	// Check that the input image is allocated
	allocationCheck(in);
	// Check that it has a size > 0 (division by getDim later on)
	if (in.size() <= 0)
	{
		throw std::invalid_argument("in has invalid size");
	}
	// Compute voxel dimension of the output
	numeric_array<t_voxsize,3> vdim = in.vdim() * div<numeric_array<t_voxsize,3> >(dim, in.dim());
	// Chech whether output image is allocated
	if (!isAllocated(out))
	{
		// if not, allocate output
		out.init(dim, vdim);
	}
	else
	{
		// If it is allocated, check that this is consistent with parameter 'dim'
		if (out.dim() != dim)
		{
			throw std::invalid_argument("out and dim arguments are inconsistent");
		}
		// Set new voxel size
		out.set_vdim(vdim);
	}
	// resample image
	resample<Interpolation, Extrapolation>(in, out, getScalingBetween(out, in));
}



/// Resample an image to match a given voxel size
template < typename Interpolation, typename Extrapolation, typename TImageIn, typename TImageOut >
void
resample
(
  const TImageIn	&in,                     ///< the input image
  TImageOut		&out,                        ///< the result of interpolation (reallocated -- no need for initialization)
  const numeric_array<t_voxsize,3> &vDim   ///< the new voxel size to match
)
{
	typedef typename TImageIn::value_type	TIn;
	typedef typename TImageOut::value_type	TOut;

	allocationCheck(in);

	typename Iterator<TImageOut>::Volumetric iOut(out);
	for (; !iOut.isAtEnd(); ++iOut)
	{/*
		*iOut = castValue<TIn, TOut>(
			Interpolation::compute(in,
			iOut.getX()*a.matrix().getXX() + a.transl().getX(),
			iOut.getY()*a.matrix().getYY() + a.transl().getY(),
			iOut.getZ()*a.matrix().getZZ() + a.transl().getZ()));
	*/
	}
}

/*
// NB: Out is reallocated, so there is no need to allocate it first
// TODO: actually I am not sure this is really good... say we want to resample inside
// an existing image?
template < typename Interpolation, typename TImageIn, typename TImageOut>
typename enable_if_c<is_Image<TImageIn>::value && is_Image<TImageOut>::value>::type
resample(const ConstPtr<TImageIn> &in, Ptr<TImageOut> &out, const Vector<int,3> &newDim)
{

	// Check whether input image is allocated
	allocationCheck(in);

	// Check whether new dimension is correct
	if (newDim.getX() < 1 ||
		newDim.getY() < 1 ||
		newDim.getZ() < 1)
	{
		throw std::invalid_argument("Invalid Size");
	}

	{
		Vector<t_voxsize,3> newVDim(EXPAND_VECTOR(in->getVDim()));
		newVDim *= in->getDim();
		//mul(newVDim, in->getDim());
		newVDim /= newDim;
		//div(newVDim, newDim);
		out = new TImageOut(newDim, newVDim);
	}


	// If input and output size are the same, simply convert image
	// NB: this is not good, because in can be a coefficient image
	// (e.g. cubic spline coefficients) and then in and out
	// would be different even though the size is the same
	/ *
	if (in->getDim() == out.getDim())
	{
		convert(in, out);
	}
	* /

	// Compute the affine transform between input and output
	AffineMap<double> a;
	computeAffineTransformBetween(getBoundingBox(out), getBoundingBox(in), a);

	// Reinterpolate image
	//resampleAffine<Interpolation>(in, out, a);
}
*/

/*

// Same thing as resample, exept that input image is smoothed if 
// downsampled
// Take much more memmory, of course.

//template < typename Interpolation, typename Tout = typename Interpolation::value_type >
template < typename Interpolation, typename Tout >
Ptr<Image<Tout> > resampleSmooth(ConstPtr<Image<typename Interpolation::value_type> > &in, Vector<int,3> &newDim)
{
	Ptr<Image<double> > in2 = Image<double>::New(param(in));
	convert(in, in2);

	Ptr<Image<double> > temp = Image<double>::New(param(in));
	Ptr<Image<double> > temp2 = in2;

	double sigma;

	// NB: right now, the cutoff is done with gaussian filtering
	// A more correct version would use a (yet unimplemented) FFT
	// to do the actual bandwith cutoff.

	const double SIGMA_FACTOR = 0.3;

	if (newDim.getX() < in2->getX())
	{
		sigma = sqrt(SQUARE(in2->getX()/newDim.getX())-1)*SIGMA_FACTOR;
		recursiveFilteringAlongSingleAxis(temp2, temp, sigma, 0, 0);
		std::swap(temp, temp2);
	}
	if (newDim.getY() < in2->getY())
	{
		sigma = sqrt(SQUARE(in2->getY()/newDim.getY())-1)*SIGMA_FACTOR;
		recursiveFilteringAlongSingleAxis(temp2, temp, sigma, 1, 0);
		std::swap(temp, temp2);
	}
	if (newDim.getZ() < in2->getZ())
	{
		sigma = sqrt(SQUARE(in2->getZ()/newDim.getZ())-1)*SIGMA_FACTOR;
		recursiveFilteringAlongSingleAxis(temp2, temp, sigma, 2, 0);
		std::swap(temp, temp2);
	}

	return resample<Interpolation, double>(temp2, newDim);
}

*/


} // namespace


#endif

