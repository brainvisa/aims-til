#ifndef TIL_IMAGETOOLS_H
#define TIL_IMAGETOOLS_H

/// \file
/// Small miscellaneous utility functions for images.

// includes from STL
#include <limits>
#include <cmath>
#include <cstdlib>		// rand
#include <typeinfo>		// typeid
#include <vector>

// includes from BOOST
#include "boost/utility/enable_if.hpp"

// includes from TIL library
#include "til/til_common.h"
#include "til/ConditionalIterator.h"
#include "til/image_common.h"
#include "til/labels.h"
#include "til/miscTools.h"
#include "til/numeric_array.h"
#include "til/proba_distributions.h"
#include "til/Range.h"
#include "til/SymMatrix3.h"


// Namespace 
namespace til
{


/// Check whether smart pointer and image are allocated
template < class TImage >
INLINE bool isAllocated(const TImage &im)
{ return im.isAllocated(); }


/// Check whether smart pointer and image are allocated
template < class TImage >
INLINE void allocationCheck(const TImage &im)
{
	if (!isAllocated(im))
	{
		throw std::invalid_argument("Unallocated image");
	}
}

/// Same thing as allocation check, for two images successively.
/// NB: for images that additionnally must have the same size, use similarityCheck
/// instead
template < class TImage1, class TImage2 >
INLINE void allocationCheck(const TImage1 &im1, const TImage2 &im2)
{
	allocationCheck(im1); allocationCheck(im2);
}

template < class TImage1, class TImage2 >
INLINE bool areSimilar(const TImage1 &im1,
					   const TImage2 &im2)
{
	return 
		(im1.dim()[0]  == im2.dim()[0]) &&
		(im1.dim()[1]  == im2.dim()[1]) &&
		(im1.dim()[2]  == im2.dim()[2]) &&
		(im1.vdim()[0] == im2.vdim()[0]) &&
		(im1.vdim()[1] == im2.vdim()[1]) &&
		(im1.vdim()[2] == im2.vdim()[2]);
}

/// Check whether both images are allocated and have the same
/// size and voxel size.
template < class TImage1, class TImage2 >
void similarityCheck(const TImage1 &im1,
					 const TImage2 &im2)
{
	allocationCheck(im1);
	allocationCheck(im2);

	if (!areSimilar(im1, im2))
	{
		throw std::invalid_argument("Incompatible images");
	}
}


/// Check that inputs do not point to the same image
template < class TImage1, class TImage2 >
bool areSameObject(const TImage1 &, const TImage2 &)
{ return false; }

template < class TImage >
typename boost::enable_if<is_Image<TImage>, bool>::type
areSameObject
(
 const TImage &im1,
 const TImage &im2
 )
{ return im1==im2; }


template < class TImage1, class TImage2 >
typename boost::enable_if_c<
 is_Image<TImage1>::value &&
 is_Image<TImage2>::value
>::type
notSameImageCheck(const TImage1 &im1,
					   const TImage2 &im2)
{
	allocationCheck(im1);
	allocationCheck(im2);
	if (areSameObject(im1, im2))
	{
		throw std::invalid_argument("Pointers point to the same image");
	}
}


/// Print basic image information (pixel type, size, voxel size)
template < class TImage >
void printInfo(const TImage &im)
{
	std::cout << "Image Type : " << typeid(typename TImage::value_type).name() << std::endl;
	std::cout << "Image Size : " << im.dim()[0] <<" "<< im.dim()[1] <<" "<< im.dim()[2] << std::endl;
	std::cout << "Image Voxel Size : " << im.vdim()[0] <<" "<< im.vdim()[1] <<" "<< im.vdim()[2] << std::endl;
}



/// Copy one image to another
template < class TImage >
void copy(const TImage &in, TImage &out)
{
	allocationCheck(in);
	allocationCheck(out);
	out.copy(in);
}



/// Creates an image with an ellipsoid centered on the image center.
/// Input:
/// - x/y/z radius (rx, ry, rz)
/// - image voxel size (vx, vy, vz)
/// - image size (sx, sy, sz). If <=0, image size is deduced
///   automatically
/// - foreground and background intensity
template < typename TImage >
void generateEllipsoidImage
(
  TImage &im,
  double rx, double ry, double rz,
  t_voxsize vx, t_voxsize vy, t_voxsize vz,
  int sx, int sy, int sz,
  typename TImage::value_type foreground,
  typename TImage::value_type background
)
{
	
	if (rx <=0 || ry <= 0 || rz <= 0)
	{
		throw std::invalid_argument("Invalid radii");
	}

	// Generate size if it is not provided
	if (sx <= 0)
	{
		sx =2*(int)(ceil(2*rx))+1;
	}
	if (sy <= 0)
	{
		sy =2*(int)(ceil(2*ry))+1;
	}
	if (sz <= 0)
	{
		sz =2*(int)(ceil(2*rz))+1;
	}

	im.init(numeric_array<int,3>(sx, sy, sz), numeric_array<float,3>(vx, vy, vz));


	// Center coordinates

	double cx = (sx-1.0)/2;
	double cy = (sy-1.0)/2;
	double cz = (sz-1.0)/2;

	int i, j, k;
	for (k = 0; k < im.dim()[2]; ++k)
	{
		for (j = 0; j < im.dim()[1]; ++j)
		{
			for (i = 0; i < im.dim()[0]; ++i)
			{
				im(numeric_array<int,3>(i,j,k)) = 
					square(i-cx) / square(rx) + 
					square(j-cy) / square(ry) +
					square(k-cz) / square(rz) <= 1?
					foreground : background;
			}
		}
	}
}




/// Creates an image with a Gaussian around image center.
/// Input:
/// - the x/y/z standard deviation of the Gaussian
/// - the amplitude of the Gaussian
/// - the image voxel siz
/// - the image size (if <=0, image size is chosen automatically)
template < typename TImage >
void generateGaussianImage(TImage &im,
						   numeric_array<double,3> sigma,
						   typename TImage::value_type amplitude,
						   numeric_array<t_voxsize,3> voxSize,
						   numeric_array<int,3> size)
{
	typedef typename TImage::value_type value_type;

	const double SIGMA_MIN = 0.0; // 1e-3;
	
	if (sigma[0] < SIGMA_MIN || sigma[1] < SIGMA_MIN || sigma[2] < SIGMA_MIN)
	{
		throw std::invalid_argument("Standard deviation is too small");
	}

	// Generate size if it is not provided

	if (size[0] <= 0)
	{
		size[0] = 2*max(1,castValue<double, int>(3.0*ceil((double)sigma[0]))) + 1;
	}
	if (size[1] <= 0)
	{
		size[1] = 2*max(1,castValue<double, int>(3.0*ceil((double)sigma[1]))) + 1;
	}
	if (size[2] <= 0)
	{
		size[2] = 2*max(1,castValue<double, int>(3.0*ceil((double)sigma[2]))) + 1;
	}

	im.init(size, voxSize);


	// Center coordinates
	// For a 3x3 image I want the center at (1,1)
	// For a 2x2 image I want the center at (0.5, 0.5)
	// Hence the '-1'

	//numeric_array<double,3> center(size-1);
	numeric_array<double,3> center;
  convert(center, (size-1));
  center *= 0.5;

	int i, j, k;
	for (k = 0; k < im.dim()[2]; ++k)
	{
		for (j = 0; j < im.dim()[1]; ++j)
		{
			for (i = 0; i < im.dim()[0]; ++i)
			{
				im(numeric_array<int,3>(i,j,k)) = value_type( amplitude * exp(-
					square(i-center[0]) / (2*square(sigma[0])) - 
					square(j-center[1]) / (2*square(sigma[1])) -
					square(k-center[2]) / (2*square(sigma[2]))));
			}
		}
	}
}



template < typename TImage >
void generateGaussianKernel(TImage &kernel,
							double sigmaX, double sigmaY, double sigmaZ,
							t_voxsize vx, t_voxsize vy, t_voxsize vz,
							int sx, int sy, int sz)
{
	generateGaussianImage(kernel, numeric_array<double,3>(sigmaX, sigmaY, sigmaZ),
		1.0, numeric_array<t_voxsize,3>(vx, vy, vz), numeric_array<int,3>(sx, sy, sz));

	setSumToOne(kernel);
}



template <typename T1, typename TImage2>
void deduceImageSizeFromGaussianStandardDeviation(T1 sigma, const TImage2 &im, int &hx, int &hy, int &hz)
{
	T1 sigmaX, sigmaY, sigmaZ;

	sigmaX = modifySigmaAccordingToVoxelSize(sigma, im.vdim()[0], im.vdim()[1], im.vdim()[2], 0);
	sigmaY = modifySigmaAccordingToVoxelSize(sigma, im.vdim()[0], im.vdim()[1], im.vdim()[2], 1);
	sigmaZ = modifySigmaAccordingToVoxelSize(sigma, im.vdim()[0], im.vdim()[1], im.vdim()[2], 2);

	hx = castValue<double, int>(3.0*ceil(sigmaX));
	hy = castValue<double, int>(3.0*ceil(sigmaY));
	hz = castValue<double, int>(3.0*ceil(sigmaZ));
}



template < typename TImage >
void getLocalHessian(const TImage &im, SymMatrix3<typename TImage::value_type> &mat, int x, int y, int z)
{
	if (!contains(getRange(im),
                      Range<int,3>(numeric_array<int, 3>(x-1, y-1, z-1),
                                   numeric_array<int, 3>(x+1, y+1, z+1))))
	{
		throw std::invalid_argument("Point is outside image or on image border");
	}

	mat.set(0,0, (im.getUnsafeValue(x+1,y,z) + im.getUnsafeValue(x-1,y,z) - 2*im.getUnsafeValue(x,y,z)) / square(im.vdim()[0]));
	mat.set(1,1, (im.getUnsafeValue(x,y+1,z) + im.getUnsafeValue(x,y-1,z) - 2*im.getUnsafeValue(x,y,z)) / square(im.vdim()[1]));
	mat.set(2,2, (im.getUnsafeValue(x,y,z+1) + im.getUnsafeValue(x,y,z-1) - 2*im.getUnsafeValue(x,y,z)) / square(im.vdim()[2]));

	mat.set(0,1, (im.getUnsafeValue(x+1,y+1,z) + im.getUnsafeValue(x-1,y-1,z) - im.getUnsafeValue(x+1,y-1,z) - im.getUnsafeValue(x-1,y+1,z)) / ( im.vdim()[0] * im.vdim()[1] ));
	mat.set(0,2, (im.getUnsafeValue(x+1,y,z+1) + im.getUnsafeValue(x-1,y,z-1) - im.getUnsafeValue(x+1,y,z-1) - im.getUnsafeValue(x-1,y,z+1)) / ( im.vdim()[0] * im.vdim()[2] ));
	mat.set(1,2, (im.getUnsafeValue(x,y+1,z+1) + im.getUnsafeValue(z,y-1,z-1) - im.getUnsafeValue(x,y+1,z-1) - im.getUnsafeValue(x,y-1,z+1)) / ( im.vdim()[1] * im.vdim()[2] ));
}


/// Get coordinates of image center
template < typename TImage >
inline numeric_array<int,3> getCenter(TImage &im)
{
  return im.dim() / 2;
	//return numeric_array<int,3>(div<Point<int,3> >(im.dim(),2));// / 2);
  /*im.dim()[0] / 2,
		im.dim()[1] / 2,
		im.dim()[2] / 2);*/
}

/// Get coordinates of image center
template < typename TImage >
void getImageCenter(const TImage &im, int &cx, int &cy, int &cz)
{
	if (!(im.dim()[0]%2 && im.dim()[1]%2 && im.dim()[2]%2))
	{
		throw std::invalid_argument("Even image size");
	}

	cx = im.dim()[0] / 2;
	cy = im.dim()[1] / 2;
	cz = im.dim()[2] / 2;
}



/*
/// Check whether image range contains position (i,j,k)
template < typename TImage >
INLINE
typename boost::enable_if<is_Image<TImage>, bool>::type
contains(const TImage &im, int i, int j, int k)
{
	return (i >= 0 && j >= 0 && k >= 0 && 
		i < im.dim()[0] && 
		j < im.dim()[1] &&
		k < im.dim()[2]);	
}
*/

/// Check whether image range contains position 'pos'
template < typename TImage >
INLINE
typename boost::enable_if<is_Image<TImage>, bool>::type
contains(const TImage & im, const numeric_array<int,3> & pos)
{
	return all_greater_equal(pos, 0) && all_less(pos, im.dim());
}


/// Tests whether image range contains position 'pos + offset', where 'pos'
/// is known to lie within the image range.
/// The idea here is that the addition is not computed for the components
/// where offset is zero. When this happens frequently, e.g. for neighborhood
/// computations, this speeds up things considerably.
template < typename TImage >
INLINE
bool contains(const TImage & im, const numeric_array<int,3> & pos, const numeric_array<int,3> & offset)
{
	/*
	// NB: this code did not make any difference in term of speed. This kind
	// of optimizations seems to be carried out well by the compiler.
	int temp;
	return
		( 
		( offset.getX() == 0  || ((temp = pos.getX() + offset.getX()) >= 0 && temp < im.dim()[0])) &&
		( offset.getY() == 0  || ((temp = pos.getY() + offset.getY()) >= 0 && temp < im.dim()[1])) &&
		( offset.getZ() == 0  || ((temp = pos.getZ() + offset.getZ()) >= 0 && temp < im.dim()[2]))
		);*/

	return
		( 
		( offset[0] == 0  || (pos[0] + offset[0] >= 0 && pos[0] + offset[0] < im.dim()[0])) &&
		( offset[1] == 0  || (pos[1] + offset[1] >= 0 && pos[1] + offset[1] < im.dim()[1])) &&
		( offset[2] == 0  || (pos[2] + offset[2] >= 0 && pos[2] + offset[2] < im.dim()[2]))
		);
}

/// Get image range
template < typename TImage >
Range<int,3> getRange(const TImage &im)
{
	return Range<int,3>(im.dim() - 1);
}


/// Get image bounding box
template < typename TImage >
Box<double,3> getBoundingBox(const TImage &im)
{
	return Box<double,3>(numeric_array<double,3>(-0.5,-0.5,-0.5), im.dim() - 0.5);
}



/// Extract the line starting from (x,y,z) in the direction of axis
/// from im to bar. bar should be allocated (of length im.getDim(axis))
// TODO: these functions are so simple, they just simply shoudn't exist. Use loop.
/*
template < typename TImage >
void extractBar(const TImage &im, int x, int y, int z, ImageAxis axis, std::vector<typename TImage::value_type> &bar)
{
	typename Iterator<TImage>::ConstVolumetric iIm(im);
	typename std::vector<typename TImage::value_type>::iterator iBar = bar.begin();
	for (iIm.setPos(x,y,z); iBar != bar.end(); iIm.next(axis), ++iBar)
	{
		*iBar = *iIm;
	}
}*/
template < typename TImage >
void extractBar(const TImage &im, const numeric_array<int,3> & pos, ImageAxis axis, std::vector<typename TImage::value_type> &bar)
{
	typename Iterator<TImage>::ConstVolumetric iIm(im);
	typename std::vector<typename TImage::value_type>::iterator iBar = bar.begin();
	for (iIm.set_pos(pos); iBar != bar.end(); iIm.next(axis), ++iBar)
	{
		*iBar = *iIm;
	}
}



/// Insert the line starting from (x,y,z) in the direction of axis
/// from bar to im.
/*
template < typename TImage >
void insertBar(TImage &im, int x, int y, int z, ImageAxis axis, const std::vector<typename TImage::value_type> &bar)
{
	typename Iterator<TImage>::Volumetric iIm(im);
	typename std::vector<typename TImage::value_type>::const_iterator iBar = bar.begin();
	for (iIm.setPos(x,y,z); iBar != bar.end(); iIm.next(axis), ++iBar)
	{
		*iIm = *iBar;
	}
}
*/
template < typename TImage >
void insertBar(TImage &im, const numeric_array<int,3> & pos, ImageAxis axis, const std::vector<typename TImage::value_type> &bar)
{
	typename Iterator<TImage>::Volumetric iIm(im);
	typename std::vector<typename TImage::value_type>::const_iterator iBar = bar.begin();
	for (iIm.set_pos(pos); iBar != bar.end(); iIm.next(axis), ++iBar)
	{
		*iIm = *iBar;
	}
}

/*
template < typename TImage >
Ptr<TImage> similarAs(Ptr<TImage> &im)
{
	Ptr<TImage> res = new TImage(param(im));
	return res;
}
*/


/// Normalize the image so that the sum of its element is equal to one
// TODO: clearly the kind of thing that should be easily done in a few lines
// without using a function: how?
template < typename TImage >
void setSumToOne(TImage &im)
{
	typedef typename TImage::value_type value_type;

	if (std::numeric_limits<value_type>::is_integer)
	{
		throw std::domain_error("Defined for floating numbers only");
	}

	value_type nrm = sum(im);

	const value_type EPSILON = 128 * std::numeric_limits<value_type>::epsilon();

	if (std::abs(nrm) <= EPSILON)
	{
		throw std::underflow_error("kernel norm too small");
	}

	mul(im, value_type(1.0/nrm));
}



/// Converts a contiguous memory buffer into a per-slice buffer
template < typename T >
T** simple2DoublePointer(T *data, int x, int y, int z)
{
	if ( x<0 || y<0 || z<0)
	{
		throw std::invalid_argument("Image size < 0");
	}
	
	T** res = new T*[z];
	
	for (int i=0; i < z; ++i)
	{
		res[i] = data + i * x * y;
	}
	
	return res;
}


template < typename TImage >
typename Iterator<TImage>::Linear
itlin(TImage &im)
{
	return typename Iterator<TImage>::Linear(im);
}
// NB: BoolFunctor is the first template parameter so that it can be given
// explicitely without having to give also the image type
template < typename BoolFunctor, typename TImage >
ConditionalIterator<typename Iterator<TImage>::Linear, BoolFunctor>
itlin(TImage &im, const BoolFunctor &boolFunctor = BoolFunctor())
{
	return ConditionalIterator<typename Iterator<TImage>::Linear, BoolFunctor>(im, boolFunctor);
}

template < typename TImage >
typename Iterator<TImage>::ConstLinear
itlin_const(const TImage &im)
{
	return typename Iterator<TImage>::ConstLinear(im);
}

template < typename BoolFunctor, typename TImage >
ConditionalIterator<typename Iterator<TImage>::ConstLinear, BoolFunctor>
itlin_const(const TImage &im, const BoolFunctor &boolFunctor = BoolFunctor())
{
	return ConditionalIterator<typename Iterator<TImage>::ConstLinear, BoolFunctor>(im, boolFunctor);
}

template < typename TImage >
typename Iterator<TImage>::Volumetric
itvol(TImage &im)
{
	return typename Iterator<TImage>::Volumetric(im);
}

template < typename BoolFunctor, typename TImage >
ConditionalIterator<typename Iterator<TImage>::Volumetric, BoolFunctor>
itvol(TImage &im, const BoolFunctor &boolFunctor = BoolFunctor())
{
	return ConditionalIterator<typename Iterator<TImage>::Volumetric, BoolFunctor>(im, boolFunctor);
}

template < typename TImage >
typename Iterator<TImage>::ConstVolumetric
itvol_const(const TImage &im)
{
	return typename Iterator<TImage>::ConstVolumetric(im);
}

template < typename BoolFunctor, typename TImage >
ConditionalIterator<typename Iterator<TImage>::ConstVolumetric, BoolFunctor>
itvol_const(const TImage &im, const BoolFunctor &boolFunctor = BoolFunctor())
{
	return ConditionalIterator<typename Iterator<TImage>::ConstVolumetric, BoolFunctor>(im, boolFunctor);
}


/// Fills the image with random numbers uniformly distributed
/// on [min, max]
template < typename TImage >
void rand
(
  TImage & im,
  typename TImage::value_type min,
  typename TImage::value_type max
)
{
	typedef typename TImage::value_type value_type;
	
	// Checking whether max < min
	if (max < min)
	{
		warning("til::rand: input max < min, setting max := min");
		max = min;
	}
	
	typename Iterator<TImage>::Linear iIm(im);
	for (; !iIm.isAtEnd(); ++iIm)
	{
		*iIm = rand(min, max);
	}
}

/// Set a value to an image slice
// TODO: clearly a task that should be easy to do in a few lines without
// having to write a function: how?
template < typename TImage >
void setSlice(TImage &im,                         ///< Image to be modified
              int sliceNumber,                    ///< Number of the slice
              ImageAxis axis,                     ///< Orientation of slice normal
              typename TImage::value_type value   ///< New value
              )
{
	Range<int,3> range = getRange(im);
  // collapse range on current axis
  range.set_bounds(axis, sliceNumber, sliceNumber);

	typename Iterator<TImage>::Volumetric iIm(im, range);
	for (; !iIm.isAtEnd(); ++iIm)
	{
		*iIm = value;
	}
}


/// Clear image borders
template < typename TImage >
void clearBorder
(
  TImage &im,
  typename TImage::value_type background = 0
)
{
	typename Iterator<TImage>::Volumetric iIm(im);
	
	for (; !iIm.isAtEnd(); ++iIm)
	{
		if (
      (iIm.pos()[0] == 0) ||
			(iIm.pos()[1] == 0) ||
			(iIm.pos()[2] == 0) ||
			(iIm.pos()[0] == im.dim()[0] - 1) ||
			(iIm.pos()[1] == im.dim()[1] - 1) ||
			(iIm.pos()[2] == im.dim()[2] - 1))
    {
			*iIm = background;
		}
	}
}


template <int dx, int dy, int dz, class VolumetricImageIterator>
bool containsNeighbor(const VolumetricImageIterator & iIm)
{
	return (
		(dx >= 0 || iIm.pos()[0] > 1+dx) &&
		(dx <= 0 || iIm.pos()[0] < iIm.image().dim()[0]-dx) &&

		(dy >= 0 || iIm.pos()[1] > 1+dy) &&
		(dy <= 0 || iIm.pos()[1] < iIm.image().dim()[1]-dy) &&

		(dz >= 0 || iIm.pos()[2] > 1+dz) &&
		(dz <= 0 || iIm.pos()[2] < iIm.image().dim()[2]-dz));
}

/// Get the world coordinates of an image point
template < typename TImage, typename T >
INLINE
typename boost::enable_if< is_Image<TImage>, numeric_array<double,3> >::type
vol2world(const TImage & im, const numeric_array<T,3> & volumePos)
{
	//return (volumePos - im.getOrigin()) * im.getVDim();
	//return volumePos * im.vdim();
  return mul<numeric_array<double,3> >(volumePos, im.vdim());
}

template < typename TImage >
INLINE
typename boost::enable_if< is_Image<TImage>, numeric_array<double,3> >::type
world2vol(const TImage & im, const numeric_array<double,3> & worldPos)
{
  // En fait c'est pas si idiot que ca de passer par le data, parce que ca nous
  // indique qu'effectivement on fait qqch de non-standard, et c'est vrai:
  // c'est une approx rapide de l'application d'une matrice que l'on fait.
  return div<numeric_array<double,3> >(worldPos, im.vdim());
	//return worldPos / im.vdim();// + im.getOrigin();
}


} // namespace

#endif

