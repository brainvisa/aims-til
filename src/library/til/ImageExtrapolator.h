#ifndef TIL_IMAGEEXTRAPOLATOR_H
#define TIL_IMAGEEXTRAPOLATOR_H

// includes from TIL library
#include "til/til_common.h"
#include "til/imageTools.h"
#include "til/labels.h"
#include "til/PointerTraits.h"


namespace til {

// TODO: actually maybe we want to have the test (im.contains) in these classes,
// first because then the user does not have to type it, and second, it makes
// possible to have an 'unsafe' extrapolation that actually does not extrapolate.
// Might be useful to use some functions (e.g. resample) where we know that
// with the parameters we give there is no extrapolation to be done.

// A macro to shorten declaration of extrapolation classes
// NB: to be undefined at the end of the declarations
/*
#define EXTRAPOL_FUNC                                           \
template < typename TImage >                                    \
INLINE                                                          \
static typename TImage::value_type                              \
getValue(const TImage &im, const Vector<int,3> & pos)           \
{ return Self::getValue(im, EXPAND_VECTOR(pos)); }              \
                                                                \
template < typename TImage >                                    \
INLINE                                                          \
static typename TImage::value_type                              \
getValue(const TImage &im, int i, int j, int k)                 \
{ if (contains(im,i,j,k)) { return im.getUnsafeValue(i,j,k); }  \
  else { return Self::getExtrapolatedValue(im, i, j, k); } }		\
                                                                \
template < typename TImage >                                    \
INLINE                                                          \
static typename TImage::value_type                              \
getExtrapolatedValue(const TImage &im, const Vector<int,3> & pos)	\
{ return Self::getExtrapolatedValue(im, EXPAND_VECTOR(pos)); }		\
																	\
template < typename TImage >										\
INLINE																\
static typename TImage::value_type										\
getExtrapolatedValue(const TImage &im, int i, int j, int k)			\
																	\
*/

#define TIL_EXTRAPOL_FUNC(imname, posname)                          \
template < typename TImage >                                        \
INLINE                                                              \
static typename TImage::value_type                                  \
getValue(const TImage &im, numeric_array<int,3> pos)                        \
{ if (contains(im,pos)) { return im.getUnsafeValue(pos); }          \
  else { return Self::getExtrapolatedValue(im, pos); } }		        \
                                                                    \
template < typename TImage >                                        \
INLINE                                                              \
static typename TImage::value_type                                  \
getExtrapolatedValue(const TImage & imname, numeric_array<int,3> posname)   \




// NB: In all extrapolator classes, it is assumed that the (i,j,k) position
// lies outside the image range
// In other words, this classes are always called after an explicit range
// checking.


/// No image extrapolation
struct RangeCheckingExtrapolator : public ImageExtrapolator_label
{
public: // typededs

	typedef RangeCheckingExtrapolator Self;

public: // static functions

	TIL_EXTRAPOL_FUNC( , pos)
	{
		//throw std::out_of_range(mySprintf("Coordinates out of image range: (%i, %i, %i)", i, j, k));
    throw std::out_of_range(mySprintf("Coordinates out of image range: (%i, %i, %i)", pos[0], pos[1], pos[2]));	}
};


/// Extrapolate image with default value of image type (zero for numeric types)
class ZeroExtrapolator : public ImageExtrapolator_label
{
public: // typedefs

	typedef ZeroExtrapolator Self;

public: // static functions

	TIL_EXTRAPOL_FUNC( , )
	{
		return typename TImage::value_type();
	}
};


/// Image extrapolation using value of nearest image point
class NearestExtrapolator : public ImageExtrapolator_label
{
public: // typededs

	typedef NearestExtrapolator Self;

public: // static functions

	TIL_EXTRAPOL_FUNC(im, pos)
	{
		if (pos[0]<0) pos[0]=0; else if (pos[0] >= im.dim()[0]) pos[0] = im.dim()[0]-1;
		if (pos[1]<0) pos[1]=0; else if (pos[1] >= im.dim()[1]) pos[1] = im.dim()[1]-1; 
		if (pos[2]<0) pos[2]=0; else if (pos[2] >= im.dim()[2]) pos[2] = im.dim()[2]-1;
		/*
    if (i<0) i=0; else if (i >= im.dim()[0]) i = im.dim()[0]-1;
		if (j<0) j=0; else if (j >= im.dim()[1]) j = im.dim()[1]-1; 
		if (k<0) k=0; else if (k >= im.dim()[2]) k = im.dim()[2]-1;*/
		return im.getUnsafeValue(pos);
	}
};


/// Extrapolation using (thin) mirror symmetry.
/// The thin means, the value at the boundary is not repeated, i.e.
/// the value at position -1 is equal to the value at position 1, not 0.
class MirrorExtrapolator : public ImageExtrapolator_label
{
public: // typededs

	typedef MirrorExtrapolator Self;

public: // static functions

	TIL_EXTRAPOL_FUNC(im, pos)
	{
		Self::mirror(pos[0], im.dim()[0]-1);
		Self::mirror(pos[1], im.dim()[1]-1);
		Self::mirror(pos[2], im.dim()[2]-1);
		return im.getUnsafeValue(pos);
	}

private: // static functions

	static void mirror(int &i, int n)
	{
		i = std::abs(i) % (2*n);
		if (i > n) i = 2*n - i;
	}
};



/// Image extrapolation using cyclic symmetry.
class CyclicExtrapolator : public ImageExtrapolator_label
{
public: // typededs

	typedef CyclicExtrapolator Self;

public: // static functions

	TIL_EXTRAPOL_FUNC(im, pos)
	{
		Self::cyclic(pos[0], im.dim()[0]);
		Self::cyclic(pos[1], im.dim()[1]);
		Self::cyclic(pos[2], im.dim()[2]);
		return im.getUnsafeValue(pos);
	}

private: // static functions

	static void cyclic(int &i, int n)
	{
		i = i%n;
		if (i<0) i+=n;
	}
};

/// This extrapolator actually does NOT extrapolate at all and does
/// not do any range checking.
/// It can be used for efficiency reasons in algorithms that would normally 
/// require extrapolation (e.g. resampling or sub-image extraction), but where 
/// we know that, given our input parameters, no extrapolation will ever be needed.
class UnsafeNoExtrapolator : public ImageExtrapolator_label
{
public: // typededs

	typedef CyclicExtrapolator Self;

public: // static functions

	template < typename TImage >
	inline
	static typename TImage::value_type
	getValue(const TImage & im, const numeric_array<int,3> & pos)
	{ return im.getUnsafeValue(pos); }

	template < typename TImage >										
	inline
	static typename TImage::value_type
	getExtrapolatedValue(const TImage &im, const numeric_array<int,3> & pos)
  { return im.getUnsafeValue(pos); }
};


#undef TIL_EXTRAPOL_FUNC

/// Enable image extrapolation by overloading the getValue function.
// still not clear whether it should inherit from image or not.
// There are many possibilities to get the extrapolated value of an image.
// One could for example code in each image class a getValue<Extrapolator> member.
// But then, do we bring really something as opposed to using directly the extrapolation
// classes? Not really, because in both cases we assume that the current code and
// the person who writes it are aware that at this point in time an Extrapolation
// is needed and used.
// This class aim at hiding the extrapolation process, by templating on it, so
// it brings another layer. It means that we can pass along this object to a function,
// and this function will do extrapolation without even knowing it.
// Initially this was motivated for image resampling. When we do image resampling
// we have to do interpolation of potentially extrapolated value. But of course,
// we should not write an interpolation that is aware of extrapolation, because both
// problems are distinct. The easiest way would be to proceed in two steps, first
// extract potentially extrapolated values, and then interpolate them, but the
// extraction of values actually depends on the interpolator, so this job has to be
// done by it.
// TODO: maybe a mechanism to inhibit from templating on itself
template < class TImageExtrapolator, class TImage >
class ExtrapolableImage : public TImage, public ExtrapolableImage_label
{
public: // typedefs
	typedef typename TImage::value_type value_type;

public: // functions

	value_type getValue(const numeric_array<int,3> & p) const
	{
	  if (this->contains(p))
	  {		
      return this->getUnsafeValue(p);
		}
		else
		{
			return TImageExtrapolator::getExtrapolatedValue(*this, p);
	  }
	}

private: // rules

	typedef typename enable_if<is_Image<TImage> >::type rule1;
	typedef typename disable_if<is_ExtrapolableImage<TImage> >::type rule2;
};

template < class TImageExtrapolator, class TImage >
ExtrapolableImage<TImageExtrapolator, TImage>
extrapolable(const TImage &im)
{
	typedef ExtrapolableImage<TImageExtrapolator, TImage> EImage;
	EImage eim;
	eim.shallowCopy(im);
	return eim;
}
/*
template < class TImageExtrapolator, class TImage >
ConstPtr<ExtrapolableImage<TImageExtrapolator, TImage> >
extrapolable(const ConstPtr<TImage> &im)
{
	typedef ExtrapolableImage<TImageExtrapolator, TImage> EImage;
	return ConstPtr<EImage>(static_cast<const EImage *>((const TImage *)im));
}
*/


} // namespace

#endif

