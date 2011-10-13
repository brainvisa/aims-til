#ifndef TIL_IMAGEARITH_H
#define TIL_IMAGEARITH_H

// includes from STL
#include <limits>

// includes from BOOST
#include <boost/type_traits.hpp>

// includes from TIL library
#include "til/til_common.h"
#include "til/imageIterator.h"
#include "til/imageTools.h"
#include "til/labels.h"

// Namespace 
namespace til {


// TODO: remove this: use til::fill instead
template < class TImage >
void set(TImage &im, typename TImage::value_type value)
{
    typename Iterator<TImage>::Linear iIm(im);
	for (; !iIm.isAtEnd(); ++iIm)
	{
		*iIm = value;
	}	
}

template < class TImage >
void neg(TImage & im)
{	
  typename Iterator<TImage>::Linear iIm(im);
	for (; !iIm.isAtEnd(); ++iIm)
	{
		*iIm = -*iIm;
	}
}


/// Multiply an image with another object.
/// Note that the type of the multiplier is not necessary of the type of the image.
/// This allows e.g. to multiply an image of vectors with either another vector
/// or a scalar or a matrix. Whichever works with the operator*= on the type of pixels
/// of the image. On the other hand, it allows stupid behavior like multiplying an
/// image of doubles with an integer without casting first the integer.

/*
template < class TImage1, class Unknown >
void mul(const Ptr<TImage1> &im1, Unknown u)
{
	mul(im1, u, int_to_type< IS_SMARTPTR_TYPE(Unknown) >());
}


template < class TImage1, class TImage2 >
void mul(const Ptr<TImage1> &im1, const ConstPtr<TImage2> &im2, int_to_type<true>)
{
	similarityCheck(im1, im2);

	typename Iterator<TImage1>::Linear iIm1(im1);
	typename Iterator<TImage2>::ConstLinear iIm2(im2);

	for (; !iIm1.isAtEnd(); ++iIm1, ++iIm2)
	{
		*iIm1 *= *iIm2;
	}
}

template < class TImage, typename TMultiplier >
void mul(const Ptr<TImage> &im, TMultiplier constant, int_to_type<false>)
{
	typename Iterator<TImage>::Linear iIm(im);
	for (; !iIm.isAtEnd(); ++iIm)
	{
		*iIm *= constant;
	}
}
*/

// NB: I chose not to ignore warnings here, as this information might still be usefull to
// the programmer.
//#pragma warning ( push )
//#pragma warning ( disable : 4244 ) // conversion from X to Y, possible loss of data

namespace detail
{
	template < class TImage1, class TImage2 >
	void
	mul_ii(TImage1 &im1, const TImage2 &im2)
	{
		similarityCheck(im1, im2);

		typename Iterator<TImage1>::Linear iIm1(im1);
		typename Iterator<TImage2>::ConstLinear iIm2(im2);

		for (; !iIm1.isAtEnd(); ++iIm1, ++iIm2)
		{
			*iIm1 *= *iIm2;
		}
	}
}

/// Multiply two images point-wise.
/// Note that we don't require here the precision to be the same for both images, as we
/// don't want to force the user to allocate an image to the appropriate precision.
/// TODO: replace all this crap by call to loop and functors.
template < class TImage1, class TImage2 >
inline
typename enable_if_c<
  is_Image<TImage1>::value &&
  is_Image<TImage2>::value
>::type
mul
(
 TImage1 &im1,			///< The first image, overwritten with the result
 const TImage2 &im2		///< The second image
)
{
	detail::mul_ii(im1, im2);
}


namespace detail{

	template < class TImage, typename TMultiplier >
	void mul_iv(TImage &im, const TMultiplier & constant)
	{
		typename Iterator<TImage>::Linear iIm(im);
		//int c = 0;
		for (; !iIm.isAtEnd(); ++iIm)
		{
			*iIm *= constant;
			/*
			if ((++c % int(im.size()/40.0)) == 0) 
			{
				std::cout << c << std::flush;
				if (c > im.size())
				{
					std::cout << "weird!!" << std::endl;
					return;
				}
			}
			*/
		}
		//std::cout << std::endl;
	}
}

/// Multiply an image with a constant.
/// Note that the constant does not necessarily has the same type as the image. For example,
/// one could imagine to multiply a tensor image with a vector. But the precision of both
/// types have to be the same.
template < class TImage, typename TMultiplier >
typename disable_if_c<
  // secong argument should not be an image
  is_Image<TMultiplier>::value ||
  // the constant should have the same precision as the image
  !boost::is_same<typename precision<TImage>::type, typename precision<TMultiplier>::type>::value
>::type
mul
(
 TImage &im,					///< The image, overwrittent with the result
 const TMultiplier & constant	///< A constant.
)
{
	detail::mul_iv(im, constant);
}



// AND is a reserved word

/*
template < class T >
void and(Ptr<Image<T> > &im1, Ptr<Image<T> > &im2)
{

	similarityCheck(im1, im2);

	T *pIm1 = im1->getPointer();
	const T *pIm2 = im2.getConstPointer();

	for (int n=0; n<im1->size(); ++n)
	{
		if (!(*pIm2)) *pIm1 = 0;
		++pIm1;
		++pIm2;
	}
}
*/


template < class TImage >
void div(TImage & im1, const TImage & im2)
{
    typedef typename TImage::value_type value_type;
    
	const value_type EPSILON = 16*std::numeric_limits<value_type>::epsilon();

	similarityCheck(im1, im2);

	typename Iterator<TImage>::Linear		iIm1(im1);
	typename Iterator<TImage>::ConstLinear	iIm2(im2);

	for (; !iIm1.isAtEnd(); ++iIm1, ++iIm2);
	{
	    if (abs(*iIm2) <= EPSILON)
	    {
			throw std::runtime_error("Division by zero");
	    }
		*iIm1 /= *iIm2;
	}
}


// Deliberately left unimplemented
// use mul(im, 1/scalar) instead

template < class TImage >
void div(TImage &im, typename TImage::value_type scalar);



template < class TImage1, class TImage2 >
typename enable_if_c<
  is_Image<TImage1>::value &&
  is_Image<TImage2>::value
>::type
sub(TImage1 & im1, const TImage2 & im2)
{
	typedef typename TImage1::value_type TPixel1;
	typedef typename TImage2::value_type TPixel2;

	similarityCheck(im1, im2);

	typename Iterator<TImage1>::Linear		iIm1(im1);
	typename Iterator<TImage2>::ConstLinear	iIm2(im2);

	for (; !iIm1.isAtEnd(); ++iIm1, ++iIm2)
	{
		*iIm1 -= castValue<TPixel2, TPixel1>(*iIm2);
	}
}


template < class TImage, typename T >
typename enable_if_c<
 is_Image<TImage>::value &&
 !is_Image<T>::value
>::type
sub(TImage & im, T scalar)
{
    typename Iterator<TImage>::Linear iIm(im);
	for (; !iIm.isAtEnd(); ++iIm)
	{
		*iIm -= scalar;
	}
}

class Add
{
public:
	template < typename T1, typename T2 >
	void operator() (const T1 &a, const T2 &b) const
	{
		add(a, b);
	}
};

template < class TImage1, class TImage2 >
typename enable_if_c<
  is_Image<TImage1>::value &&
  is_Image<TImage2>::value
>::type
add(TImage1 & im1, const TImage2 & im2)
{
	typedef typename TImage1::value_type TPixel1;
	typedef typename TImage2::value_type TPixel2;

	similarityCheck(im1, im2);

	typename Iterator<TImage1>::Linear iIm1(im1);
	typename Iterator<TImage2>::ConstLinear iIm2(im2);

	for (; !iIm1.isAtEnd(); ++iIm1, ++iIm2)
	{
		*iIm1 += castValue<TPixel2, TPixel1>(*iIm2);
	}
}

template < class TImage, typename T >
typename enable_if_c<
 is_Image<TImage>::value &&
 !is_Image<T>::value
>::type
add(TImage & im, T scalar)
{
    typename Iterator<TImage>::Linear iIm(im);
	for (; !iIm.isAtEnd(); ++iIm)
	{
		*iIm += scalar;
	}
}



template < class TImage >
void square(const TImage & in, TImage & out)
{
	similarityCheck(in, out);

	typename Iterator<TImage>::ConstLinear iIn(in);
	typename Iterator<TImage>::Linear iOut(out);

	for (; !iIn.isAtEnd(); ++iIn, ++iOut)
	{
		*iOut = square(*iIn);
	}
}



template < class TImage >
void sqrt(const TImage &in, TImage &out)
{
	similarityCheck(in, out);

	typename Iterator<TImage>::ConstLinear iIn(in);
	typename Iterator<TImage>::Linear iOut(out);
	
	for (; !iIn.isAtEnd(); ++iIn, ++iOut)
	{
		*iOut = sqrt(*iIn);
	}
}

namespace functor
{

} // namespace functor

} // namespace til


#endif

