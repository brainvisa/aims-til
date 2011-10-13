#ifndef TIL_BASICSTATS_H
#define TIL_BASICSTATS_H

// includes from BOOST
#include "boost/utility/enable_if.hpp"

// includes from TIL library
#include "til/til_common.h"
#include "til/labels.h"

// Namespace 
namespace til {

typedef double t_bsprec;


/// Computes the sum of the intensities of the input image
template < class TImage >
t_bsprec sum(const TImage &im)
{
	t_bsprec res = t_bsprec(0);
	typename Iterator<TImage>::ConstLinear iIm(im);
	
	for (; !iIm.isAtEnd(); ++iIm)
	{
		res += t_bsprec(*iIm);
	}
	return res;
}


/// Returns the mean of the intensities of the input image
template < class TImage >
t_bsprec mean(const TImage &im)
{
	return sum(im) / im.size();
}



/// Returns the estimated variance of the intensities of the input image.
/// (Normalization by (size - 1))
template < class TImage >
t_bsprec varEst(const TImage &im)
{
	t_bsprec meanIm = mean(im);
	t_bsprec res = t_bsprec(0);

	typename Iterator<TImage>::ConstLinear iIm(im);

	for (; !iIm.isAtEnd(); ++iIm)
	{
		res += square(*iIm-meanIm);
	}
	
	return res / (im.size()-1);
}



/// Returns the variance of the intensities of the input image
/// (Normalization by size)
template < class TImage >
typename boost::enable_if<is_Image<TImage>, t_bsprec>::type
var(const TImage &im, t_bsprec meanIm)
{
	t_bsprec res = t_bsprec(0);

	typename Iterator<TImage>::ConstLinear iIm(im);
	
	for (; !iIm.isAtEnd(); ++iIm)
	{
		res += square(*iIm-meanIm);
	}

	return res / im.size();
}


// Returns the variance of the intensities
// (Normalization by size)
template < class TImage >
INLINE 
typename boost::enable_if<is_Image<TImage>, t_bsprec>::type
var(const TImage &im)
{
	return var(im, mean(im));
}



/*
template < class TImage >
typename TImage::value_type max(const ConstPtr<TImage> &im)
{
	typename TImage::ConstLinearIterator iIm(im);
	typename TImage::value_type res = *iIm;

    for (;!iIm.isAtEnd(); ++iIm)
    {
        if (*iIm > res) res = *iIm;   
    }
    return res;
}
*/


/// Returns the maximum intensity of the input image

//(NB: all the following -> solved switching to .NET)

// NB : returns a double instead of a 'typename TImage::value_type'
// because of a VC++6 bug 
// An other work around would be to use a class, since VC can handle
// typenames in return types for class members. However, it is not
// convenient for two reasons:
// (1) no class name overloading
// (2) no template argument deduction, leading to very heavy
//     expressions, like:    maxIm<Image<ushort> >::value(im)
//     (as opposed to a simple max(im)!)

template < class TImage >
typename TImage::value_type max(const TImage &im)
{
	typename Iterator<TImage>::ConstLinear iIm(im);
	typename TImage::value_type res = *iIm;

    for (;!iIm.isAtEnd(); ++iIm)
    {
        if (*iIm > res) res = *iIm;   
    }
    return res;
}


// Returns the minimum intensity of the input image
template < class TImage >
typename boost::enable_if<is_Image<TImage>, typename TImage::value_type>::type
min(const TImage &im)
{
	typename Iterator<TImage>::ConstLinear iIm(im);
	typename TImage::value_type res = *iIm;

    for (; !iIm.isAtEnd(); ++iIm)
    {
        if (*iIm < res) res = *iIm;   
    }
    return res;
}

} // namespace


#endif

