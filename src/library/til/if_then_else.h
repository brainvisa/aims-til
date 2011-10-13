#ifndef TIL_IF_THEN_ELSE_H
#define TIL_IF_THEN_ELSE_H

// Local includes
#include "til/til_common.h"


///////////////////////////////////////////////////////////////////
//
// if_then_else_old functions, to be used with Pixel Test and
// Pixel Action classes
//
///////////////////////////////////////////////////////////////////


// Namespace 

namespace til {




template < class TImage1, class TImage2, class PixelTest, class PixelAction1, class PixelAction2 >
void if_then_else_old(const TImage1 &in,
				  TImage2 &out,
				  const PixelTest &test,
				  PixelAction1 &actionIf,
				  PixelAction2 &actionElse)
{
	typename Iterator<TImage1>::ConstVolumetric iIn(in);
	typename Iterator<TImage2>::Volumetric iOut(out);
	for (; !iIn.isAtEnd(); ++iIn, ++iOut, test.next(), actionIf.next(), actionElse.next())
	{
		if (test.compute(iIn))
		{
			actionIf.apply(iOut);
		}
		else
		{
			actionElse.apply(iOut);
		}
	}
}


template < class TImage, class PixelTest, class PixelAction1, class PixelAction2 >
void if_then_else_old(TImage &im,
				  const PixelTest &test,
				  PixelAction1 &actionIf,
				  PixelAction2 &actionElse)
{
	typename Iterator<TImage>::Volumetric iIm(im);
	for (; !iIm.isAtEnd(); ++iIm, test.next(), actionIf.next(), actionElse.next())
	{
		if (test.compute(iIm))
		{
			actionIf.apply(iIm);
		}
		else
		{
			actionElse.apply(iIm);
		}
	}
}



template < class TImage1, class TImage2, class PixelTest, class PixelAction >
void if_then(const TImage1 &in,
			 TImage2 &out,
			 const PixelTest &test,
			 PixelAction &action)
{
	typename Iterator<TImage1>::ConstVolumetric iIn(in);
	typename Iterator<TImage2>::Volumetric iOut(out);
	for (; !iIn.isAtEnd(); ++iIn, ++iOut, test.next(), action.next())
	{
		if (test.compute(iIn))
		{
			action.apply(iOut);
		}
	}
}


template < class TImage, class PixelTest, class PixelAction >
void if_then(TImage &im,
             const PixelTest &test,
             PixelAction &action)
{
	typename Iterator<TImage>::Volumetric iIm(im);
	for (; !iIm.isAtEnd(); ++iIm, test.next(), action.next())
	{
		if (test.compute(iIm))
		{
			action.apply(iIm);
		}
	}
}

template < class TImage, class PixelTest, class PixelAction >
void if_then_inInterior(TImage &im,
			 const PixelTest &test,
			 PixelAction &action)
{
	// If image has not interior: do nothing & quit!

	if ((im.dim()[0] < 3) || 
		(im.dim()[1] < 3) || 
		(im.dim()[2] < 3))
	{
		return;
	}


	// The range excluding image border

	Range<int,3> range(1,1,1,im.dim()[0]-2,im.dim()[1]-2,im.dim()[2]-2);

	typename Iterator<TImage>::Volumetric iIm(im, range);
	for (; !iIm.isAtEnd(); ++iIm, test.next(), action.next())
	{
		if (test.compute(iIm))
		{
			action.apply(iIm);
		}
	}
}


} // namespace

#endif

