#ifndef TIL_CONVERT_H
#define TIL_CONVERT_H

// includes from STL
#include <limits>

// includes from TIL library
#include "til/til_common.h"
#include "til/imageBasicStats.h"
#include "til/miscTools.h"

// Namespace 
namespace til
{
  
  /// Convert an image into the type of the other image
  /// while scaling its intensity range to span through the
  /// whole possible range of TOut.
  /// Typically to convert floats to integers
  /// out image should have been previously allocated
  // NB: The fact that out has to be allocated is not necessary. We 
  // could check that out is allocated, and allocate it if needed. 
  // However, as a general philosophy, I think it is better that 
  // no allocation is done on images without the user knowing it,
  // because memory size can be huge, and the user has to be in
  // control and aware of everything on this level.
  template < class TImageIn, class TImageOut >
  void convertScale
  (
   TImageIn const & in,  ///< the input image
   TImageOut      & out  ///< the output image, previously allocated
   )
  {
  	typedef typename TImageIn::value_type TPixelIn;
  	typedef typename TImageOut::value_type TPixelOut;
  
  	const TPixelIn EPSILON = 32 * std::numeric_limits<TPixelIn>::epsilon();
  
  	// Check that images have same size
  	similarityCheck(in, out);
  
      // Test if in has too little variation to be rescaled
      // If so, output image is set to zero
  	double maxIm = (double) max(in);
  	double minIm = (double) min(in);
  	if (std::abs(maxIm - minIm) <= EPSILON)
  	{
  		out.reset();
  		return;
  	}
  
      // Compute the parameters of the affine rescaling of the intensity
  	// TODO: choose automatically type of coefficients a and b
  	// using traits (e.g. double is not enough is TIn is of type long double
  	TPixelIn t1max = std::numeric_limits<TPixelIn>::max();
  	TPixelIn t1min = std::numeric_limits<TPixelIn>::min();
  	
  	double a = (t1max - t1min) / (maxIm - minIm);
  	double b = (t1min * maxIm - t1max*minIm) / (maxIm - minIm);
  
  	// Apply the intensity affine transform from in to out
  	typename Iterator<TImageIn>::ConstLinear iIn(in);
  	typename Iterator<TImageOut>::Linear iOut(out);
  	do
  	{
  		*iOut = castValue<TPixelIn,TPixelOut>(*iIn*a + b);
  	}
  	while (iIn.next() && iOut.next());
  }
  
  
  
  /// Convert an image into the type of the other image.
  /// out image should have been previously allocated
  // TODO: I renammed this into convert_im because of a conflict with a general convert,
  // that will eventually take over when the Grand Plan is finished.
  template < class TImageIn, class TImageOut >
  void convert_im
  (
   TImageIn const & in, ///< The input image
   TImageOut      & out ///< The output image, previously allocated
   )
  {
  	typedef typename TImageIn::value_type TPixelIn;
  	typedef typename TImageOut::value_type TPixelOut;
  
  	// Check that image have same size
  	similarityCheck(in, out);
  
  	// If the maximum value possible of new type is less than the
  	// maximum value of old type...
  	if (std::numeric_limits<TPixelOut>::max() < std::numeric_limits<TPixelIn>::max())
  	{
  		// ...check whether maximum input image value is greater than new type max
  		if (max(in) > std::numeric_limits<TPixelOut>::max())
  		{
  			// Throw an exception: The user will have to think of an
  			// alternative
  			throw std::overflow_error("Output type cannot handle input image range");
  		}
  	}
  
  	// Same thing with minimums:
  	// If the minimum value possible of new type is greater than the
  	// minimum value of old type...
  	if (type_min<TPixelOut>() > type_min<TPixelIn>())
  	{
  		// ...check whether minimum image value is lesser than new type min
  		if (min(in) < type_min<TPixelOut>())
  		{
  			// Throw an exception: The user will have to think of an
  			// alternative
  			throw std::overflow_error("Output type cannot handle input image range");
  		}
  	}
  
  	// Do the conversion
  	typename Iterator<TImageIn>::ConstLinear iIn(in);
  	typename Iterator<TImageOut>::Linear iOut(out);
  	for (; !iIn.isAtEnd(); ++iIn, ++iOut)
  	{
  		*iOut = castValue<TPixelIn,TPixelOut>(*iIn);
  	}
  }
  
} // namespace


#endif

