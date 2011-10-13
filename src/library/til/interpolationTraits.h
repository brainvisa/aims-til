#ifndef TIL_INTERPOLATION_TRAITS_H
#define TIL_INTERPOLATION_TRAITS_H

// includes from TIL library
#include "til/CatmullRomInterpolation.h"
#include "til/CubicSplineInterpolation.h"
#include "til/KeysInterpolation.h"
#include "til/LagrangeInterpolation.h"
#include "til/TrilinearInterpolation.h"

namespace til
{

  // TODO: shouldn'it be interpolation_traits::sample_size, rather than interpolationssamplesize::value?... pros/cons?

	// Output precision of interpolation
	// TODO: use traits on value_type rather than fixed precision...
	typedef double t_interp;

	/// Traits for interpolators representing the number of sample points
	/// each of these schemes use
	template < typename Interpolator >
	struct InterpolationSampleSize {};


	// A macro to define specialization of interpolationSampleSize.
	// NB: to be undefined at the end of the definitions
#define DEFINE_INTSAMPLESIZE_SPECIALIZATION(interpolator, sampleSize)   \
  template <>                                                           \
  struct InterpolationSampleSize< interpolator >                        \
  {                                                                     \
    static const int value = (sampleSize);                              \
  };                                                                    \

	//DEFINE_INTSAMPLESIZE_SPECIALIZATION(CatmullRomInterpolation, 4);
	//DEFINE_INTSAMPLESIZE_SPECIALIZATION(CubicSplineInterpolation, 4);
	//DEFINE_INTSAMPLESIZE_SPECIALIZATION(KeysInterpolation, 6);
	//DEFINE_INTSAMPLESIZE_SPECIALIZATION(Lagrange4Interpolation, 4);
	//DEFINE_INTSAMPLESIZE_SPECIALIZATION(LinearInterpolation, 2);

  template < typename T >
  struct InterpolationSampleSize<CubicSplineInterpolation<T> >
  { static const int value = 4; };

  template < typename T >
  struct InterpolationSampleSize<CatmullRomInterpolation<T> >
  { static const int value = 4; };

  template < typename T >
  struct InterpolationSampleSize<KeysInterpolation<T> >
  { static const int value = 6; };

  template < typename T >
  struct InterpolationSampleSize<Lagrange4Interpolation<T> >
  { static const int value = 4; };

  template < typename T >
  struct InterpolationSampleSize<LinearInterpolation<T> >
  { static const int value = 2; };

#undef DEFINE_INTSAMPLESIZE_SPECIALIZATION




} // namespace til


#endif

