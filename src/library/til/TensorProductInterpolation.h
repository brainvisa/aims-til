#ifndef TIL_TENSOR_PRODUCT_INTERPOLATION_H
#define TIL_TENSOR_PRODUCT_INTERPOLATION_H

// include from STL library
#include <cassert>

// includes from TIL library
#include "til/til_common.h"
#include "til/til_declarations.h"
#include "til/interpolationTraits.h"


// The TensorProductInterpolation class is used to interpolate images
// It builds up a tensor interpolation scheme starting from a standard
// 1-D interpolation, by interpolating intensity values in the X
// direction, then the Y, then the Z.

// Example: To use trilinear interpolation:
// typedef TensorProductInterpolation<MyImage, LinearInterpolation> TrilinearInterpolation;


namespace til {

template < typename Interpolation, int TNNeighbors = InterpolationSampleSize<Interpolation>::value >
class TensorProductInterpolation;

template < typename Interpolation>
class TensorProductInterpolation< Interpolation, 2>
{
public: // typedefs

	typedef t_interp ReturnType;
	typedef TensorProductInterpolation< Interpolation, 2> Self;

public:

	template < class Extrapolator, class TImage >
	INLINE static t_interp computeUnsafe(const TImage &im, const numeric_array<t_interp,3> &pos)
	{ return Self::template computeUnsafe<Extrapolator>(im, EXPAND_VECTOR(pos)); }

	template < class Extrapolator, class TImage >
	INLINE static t_interp compute(const TImage &im, const numeric_array<t_interp,3> &pos)
	{ return Self::template compute<Extrapolator>(im, EXPAND_VECTOR(pos)); }

	template < class Extrapolator, class TImage >
	INLINE static t_interp computeUnsafe(const TImage &im, t_interp x, t_interp y, t_interp z)
	{
		int i1, j1, k1;
		int i2, j2, k2;
		
		t_interp dx, dy, dz;

		dx = x - (i1 = floor(x)); i2 = i1+1;
		dy = y - (j1 = floor(y)); j2 = j1+1;
		dz = z - (k1 = floor(z)); k2 = k1+1;

		typename Iterator<TImage>::ConstVolumetric iIm(im);
		iIm.setPos(i1, j1, k1);

		assert(dx>=0); assert(dx<1);
		assert(dy>=0); assert(dy<1);
		assert(dz>=0); assert(dz<1);

		return 
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			*iIm,
			iIm.template getUnsafeValue<Extrapolator, +1, 0, 0>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getUnsafeValue<Extrapolator,  0,+1, 0>(),
			iIm.template getUnsafeValue<Extrapolator, +1,+1, 0>(), dx), dy),
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			iIm.template getUnsafeValue<Extrapolator,  0, 0,+1>(),
			iIm.template getUnsafeValue<Extrapolator, +1, 0,+1>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getUnsafeValue<Extrapolator,  0,+1,+1>(),
			iIm.template getUnsafeValue<Extrapolator, +1,+1,+1>(), dx), dy), dz);
	}

	template < class Extrapolator, class TImage >
	INLINE static t_interp compute(const TImage &im, t_interp x, t_interp y, t_interp z)
	{
		int i1, j1, k1;
		int i2, j2, k2;
		
		t_interp dx, dy, dz;

		dx = x - (i1 = floor(x)); i2 = i1+1;
		dy = y - (j1 = floor(y)); j2 = j1+1;
		dz = z - (k1 = floor(z)); k2 = k1+1;

		assert(dx>=0); assert(dx<1);
		assert(dy>=0); assert(dy<1);
		assert(dz>=0); assert(dz<1);

		typename Iterator<TImage>::ConstVolumetric iIm(im);
		iIm.setPos(i1, j1, k1);


/*		if (contains(getRange(im), Range(i1, j1, k1, i4, j4, k4)))
		{
			return <t_interp>Unsafe(im, x, y, z);
		}
		else*/
		{
		return
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			*iIm,
			iIm.template getValue<Extrapolator,+1, 0, 0>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getValue<Extrapolator, 0,+1, 0>(),
			iIm.template getValue<Extrapolator,+1,+1, 0>(), dx), dy),
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			iIm.template getValue<Extrapolator, 0, 0,+1>(),
			iIm.template getValue<Extrapolator,+1, 0,+1>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getValue<Extrapolator, 0,+1,+1>(),
			iIm.template getValue<Extrapolator,+1,+1,+1>(), dx), dy), dz);
		}
	}

};

template < typename Interpolation>
class TensorProductInterpolation<Interpolation, 4>
{
public: // typedefs

	typedef t_interp ReturnType;
	typedef TensorProductInterpolation<Interpolation, 4> Self;

public:

	template < class Extrapolator, class TImage >
	INLINE static t_interp computeUnsafe(const TImage &im, const numeric_array<t_interp,3> &pos)
	{ return Self::template computeUnsafe<Extrapolator>(im, EXPAND_VECTOR(pos)); }
	template < class Extrapolator, class TImage >
	INLINE static t_interp compute(const TImage &im, const numeric_array<t_interp,3> &pos)
	{ return Self::template compute<Extrapolator>(im, EXPAND_VECTOR(pos)); }

	template < class Extrapolator, class TImage >
	INLINE static t_interp computeUnsafe(const TImage &im, t_interp x, t_interp y, t_interp z)
	{
		int i1, j1, k1;
		int i2, j2, k2;
		int i3, j3, k3;
		int i4, j4, k4;
		
		t_interp dx, dy, dz;

		dx = x - (i2 = floor(x)); Self::span(i2, i1, i3, i4);
		dy = y - (j2 = floor(y)); Self::span(j2, j1, j3, j4);
		dz = z - (k2 = floor(z)); Self::span(k2, k1, k3, k4);		

		assert(dx>=0); assert(dx<1);
		assert(dy>=0); assert(dy<1);
		assert(dz>=0); assert(dz<1);

		typename Iterator<TImage>::ConstVolumetric iIm(im);
		iIm.setPos(i2, j2, k2);


		return Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			iIm.template getUnsafeValue<Extrapolator,-1,-1,-1>(),
			iIm.template getUnsafeValue<Extrapolator, 0,-1,-1>(),
			iIm.template getUnsafeValue<Extrapolator,+1,-1,-1>(),
			iIm.template getUnsafeValue<Extrapolator,+2,-1,-1>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getUnsafeValue<Extrapolator,-1, 0,-1>(),
			iIm.template getUnsafeValue<Extrapolator, 0, 0,-1>(),
			iIm.template getUnsafeValue<Extrapolator,+1, 0,-1>(),
			iIm.template getUnsafeValue<Extrapolator,+2, 0,-1>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getUnsafeValue<Extrapolator,-1,+1,-1>(),
			iIm.template getUnsafeValue<Extrapolator, 0,+1,-1>(),
			iIm.template getUnsafeValue<Extrapolator,+1,+1,-1>(),
			iIm.template getUnsafeValue<Extrapolator,+2,+1,-1>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getUnsafeValue<Extrapolator,-1,+2,-1>(),
			iIm.template getUnsafeValue<Extrapolator, 0,+2,-1>(),
			iIm.template getUnsafeValue<Extrapolator,+1,+2,-1>(),
			iIm.template getUnsafeValue<Extrapolator,+2,+2,-1>(), dx), dy),
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			iIm.template getUnsafeValue<Extrapolator,-1,-1, 0>(),
			iIm.template getUnsafeValue<Extrapolator, 0,-1, 0>(),
			iIm.template getUnsafeValue<Extrapolator,+1,-1, 0>(),
			iIm.template getUnsafeValue<Extrapolator,+2,-1, 0>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getUnsafeValue<Extrapolator,-1, 0, 0>(),
			*iIm,
			iIm.template getUnsafeValue<Extrapolator,+1, 0, 0>(),
			iIm.template getUnsafeValue<Extrapolator,+2, 0, 0>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getUnsafeValue<Extrapolator,-1,+1, 0>(),
			iIm.template getUnsafeValue<Extrapolator, 0,+1, 0>(),
			iIm.template getUnsafeValue<Extrapolator,+1,+1, 0>(),
			iIm.template getUnsafeValue<Extrapolator,+2,+1, 0>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getUnsafeValue<Extrapolator,-1,+2, 0>(),
			iIm.template getUnsafeValue<Extrapolator, 0,+2, 0>(),
			iIm.template getUnsafeValue<Extrapolator,+1,+2, 0>(),
			iIm.template getUnsafeValue<Extrapolator,+2,+2, 0>(), dx), dy),
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			iIm.template getUnsafeValue<Extrapolator,-1,-1,+1>(),
			iIm.template getUnsafeValue<Extrapolator, 0,-1,+1>(),
			iIm.template getUnsafeValue<Extrapolator,+1,-1,+1>(),
			iIm.template getUnsafeValue<Extrapolator,+2,-1,+1>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getUnsafeValue<Extrapolator,-1, 0,+1>(),
			iIm.template getUnsafeValue<Extrapolator, 0, 0,+1>(),
			iIm.template getUnsafeValue<Extrapolator,+1, 0,+1>(),
			iIm.template getUnsafeValue<Extrapolator,+2, 0,+1>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getUnsafeValue<Extrapolator,-1,+1,+1>(),
			iIm.template getUnsafeValue<Extrapolator, 0,+1,+1>(),
			iIm.template getUnsafeValue<Extrapolator,+1,+1,+1>(),
			iIm.template getUnsafeValue<Extrapolator,+2,+1,+1>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getUnsafeValue<Extrapolator,-1,+2,+1>(),
			iIm.template getUnsafeValue<Extrapolator, 0,+2,+1>(),
			iIm.template getUnsafeValue<Extrapolator,+1,+2,+1>(),
			iIm.template getUnsafeValue<Extrapolator,+2,+2,+1>(), dx), dy),
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			iIm.template getUnsafeValue<Extrapolator,-1,-1,+2>(),
			iIm.template getUnsafeValue<Extrapolator, 0,-1,+2>(),
			iIm.template getUnsafeValue<Extrapolator,+1,-1,+2>(),
			iIm.template getUnsafeValue<Extrapolator,+2,-1,+2>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getUnsafeValue<Extrapolator,-1, 0,+2>(),
			iIm.template getUnsafeValue<Extrapolator, 0, 0,+2>(),
			iIm.template getUnsafeValue<Extrapolator,+1, 0,+2>(),
			iIm.template getUnsafeValue<Extrapolator,+2, 0,+2>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getUnsafeValue<Extrapolator,-1,+1,+2>(),
			iIm.template getUnsafeValue<Extrapolator, 0,+1,+2>(),
			iIm.template getUnsafeValue<Extrapolator,+1,+1,+2>(),
			iIm.template getUnsafeValue<Extrapolator,+2,+1,+2>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getUnsafeValue<Extrapolator,-1,+2,+2>(),
			iIm.template getUnsafeValue<Extrapolator, 0,+2,+2>(),
			iIm.template getUnsafeValue<Extrapolator,+1,+2,+2>(),
			iIm.template getUnsafeValue<Extrapolator,+2,+2,+2>(), dx), dy), dz);
	}

	template < class Extrapolator, class TImage >
	INLINE static t_interp compute(const TImage &im, t_interp x, t_interp y, t_interp z)
	{
		int i1, j1, k1;
		int i2, j2, k2;
		int i3, j3, k3;
		int i4, j4, k4;
		
		t_interp dx, dy, dz;

		dx = x - (i2 = (int)floor(x)); Self::span(i2, i1, i3, i4);
		dy = y - (j2 = (int)floor(y)); Self::span(j2, j1, j3, j4);
		dz = z - (k2 = (int)floor(z)); Self::span(k2, k1, k3, k4);		

		assert(dx>=0); assert(dx<1);
		assert(dy>=0); assert(dy<1);
		assert(dz>=0); assert(dz<1);

		typename Iterator<TImage>::ConstVolumetric iIm(im);
		iIm.set_pos(numeric_array<int,3>(i2, j2, k2));
		
		/*
		std::cout << "i" <<" "<< i1 <<" "<< i2 <<" "<< i3 <<" "<< i4 << std::endl;
		std::cout << "j" <<" "<< j1 <<" "<< j2 <<" "<< j3 <<" "<< j4 << std::endl;
		std::cout << "k" <<" "<< k1 <<" "<< k2 <<" "<< k3 <<" "<< k4 << std::endl;
*/
/*		if (contains(getRange(im), Range(i1, j1, k1, i4, j4, k4)))
		{
			return Self::compute<t_interp>Unsafe(im, x, y, z);
		}
		else*/
		{
		return
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			iIm.template getValue<Extrapolator,-1,-1,-1>(),
			iIm.template getValue<Extrapolator, 0,-1,-1>(),
			iIm.template getValue<Extrapolator,+1,-1,-1>(),
			iIm.template getValue<Extrapolator,+2,-1,-1>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getValue<Extrapolator,-1, 0,-1>(),
			iIm.template getValue<Extrapolator, 0, 0,-1>(),
			iIm.template getValue<Extrapolator,+1, 0,-1>(),
			iIm.template getValue<Extrapolator,+2, 0,-1>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getValue<Extrapolator,-1,+1,-1>(),
			iIm.template getValue<Extrapolator, 0,+1,-1>(),
			iIm.template getValue<Extrapolator,+1,+1,-1>(),
			iIm.template getValue<Extrapolator,+2,+1,-1>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getValue<Extrapolator,-1,+2,-1>(),
			iIm.template getValue<Extrapolator, 0,+2,-1>(),
			iIm.template getValue<Extrapolator,+1,+2,-1>(),
			iIm.template getValue<Extrapolator,+2,+2,-1>(), dx), dy),
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			iIm.template getValue<Extrapolator,-1,-1, 0>(),
			iIm.template getValue<Extrapolator, 0,-1, 0>(),
			iIm.template getValue<Extrapolator,+1,-1, 0>(),
			iIm.template getValue<Extrapolator,+2,-1, 0>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getValue<Extrapolator,-1, 0, 0>(),
			*iIm,
			iIm.template getValue<Extrapolator,+1, 0, 0>(),
			iIm.template getValue<Extrapolator,+2, 0, 0>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getValue<Extrapolator,-1,+1, 0>(),
			iIm.template getValue<Extrapolator, 0,+1, 0>(),
			iIm.template getValue<Extrapolator,+1,+1, 0>(),
			iIm.template getValue<Extrapolator,+2,+1, 0>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getValue<Extrapolator,-1,+2, 0>(),
			iIm.template getValue<Extrapolator, 0,+2, 0>(),
			iIm.template getValue<Extrapolator,+1,+2, 0>(),
			iIm.template getValue<Extrapolator,+2,+2, 0>(), dx), dy),
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			iIm.template getValue<Extrapolator,-1,-1,+1>(),
			iIm.template getValue<Extrapolator, 0,-1,+1>(),
			iIm.template getValue<Extrapolator,+1,-1,+1>(),
			iIm.template getValue<Extrapolator,+2,-1,+1>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getValue<Extrapolator,-1, 0,+1>(),
			iIm.template getValue<Extrapolator, 0, 0,+1>(),
			iIm.template getValue<Extrapolator,+1, 0,+1>(),
			iIm.template getValue<Extrapolator,+2, 0,+1>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getValue<Extrapolator,-1,+1,+1>(),
			iIm.template getValue<Extrapolator, 0,+1,+1>(),
			iIm.template getValue<Extrapolator,+1,+1,+1>(),
			iIm.template getValue<Extrapolator,+2,+1,+1>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getValue<Extrapolator,-1,+2,+1>(),
			iIm.template getValue<Extrapolator, 0,+2,+1>(),
			iIm.template getValue<Extrapolator,+1,+2,+1>(),
			iIm.template getValue<Extrapolator,+2,+2,+1>(), dx), dy),
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			iIm.template getValue<Extrapolator,-1,-1,+2>(),
			iIm.template getValue<Extrapolator, 0,-1,+2>(),
			iIm.template getValue<Extrapolator,+1,-1,+2>(),
			iIm.template getValue<Extrapolator,+2,-1,+2>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getValue<Extrapolator,-1, 0,+2>(),
			iIm.template getValue<Extrapolator, 0, 0,+2>(),
			iIm.template getValue<Extrapolator,+1, 0,+2>(),
			iIm.template getValue<Extrapolator,+2, 0,+2>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getValue<Extrapolator,-1,+1,+2>(),
			iIm.template getValue<Extrapolator, 0,+1,+2>(),
			iIm.template getValue<Extrapolator,+1,+1,+2>(),
			iIm.template getValue<Extrapolator,+2,+1,+2>(), dx),
			Interpolation::template compute<t_interp>(
			iIm.template getValue<Extrapolator,-1,+2,+2>(),
			iIm.template getValue<Extrapolator, 0,+2,+2>(),
			iIm.template getValue<Extrapolator,+1,+2,+2>(),
			iIm.template getValue<Extrapolator,+2,+2,+2>(), dx), dy), dz);
		}
	}

private:

	INLINE static void span(int i2, int &i1, int &i3, int &i4)
	{
		i1 = i2-1;
		i4 = (i3 = i2 + 1) + 1;
	}
};

template < typename Interpolation>
class TensorProductInterpolation<Interpolation, 6>
{
public: // typedefs

	typedef t_interp ReturnType;
	typedef TensorProductInterpolation<Interpolation, 6> Self;

public:

	template < class TImage >
	INLINE static t_interp computeUnsafe(const TImage &im, const numeric_array<t_interp,3> &pos)
	{ return Self::computeUnsafe(im, EXPAND_VECTOR(pos)); }
	template < class TImage >
	INLINE static t_interp compute(const TImage &im, const numeric_array<t_interp,3> &pos)
	{ return Self::compute(im, EXPAND_VECTOR(pos)); }

	template < class TImage >
	INLINE static t_interp computeUnsafe(const TImage &im, t_interp x, t_interp y, t_interp z)
	{
		int i1, j1, k1;
		int i2, j2, k2;
		int i3, j3, k3;
		int i4, j4, k4;
		int i5, j5, k5;
		int i6, j6, k6;
		
		t_interp dx, dy, dz;

		dx = x - (i3 = floor(x)); Self::span(i3, i1, i2, i4, i5, i6);
		dy = y - (j3 = floor(y)); Self::span(j3, j1, j2, j4, j5, j6);
		dz = z - (k3 = floor(z)); Self::span(k3, k1, k2, k4, k5, k6);		

		assert(dx>=0); assert(dx<1);
		assert(dy>=0); assert(dy<1);
		assert(dz>=0); assert(dz<1);

		return Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j1,k1),
			im.getUnsafeValue(i2,j1,k1),
			im.getUnsafeValue(i3,j1,k1),
			im.getUnsafeValue(i4,j1,k1),
			im.getUnsafeValue(i5,j1,k1),
			im.getUnsafeValue(i6,j1,k1), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j2,k1),
			im.getUnsafeValue(i2,j2,k1),
			im.getUnsafeValue(i3,j2,k1),
			im.getUnsafeValue(i4,j2,k1),
			im.getUnsafeValue(i5,j2,k1),
			im.getUnsafeValue(i6,j2,k1), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j3,k1),
			im.getUnsafeValue(i2,j3,k1),
			im.getUnsafeValue(i3,j3,k1),
			im.getUnsafeValue(i4,j3,k1),
			im.getUnsafeValue(i5,j3,k1),
			im.getUnsafeValue(i6,j3,k1), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j4,k1),
			im.getUnsafeValue(i2,j4,k1),
			im.getUnsafeValue(i3,j4,k1),
			im.getUnsafeValue(i4,j4,k1),
			im.getUnsafeValue(i5,j4,k1),
			im.getUnsafeValue(i6,j4,k1), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j5,k1),
			im.getUnsafeValue(i2,j5,k1),
			im.getUnsafeValue(i3,j5,k1),
			im.getUnsafeValue(i4,j5,k1),
			im.getUnsafeValue(i5,j5,k1),
			im.getUnsafeValue(i6,j5,k1), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j6,k1),
			im.getUnsafeValue(i2,j6,k1),
			im.getUnsafeValue(i3,j6,k1),
			im.getUnsafeValue(i4,j6,k1),
			im.getUnsafeValue(i5,j6,k1),
			im.getUnsafeValue(i6,j6,k1), dx), dy),
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j1,k2),
			im.getUnsafeValue(i2,j1,k2),
			im.getUnsafeValue(i3,j1,k2),
			im.getUnsafeValue(i4,j1,k2),
			im.getUnsafeValue(i5,j1,k2),
			im.getUnsafeValue(i6,j1,k2), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j2,k2),
			im.getUnsafeValue(i2,j2,k2),
			im.getUnsafeValue(i3,j2,k2),
			im.getUnsafeValue(i4,j2,k2),
			im.getUnsafeValue(i5,j2,k2),
			im.getUnsafeValue(i6,j2,k2), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j3,k2),
			im.getUnsafeValue(i2,j3,k2),
			im.getUnsafeValue(i3,j3,k2),
			im.getUnsafeValue(i4,j3,k2),
			im.getUnsafeValue(i5,j3,k2),
			im.getUnsafeValue(i6,j3,k2), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j4,k2),
			im.getUnsafeValue(i2,j4,k2),
			im.getUnsafeValue(i3,j4,k2),
			im.getUnsafeValue(i4,j4,k2),
			im.getUnsafeValue(i5,j4,k2),
			im.getUnsafeValue(i6,j4,k2), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j5,k2),
			im.getUnsafeValue(i2,j5,k2),
			im.getUnsafeValue(i3,j5,k2),
			im.getUnsafeValue(i4,j5,k2),
			im.getUnsafeValue(i5,j5,k2),
			im.getUnsafeValue(i6,j5,k2), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j6,k2),
			im.getUnsafeValue(i2,j6,k2),
			im.getUnsafeValue(i3,j6,k2),
			im.getUnsafeValue(i4,j6,k2),
			im.getUnsafeValue(i5,j6,k2),
			im.getUnsafeValue(i6,j6,k2), dx), dy),
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j1,k3),
			im.getUnsafeValue(i2,j1,k3),
			im.getUnsafeValue(i3,j1,k3),
			im.getUnsafeValue(i4,j1,k3),
			im.getUnsafeValue(i5,j1,k3),
			im.getUnsafeValue(i6,j1,k3), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j2,k3),
			im.getUnsafeValue(i2,j2,k3),
			im.getUnsafeValue(i3,j2,k3),
			im.getUnsafeValue(i4,j2,k3),
			im.getUnsafeValue(i5,j2,k3),
			im.getUnsafeValue(i6,j2,k3), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j3,k3),
			im.getUnsafeValue(i2,j3,k3),
			im.getUnsafeValue(i3,j3,k3),
			im.getUnsafeValue(i4,j3,k3),
			im.getUnsafeValue(i5,j3,k3),
			im.getUnsafeValue(i6,j3,k3), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j4,k3),
			im.getUnsafeValue(i2,j4,k3),
			im.getUnsafeValue(i3,j4,k3),
			im.getUnsafeValue(i4,j4,k3),
			im.getUnsafeValue(i5,j4,k3),
			im.getUnsafeValue(i6,j4,k3), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j5,k3),
			im.getUnsafeValue(i2,j5,k3),
			im.getUnsafeValue(i3,j5,k3),
			im.getUnsafeValue(i4,j5,k3),
			im.getUnsafeValue(i5,j5,k3),
			im.getUnsafeValue(i6,j5,k3), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j6,k3),
			im.getUnsafeValue(i2,j6,k3),
			im.getUnsafeValue(i3,j6,k3),
			im.getUnsafeValue(i4,j6,k3),
			im.getUnsafeValue(i5,j6,k3),
			im.getUnsafeValue(i6,j6,k3), dx), dy),
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j1,k4),
			im.getUnsafeValue(i2,j1,k4),
			im.getUnsafeValue(i3,j1,k4),
			im.getUnsafeValue(i4,j1,k4),
			im.getUnsafeValue(i5,j1,k4),
			im.getUnsafeValue(i6,j1,k4), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j2,k4),
			im.getUnsafeValue(i2,j2,k4),
			im.getUnsafeValue(i3,j2,k4),
			im.getUnsafeValue(i4,j2,k4),
			im.getUnsafeValue(i5,j2,k4),
			im.getUnsafeValue(i6,j2,k4), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j3,k4),
			im.getUnsafeValue(i2,j3,k4),
			im.getUnsafeValue(i3,j3,k4),
			im.getUnsafeValue(i4,j3,k4),
			im.getUnsafeValue(i5,j3,k4),
			im.getUnsafeValue(i6,j3,k4), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j4,k4),
			im.getUnsafeValue(i2,j4,k4),
			im.getUnsafeValue(i3,j4,k4),
			im.getUnsafeValue(i4,j4,k4),
			im.getUnsafeValue(i5,j4,k4),
			im.getUnsafeValue(i6,j4,k4), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j5,k4),
			im.getUnsafeValue(i2,j5,k4),
			im.getUnsafeValue(i3,j5,k4),
			im.getUnsafeValue(i4,j5,k4),
			im.getUnsafeValue(i5,j5,k4),
			im.getUnsafeValue(i6,j5,k4), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j6,k4),
			im.getUnsafeValue(i2,j6,k4),
			im.getUnsafeValue(i3,j6,k4),
			im.getUnsafeValue(i4,j6,k4),
			im.getUnsafeValue(i5,j6,k4),
			im.getUnsafeValue(i6,j6,k4), dx), dy),
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j1,k5),
			im.getUnsafeValue(i2,j1,k5),
			im.getUnsafeValue(i3,j1,k5),
			im.getUnsafeValue(i4,j1,k5),
			im.getUnsafeValue(i5,j1,k5),
			im.getUnsafeValue(i6,j1,k5), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j2,k5),
			im.getUnsafeValue(i2,j2,k5),
			im.getUnsafeValue(i3,j2,k5),
			im.getUnsafeValue(i4,j2,k5),
			im.getUnsafeValue(i5,j2,k5),
			im.getUnsafeValue(i6,j2,k5), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j3,k5),
			im.getUnsafeValue(i2,j3,k5),
			im.getUnsafeValue(i3,j3,k5),
			im.getUnsafeValue(i4,j3,k5),
			im.getUnsafeValue(i5,j3,k5),
			im.getUnsafeValue(i6,j3,k5), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j4,k5),
			im.getUnsafeValue(i2,j4,k5),
			im.getUnsafeValue(i3,j4,k5),
			im.getUnsafeValue(i4,j4,k5),
			im.getUnsafeValue(i5,j4,k5),
			im.getUnsafeValue(i6,j4,k5), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j5,k5),
			im.getUnsafeValue(i2,j5,k5),
			im.getUnsafeValue(i3,j5,k5),
			im.getUnsafeValue(i4,j5,k5),
			im.getUnsafeValue(i5,j5,k5),
			im.getUnsafeValue(i6,j5,k5), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j6,k5),
			im.getUnsafeValue(i2,j6,k5),
			im.getUnsafeValue(i3,j6,k5),
			im.getUnsafeValue(i4,j6,k5),
			im.getUnsafeValue(i5,j6,k5),
			im.getUnsafeValue(i6,j6,k5), dx), dy),
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j1,k6),
			im.getUnsafeValue(i2,j1,k6),
			im.getUnsafeValue(i3,j1,k6),
			im.getUnsafeValue(i4,j1,k6),
			im.getUnsafeValue(i5,j1,k6),
			im.getUnsafeValue(i6,j1,k6), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j2,k6),
			im.getUnsafeValue(i2,j2,k6),
			im.getUnsafeValue(i3,j2,k6),
			im.getUnsafeValue(i4,j2,k6),
			im.getUnsafeValue(i5,j2,k6),
			im.getUnsafeValue(i6,j2,k6), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j3,k6),
			im.getUnsafeValue(i2,j3,k6),
			im.getUnsafeValue(i3,j3,k6),
			im.getUnsafeValue(i4,j3,k6),
			im.getUnsafeValue(i5,j3,k6),
			im.getUnsafeValue(i6,j3,k6), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j4,k6),
			im.getUnsafeValue(i2,j4,k6),
			im.getUnsafeValue(i3,j4,k6),
			im.getUnsafeValue(i4,j4,k6),
			im.getUnsafeValue(i5,j4,k6),
			im.getUnsafeValue(i6,j4,k6), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j5,k6),
			im.getUnsafeValue(i2,j5,k6),
			im.getUnsafeValue(i3,j5,k6),
			im.getUnsafeValue(i4,j5,k6),
			im.getUnsafeValue(i5,j5,k6),
			im.getUnsafeValue(i6,j5,k6), dx),
			Interpolation::template compute<t_interp>(
			im.getUnsafeValue(i1,j6,k6),
			im.getUnsafeValue(i2,j6,k6),
			im.getUnsafeValue(i3,j6,k6),
			im.getUnsafeValue(i4,j6,k6),
			im.getUnsafeValue(i5,j6,k6),
			im.getUnsafeValue(i6,j6,k6), dx), dy), dz);
	}

	template < class TImage >
	INLINE static t_interp compute(const TImage &im, t_interp x, t_interp y, t_interp z)
	{
		int i1, j1, k1;
		int i2, j2, k2;
		int i3, j3, k3;
		int i4, j4, k4;
		int i5, j5, k5;
		int i6, j6, k6;
		
		t_interp dx, dy, dz;

		dx = x - (i3 = floor(x)); Self::span(i3, i1, i2, i4, i5, i6);
		dy = y - (j3 = floor(y)); Self::span(j3, j1, j2, j4, j5, j6);
		dz = z - (k3 = floor(z)); Self::span(k3, k1, k2, k4, k5, k6);		

		assert(dx>=0); assert(dx<1);
		assert(dy>=0); assert(dy<1);
		assert(dz>=0); assert(dz<1);

		
		/*
		std::cout << "i" <<" "<< i1 <<" "<< i2 <<" "<< i3 <<" "<< i4 << std::endl;
		std::cout << "j" <<" "<< j1 <<" "<< j2 <<" "<< j3 <<" "<< j4 << std::endl;
		std::cout << "k" <<" "<< k1 <<" "<< k2 <<" "<< k3 <<" "<< k4 << std::endl;
*/
/*		if (contains(getRange(im), Range(i1, j1, k1, i4, j4, k4)))
		{
			return Self::compute<t_interp>Unsafe(im, x, y, z);
		}
		else*/
		{

    // TODO: obviously, some kind of iteration/loop/recursion is needed here.

		return Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j1,k1),
			im.getValue(i2,j1,k1),
			im.getValue(i3,j1,k1),
			im.getValue(i4,j1,k1),
			im.getValue(i5,j1,k1),
			im.getValue(i6,j1,k1), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j2,k1),
			im.getValue(i2,j2,k1),
			im.getValue(i3,j2,k1),
			im.getValue(i4,j2,k1),
			im.getValue(i5,j2,k1),
			im.getValue(i6,j2,k1), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j3,k1),
			im.getValue(i2,j3,k1),
			im.getValue(i3,j3,k1),
			im.getValue(i4,j3,k1),
			im.getValue(i5,j3,k1),
			im.getValue(i6,j3,k1), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j4,k1),
			im.getValue(i2,j4,k1),
			im.getValue(i3,j4,k1),
			im.getValue(i4,j4,k1),
			im.getValue(i5,j4,k1),
			im.getValue(i6,j4,k1), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j5,k1),
			im.getValue(i2,j5,k1),
			im.getValue(i3,j5,k1),
			im.getValue(i4,j5,k1),
			im.getValue(i5,j5,k1),
			im.getValue(i6,j5,k1), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j6,k1),
			im.getValue(i2,j6,k1),
			im.getValue(i3,j6,k1),
			im.getValue(i4,j6,k1),
			im.getValue(i5,j6,k1),
			im.getValue(i6,j6,k1), dx), dy),
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j1,k2),
			im.getValue(i2,j1,k2),
			im.getValue(i3,j1,k2),
			im.getValue(i4,j1,k2),
			im.getValue(i5,j1,k2),
			im.getValue(i6,j1,k2), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j2,k2),
			im.getValue(i2,j2,k2),
			im.getValue(i3,j2,k2),
			im.getValue(i4,j2,k2),
			im.getValue(i5,j2,k2),
			im.getValue(i6,j2,k2), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j3,k2),
			im.getValue(i2,j3,k2),
			im.getValue(i3,j3,k2),
			im.getValue(i4,j3,k2),
			im.getValue(i5,j3,k2),
			im.getValue(i6,j3,k2), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j4,k2),
			im.getValue(i2,j4,k2),
			im.getValue(i3,j4,k2),
			im.getValue(i4,j4,k2),
			im.getValue(i5,j4,k2),
			im.getValue(i6,j4,k2), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j5,k2),
			im.getValue(i2,j5,k2),
			im.getValue(i3,j5,k2),
			im.getValue(i4,j5,k2),
			im.getValue(i5,j5,k2),
			im.getValue(i6,j5,k2), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j6,k2),
			im.getValue(i2,j6,k2),
			im.getValue(i3,j6,k2),
			im.getValue(i4,j6,k2),
			im.getValue(i5,j6,k2),
			im.getValue(i6,j6,k2), dx), dy),
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j1,k3),
			im.getValue(i2,j1,k3),
			im.getValue(i3,j1,k3),
			im.getValue(i4,j1,k3),
			im.getValue(i5,j1,k3),
			im.getValue(i6,j1,k3), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j2,k3),
			im.getValue(i2,j2,k3),
			im.getValue(i3,j2,k3),
			im.getValue(i4,j2,k3),
			im.getValue(i5,j2,k3),
			im.getValue(i6,j2,k3), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j3,k3),
			im.getValue(i2,j3,k3),
			im.getValue(i3,j3,k3),
			im.getValue(i4,j3,k3),
			im.getValue(i5,j3,k3),
			im.getValue(i6,j3,k3), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j4,k3),
			im.getValue(i2,j4,k3),
			im.getValue(i3,j4,k3),
			im.getValue(i4,j4,k3),
			im.getValue(i5,j4,k3),
			im.getValue(i6,j4,k3), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j5,k3),
			im.getValue(i2,j5,k3),
			im.getValue(i3,j5,k3),
			im.getValue(i4,j5,k3),
			im.getValue(i5,j5,k3),
			im.getValue(i6,j5,k3), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j6,k3),
			im.getValue(i2,j6,k3),
			im.getValue(i3,j6,k3),
			im.getValue(i4,j6,k3),
			im.getValue(i5,j6,k3),
			im.getValue(i6,j6,k3), dx), dy),
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j1,k4),
			im.getValue(i2,j1,k4),
			im.getValue(i3,j1,k4),
			im.getValue(i4,j1,k4),
			im.getValue(i5,j1,k4),
			im.getValue(i6,j1,k4), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j2,k4),
			im.getValue(i2,j2,k4),
			im.getValue(i3,j2,k4),
			im.getValue(i4,j2,k4),
			im.getValue(i5,j2,k4),
			im.getValue(i6,j2,k4), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j3,k4),
			im.getValue(i2,j3,k4),
			im.getValue(i3,j3,k4),
			im.getValue(i4,j3,k4),
			im.getValue(i5,j3,k4),
			im.getValue(i6,j3,k4), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j4,k4),
			im.getValue(i2,j4,k4),
			im.getValue(i3,j4,k4),
			im.getValue(i4,j4,k4),
			im.getValue(i5,j4,k4),
			im.getValue(i6,j4,k4), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j5,k4),
			im.getValue(i2,j5,k4),
			im.getValue(i3,j5,k4),
			im.getValue(i4,j5,k4),
			im.getValue(i5,j5,k4),
			im.getValue(i6,j5,k4), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j6,k4),
			im.getValue(i2,j6,k4),
			im.getValue(i3,j6,k4),
			im.getValue(i4,j6,k4),
			im.getValue(i5,j6,k4),
			im.getValue(i6,j6,k4), dx), dy),
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j1,k5),
			im.getValue(i2,j1,k5),
			im.getValue(i3,j1,k5),
			im.getValue(i4,j1,k5),
			im.getValue(i5,j1,k5),
			im.getValue(i6,j1,k5), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j2,k5),
			im.getValue(i2,j2,k5),
			im.getValue(i3,j2,k5),
			im.getValue(i4,j2,k5),
			im.getValue(i5,j2,k5),
			im.getValue(i6,j2,k5), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j3,k5),
			im.getValue(i2,j3,k5),
			im.getValue(i3,j3,k5),
			im.getValue(i4,j3,k5),
			im.getValue(i5,j3,k5),
			im.getValue(i6,j3,k5), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j4,k5),
			im.getValue(i2,j4,k5),
			im.getValue(i3,j4,k5),
			im.getValue(i4,j4,k5),
			im.getValue(i5,j4,k5),
			im.getValue(i6,j4,k5), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j5,k5),
			im.getValue(i2,j5,k5),
			im.getValue(i3,j5,k5),
			im.getValue(i4,j5,k5),
			im.getValue(i5,j5,k5),
			im.getValue(i6,j5,k5), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j6,k5),
			im.getValue(i2,j6,k5),
			im.getValue(i3,j6,k5),
			im.getValue(i4,j6,k5),
			im.getValue(i5,j6,k5),
			im.getValue(i6,j6,k5), dx), dy),
			Interpolation::template compute<t_interp>(
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j1,k6),
			im.getValue(i2,j1,k6),
			im.getValue(i3,j1,k6),
			im.getValue(i4,j1,k6),
			im.getValue(i5,j1,k6),
			im.getValue(i6,j1,k6), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j2,k6),
			im.getValue(i2,j2,k6),
			im.getValue(i3,j2,k6),
			im.getValue(i4,j2,k6),
			im.getValue(i5,j2,k6),
			im.getValue(i6,j2,k6), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j3,k6),
			im.getValue(i2,j3,k6),
			im.getValue(i3,j3,k6),
			im.getValue(i4,j3,k6),
			im.getValue(i5,j3,k6),
			im.getValue(i6,j3,k6), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j4,k6),
			im.getValue(i2,j4,k6),
			im.getValue(i3,j4,k6),
			im.getValue(i4,j4,k6),
			im.getValue(i5,j4,k6),
			im.getValue(i6,j4,k6), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j5,k6),
			im.getValue(i2,j5,k6),
			im.getValue(i3,j5,k6),
			im.getValue(i4,j5,k6),
			im.getValue(i5,j5,k6),
			im.getValue(i6,j5,k6), dx),
			Interpolation::template compute<t_interp>(
			im.getValue(i1,j6,k6),
			im.getValue(i2,j6,k6),
			im.getValue(i3,j6,k6),
			im.getValue(i4,j6,k6),
			im.getValue(i5,j6,k6),
			im.getValue(i6,j6,k6), dx), dy), dz);
		}
	}

private:

	INLINE static void span(int i3, int &i1, int &i2, int &i4, int &i5, int &i6)
	{
		i1 = (i2 = i3 - 1) - 1;
		i6 = (i5 = (i4 = i3 + 1) + 1) +1;
	}
};



} // namespace

#endif

