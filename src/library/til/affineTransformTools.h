#ifndef TIL_AFFINETOOLS_H
#define TIL_AFFINETOOLS_H

/// \file Belongs to affineTransform package -- do not include directly, include til/AffineTransform.h instead

// includes from STL
#include <cassert>
#include <iostream>


// Namespace 
namespace til
{
	//namespace linalg {

/// Compute the affine transform that would send a Box on another Box.
template < typename T >
void computeAffineTransformBetween(const Box<T,3> & from, const Box<T,3> & to, Affine<T> & a)
{
  assert(all_greater_equal(size(from), 128*std::numeric_limits<T>::epsilon()));

	numeric_array<T,3> diag = size(to) / size(from);
  // TODO: create a function that initializes a matrix to a diagonal from a vector
	a.reset();
	a.getMatrix()(0,0) = diag[0];
	a.getMatrix()(1,1) = diag[1];
	a.getMatrix()(2,2) = diag[2];

	apply(a.getMatrix(), from.min_bounds(), a.getTransl());
	
	//print(a.transl());

	a.getTransl() -= to.min_bounds();
	neg(a.getTransl());
}

// TODO: is it really different from out = a*in?
template < typename T1, typename T2 >
INLINE
void apply
(
  const Affine<T1> & a,
  const numeric_array<T2,3> & in,
  numeric_array<typename combine<T1,T2>::type,3> & out
)
{
	apply(a.getMatrix(), in, out);
	out += a.getTransl();
}


template < typename T >
void print(const Affine<T> & a)
{
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			std::cout << a.getMatrix()(i,j) << " ";
		}
		std::cout << a.getTransl()[i] << std::endl;
	}
}
	//} // namespace linalg
} // namespace til

#endif

