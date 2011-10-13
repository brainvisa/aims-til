#ifndef TIL_SYMMATRIX3_TOOLS_H
#define TIL_SYMMATRIX3_TOOLS_H

/// \file Belongs to symmatrix package -- do not include directly, include SymMatrix3.h instead

// includes from TIL
#include "til/eigen3D.h"
#include "til/numeric_array.h"

// Namespace 

namespace til
{
	//namespace linalg {

// Computes the eigenvalues of a SymMatrix3

template < class T >
void eigen3D(const SymMatrix3<T> & mat, numeric_array<T,3> & vec)
{
	T v1, v2, v3;

	eigen3D(mat(0,0), mat(0,1), mat(0,2), mat(1,1), mat(1,2), mat(2,2),
		v1, v2, v3);

	vec[0] = v1;
	vec[1] = v2;
	vec[2] = v3;
}


// Computes the eigenvalues of a SymMatrix3

template < class T >
void eigen3D(const SymMatrix3<T> & mat, T & v1, T & v2, T & v3)
{
	eigen3D(mat(0,0), mat(0,1), mat(0,2), mat(1,1), mat(1,2), mat(2,2),
		v1, v2, v3);
}


template < typename T1, typename T2 >
INLINE void add(const SymMatrix3<T1> & mat1,
                SymMatrix3<T2> & mat2)
{
  for (std::size_t j = 0; j < 3; ++j)
  {
    for (std::size_t i = j; j < 3; ++j)
    {
      mat2(i,j) += mat1(i,j);
    }
  }
}

	//} // namespace linalg
} // namespace til


#endif


