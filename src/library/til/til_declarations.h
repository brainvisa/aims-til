#ifndef TIL_DECLARATIONS_H
#define TIL_DECLARATIONS_H

// includes from STL
#include <cstddef>  // std::size_t

/// \file
/// This file contains forward declarations of classes defined in the
/// TIL library. This is used in particular to speed up computation time

namespace til
{
	template < typename T >                    class Affine;
	template < typename T >                    class AffineMap;
  template < typename T, std::size_t D >     class Box;
  template < class TClass, typename TNewPrecision > 
                                             class ChangePrecision;
  template < typename TContainer >           class const_cyclic_iterator;
  template < typename TContainer >           class cyclic_iterator;
	template < typename T >					           class ImageC;
	template < typename T >					           class ImageNC;
	template < typename T >                    class ImageRLE;
	template < typename T >                    class Matrix3;
  template < typename T, typename TAccumulation >
                                             class MeanAccumulator;
	template < typename TParam >               class Mesh;
  template < class TMesh >                   struct MeshAttributes;
  template < typename TParam >               class MeshFaceCollection;
  template < typename T, std::size_t N >     class NaryTree;
  template < typename T, std::size_t D >     class numeric_array;
  //template < typename T, std::size_t D, typename TStorage > class Point;
	//template < typename T >                    class Point3D;
	template < typename T >                    class Scaling;
	template < typename T >                    class ScalingMap;
  template < typename T, typename BaselinePolicy >
                                             class sparse_vector;
	template < typename T >                    class SymMatrix3;
	//template < typename T, std::size_t D, typename TStorage > class Vector;

  template < class TMesh >                   struct  MeshTraits;

  namespace detail
  {
    template < class TMesh, typename TParam >      class   AddNormalAttribute;
    template < class TMesh, typename TParam >      class   AddNeighborIndexAttribute;
  }

}

#endif

