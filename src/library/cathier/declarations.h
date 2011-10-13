#ifndef _DECLARATIONS_H_
#define _DECLARATIONS_H_

// External declarations
//template <int D, typename T> class AimsTimeSurfaceFaceCollection;
//template <int D, typename T> class AimsSurfaceFaceCollection;

// TIL declarations
namespace til
{
  template < class TClass, typename TNewPrecision > 
                                                class ChangePrecision;
  template < typename TContainer >              class const_cyclic_iterator;
  template < typename TContainer >              class cyclic_iterator;
  template < typename TParam >                  class Mesh;
  template < class TMesh >                      class MeshAttributes;
  template < typename TParam >                  class MeshFaceCollection;
  template < typename T, std::size_t N >        class NaryTree;
  template < typename T, std::size_t D>         class Point;
  template < typename T >                       class sparse_vector;
  template < typename T, std::size_t D >        class Vector;
  

  template < class TMesh >                struct  MeshTraits;

  namespace detail
  {
    template < class TMesh, typename TParam >      class   AddNormalAttribute;
    template < class TMesh, typename TParam >      class   AddNeighborIndexAttribute;
  }

} // namespace til

#endif //_DECLARATIONS_H_
