#ifndef TIL_ISTRAITS_H_
#define TIL_ISTRAITS_H_

// includes from TIL
#include "til/ext_declarations.h"
#include "til/til_declarations.h"


namespace til {

#define TIL_DECLARE_IS(typeM)  TIL_DECLARE_IS_VALUE_TEMPLATE(typeM, 0, typename T)
#define TIL_DECLARE_IS_VALUE(typeM, valueM)  TIL_DECLARE_IS_VALUE_TEMPLATE(typeM, valueM, typename T)
// The following macro does not rely on TIL_DECLARE_IS_VALUE_TEMPLATE because 
// TIL_COMMA, sometimes used in the tplM argument, cannot be carried over: it is
// interpreted straight away, and then TIL_DECLARE_IS_VALUE_TEMPLATE receives more
// arguments than it should.
#define TIL_DECLARE_IS_TEMPLATE(typeM, tplM)              \
template < tplM >                                         \
struct is_##typeM                                         \
{                                                         \
  enum { value = 0 };                                     \
};
#define TIL_DECLARE_IS_VALUE_TEMPLATE(typeM, valM, tplM)  \
template < tplM >                                         \
struct is_##typeM                                         \
{                                                         \
  enum { value = valM };                                  \
};

#define TIL_DECLARE_IS_SPEC(typeM, cnameM) TIL_DECLARE_IS_SPEC_T(typeM, , cnameM)
#define TIL_DECLARE_IS_SPEC_T(typeM, templateArgsM, cnameM)   \
template < templateArgsM >                                    \
struct is_##typeM < cnameM  >                                 \
{                                                             \
  enum { value = 1 };                                         \
};

#define TIL_COMMA ,

//////////////////// is_3DPoint ////////////////////

/// True if type implements a 3D Point
TIL_DECLARE_IS(3DPoint);

//TIL_DECLARE_IS_SPEC_T(3DPoint, typename T TIL_COMMA typename TStorage, Point<T TIL_COMMA 3 TIL_COMMA TStorage>);
//TIL_DECLARE_IS_SPEC_T(3DPoint, typename T, AimsVector<T TIL_COMMA 3>);
TIL_DECLARE_IS_SPEC_T(3DPoint, typename T, boost::array<T TIL_COMMA 3>);

//////////////////// is_3DVector ////////////////////

/// True if type implements a 3D Vector
// NB: by default, it takes the value of is_3DPoint, because a 3D Point
// is also a 3D vector. So things that were already specified for 3D Points
// do not need to be added again.
TIL_DECLARE_IS_VALUE(3DVector, is_3DPoint<T>::value);
//TIL_DECLARE_IS_SPEC_T(3DVector, typename T TIL_COMMA typename TStorage, Vector<T TIL_COMMA 3 TIL_COMMA TStorage>);


//////////////////// is_BoostArray_N ////////////////////

/// True if type is boost::array of size N
TIL_DECLARE_IS_TEMPLATE(BoostArray_N, typename T TIL_COMMA std::size_t N);
TIL_DECLARE_IS_SPEC_T(BoostArray_N, typename T TIL_COMMA std::size_t N, boost::array<T TIL_COMMA N> TIL_COMMA N);

//////////////////// is_mesh ////////////////////

/// True if type implements a mesh
TIL_DECLARE_IS(mesh);
TIL_DECLARE_IS_SPEC_T(mesh, typename TParam, Mesh<TParam>);
//TIL_DECLARE_IS_SPEC_T(Mesh, class TMesh TIL_COMMA typename TParam, detail::AddNeighborIndexAttribute<TMesh TIL_COMMA TParam>);

//////////////////// is_container ////////////////////

/// True if type is a (STL) container
TIL_DECLARE_IS(container);
TIL_DECLARE_IS_SPEC_T(container, typename T TIL_COMMA typename A, std::vector<T TIL_COMMA A>);
TIL_DECLARE_IS_SPEC_T(container, typename T TIL_COMMA typename A, std::list<T TIL_COMMA A>);
TIL_DECLARE_IS_SPEC_T(container, typename K TIL_COMMA typename P TIL_COMMA typename A, std::set<K TIL_COMMA P TIL_COMMA A>);
TIL_DECLARE_IS_SPEC_T(container, typename K TIL_COMMA typename T TIL_COMMA typename P TIL_COMMA typename A, std::map<K TIL_COMMA T TIL_COMMA P TIL_COMMA A>);
TIL_DECLARE_IS_SPEC_T(container, typename T TIL_COMMA typename BaselinePolicy, til::sparse_vector<T TIL_COMMA BaselinePolicy>);

//////////////////// is_map ////////////////////

/// True if type is a map container
TIL_DECLARE_IS(map);
TIL_DECLARE_IS_SPEC_T(map, typename TKey TIL_COMMA typename TValue TIL_COMMA typename TCompare TIL_COMMA typename TAlloc, std::map<TKey TIL_COMMA TValue TIL_COMMA TCompare TIL_COMMA TAlloc>);
TIL_DECLARE_IS_SPEC_T(map, typename TKey TIL_COMMA typename TValue TIL_COMMA typename TCompare TIL_COMMA typename TAlloc, std::multimap<TKey TIL_COMMA TValue TIL_COMMA TCompare TIL_COMMA TAlloc>);

//////////////////// is_pointer ////////////////////

TIL_DECLARE_IS(pointer)
TIL_DECLARE_IS_SPEC_T(pointer, typename T, T*);
TIL_DECLARE_IS_SPEC_T(pointer, typename T, boost::shared_ptr<T>);


// undef all local macros
#undef TIL_DECLARE_IS
#undef TIL_DECLARE_IS_VALUE
#undef TIL_DECLARE_IS_TEMPLATE
#undef TIL_DECLARE_IS_VALUE_TEMPLATE
// To allow others to use it
//#undef TIL_DECLARE_IS_SPEC_T
#undef TIL_DECLARE_IS_SPEC
#undef TIL_COMMA

} // namespace til

#endif //_ISTRAITS_H_

