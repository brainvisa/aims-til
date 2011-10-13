#ifndef TIL_LABEL_CLASSES_H
#define TIL_LABEL_CLASSES_H

/// \file Defines empty classes that serves as labels.
/// Labels may have multiple purpose. They can be used
/// to detect if a template argument is indeed of the expected type.
/// They can also help collecting similar object together and create
/// a hierarchy that help user browse through the library API (e.g. when
/// using an automatically generated documentation).
/// Basically, it replaces the standard inheritance scheme that is avoided
/// in template programming for efficiency reasons, replaces namespaces
/// to avoid the user having to jongle with potentially dozens of namespaces,
/// and help the compiler having readable compiler errors when substitution
/// fails.


// includes from TIL library
#include "til/templateTools.h"


namespace til
{

	// Predeclarations
	template < typename T > struct PointerTraits;

	// TODO: isn't it better to just have traits somewhere?
	// Especially since inheritance here is never used...


// This macro declares both a label class and the associated test class
// to test whether a class has this label or not.
// NB: to be undef at the end of the file
#define DEFINE_LABEL_CLASS(name) DEFINE_LABEL_CLASS_WITH_INHERITANCE(name, )

#define DEFINE_LABEL_CLASS_WITH_INHERITANCE(name, inherit)                        \
DEFINE_LABEL_CLASS_FULL(name##_label, Is##name)                                   \
DEFINE_TEST_CLASS_FULL(is_##name, Is##name, inherit)                              \

#define DEFINE_LABEL_CLASS_FULL(labelname, memtypename)                           \
class labelname { public: typedef void memtypename; };                            \


// Same as above, with possible inheritance
// NB: to be undef at the end of the file
// NB: using an enum instead of a static const bool brings in somw problems with template
// argument deduction, especially in CPPUNIT macros.
// NB: I used protected instead of private just to have gcc shut up that "all member functions 
// in class X are private"
#ifdef _MSC_VER
#define DEFINE_TEST_CLASS_FULL(testname, memtypename, inherit)                    \
template <typename T>                                                             \
class testname inherit                                                            \
{                                                                                 \
protected:                                                                        \
  typedef char Yes;                                                               \
  struct No { char a[2]; };                                                       \
  template <typename U> static Yes testfunc (typename U::memtypename const*);     \
  template <typename U> static No  testfunc (...);                                \
public:                                                                           \
  static const bool value = (sizeof(testfunc <T>(0)) == sizeof(Yes));			        \
};
#else
#define DEFINE_TEST_CLASS_FULL(testname, memtypename, inherit)                    \
template <typename T>                                                             \
class testname inherit                                                            \
{                                                                                 \
protected:                                                                        \
  typedef char Yes;                                                               \
  struct No { char a[2]; };                                                       \
  template <typename U> static Yes testfunc (typename U:: memtypename const*);    \
  template <typename U> static No  testfunc (...);                                \
public:                                                                           \
  enum {                                                                          \
    value = (sizeof(testname <T>::template testfunc <T>(0)) == sizeof(Yes))       \
  };                                                                              \
};
#endif

  DEFINE_TEST_CLASS_FULL(has_result_type, result_type, );

template < typename T >
class is_ImagePointer
{
protected:
	typedef char Yes;
	struct No { char a[2]; };
	template <typename U> static Yes isImagePointer(typename PointerTraits<U>::DataType::IsImagePointer const*);
	template <typename U> static No  isImagePointer(...);
public:
	enum {
	 	value = (sizeof(is_ImagePointer<T>::template isImagePointer<T>(0)) == sizeof(Yes))
	};
};
/*
template < typename T >
class is_detemplated_functor
{
private:
	typedef char Yes;
	struct No { char a[2]; };
	template
	template <typename U> static Yes zztest(typename U::TypeStruct<double>::type const*);
	template <typename U> static Yes zztest(typename U::TypeStruct<double,double>::type const*);
	template <typename U> static No  zztest(...);
public:
	/ *
	enum {
	 	value = (sizeof(is_detemplated_functor<T>::template zztest<T>(0)) == sizeof(Yes))
	};
	* /
	static const bool value = (sizeof(zztest<T>(0)) == sizeof(Yes));
};
*/

	// Definitions begin here

	/// Label for extrapolable images
	DEFINE_LABEL_CLASS(ExtrapolableImage);
	/// Label for all image classes
	DEFINE_LABEL_CLASS(Image);
	/// Label for image extrapolators
	DEFINE_LABEL_CLASS(ImageExtrapolator);
	/// Label for image functors
	DEFINE_LABEL_CLASS(ImageFunctor);
	/// Label for all image iterator classes
	DEFINE_LABEL_CLASS(ImageIterator);
	/// Label for interpolator of a sequence of number
	DEFINE_LABEL_CLASS(Interpolator);
	/// Label for coordinate transforms
	DEFINE_LABEL_CLASS(Mapping);
	/// Label for TIL smart pointers
	//DEFINE_LABEL_CLASS(SmartPtr);
	/// Label for all volumetric image iterator classes
	DEFINE_LABEL_CLASS_WITH_INHERITANCE(VolumetricImageIterator, : public ImageIterator_label);
	/// Label for neighborhoods
	//DEFINE_LABEL_CLASS(Neighborhood);
	DEFINE_LABEL_CLASS(detemplated_functor);
  /// Label for multidimensional minimization algorithms
  DEFINE_LABEL_CLASS(minimization_algorithm_nd);
	

class Neighborhood_label { typedef void IsNeighborhood; };
template <typename T>
class is_Neighborhood
{
protected:
	typedef char Yes;
	struct No { char a[2]; };
	template <typename U> static Yes isNeighborhood (typename U::IsNeighborhood const*);
	template <typename U> static No  isNeighborhood (...);
public:
	//static const bool value = (sizeof(isNeighborhood <T>(0)) == 1);
	enum {
#ifdef _MSC_VER
		value = (sizeof(isNeighborhood <T>(0)) == 1)
#else
		value = (sizeof(is_Neighborhood<T>::template isNeighborhood<T>(0)) == sizeof(Yes))
#endif

	};
};


#undef DEFINE_LABEL_CLASS
#undef DEFINE_LABEL_CLASS_WITH_INHERITANCE
} // namespace

#endif

