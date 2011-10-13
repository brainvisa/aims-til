#ifndef TIL_POINTER_TRAITS_H
#define TIL_POINTER_TRAITS_H

// includes from TIL library
#include "til/til_declarations.h"
#include "til/ImageBase.h"

namespace til {

	/// Pointer traits
	template < typename Pointer >
	struct PointerTraits {};

	/*
	/// Traits for ConstPtr
	template < typename T >
	struct PointerTraits<ConstPtr<T> >
	{
		typedef T DataType;
		static const bool is_const = true;
	};

	/// Traits for Ptr
	template < typename T >
	struct PointerTraits<Ptr<T> >
	{
		typedef T DataType;
		static const bool is_const = false;
	};
	*/
	/*
	/// Traits for linear image iterators
	template < typename TImage >
	struct PointerTraits<typename Iterator<TImage>::Linear>
	{
		typedef typename Iterator<TImage>::Linear::value_type TData;
		static const bool is_const = false;
	};

	/// Traits for const linear image iterators
	template < typename TImage >
	struct PointerTraits<typename Iterator<TImage>::ConstLinear>
	{
		typedef typename Iterator<TImage>::ConstLinear::value_type TData;
		static const bool is_const = true;
	};

	/// Traits for volumetric image iterators
	template < typename TImage >
	struct PointerTraits<typename Iterator<TImage>::Volumetric>
	{
		typedef typename Iterator<TImage>::Volumetric::value_type TData;
		static const bool is_const = false;
	};

	/// Traits for const volumetric image iterators
	template < typename TImage>
	struct PointerTraits<typename Iterator<TImage>::ConstVolumetric>
	{
		typedef typename Iterator<TImage>::ConstVolumetric::value_type TData;
		static const bool is_const = true;
	};
	*/

	/// Traits for const C pointerss
	template < typename T >
	struct PointerTraits<const T*>
	{
		typedef T DataType;
		static const bool is_const = true;
	};

	/// Traits for C pointers
	template < typename T >
	struct PointerTraits<T*>
	{
		typedef T DataType;
		static const bool is_const = false;
	};
}

#endif


