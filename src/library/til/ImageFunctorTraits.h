#ifndef TIL_IMAGE_FUNCTOR_TRAITS_H
#define TIL_IMAGE_FUNCTOR_TRAITS_H

// include from TIL library
#include "til/til_common.h"
#include "til/NeighborhoodConfigurations.h"
#include "til/TExprMacros.h"

namespace til		{
namespace expr		{
namespace functor	{

	/// Traits for image functors.
	/// Image functors are functors designed for images only. I.e. they
	/// concern things less trivial that simple intensity transformation.
	/// For example, they could use the value of their neighbors.
	/// Therefore, the first difference between a numerical functor
	/// and an image functor is that an image functor takes image
	/// iterators as an input.
	/// Still -- image iterators are something more general and independant
	/// of template expressions, which is why they are separated here and
	/// need these traits to embed them in template expressions.
	/// Actually, again, we need these traits only to get the type of the
	/// returned value. In other words, to compensate for the lack
	/// of typeof...
	// TODO: do we really need to distinguish FunctorTraits and
	// ImageFUnctorTraits?

	//template < typename ImageFunctor > struct ImageFunctorTraits;
	
	template < class Functor > struct FunctorTraits;

	template < class _TPixel, class _TNeighborhood >
	struct FunctorTraits<IsIsolated<_TPixel, _TNeighborhood> >
	{
		template < typename T > struct TypeStruct
		{
			typedef bool Type;
		};
	};

} // namespace functor
} // namespace expr
} // namespace til


#endif

