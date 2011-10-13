#ifndef TIL_DYNIMAGEWRAPPERS_H
#define TIL_DYNIMAGEWRAPPERS_H

// include from TIL library
#include "til/til_common.h"
#include "til/DynImage.h"
#include "til/imageArith.h"
#include "til/templateTools.h"

// ignore warnings
//#pragma warning ( push )
//#pragma warning ( disable: C4244 ) // conversion from X to Y, possible loss of data

namespace til
{
	// predeclarations
	//template < class T > class Ptr;

	/*
	template < class TTypeCollection, typename T >
	void
	mul2(Dyn<TTypeCollection> dyn, T value)
	{
		template < typename T >
		struct Temp
		{
			void operator()(){ if (T object = dynamic_cast<T>(m_object)) mul(object}
		}
		mpl::for_each<TTypeCollection>(
	}
	*/
/*

#define DYNIM_FUNC_WRAP_DYNIM_VALUE(funcname)                         \
template < template <typename> class TImage , typename T>             \
typename disable_if<is_Image<T> >::type                               \
funcname (DynImage<TImage> & dynim, T value)                          \
{                                                                     \
	TIL_DYNIM_BLOCKS( detail::funcname##_iv (pim->image(), static_cast<typename change_precision<T,PixelType>::type>(value)); return;) 	\
}

  DYNIM_FUNC_WRAP_DYNIM_VALUE(add);
  DYNIM_FUNC_WRAP_DYNIM_VALUE(sub);
  DYNIM_FUNC_WRAP_DYNIM_VALUE(mul);
  DYNIM_FUNC_WRAP_DYNIM_VALUE(div);

#undef DYNIM_FUNC_WRAP_DYNIM_VALUE


#define DYNIM_FUNC_WRAP_2DYNIM(function)												\
template < template <typename> class TImage>											\
void function (DynImage<TImage> & dynim, DynImage<TImage> & dynim2)						\
{																						\
  TIL_DYNIM_BLOCKS( (DynImage_typed<TImage TIL_COMMA PixelType> *pim2 = dynamic_cast<DynImage_typed<TImage TIL_COMMA PixelType>*>(&dynim2) ) ; detail::function##_ii (pim->image(), pim2->image()); return;)					\
}

// TODO: right now, dynim is not const because the dynamic_cast assumes it is const.
// Change that (maybe using some boost traits or something).
#define DYNIM_FUNC_WRAP_IM_DYNIM(function)										\
template < template <typename> class TImage , typename TImage2>					\
typename enable_if<is_Image<TImage2> >::type									\
function (TImage2 & im2, DynImage<TImage> &dynim)								\
{																				\
  TIL_DYNIM_BLOCKS( detail::function##_ii (im2, pim->image()); return; )		\
}

#define DYNIM_FUNC_WRAP_DYNIM_IM(function)										\
template < template <typename> class TImage , typename TImage2>					\
typename enable_if<is_Image<TImage2> >::type									\
function (DynImage<TImage> & dynim, const TImage2 & im2)						\
{																				\
  TIL_DYNIM_BLOCKS( detail::function##_ii (pim->image(), im2); return; )		\
}

#define DYNIM_FUNC_WRAP_TWO_IMAGES(function) \
        DYNIM_FUNC_WRAP_IM_DYNIM(function)   \
        DYNIM_FUNC_WRAP_DYNIM_IM(function)   \
        DYNIM_FUNC_WRAP_2DYNIM(function)

#define TIL_COMMA ,
*/

/*
	template < template <typename> class TImage , typename TImage2>
void mul (TImage2 im2, DynImage<TImage> & dynim)
{                                                                 
  if(0) {}
  else if (DynImage_typed<TImage, char> *im = dynamic_cast<DynImage_typed<TImage, char >*>((DynImage<TImage>*)dynim))
  {
    typedef char PixelType;
    mul(im2, im.image()); return;
  }
  else { throw std::runtime_error("Unknown dynamic image type"); }
}
*/

// Well, we know that at this point, we can't control which type is multiplied by which.
#pragma warning ( push )
#pragma warning ( disable : 4244 ) // conversion from X to Y, possible loss of data

  /*
  DYNIM_FUNC_WRAP_TWO_IMAGES(mul);
  */
  
#pragma warning ( pop )

#undef DYNIM_FUNC_WRAP_TWO_IMAGES
#undef DYNIM_FUNC_WRAP_DYNIM_IM
#undef DYNIM_FUNC_WRAP_IM_DYNIM
#undef DYNIM_FUNC_WRAP_2DYNIM


#undef TIL_COMMA

} // namespace til

//#pragma warning ( pop )


#endif

