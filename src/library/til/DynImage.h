#ifndef TIL_DYNIMAGE_H
#define TIL_DYNIMAGE_H

// include from STL
#include <typeinfo>

// include from BOOST
#include "boost/mpl/assert.hpp"
#include "boost/mpl/equal.hpp"
#include "boost/mpl/int.hpp"
#include "boost/mpl/placeholders.hpp"
#include "boost/mpl/transform.hpp"
#include "boost/mpl/vector.hpp"

// include from TIL library
#include "til/til_common.h"
#include "til/templateTools.h"

namespace til {

	/// Base class for dynamic type images.
	/// Note that this class is templated over an image container
	/// type. I.e., this class only brings the "dynamic typing" 
	/// functionality to an existing, static-typed image container class.
	/// Note: we have to template over something at this level, because otherwise
	/// we would have to know all possible images. This means that we would have
	/// to store that somewhere. Not good for extensability, and for performance.
	/// I would love to template over a general TImage, that would accept TImage's with
	/// any number of templateing arguments, and then use traits to decide how to construct
	/// an image of a specific type. Unfortunately C++ forbid this kind of templation: we
	/// have to know how many template arguments a template parameter accepts.
	/// Note on the overall strategy taken here for dynamic images: this is a sort of
	/// "top-down" approach. The DynImage never intend to mimic the behavior of a standard
	/// image so that algorithms could use it as a fully compliant image parameter and
	/// use it transparently. As a matter of fact, DynImage cannot be passed to any
	/// normal image function because it does not implement any image feature, like
	/// indexing, iterators, or even like returning the value of one of its element!
	/// Instead, the top-down approach consists in deciding at the earliest level possible
	/// which image we are dealing with, and calling the appropriate function. This is
	/// done with high-level wrappers like 'Dyn', that transform a normal image processing
	/// functor into a dynamic image processing functor, basically by simply adding
	/// switches.
	/// The reason for this is efficiency. The earliest we know what type we are dealing
	/// with, the faster. Just imagine if we have to check for the image type everytime
	/// we access one of its element! This would be the case if we had an iterator for
	/// dynamic images
	/// The bad side of the top down approach is that things get complicated when using
	/// more than one dynamic images. There is an unavoidable combinatoric explosion.

	template < template < typename > class TImage >
	class DynImage : public Image_label
	{
	public: // constructors & destructor
		virtual ~DynImage() {};
	};

	

	/// Child class of dynamic image base class, implementing
	/// a specific type.
	template < template < typename > class TImage, typename T >
	class DynImage_typed : public DynImage<TImage>
	{
	public:
		DynImage_typed(TImage<T> & im) : m_im(im) {}
	public:
		TImage<T> & image() { return m_im; }

	private:
		TImage<T> & m_im;
	};

	/*
	template < template < typename > class TImage, typename value_type >
	DynImage<TImage> * dynImageFactory(TImage<value_type> *im)
	{
		return new DynImage_typed<TImage, value_type>(im);
	}
	*/

	/*
	template < template <typename> class TImage, typename value_type>
	DynImage_typed<TImage, value_type> *
	typecast(DynImage<TImage> * dynim)
	{
		return dynamic_cast<DynImage_typed<TImage, value_type>*>((DynImage<TImage>*)dynim));
	}
	*/

	// Note that the first template parameter is redundant and thus theoretically
	// not needed, however I could not find a C++ way to deduce it from the other :(
	template < template <typename> class TImage, class TImage0 >
	DynImage<TImage> *
	dynImageFactory(TImage0 & im)
	{
		return new DynImage_typed<TImage, typename TImage0::value_type>(im);
	}

/// Macro for code to do dynamic type resolution.
#define TIL_DYNIM_BLOCKS(funcall)  \
TIL_DYNIM_BLOCKS_BEGIN             \
TIL_DYNIM_BLOCK(funcall, uchar)    \
TIL_DYNIM_BLOCK(funcall, char)     \
TIL_DYNIM_BLOCK(funcall, ushort)   \
TIL_DYNIM_BLOCK(funcall, short)    \
TIL_DYNIM_BLOCK(funcall, int)      \
TIL_DYNIM_BLOCK(funcall, float)    \
TIL_DYNIM_BLOCK(funcall, double)   \
TIL_DYNIM_BLOCKS_END               \

#define TIL_DYNIM_BLOCKS_BEGIN if(0) {}
#define TIL_DYNIM_BLOCKS_END else { throw std::runtime_error("Unknown dynamic image type"); }
#define TIL_DYNIM_BLOCK(funcall, type)                                                                        \
else if (DynImage_typed<TImage, type > *pim = dynamic_cast<DynImage_typed<TImage, type >*>(&dynim))  \
{                                                                                                             \
  typedef type PixelType;                                                                                     \
  funcall                                                                                                     \
}                                                                                                             \



	/// Functor wrapper for dynamic type images.
	/// This functor helps reusing functors that were written for
	/// static-type images with dynamic type images. Whenever Functor was
	/// used before, simply use Dyn<Functor>. If you expect a return value,
	/// use Dyn<Functor, TReturn>. The constraint on the Functor class here is
	/// that a call to its operator() is enough to call all parts of
	/// the functor that need to process the image (e.g. there is not Functor::init(im)).
	/// If not, you will have to create your own functor wrapper.
	template < typename Functor, typename TReturn = void >
	class Dyn : public Functor
	{
	public: // constructors & destructor

		Dyn(const Functor &functor) : Functor(functor) {}

	public: // operators

		template < template < typename > class TImage >
		TReturn operator()(DynImage<TImage> * dynim)
		{
			TIL_DYNIM_BLOCKS( return m_functor(dynim.image()); )
		}
	};

	/// Conveniance function to create dynamic functor wrapper with return type
	template < typename TReturn, typename Functor >
	Dyn<Functor, TReturn> dyn(const Functor &functor)
	{
		return Dyn<Functor, TReturn>(functor);
	}

	/// Conveniance function to create dynamic functor wrapper witout return type
	template < typename Functor >
	Dyn<Functor, void> dyn(const Functor &functor)
	{
		return Dyn<Functor, void>(functor);
	}

	/*
	template < int N >
	struct MyDynCollection {};
	
	template <> struct MyDynCollection<0> { static const int N  = 0; typedef Dyn<ImageC<int> >		Type; };
	template <> struct MyDynCollection<1> { static const int N  = 1; typedef Dyn<ImageC<short> >	Type; };
	template <> struct MyDynCollection<2> { static const int N  = 2; typedef Dyn<ImageC<double> >	Type; };

	template < class TDynCollection >
	struct DynCollectionTraits {};

	template <> struct DynCollectionTraits<MyDynCollection> { static int N = 3; };
	*/

	class DynBase {};
	
	template < typename T >
	class DynTyped : public DynBase
	{
	public: // typedefs
		typedef T Type;
	public: // constructors & destructor
		DynTyped(T elem) : m_elem(elem) {}
	private: // data
		T m_elem;
	};
	
	typedef boost::mpl::vector<
		unsigned char,
		char,
		unsigned short,
		short,
		int,
		float,
		double
	> StandardNumericTypes;

	typedef boost::mpl::vector<
		ImageC<unsigned char>,
		ImageC<char>,
		ImageC<unsigned short>,
		ImageC<short>,
		ImageC<int>,
		ImageC<float>,
		ImageC<double>
	> imeu;


	typedef boost::mpl::transform<StandardNumericTypes, ImageC<boost::mpl::_1> >::type toto;

	//boost::mpl::at<toto, boost::mpl::int_<0> >::type i(8);

	BOOST_MPL_ASSERT(( boost::mpl::equal<imeu,toto> ));


	/*
	template < typename T >
	class DynamicCast
	{
		T *im = dynamic_cast<T*>()
		DynImage_typed<TImage, type > *im = dynamic_cast<DynImage_typed<TImage, type >*>((DynImage<TImage>*)dynim)

	};
	*/


	/*
	template < int TN, template <int> typename TFunctor >
	class Loop
	{
	public:
		static const int N = TN;
		static void execute()
		{
			TFunctor<N>();
			Loop<N-1>::execute();
		}
	};

	template <> class loop<0>
	{
	public:
		static void execute() {}
	};
	*/

	/*
	template <class TClass, class TFunctor>
	class ExecuteIfIsA
	{
	public:
		void operator()
		{
			if (
		}
	};
	*/


	template < class TMPLContainer >
	class Loop
	{
	public: // typedefs
		typedef Loop<TMPLContainer> Self;
		typedef typename boost::mpl::begin<TMPLContainer>	begin;
		typedef typename boost::mpl::end<TMPLContainer>		end;

	public: // static functions

		/// Execute a templated functor.
		/// The type is passed here as the template parameter of the functor class.
		/// This is a memoryless process, as a new functor will be created everytime.
		/// Therefore, this function takes no functor parameter
		template < template <typename> class TFunctor >
		static void execute()
		{
			Self::template execute<Self::begin, TFunctor>();
		}

		/// Execute a functor.
		/// The type is passed here as the template parameter of the operator() of the functor class.
		/// This method has memory of history, thus a functor can be passed as an argument
		template < class TFunctor >
		static void execute(const TFunctor &functor = TFunctor())
		{
			Self::template execute<Self::begin, TFunctor>(functor);
		}

	private:

		template < class TMPLIterator, template < typename > class TFunctor >
		static void execute()
		{
			// Stop here if we reached the end of the container
			if (boost::is_same<typename TMPLIterator::type, typename boost::mpl::end<TMPLContainer>::type>::value) return;

			// Call functor with current type as its template parameter
			TFunctor<typename TMPLIterator::type> functor;
			functor();
			
			// go to next element
			execute<typename boost::mpl::next<TMPLIterator>::type, TFunctor>(functor);
		}
		template < class TMPLIterator, class TFunctor >
		static void execute(const TFunctor &functor)
		{
			// Stop here if we reached the end of the container
			if (boost::is_same<typename TMPLIterator::type, typename Self::end::type>::value) return;

			// Call functor with current type as its template parameter
			functor.template operator()<typename TMPLIterator::type>();
			
			// go to next element
			execute<typename boost::mpl::next<TMPLIterator>::type, TFunctor>(functor);			
		}
	};

	template < typename TFunctor, typename T >
	class ExecuteIfDynCast
	{
	public:

		template < typename V >
		void operator()()
		{
			if (V object = dynamic_cast<V>(m_object))
			{
				m_functor(object);
			}
		}
	
	private: // data

		T m_object;
		TFunctor m_functor;
	};


	//boost::mpl::for_each<TypeCollection>(ExecuteIfDynCast(object, functor));

	template < class TTypeCollection, class TFunctor >
	void dyn(const TFunctor & functor)
	{
		Loop<TTypeCollection>::execute(functor);
	}


	/// A generic object.
	class Object
	{
		//virtual void operator()() = 0;

		/*
		template < template < typename, typename > class TBinaryFunctor >
		virtual TBinaryFunctor<Object, Object> makeFunctor(TBinaryFunctor f, Object i);

		template < class TFunctor, TTypeCollection >
		virtual WrapFunctor<TFunctor, TTypeCollection>
		peelFunctor(WrapFunctor<TFunctor, TTypeCollection>);
		*/
	};
	
	/// A concrete instance of a generic object, containing a real object of type T.
	template < typename T >
	class Object_typed
	{
		//virtual void operator()() { this->operator(); }

		/*
		template < template < typename, typename > class TBinaryFunctor >
		virtual TBinaryFunctor<Object, Object> makeFunctor(TBinaryFunctor f, Object i)
		{
			return i.makeFunctor(f, *this);
		}

		template < template < typename, typename > class TBinaryFunctor, typename U >
		virtual TBinaryFunctor<Object, Object> makeFunctor(TBinaryFunctor f, Object_typed<U> i)
		{
			return TBinaryFunctor<Object_typed<U>, Object_typed<T> >(f);
		}

		template < class TFunctor, TTypeCollection >
		virtual WrapFunctor<TFunctor, TTypeCollection>
		peelFunctor(WrapFunctor<TFunctor, TTypeCollection>)
		{
		}
		*/
	};

	//wrapFunctor(Mul)(i,j,k);

/*
	template < typename TFunctor, typename TTypeCollection = None >
	class WrapFunctor
	{
		operator(Object i)
		{
			i->peelFunctor(*this);	
		};

		TFunctor m_functor;
	};


	mul(Image i, Image j)
	{
		any f = FunctorFactory(Mul, i, j);
		f(i,j);
	}

	mul(Object i, Object j)
	{
		GenericBinaryFunctor f = binaryFunctorFactory(Mul, i, j);
		f(i,j);
	};

	template < typename TBinaryFunctor >
	binaryFunctorFactory(TBinaryFunctor f, Object i, Object j)
	{
		return i->makeFunctor(f, j);
	};
*/
	/*
	template < template <int> class TDynCollection, typename TFunctor >
	void dyn(const TFunctor & functor)
	{
		Loop<DynCollectionTraits<TDynCollection>::N, DynCompareAndExecute>::execute();
	}
	*/

//#undef TIL_DYNIM_BLOCKS
//#undef TIL_DYNIM_BLOCK
//#undef TIL_DYNIM_BLOCK_BEGIN
//#undef TIL_DYNIM_BLOCK_END
} //namespace

#endif


