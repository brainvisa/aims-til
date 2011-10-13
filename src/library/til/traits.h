#ifndef TIL_TRAITS_H
#define TIL_TRAITS_H

// include from STL
#include <limits>

// includes from BOOST
#include "boost/call_traits.hpp"
#include "boost/type_traits.hpp"
#include "boost/utility/enable_if.hpp"

// includes from TIL library
#include "til/til_common.h"
#include "til/til_declarations.h"
#include "til/templateTools.h"

namespace til {

#define TIL_COMMA ,

  //------------------------------------------------------------------------------

    //----------------//
   //  require_that  //
  //----------------//
  
  template < bool B > class require_that;
  template <> class require_that<true> {};

  //------------------------------------------------------------------------------

    //-------------//
   //  precision  //
  //-------------//

	/// Numerical precision of the data for storage classes.
	/// Stores the precision of the class, i.e. the format in which ultimately
	/// the data is stored. E.g. precision<ImageC<double> >::type is double.
  /// precision<std::complex<float> > is float, not complex<float>.
	/// This allows to check whether data types are compatible, independently of
	/// the object itself. Mainly for performance reasons.
	/// E.g. when multiplying a matrix image with a constant vector, this could
	/// be used to check that the precision of matrices in the image and the constant
	/// vector are the same, to avoid to do a cast for every multiplication.
  /// Note that this is completely different from Container::value_type, which
  /// merely gives the type of the containee without cascading down to the encoding type.
	template < typename T >
	struct precision { };

	/// Specialization for numerical types
#define TIL_DEFINE_PRECISION_FOR_NUMERIC_TYPES(argtype)   \
  template <>                                             \
  struct precision< argtype >                             \
  { typedef argtype type; };

	TIL_FOR_ALL_NUMERIC_TYPES(TIL_DEFINE_PRECISION_FOR_NUMERIC_TYPES);

	// If you happen to have build another base numeric class, like double_double,
	// simply add a similar precision class. For all other non-fundamental classes,
	// and this includes complex number, rational fraction, etc., the precision should
	// reflect the encoding.

#undef TIL_DEFINE_PRECISION_FOR_NUMERIC_TYPES

	// A macro to define specialization of precision trait.
	// The template argument on which recursion is done has to be named 'T'.
	// NB: to be undefined at the end of the definitions
#define DEFINE_PRECISION_RECURSIVE_SPECIALIZATION_T(argtype, targs)   \
	template < targs >                                                  \
	struct precision< argtype >                                         \
	{ typedef typename precision<T>::type type; }                       \

#define DEFINE_PRECISION_RECURSIVE_SPECIALIZATION(argtype)            \
	DEFINE_PRECISION_RECURSIVE_SPECIALIZATION_T(argtype, typename T)

	// Specialization for images
	DEFINE_PRECISION_RECURSIVE_SPECIALIZATION(ImageC<T>);
	DEFINE_PRECISION_RECURSIVE_SPECIALIZATION(ImageNC<T>);
	DEFINE_PRECISION_RECURSIVE_SPECIALIZATION(ImageRLE<T>);

	// Specialization for vectors
  //DEFINE_PRECISION_RECURSIVE_SPECIALIZATION_T(Vector<T TIL_COMMA D TIL_COMMA TStorage>, typename T TIL_COMMA std::size_t D TIL_COMMA typename TStorage);

  // Specialization for points
  //DEFINE_PRECISION_RECURSIVE_SPECIALIZATION_T(Point<T TIL_COMMA D TIL_COMMA TStorage>, typename T TIL_COMMA std::size_t D TIL_COMMA typename TStorage);

  // Specialization for arrays
  DEFINE_PRECISION_RECURSIVE_SPECIALIZATION_T(numeric_array<T TIL_COMMA D>, typename T TIL_COMMA std::size_t D);
  DEFINE_PRECISION_RECURSIVE_SPECIALIZATION_T(sparse_vector<T TIL_COMMA BaselinePolicy>, typename T TIL_COMMA typename BaselinePolicy);

	// Specialization for matrix
	DEFINE_PRECISION_RECURSIVE_SPECIALIZATION(Matrix3<T>);
	DEFINE_PRECISION_RECURSIVE_SPECIALIZATION(SymMatrix3<T>);

#undef DEFINE_PRECISION_RECURSIVE_SPECIALIZATION

  //------------------------------------------------------------------------------

    //--------------------//
   //  change_precision  //
  //--------------------//

	/// Changing the numerical precision of an object.
	/// This allows to change the precision to a desired level of accuracy/compatibility.
	template < typename T, typename TNewPrecision >
	struct change_precision {};

	// A macro to define change_precision for numerical types.
	// NB: to be undef'd before the end of this file.
#define TIL_DEFINE_CHANGE_PRECISION_FOR_NUMERIC_TYPES(argtype)    \
	template < typename TNewPrecision >                             \
	struct change_precision< argtype, TNewPrecision >               \
	{                                                               \
		typedef TNewPrecision type;                                   \
	private:                                                        \
	  typedef typename enable_if_c<std::numeric_limits<TNewPrecision>::is_specialized>::type CheckIsNumeric;	\
	};                                                              \

	TIL_FOR_ALL_NUMERIC_TYPES(TIL_DEFINE_CHANGE_PRECISION_FOR_NUMERIC_TYPES);

#undef TIL_DEFINE_CHANGE_PRECISION_FOR_NUMERIC_TYPES


	// A macro to define change_precision recursively for non-numerical types.
	// NB: to be undef'd before the end of this file.
#define TIL_DEFINE_CHANGE_PRECISION_RECURSIVE_SPECIALIZATION_T(argtype, targs)  \
	template < targs , typename TNewPrecision>                                    \
	struct change_precision< argtype , TNewPrecision>                             \
	{                                                                             \
		typedef typename change_precision<T, TNewPrecision>::type type;             \
	private:                                                                      \
		typedef typename enable_if_c<std::numeric_limits<TNewPrecision>::is_specialized>::type CheckIsNumeric;	\
	};                                                                            \


#define TIL_DEFINE_CHANGE_PRECISION_RECURSIVE_SPECIALIZATION(argtype)			\
	TIL_DEFINE_CHANGE_PRECISION_RECURSIVE_SPECIALIZATION_T(argtype, typename T)


	// Specialization for images
	TIL_DEFINE_CHANGE_PRECISION_RECURSIVE_SPECIALIZATION(ImageC<T>);
	TIL_DEFINE_CHANGE_PRECISION_RECURSIVE_SPECIALIZATION(ImageNC<T>);
	TIL_DEFINE_CHANGE_PRECISION_RECURSIVE_SPECIALIZATION(ImageRLE<T>);

	// Specialization for vectors
  //TIL_DEFINE_CHANGE_PRECISION_RECURSIVE_SPECIALIZATION_T(Vector<T TIL_COMMA D TIL_COMMA TStorage>, typename T TIL_COMMA std::size_t D TIL_COMMA typename TStorage);

	// Specialization for matrix
	TIL_DEFINE_CHANGE_PRECISION_RECURSIVE_SPECIALIZATION(Matrix3<T>);
	TIL_DEFINE_CHANGE_PRECISION_RECURSIVE_SPECIALIZATION(SymMatrix3<T>);



  //------------------------------------------------------------------------------

    //-----------//
   //  combine  //
  //-----------//

	/// Defines the return type of an operation between numbers of different types.
	template < typename T1, typename T2 >
	struct combine
	{
		typedef typename type_if<
			// We choose T1 if T1 is floating point and T2 is not...
			(!std::numeric_limits<T1>::is_integer && std::numeric_limits<T2>::is_integer) ||
			// or, if not (T2 is floating and T1 is not) and...
			((!std::numeric_limits<T1>::is_integer || std::numeric_limits<T2>::is_integer) &&
			// if encoding of T1 is larger than T2
			(sizeof(T1) >= sizeof(T2))),
			T1, T2 >::type type;
		// NB: Note that we are in deep trouble if T1 and T2 have the same size and
		// are both either integer or floating point, because then when arbitrarily
		// chose T1... 
		// TODO: How is it done in the compiler?
	};

  //------------------------------------------------------------------------------

  /*
  template < typename TFunctor, typename T1, typename T2,  >
  struct result_of
  {
  };
  */

  //------------------------------------------------------------------------------
  
    //--------------//
   //  value_type  //
  //--------------//

	/// Externalization of the standard value_type member typedef.
	/// The value_type is actually a bit schizophrenic of a type. For containers,
	/// it is the type of the containees. For pointers, it is the type of the
	/// pointees. For iterators, it is the type of the iteratees.
	template < typename T >
	struct value_type_of
	{
		typedef typename T::value_type type;
	};

	/// Specialization for pointers.
	/// Note that this specialization is correct for all interpretations of pointers
	/// as a container, a pointer, or an iterator.
	template < typename T >
	struct value_type_of<T*>
	{
		typedef T type;
	};

  /// Specialization for references.
  /// I added this because sometimes, a functor takes a reference as an input (obviously when the input is large,
  /// such as a container). Now, should the argument_type of those functors be T&, or simply T? So far I chose to
  /// use T only, which eases those things but is not true, and the exact information of "is it T or T&" is lost
  /// when querying argument_type, which might be necessary later on. So I decided to once again change my mind and
  /// thus brought me to add this crap to ease things. The other solution would be to call value_typeof<naked_type>
  /// whereever one things either T or T& might be used. But it's kind of a nightmare if we have to think
  /// about that each time.
  template < typename T >
  struct value_type_of<T&>
  {
    typedef typename value_type_of<T>::type type;
  };

  /// Specialization for None type.
  template <>
  struct value_type_of<None>
  {
    typedef None type;
  };


  //------------------------------------------------------------------------------

    //-----------------//
   // deref           //
  //-----------------//

  /// Type return by T when unitary T::operator*() is called.
	template < typename T >
	struct deref
	{
		typedef typename T::reference type;
	};

	template < typename T >
	struct deref<T*>
	{
		typedef T & type;
	};

  //------------------------------------------------------------------------------

  /*
    //----------------//
   //  if_then_else  //
  //----------------//

  template < bool B, typename TTrue, typename TFalse >
  struct if_then_else
  {
    typedef TTrue type;
  };

  template < typename TTrue, typename TFalse >
  struct if_then_else<false, TTrue, TFalse >
  {
    typedef TFalse type;
  };
  */
  
  //------------------------------------------------------------------------------

  // Actually, can't work, because sometimes we want the return value like Vactor<T,3>, that doesn't
  // fit the template template parameter format.
  template < template <typename> class TReturn >
  struct unary_detemplated_functor
  {
    template < typename T >
    struct result_type
    {
      typedef TReturn<T> type;
    };
  };


  //------------------------------------------------------------------------------

    //------------//
   //  range_of  //
  //------------//

  template < typename T >
  struct range_of
  { typedef typename T::range type; };

  template < typename T >
  struct range_of<const T>
  { typedef typename T::const_range type; };

    //------------------//
   //  const_range_of  //
  //------------------//

  template < typename T >
  struct const_range_of
  { typedef typename T::const_range type; };


  //------------------------------------------------------------------------------

    //----------------//
   //  reference_of  //
  //----------------//

  /// Returns T::reference or T::const_reference, depending on the constness of T.
  template < typename T >
  struct reference_of
    : public type_if<boost::is_const<T>::value, typename T::const_reference, typename T::reference>
  {};

  //------------------------------------------------------------------------------

    //---------//
   //  opti   // (was: param)
  //---------//

  /// Return T is T is stateless^H^H^H^Hempty, otherwise boost::call_traits<T>::param_type
  template < typename T >
  struct opti
    //: public type_if<boost::is_stateless<T>::value, T, typename boost::call_traits<T>::param_type>
    : public type_if<boost::is_empty<T>::value, T, typename boost::call_traits<T>::param_type>
  {};


  //------------------------------------------------------------------------------
  
    //-------------//
   //  is_signed  //
  //-------------//

  template < typename T >
  struct is_signed
    : public type_if<std::numeric_limits<T>::is_signed, true_type, false_type>::type
  {};

  //------------------------------------------------------------------------------
  
    //--------------//
   //  is_integer  //
  //--------------//
  
  template < typename T >
  struct is_integer
    : public type_if<std::numeric_limits<T>::is_integer, true_type, false_type>::type
  {};

  //------------------------------------------------------------------------------
  
    //---------------------//
   //  is_floating_point  //
  //---------------------//

  // NB: small subtelty: we have to check that T is not an integer, but also that T is a number!
  template < typename T >
  struct is_floating_point
    : public type_if<std::numeric_limits<T>::is_specialized && !std::numeric_limits<T>::is_integer, true_type, false_type>::type
  {};

  //------------------------------------------------------------------------------

    //-------------------//
   //  accumulation_of  //
  //-------------------//

  template < typename T >
  struct accumulation_of {};
  
  template <>
  struct accumulation_of<float> { typedef double type; };

  //------------------------------------------------------------------------------

//#undef TIL_COMMA

} // namespace



#endif
