#ifndef TIL_TEMPLATE_TOOLS_H
#define TIL_TEMPLATE_TOOLS_H

/// \file
/// Collects template tools used for library implementation.
/// Do not have any interest whatsoever to a library user.

// includes from BOOST library
#include "boost/type_traits.hpp"

// includes from TIL library
#include "til/til_common.h"

/// Apply a given macro for numeric types
#define TIL_FOR_ALL_NUMERIC_TYPES(macro)    \
  macro(bool)                               \
  macro(unsigned char)                      \
  macro(char)                               \
  macro(unsigned short)                     \
  macro(short)                              \
  macro(int)                                \
  macro(unsigned int)                       \
  macro(long)                               \
  macro(unsigned long)                      \
  macro(float)                              \
  macro(double)                             \
  macro(long double)                        \


namespace til
{

  //------------------------------------------------------------------------------------------------

    //-------------------------//
   //  true_type, false_type  //
  //-------------------------//
  
  struct true_type  { static bool const value = true; };
  struct false_type { static bool const value = false; };
  
  //------------------------------------------------------------------------------------------------

    //-------------//
   //  bool_type  //
  //-------------//

  template < bool B >
  struct bool_type : public false_type {};
  template < >
  struct bool_type<true> : public true_type {};

  //------------------------------------------------------------------------------------------------
  
    //-----------//
   //  type_if  //
  //-----------//

	/// Template to choose between one type or the other depending on some condition
	template < bool B, typename TypeIfTrue, typename TypeIfFalse >
	struct type_if { typedef TypeIfTrue type; };
	template < typename TypeIfTrue, typename TypeIfFalse >
	struct type_if<false, TypeIfTrue, TypeIfFalse> { typedef TypeIfFalse type; };


  //------------------------------------------------------------------------------------------------

    //-----------//
   //  is_same  //
  //-----------//

  template < typename T1, typename T2 >
  struct is_same : public false_type {};
  template < typename T >
  struct is_same<T,T> : public true_type {};
  

  //------------------------------------------------------------------------------------------------

    //--------------//
   //  naked_type  //
  //--------------//

	template < typename T >
	struct naked_type : public boost::remove_const<typename boost::remove_reference<T>::type> {};

  //------------------------------------------------------------------------------------------------

	/// Test whether two types are basically the same as far as rvalue goes, i.e are the
	/// same modulo possible const and reference differences.
	template < typename T1, typename T2 >
	struct is_same_rvalue_type
		: public is_same
      <
			  typename naked_type<T1>::type,
			  typename naked_type<T2>::type
		  >			   
	{};

  //------------------------------------------------------------------------------------------------

    //-----------------//
   //  is_same_naked  //
  //-----------------//

  template < typename T1, typename T2 >
  struct is_same_naked : public boost::is_same<typename naked_type<T1>::type, typename naked_type<T2>::type> {};

  //------------------------------------------------------------------------------------------------


    //-------------------------//
   //  enable_if, disable_if  //
  //-------------------------//


	/// Template tool to enable a template specialization under some condition
	template <bool B, class T = void>
	struct enable_if_c {
		typedef T type;
	};
	template <class T>
	struct enable_if_c<false, T> {};
	template <class Cond, class T = void>
	struct enable_if : public enable_if_c<Cond::value, T> {};

	template <bool B, class T = void>
	struct disable_if_c {
		typedef T type;
	};
	template <class T>
	struct disable_if_c<true, T> {};
	template <class Cond, class T = void>
	struct disable_if : public disable_if_c<Cond::value, T> {};


  //------------------------------------------------------------------------------------------------

  /*
    //---------------//
   //  int_to_type  //
  //---------------//

	/// A template tool to overcome the lack of partial specialization for functions
	template <int N>
	struct int_to_type
	{
		static int const value = N;
	};
  */

  //------------------------------------------------------------------------------------------------

    //--------//
   //  prec  //
  //--------//

  /// A dummy class used to pass a precision type for computations to some functions.
  template < typename T >
  struct prec { typedef T type; };

  //------------------------------------------------------------------------------------------------

    //--------//
   // check  //
  //--------//

  /// Fails to compile if b is false.
  /// Offers a handy short-cut in in-class checks to the "typename enable_if_c<>::type".
  
  template < bool b >
  struct check { typedef typename enable_if_c<b>::type type; };

} // namespace til


#endif
