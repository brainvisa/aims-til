#ifndef TIL_TEXPR_NUMERICAL_FUNCTORS_H
#define TIL_TEXPR_NUMERICAL_FUNCTORS_H

/// \file Detemplated functors for template expressions. Belongs to TExpr package.
/// Do not include directly, include til/TExpr.h instead.

// includes from STL
#include <functional>
#include <cmath>

// includes from BOOST
#include "boost/type_traits.hpp"
#include "boost/call_traits.hpp"

// includes from TIL
#include "til/cat2type.h"
#include "til/functors.h"
#include "til/labels.h"


namespace til { namespace expr { namespace functor {

  //--------------------------------------------------------------------------------------------
  
    //-----------------------//
   //  DetemplateOperator1  //
  //-----------------------//

  /// A class to "detemplate" a stateless operator with a single template argument.
	/// I.e. to transform a F<T> into a non-templated F that has an operator()<T>.
	/// Detemplation is necessary for template expressions: functors involved in
	/// template expressions are not type-specific.
	/// Note that to be recognized as a Genuine Detemplated Functor, a class
	/// must (1), have a templated TypeStruct, (2), have (inherit from) the label 
  /// detemplatedFunctor.
	/// Actually, (2) comes from the fact that I don't know any mean to
	/// automatically check for (1)  :(
	// NB: the use of "operator" is to emphasize this functor should be
	// stateless.
	template < template <typename> class TOperator >
	struct DetemplateOperator1 : public detemplated_functor_label
	{
		template < typename T, typename TUnused = T > struct TypeStruct
		{
			typedef typename TOperator<T>::result_type Type;
		};

		template < typename T >
		//typename boost::enable_if<
		//	boost::is_same<typename TOperator<T>::argument_type, T>
		//	,
			typename TOperator<T>::result_type
		//>::type
		operator()(const T & x) const
		{
			return TOperator<T>()(x);
		}

		
		template < typename T >
		//typename boost::enable_if_c<
			//boost::is_same<typename TOperator<T>::first_argument_type, T>::value
			//&& boost::is_same<typename TOperator<T>::second_argument_type, const T &>::value
			//,
			typename TOperator<T>::result_type
		//>::type
		//operator()(const T & x, const T & y) const
		operator()(const T & x, const T & y) const
		{
			return TOperator<T>()(x,y);
		}
		
		/*
		/// This is for assignments
		template < typename T >
		//typename boost::enable_if_c<
			//boost::is_same<typename TOperator<T>::first_argument_type, T&>::value
			//&& boost::is_same<typename TOperator<T>::second_argument_type, T>::value
			//,
			typename TOperator<T>::result_type
		//>::type
		operator()(T & x, const T & y) const
		{
			return TOperator<T>()(x,y);
		}
		*/
	};

  //--------------------------------------------------------------------------------------------
  
    //-----------------------------//
   //  DetemplateAssignOperator1  //
  //-----------------------------//

  /// Detemplation of assignment functors taking one template parameter.
	template < template <typename> class TOperator >
	struct DetemplateAssignOperator1 : public detemplated_functor_label
	{
		template < typename T, typename TUnused = T > struct TypeStruct
		{
			typedef typename TOperator<T>::result_type Type;
		};

		template < typename T >
		//typename boost::enable_if_c<
			//boost::is_same<typename TOperator<T>::first_argument_type, T&>::value
			//&& boost::is_same<typename TOperator<T>::second_argument_type, T>::value
			//,
			typename TOperator<T>::result_type
		//>::type
		operator()(T & x, const T & y) const
		{
			return TOperator<T>()(x,y);
		}
		
	};

  //--------------------------------------------------------------------------------------------
  
    //-----------------------//
   //  DetemplateOperator2  //
  //-----------------------//


	/// A class to "detemplate" a pure functor with two template arguments.
	/// i.e. to transform a F<T> into a F that has an operator()<T>.
	template < template <typename, typename> class TOperator >
	struct DetemplateOperator2 : public detemplated_functor_label
	{
		template < typename T1, typename T2 > struct TypeStruct
		{
			typedef typename TOperator<T1,T2>::result_type Type;
		};

		template < typename T1, typename T2 >
		//typename boost::enable_if_c<
		//boost::is_same<typename TOperator<T1,T2>::first_argument_type, T1>::value &&
		//boost::is_same<typename TOperator<T1,T2>::second_argument_type, T2>::value
		//	,
			typename TOperator<T1,T2>::result_type
		//>::type
		operator()(const T1 & x, const T2 & y) const
		{
			return TOperator<T1,T2>()(x,y);
		}

		template < typename T1, typename T2 >
		typename TOperator<T1,T2>::result_type
		operator()(T1 & x, const T2 & y) const
		{
			return TOperator<T1,T2>()(x,y);
		}
	};


  //--------------------------------------------------------------------------------------------
  
    //----------------//
   //  DefaultThird  //
  //----------------//


  /// A class to replace the third (result) template parameter of a functor
  /// by a default value. Necessary for detemplation of these functors...
  /// C++  >:(
  template < template <typename,typename,typename> class T >
  struct DefaultThird
  {
    template < typename T1, typename T2 >
    struct type : public T<T1,T2,typename combine<T1,T2>::type> {};
  };


  //--------------------------------------------------------------------------------------------

    //------------------------------------//
   //  Detemplated functor declarations  //
  //------------------------------------//

  typedef DetemplateOperator1<std::negate>          Negate;
  typedef DetemplateOperator1<til::functor::Abs>		Abs;
	typedef DetemplateOperator1<til::functor::Sqrt>		Sqrt;

  /*
  typedef DetemplateOperator2<til::functor::Add>		Plus;
  typedef DetemplateOperator2<til::functor::Sub>		Minus;
  typedef DetemplateOperator2<til::functor::Mul>		Multiplies;
  typedef DetemplateOperator2<til::functor::Div>		Divides;
  */

  typedef DetemplateOperator2<DefaultThird<til::functor::Add>::type>		Plus;
  typedef DetemplateOperator2<DefaultThird<til::functor::Sub>::type>		Minus;
  typedef DetemplateOperator2<DefaultThird<til::functor::Mul>::type>		Multiplies;
  typedef DetemplateOperator2<DefaultThird<til::functor::Div>::type>		Divides;

  // NB: I think this is to be used only when the functor is itself a template expression,
  // otherwise bind should be used (obviously, since a(expr) is likely to fail).
  // Right now the functor is passed by value, it might be a problem, I dunno.
  typedef DetemplateOperator2<til::functor::Call>   Call;

	//typedef DetemplateOperator2<til::functor::AddTo>	AddTo;
	//typedef DetemplateOperator2<til::functor::SubTo>	SubTo;
	//typedef DetemplateOperator2<til::functor::MulTo>	MulTo;
	//typedef DetemplateOperator2<til::functor::DivTo>	DivTo;
	typedef DetemplateOperator2<til::functor::AddTo>	AddTo;
	typedef DetemplateOperator2<til::functor::SubTo>	SubTo;
	typedef DetemplateOperator2<til::functor::MulTo>	MulTo;
	typedef DetemplateOperator2<til::functor::DivTo>	DivTo;

  typedef DetemplateOperator2<til::functor::CastTo> CastTo;

	typedef DetemplateOperator1<std::equal_to>			   Equal_To;
	typedef DetemplateOperator1<std::not_equal_to>		 Not_Equal_To;
	typedef DetemplateOperator1<std::greater>			     Greater;
	typedef DetemplateOperator1<std::less>				     Less;
	typedef DetemplateOperator1<std::greater_equal>		 Greater_Equal;
	typedef DetemplateOperator1<std::less_equal>		   Less_Equal;

	//typedef DetemplateOperator1<til::functor::Assign>	Assign;
	typedef DetemplateAssignOperator1<til::functor::Assign>	Assign;
	typedef DetemplateOperator1<til::functor::Deref>	Deref;



  /// Cast functor
	template < typename TTo >
	struct Cast
	{
		template < typename TFrom >
		struct TypeStruct
		{
			typedef TTo Type;
		};

		template < typename TFrom > TTo operator()(TFrom a) const
		{
			return TTo(a);
		}
	};


  //--------------------------------------------------------------------------------------------

    //--------//
   //  Wrap  //
  //--------//
	
  namespace detail
  {
    template < typename T, bool b > 
    struct WrapReturnValue
    { typedef typename T::result_type type; };

    template < typename T >
    struct WrapReturnValue<T,false>
    { typedef void type; };        
  }

  /// Wrap a standard functor inside a detemplated functor.
  /// The operation at hand here is different from detemplation.
  /// Detemplation is a mean to reuse pure templated predicates as genuine
  /// template expression functors. The resulting detemplated functors are
  /// defined online, inside template expressions.
  /// Functor wrapping, on the other hand, deals with functors (usually, non-stateless
  /// functors) that are defined beforehand and that we wish to use in a template
  /// expression.
  // NB: Historical notes: Wrap inherits publicly from TFunctor to solve two problems: (1) the fact that
  // operator() cannot be templated because if TFunctor requires a non-const reference argument, it better
  // be that way; (2), because of the first point, to avoid having to have different Wrap functions for
  // different number of arguments (one for unary, one for binary, etc.).
  template < typename TFunctor >
  class Wrap : public TFunctor
  {
  public: // typedefs
    template < typename T > struct TypeStruct
    {
      typedef typename detail::WrapReturnValue<TFunctor,has_result_type<TFunctor>::value>::type Type;
      //typedef typename TFunctor::result_type Type;
    };
  public: // typedefs
    Wrap(TFunctor & f) : TFunctor(f) {}
  };

}}} // namespace til::expr::functor


#endif

