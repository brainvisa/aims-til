#ifndef TIL_TEXPR_H
#define TIL_TEXPR_H

/// \file
/// This file contains all the material a library user should need to
/// use template expressions.
/// This mainly comprises placeholders (_1, _2, _3) and the loop functions

// NB: I think the reason to rewrite placeholders was that standard placeholders would
// be a bit cumbersome to use with image iterators. More precisely, when we pass an image
// iterator to a functor like _1 += 10, we would like _1 to be actually interpreted as
// *(iterator), i.e. the value of this point. But we also would like to write more advanced
// stuff by binding the iterator to some image measure, say bind(isIsolatedPoint, _1), to test
// whether current pixel is isolated. Now, what is passed to this function is not
// the value of the pixel anymore, but truely the iterator, that has access to all of its
// neighbors.


// includes from STL
#include <algorithm>
#include <functional>

// includes from TIL library
#include "til/til_common.h"
#include "til/templateTools.h"
#include "til/TExprBasicFunctors.h"
#include "til/TExprConcatenation.h"
#include "til/TExprPlaceHolders.h"
#include "til/TExprMacros.h"


namespace til {

  /// namespace for template expressions.
	namespace expr {


  //---------------------------------------------------------------------------------------------------

    //-----------------//
   //  TExprConstant  //
  //-----------------//
	
  	
  	/// A template expression class for constant values.
    // TODO: I wonder if declaring a 'const T' brings any possible compiler optimization here... 
    // On the other hand, it might bring us some painful headaches when composed with other TExpr's?
    // Well, let's keep it till it breaks :/    
  	template < typename T >
  	class TExprConstant
  	{
  	public: // typedefs
  
  		EXPR_RESULT_TYPE(const T);
  
  	public: // constuctors
    
  		TExprConstant(T value) : m_value(value) {}
  
  	public: // operators
  		
      // I am not using the macros here to avoid painful warnings on the lines 'unused variable i1'
      template < class Iterator1 >
      typename TypeStruct<Iterator1>::Type
      operator()(Iterator1 &)
      { return m_value; }
      
      template < class Iterator1, class Iterator2>
      typename TypeStruct<Iterator1, Iterator2>::Type
      operator()(Iterator1 &, Iterator2 &)
      { return m_value; }
      
      template < class Iterator1, class Iterator2, class Iterator3 >
      typename TypeStruct<Iterator1, Iterator2, Iterator3>::Type
      operator()(Iterator1 &, Iterator2 &, Iterator3 &)
      { return m_value; }
      
      /*
  		EXPRFUNC_1ARG(operator(), return m_value;);
  		EXPRFUNC_2ARG(operator(), return m_value;);
  		EXPRFUNC_3ARG(operator(), return m_value;);
      */
      
  	private: // data
  
  		const T m_value;
  	};
  

  //---------------------------------------------------------------------------------------------------

    //---------//
   //  TExpr  //
  //---------//



#define TIL_DEFINE_TEXPR_ARITHMETIC_OPERATOR(name, functor)								      \
TIL_DEFINE_TEXPR_ARITHMETIC_OPERATOR_EXPR(name, functor)                        \
TIL_DEFINE_TEXPR_ARITHMETIC_OPERATOR_VALUE(name, functor)                       \

// NB: to be undefined after the definition of TExpr
#define TIL_DEFINE_TEXPR_ARITHMETIC_OPERATOR_EXPR(name, functor)								  \
template < typename Expr2 >															                          \
TExpr<TExprBinaryOperator_NoRes<Expr, Expr2, functor> >									          \
name (const TExpr<Expr2> &e) const													                      \
{																					                                        \
	typedef TExprBinaryOperator_NoRes< Expr, Expr2, functor > TExprRet;					    \
	return TExpr<TExprRet>(TExprRet(this->getExpr(), e.getExpr(), functor ()));		  \
}																					                                        \


#define TIL_DEFINE_TEXPR_ARITHMETIC_OPERATOR_VALUE(name, functor)								            \
template < typename T >															                                        \
TExpr<TExprBinaryOperator_NoRes<Expr, TExprConstant<T>, functor> >							            \
name (const T & value) const													                                      \
{																					                                                  \
	typedef TExprBinaryOperator_NoRes< Expr, TExprConstant<T>, functor > TExprRet;					  \
	return TExpr<TExprRet>(TExprRet(this->getExpr(), TExprConstant<T>(value), functor ()));		\
}																					                                                  \


  	/// A wrapper class of a template expression.
  	/// The reason for this is that we want all template expressions to be
  	/// variations of the same template class, so that manipulation with
  	/// operators is made easier, if not just plain possible.
  	template < class Expr >
  	class TExpr
  	{
  	public: // typedefs
  	
  		typedef void IsATExpr;
  
  		EXPR_RESULT_TYPE(typename Expr::template TypeStruct<Iterator1 TIL_COMMA Iterator2 TIL_COMMA Iterator3>::Type);
  
  	public: // constructors & destructor
  
  		explicit TExpr() : m_e() {}
  		TExpr(const Expr &e) : m_e(e) {}
  		
  	public: // set & get
  		const Expr & getExpr() const { return m_e; }
  
  	public: // functions
  
  		EXPRFUNC_1ARG(operator(), return m_e(i1););
  		EXPRFUNC_2ARG(operator(), return m_e(i1, i2););
  		EXPRFUNC_3ARG(operator(), return m_e(i1, i2, i3););
  
  	public: // operators
  
  		/// Assignement to one expression
  		template < typename Expr2 >
  			//TExpr< TExprAssign< Expr, Expr2 > >
  			TExpr< TExprBinaryOperator_NoRes < Expr, Expr2, functor::Assign> > 
  			operator=(const TExpr<Expr2> &e2) const
  		{
  			typedef TExprBinaryOperator_NoRes < Expr, Expr2, functor::Assign> TExprRet;
  			return TExpr<TExprRet>(TExprRet(this->getExpr(), e2.getExpr(), functor::Assign()));
  		}
  
  		/// Assignement to a constant
  		template < typename T >
  		TExpr< TExprBinaryOperator_NoRes < Expr, TExprConstant<T>, functor::Assign > >
  		operator=(const T &value) const
  		{
  			typedef TExprBinaryOperator_NoRes < Expr, TExprConstant<T>, functor::Assign > TExprRet;
  			return TExpr<TExprRet>(TExprRet(this->getExpr(), TExprConstant<T>(value), functor::Assign()));
  		}
  
      template < typename Expr2 >
      TExpr < TExprBinaryOperator<Expr, Expr2, functor::Call> >
      operator()(const TExpr<Expr2> & e2) const
      {
        typedef TExprBinaryOperator<Expr, Expr2, functor::Call> TExprRet;
        return TExpr<TExprRet>(TExprRet(this->getExpr(), e2.getExpr(), functor::Call()));
      }
  		
  	public: // self-arithmetic operators
  		
  		TIL_DEFINE_TEXPR_ARITHMETIC_OPERATOR(operator+=, til::expr::functor::AddTo);
  		TIL_DEFINE_TEXPR_ARITHMETIC_OPERATOR(operator-=, til::expr::functor::SubTo);
  		TIL_DEFINE_TEXPR_ARITHMETIC_OPERATOR(operator*=, til::expr::functor::MulTo);
  		TIL_DEFINE_TEXPR_ARITHMETIC_OPERATOR(operator/=, til::expr::functor::DivTo);
  
  	private: // data
  		Expr m_e;
  	};

#undef TIL_DEFINE_TEXPR_ARITHMETIC_OPERATOR
#undef TIL_DEFINE_TEXPR_ARITHMETIC_OPERATOR_EXPR
#undef TIL_DEFINE_TEXPR_ARITHMETIC_OPERATOR_VALUE
  

  //---------------------------------------------------------------------------------------------------

      //--------------------------//
     //  Placeholders instances  //
    //--------------------------//
    
    // TODO: is it really bad to have const instead of singletons? I think
    // in case of a stateless operator, this does not really matter, so simpler is 
    // probably better.
    /// Placeholder for the first argument.
	  const TExpr<expr::FirstArgument>	_1 = TExpr<FirstArgument>();
    /// Placeholder for the second argument.
	  const TExpr<expr::SecondArgument> _2 = TExpr<SecondArgument>();
    /// Placeholder for the third argument.
	  const TExpr<expr::ThirdArgument>	_3 = TExpr<ThirdArgument>();

  
  //---------------------------------------------------------------------------------------------------

      //--------------//
     //  TExprValue  //
    //--------------//
  
  	/// Represent a constant in template expression.
  	/// Unfortunately, as for now, only integer can be templated, so this
  	/// constant has to be an integer
  	// depraciated
  	template < int N, typename T >
  	class TExprValue
  	{
  	public: // typdefs
  		//typedef const T Type;
  
  		template < typename Iterator1, typename Iterator2 = Iterator1, typename Iterator3 = Iterator1 >
  		struct TypeStruct
  		{
  			typedef const T Type;
  		};
  		
  
  	public:
  	 
      /*
  		EXPRFUNC_1ARG(operator(), return N; );
  		EXPRFUNC_2ARG(operator(), return N; );
  		EXPRFUNC_3ARG(operator(), return N; );
      */
      EXPRFUNC_1ARG_ARG(operator(), return N;, );
      EXPRFUNC_2ARG_ARG(operator(), return N;, , );
      EXPRFUNC_3ARG_ARG(operator(), return N;, , , );

  	};
  	
  //---------------------------------------------------------------------------------------------------

    //---------------//
   //  TExprLValue  //
  //---------------//

    /// A template expression class for left-values.
  	template < typename T >
  	class TExprLValue
  	{
  	public: // typedefs
  		EXPR_RESULT_TYPE( T & );
  
  	public: // constructors & destructor
  		TExprLValue(T & lvalue) : m_pvalue(lvalue) {}
  
  	public: // functions
  
      /*
  		EXPRFUNC_1ARG( operator(), return m_pvalue; );
  		EXPRFUNC_2ARG( operator(), return m_pvalue; );
  		EXPRFUNC_3ARG( operator(), return m_pvalue; );
      */
      
      EXPRFUNC_1ARG_ARG( operator(), return m_pvalue;, );
      EXPRFUNC_2ARG_ARG( operator(), return m_pvalue;, , );
      EXPRFUNC_3ARG_ARG( operator(), return m_pvalue;, , , );
  
  	private: // data
  
  		T & m_pvalue;
  	};
  
  	/*
  	template < typename T >
  	class TExprLValue
  	{
  	public: // typedefs
  		EXPR_RESULT_TYPE( T );
  
  	public: // constructors & destructor
  		TExprLValue(T & lvalue) : m_lvalue(lvalue) {}
  
  	public: // functions
  
  		EXPRFUNC_1ARG( operator(), return m_lvalue; );
  		EXPRFUNC_2ARG( operator(), return m_lvalue; );
  		EXPRFUNC_3ARG( operator(), return m_lvalue; );
  
  		EXPRFUNC_1ARG_RET( getLValue, return m_lvalue; , T & );
  		EXPRFUNC_2ARG_RET( getLValue, return m_lvalue; , T & );
  		EXPRFUNC_3ARG_RET( getLValue, return m_lvalue; , T & );
  
  
  	private: // data
  
  		T &m_lvalue;
  	};
  */


  //---------------------------------------------------------------------------------------------------

    //----------------------//
   //  TExprFunctorHelper  //
  //----------------------//

    // NB: this is not technically a TExpr (it has no TypeStruct and its operator() are non standard).
    // It is just an overwrap over functors to redefine their operator(). Some kind of automatic
    // binding.
    template < typename TFunctor >
    class TExprFunctorHelper
    {
    public:
      TExprFunctorHelper(TFunctor f) : m_functor(f) {}

    public:
      template < typename Expr1 >
      TExpr<TExprUnaryOperator<Expr1, TFunctor> >
      operator()(TExpr<Expr1> e)
      {
        typedef TExprUnaryOperator<Expr1, TFunctor> TExprRet;
        return TExpr<TExprRet>(TExprRet(e.getExpr(), m_functor));
      }

      template < typename Expr1, typename Expr2 >
      TExpr<TExprBinaryOperator<Expr1, Expr2, TFunctor> >
      operator()(TExpr<Expr1> e1, TExpr<Expr2> e2)
      {
        typedef TExprBinaryOperator<Expr1, Expr2, TFunctor> TExprRet;
        return TExpr<TExprRet>(TExprRet(e1.getExpr(), e2.getExpr(), m_functor));
      }
      
    private:
      TFunctor m_functor;
    };
  
	} // namespace expr


  //---------------------------------------------------------------------------------------------------

    //------------------//
   //  loop functions  //
  //------------------//

  template < typename Functor, typename TIterator >
  void floop(Functor & functor, TIterator iIm)
  { for (; !iIm.isAtEnd(); ++iIm) functor(iIm); }

  template < typename Functor, typename TIterator >
  void floop(Functor functor, TIterator iIm)
  { for (; !iIm.isAtEnd(); ++iIm) functor(iIm); }

  template < typename Functor, typename TIterator >
  void floop2(TIterator iIm)
  { for (; !iIm.isAtEnd(); ++iIm) Functor::compute(*iIm); }

  template < typename Expr, typename TIterator1 >
  //typename boost::enable_if<is_ImageIterator<TIterator1> >::type
  void
  loop(expr::TExpr<Expr> & expr, TIterator1 iIm1)
	{
		for (; !iIm1.isAtEnd(); ++iIm1)
		{
			expr(iIm1);
		}
	}

  template < typename Expr, typename TIterator1 >
  //typename boost::enable_if<is_ImageIterator<TIterator1> >::type
  void
  loop(expr::TExpr<Expr> expr, TIterator1 iIm1)
  {
    for (; !iIm1.isAtEnd(); ++iIm1)
    {
      expr(iIm1);
    }
  }

	namespace detail
	{
		template < typename Expr, typename TIterator >
		void
		loop_c(expr::TExpr<Expr> expr, TIterator start, const TIterator & end)
		{
			for (; start != end; ++start)
			{
				expr(start);
			}
		}

		template < typename Expr, typename TContainer >
		void
		loop_x(expr::TExpr<Expr> expr, TContainer & c)
		{
			loop_c(expr, c.begin(), c.end());
		}

		
		template < typename Expr, typename TIterator1, typename TIterator2, typename TIterator3 >
		void
		loop_xxx(expr::TExpr<Expr> expr, TIterator1 start1, const TIterator1 end1, TIterator2 start2, TIterator3 start3)
		{
			for (; start1 != end1; ++start1, ++start2, ++start3)
			{
				expr(start1, start2, start3);
			}
		}

		template < typename Expr, typename TContainer1, typename TContainer2, typename TContainer3 >
		void
		loop_xxx(expr::TExpr<Expr> expr, TContainer1 & c1, TContainer2 & c2, TContainer3 & c3)
		{
			loop_xxx(expr, c1.begin(), c1.end(), c2.begin(), c3.begin());
		}

		template < typename Expr, typename TIterator1, typename TIterator2 >
		void
		loop_xx(expr::TExpr<Expr> expr, TIterator1 start1, const TIterator1 end1, TIterator2 start2)
		{
			for (; start1 != end1; ++start1, ++start2)
			{
				expr(start1, start2);
			}
		}

		template < typename Expr, typename TContainer1, typename TContainer2 >
		void
		loop_xx(expr::TExpr<Expr> expr, TContainer1 & c1, TContainer2 & c2)
		{
			loop_xx(expr, c1.begin(), c1.end(), c2.begin());
		}
	}

  template < int N, typename Expr, typename TIterator1, typename TIterator2 >
  void
  depth_loop(expr::TExpr<Expr> expr, TIterator1 begin1, TIterator1 end1, TIterator2 begin2)
  {
    if (N==1)
    {
      detail::loop_xx(expr, begin1, end1, begin2);
    }
    else
    {
      for (; begin1 != end1; ++begin1, ++begin2)
      {
        depth_loop<N-1>(begin1->begin(), begin1->end(), begin2->begin());
      }
    }
  }

  template < int N, typename Expr, typename TContainer1, typename TContainer2 >
  void
  depth_loop(expr::TExpr<Expr> expr, TContainer1 & c1, TContainer2 & c2)
  {
    depth_loop<N>(expr, c1.begin(), c1.end(), c2.begin());
  }

	template < typename Functor, typename TIterator >
		void
		loopboost(Functor &functor, TIterator iIm)
	{
		for (; !iIm.isAtEnd(); ++iIm) { functor(iIm); }
	}


  // TODO: Actually, all these test, that we might want to keep for example for performance comparison,
  // could be put in a special namespace, perhaps with special compilation flags.
  template < typename Expr, typename TIterator1 >
  void
  loop2(expr::TExpr<Expr> & expr, TIterator1 iIm1)
  {
    do
    {
      expr(iIm1);
    }
    while (iIm1.next());
  }

  template < typename Expr, typename TIterator1 >
  void
  loop2(expr::TExpr<Expr> expr, TIterator1 iIm1)
  {
    do
    {
      expr(iIm1);
    }
    while (iIm1.next());
  }

  /*
  template < typename Expr, typename TIterator1, typename TIterator2 >
  typename enable_if_c<
    is_ImageIterator<TIterator1>::value &&
    is_ImageIterator<TIterator2>::value
  >::type
  loop(expr::TExpr<Expr> & expr, TIterator1 iIm1, TIterator2 iIm2)
  { for (; !iIm1.isAtEnd(); ++iIm1, ++iIm2) expr(iIm1, iIm2); }
  */

  template < typename Expr, typename TIterator1, typename TIterator2 >
  typename enable_if_c<
    is_ImageIterator<TIterator1>::value &&
    is_ImageIterator<TIterator2>::value,
    expr::TExpr<Expr>
  >::type
  loop(expr::TExpr<Expr> expr, TIterator1 iIm1, TIterator2 iIm2)
  {
    for (; !iIm1.isAtEnd(); ++iIm1, ++iIm2) expr(iIm1, iIm2);
    return expr;
  }
  	
  template < typename TFunctor, typename TRange >
  inline void loop_r(TFunctor f, TRange r)
  {
    for (; r.ok(); ++r)
    {
      f(*r);
    }
  }

  /*
	template < typename Expr, typename TImage1, typename TImage2 >
	typename enable_if_c<is_Image<TImage1>::value && is_Image<TImage2>::value>::type
	loop(expr::TExpr<Expr> &expr, TImage1 &im1, TImage2 &im2)
	{
		loop(expr, itlin(im1), itlin(im2));
	}
  */


	template < typename Expr, typename TIterator1, typename TIterator2, typename TIterator3 >
	typename enable_if_c<
		is_ImageIterator<TIterator1>::value &&
		is_ImageIterator<TIterator2>::value &&
		is_ImageIterator<TIterator3>::value,
    expr::TExpr<Expr>
	>::type
	loop(expr::TExpr<Expr> expr, TIterator1 iIm1, TIterator2 iIm2, TIterator3 iIm3)
	{
		for (; !iIm1.isAtEnd(); ++iIm1, ++iIm2, ++iIm3) expr(iIm1, iIm2, iIm3);
    return expr;
	}

  /*
	template < typename Expr, typename TImage1, typename TImage2, typename TImage3 >
		typename enable_if_c<is_Image<TImage1>::value && is_Image<TImage2>::value && is_Image<TImage3>::value>::type
	loop_iii(expr::TExpr<Expr> &expr, TImage1 &im1, TImage2 &im2, TImage3 &im3)
	{
		loop(expr, itlin(im1), itlin(im2), itlin(im3));
	}
  */


  //---------------------------------------------------------------------------------------------------
  
    //---------//
   //  check  //  
  //---------//

  // use std::equal instead
  /*
  template < typename TPredicate, typename TIterator, typename TIterator2 >
  inline bool check(TIterator begin, TIterator end, TIterator2 begin2, TPredicate f)
  {
    for (; begin != end; ++begin, ++begin2)
    {
      if (!f(*begin, *begin2)) return false;
      return true;
    }
  }
  
  template < typename TPredicate, typename TRange1, typename TRange2 >
  inline bool check(TPredicate f, TRange1 r1, TRange2 r2)
  {
    for (; r1.ok(); ++r1, ++r2)
    {
      if (!f(*r1,*r2)) return false;
      return true;
    }
  }
  */

  /*
  template < typename TPredicate, typename TContainer1, typename TContainer2 >
  inline bool check(TPredicate f, const TContainer1 & x, const TContainer2 & y)
  {
    typedef typename oprange<TContainer1,TContainer2>::first_range Range1;
    typedef typename oprange<TContainer1,TContainer2>::second_range Range2;
    check(f, whole_range<Range1>(x), whole_range<Range2>(y));
  }
  */

	/*
	template < typename ImageFunctor, typename Expr>
	typename enable_if<is_ImageFunctor<ImageFunctor>,
	expr::TExpr<expr::TExprImageFunctor<ImageFunctor, Expr> >
	>::type
	bind(const ImageFunctor &functor, const expr::TExpr<Expr> &e)
	{
		typedef expr::TExprImageFunctor<ImageFunctor, Expr> TExprRet;
		return expr::TExpr<TExprRet>(TExprRet(functor, e.getExpr()));
	}
	*/


  /*
	/// Bind a detemplated functor.
	template < typename TFunctor, typename Expr>
	typename boost::enable_if<
		is_detemplated_functor<TFunctor>
		,
		expr::TExpr<expr::TExprUnaryOperator<Expr, TFunctor> >
	>::type
	//bind(TFunctor & functor, const expr::TExpr<Expr> & e)
  bind(TFunctor & functor, expr::TExpr<Expr> e)
	{
		typedef expr::TExprUnaryOperator<Expr, TFunctor> TExprRet;
		return expr::TExpr<TExprRet>(TExprRet(e.getExpr(), functor));
	}

	/// Bind a standard functor.
	// NB: e is const, (1) because it can: we make copy and return a different functor
	// anyway, (2) because it must: bind should accept const texpr such as placeholders.
	template < typename TFunctor, typename Expr>
	typename boost::disable_if<
		is_detemplated_functor<TFunctor>
		,
		expr::TExpr<expr::TExprUnaryOperator<Expr, expr::functor::Wrap<TFunctor> > >
	>::type
	bind(TFunctor & functor, expr::TExpr<Expr> e)
  //bind(TFunctor & functor, const expr::TExpr<Expr> & e)
	{
		typedef expr::TExprUnaryOperator<Expr, expr::functor::Wrap<TFunctor> > TExprRet;
		return expr::TExpr<TExprRet>(TExprRet(e.getExpr(), expr::functor::Wrap<TFunctor>(functor)));
	}
  */

} // namespace til


// Package includes
#include "til/TExprOperators.h"
#include "til/TExprFunctions.h"

#endif
