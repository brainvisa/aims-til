#ifndef TIL_TEXPR_CONCATENATION_H
#define TIL_TEXPR_CONCATENATION_H

// includes from STL
#include <cmath>		// sqrt

// includes from BOOST
#include "boost/type_traits.hpp"

// includes from TIL library
//#include "til/ImageFunctorTraits.h"
#include "til/TExprBasicFunctors.h"
#include "til/TExprMacros.h"
#include "til/traits.h"

namespace til
{
	namespace expr
	{
		/*
		template < typename Expr, typename TDest >
		class TExprCast
		{
		public: // typedefs
			EXPR_RESULT_TYPE( TDest );
		public: // constructors & destructor
			TExprCast(const Expr &e) : m_e(e) {}
		public: // functions
			EXPRFUNC_1ARG(operator(), return TDest(m_e(i1)); );
			EXPRFUNC_2ARG(operator(), return TDest(m_e(i1, i2)); );
			EXPRFUNC_3ARG(operator(), return TDest(m_e(i1, i2, i3)); );
		private: // data
			Expr m_e;
		};
		*/

		/*
		template < typename UnaryImageFunctor, typename Expr >
		class TExprImageFunctor
		{
		public: // typedefs
			EXPR_RESULT_TYPE(typename functor::ImageFunctorTraits<UnaryImageFunctor>::TypeStruct<typename Expr::IteratorStruct<Iterator1 TIL_COMMA Iterator2 TIL_COMMA Iterator3>::Type >::Type );

		public: // constructors & destructor
			TExprImageFunctor(const UnaryImageFunctor &functor, const Expr &e) : m_functor(functor), m_e(e) {}
		
		public: // functions

			EXPRFUNC_1ARG(operator(), return m_functor(m_e(i1)); );
			EXPRFUNC_2ARG(operator(), return m_functor(m_e(i1, i2)); );
			EXPRFUNC_3ARG(operator(), return m_functor(m_e(i1, i2, i3)); );

		private: // data
			UnaryImageFunctor m_functor;
			Expr m_e;
		};
		*/

		/// Apply a unary numerical functor to a template expression
		template < typename Expr, typename UnaryFunctor >
		class TExprUnaryOperator
		{
		public: // typedefs
			EXPR_RESULT_TYPE(
typename UnaryFunctor::template TypeStruct<
  typename naked_type<
//  typename boost::remove_const<
    typename Expr::template TypeStruct<Iterator1 TIL_COMMA Iterator2 TIL_COMMA Iterator3>::Type
  >::type
>::Type );

		public: // constructors & destructor
			TExprUnaryOperator(const Expr &e, const UnaryFunctor &functor) : m_e(e), m_functor(functor) {}

		public: // functions
			EXPRFUNC_1ARG(operator(), return m_functor(m_e(i1)); );
			EXPRFUNC_2ARG(operator(), return m_functor(m_e(i1, i2)); );
			EXPRFUNC_3ARG(operator(), return m_functor(m_e(i1, i2, i3)); );

		private: // data

			Expr m_e;
			UnaryFunctor m_functor;
		};


		/// Apply a binary numerical functor to two template expressions
		template < typename TExpr1, typename TExpr2, typename BinaryOperator >
		class TExprBinaryOperator
		{
		public: //typedefs

			EXPR_RESULT_TYPE(
typename BinaryOperator::template TypeStruct<
    typename naked_type<
//    typename boost::remove_const<
//	  typename boost::remove_reference<
		typename TExpr1::template TypeStruct<Iterator1 TIL_COMMA Iterator2 TIL_COMMA Iterator3>::Type
	>::type 
	TIL_COMMA 
    typename naked_type<
//    typename boost::remove_const<
//    typename boost::remove_reference<
		typename TExpr2::template TypeStruct<Iterator1 TIL_COMMA Iterator2 TIL_COMMA Iterator3>::Type
	>::type
>::Type );

		public: // constructors & destructor

			TExprBinaryOperator(const TExpr1 &e1, const TExpr2 &e2, const BinaryOperator &functor)
				: m_e1(e1), m_e2(e2), m_functor(functor) {}

		public: // functions
		
			//EXPRFUNC_1ARG(operator(), return m_functor.operator()<typename boost::call_traits<typename TExpr1::TypeStruct< Iterator1 >::Type >::value_type>(m_e1(i1),         m_e2(i1)););
			EXPRFUNC_1ARG(operator(), return m_functor(m_e1(i1),         m_e2(i1)););
			EXPRFUNC_2ARG(operator(), return m_functor(m_e1(i1, i2),     m_e2(i1, i2)););
			EXPRFUNC_3ARG(operator(), return m_functor(m_e1(i1, i2, i3), m_e2(i1, i2, i3)););

		private: // data

			TExpr1 m_e1;
			TExpr2 m_e2;
			BinaryOperator m_functor;
		};


    template < typename TExpr1, typename TExpr2, typename BinaryOperator >
    class TExprBinaryOperator_NoRes
    {
    public: //typedefs

      EXPR_RESULT_TYPE(
typename BinaryOperator::template TypeStruct<
    typename naked_type<
//    typename boost::remove_const<
//    typename boost::remove_reference<
    typename TExpr1::template TypeStruct<Iterator1 TIL_COMMA Iterator2 TIL_COMMA Iterator3>::Type
  >::type 
  TIL_COMMA 
    typename naked_type<
//    typename boost::remove_const<
//    typename boost::remove_reference<
    typename TExpr2::template TypeStruct<Iterator1 TIL_COMMA Iterator2 TIL_COMMA Iterator3>::Type
  >::type
>::Type );

    public: // constructors & destructor

      TExprBinaryOperator_NoRes(const TExpr1 &e1, const TExpr2 &e2, const BinaryOperator &functor)
        : m_e1(e1), m_e2(e2), m_functor(functor) {}

    public: // functions
    
      //EXPRFUNC_1ARG(operator(), return m_functor.operator()<typename boost::call_traits<typename TExpr1::TypeStruct< Iterator1 >::Type >::value_type>(m_e1(i1),         m_e2(i1)););
      EXPRFUNC_1ARG(operator(), m_functor(m_e1(i1),         m_e2(i1)););
      EXPRFUNC_2ARG(operator(), m_functor(m_e1(i1, i2),     m_e2(i1, i2)););
      EXPRFUNC_3ARG(operator(), m_functor(m_e1(i1, i2, i3), m_e2(i1, i2, i3)););

    private: // data

      TExpr1 m_e1;
      TExpr2 m_e2;
      BinaryOperator m_functor;
    };



		/*
		/// Apply a binary numerical functor to two template expressions
		template < typename TExpr1, typename TExpr2, typename BinaryToOperator >
		class TExprBinaryToOperator
		{
		public: //typedefs

			EXPR_RESULT_TYPE(typename functor::FunctorTraits<BinaryToOperator>::template TypeStruct<typename TExpr1::template TypeStruct<Iterator1 TIL_COMMA Iterator2 TIL_COMMA Iterator3>::Type TIL_COMMA typename TExpr2::template TypeStruct<Iterator1 TIL_COMMA Iterator2 TIL_COMMA Iterator3>::Type >::Type );

		public: // constructors & destructor

			TExprBinaryToOperator(const TExpr1 &e1, const TExpr2 &e2, const BinaryToOperator &functor)
				: m_e1(e1), m_e2(e2), m_functor(functor) {}

		public: // functions
		
			EXPRFUNC_1ARG(operator(), return m_functor(m_e1(i1),         m_e2(i1)););
			EXPRFUNC_2ARG(operator(), return m_functor(m_e1(i1, i2),     m_e2(i1, i2)););
			//EXPRFUNC_2ARG(operator(), if (m_e2(i1, i2) > 0 && *(m_e1.getLValue(i1, i2)) % 1000 == 0) std::cout << *(m_e1.getLValue(i1, i2)) << std::endl;; return m_functor(*(m_e1.getLValue(i1, i2)),     m_e2(i1, i2)););
			//EXPRFUNC_2ARG(operator(), if (*i1 != *i2) std::cout << (int)*i1 << " " << (int)*i2 << " " << m_e2(i1,i2) << std::endl; return m_functor(*(m_e1.getLValue(i1, i2)),     m_e2(i1, i2)););
			EXPRFUNC_3ARG(operator(), return m_functor(m_e1(i1, i2, i3), m_e2(i1, i2, i3)););

		private: // data

			TExpr1 m_e1;
			TExpr2 m_e2;
			BinaryToOperator m_functor;
		};
		*/

		
		/// If/then block using template expressions
		template < typename TExprIf, typename TExprThen >
		class TExprIfThen
		{
		public: // typedefs

			EXPR_RESULT_TYPE(void);

		public: // constructors & destructor

			TExprIfThen(const TExprIf &eIf, const TExprThen &eThen)
				: m_eIf(eIf), m_eThen(eThen) {};

		public: // functions

			EXPRFUNC_1ARG(operator(), if (m_eIf(i1))         m_eThen(i1););
			EXPRFUNC_2ARG(operator(), if (m_eIf(i1, i2))     m_eThen(i1, i2););
			EXPRFUNC_3ARG(operator(), if (m_eIf(i1, i2, i3)) m_eThen(i1, i2, i3););

		private: // data

			TExprIf		m_eIf;
			TExprThen	m_eThen;
		};
		


		/// If/then block using template expressions
		template < typename TExprIf, typename TExprThen, typename TExprElse >
		class TExprIfThenElse
		{
		public: // typedefs

			EXPR_RESULT_TYPE(typename TExprThen::template TypeStruct<Iterator1 TIL_COMMA Iterator2 TIL_COMMA Iterator3>::Type );

		public: // constructors & destructor

			TExprIfThenElse(const TExprIf &eIf, const TExprThen &eThen, const TExprElse &eElse)
				: m_eIf(eIf), m_eThen(eThen), m_eElse(eElse) {};

		public: // functions

			EXPRFUNC_1ARG(operator(), if (m_eIf(i1))         return m_eThen(i1);			return m_eElse(i1); );
			EXPRFUNC_2ARG(operator(), if (m_eIf(i1, i2))     return m_eThen(i1, i2);		return m_eElse(i1, i2); );
			EXPRFUNC_3ARG(operator(), if (m_eIf(i1, i2, i3)) return m_eThen(i1, i2, i3);	return m_eElse(i1, i2, i3); );

		private: // data

			TExprIf		m_eIf;
			TExprThen	m_eThen;
			TExprElse	m_eElse;
		};

		/*
		/// Assign one expression to the other
		template < typename TExpr1, typename TExpr2 >
		class TExprAssign
		{
		public: // typedefs

			EXPR_RESULT_TYPE(void);

		public: // constructors & destructor

			TExprAssign(const TExpr1 &e1, const TExpr2 &e2)
				: m_e1(e1), m_e2(e2) {}

		public: // functions

			//EXPRFUNC_1ARG(operator(), *(m_e1.getLValue(i1)) = m_e2(i1););
			//EXPRFUNC_2ARG(operator(), *(m_e1.getLValue(i1, i2)) = m_e2(i1, i2););
			//EXPRFUNC_3ARG(operator(), *(m_e1.getLValue(i1, i2, i3)) = m_e2(i1, i2, i3););
			EXPRFUNC_1ARG(operator(), m_e1(i1) = m_e2(i1););
			EXPRFUNC_2ARG(operator(), m_e1(i1, i2) = m_e2(i1, i2););
			EXPRFUNC_3ARG(operator(), m_e1(i1, i2, i3) = m_e2(i1, i2, i3););

		private: // data

			TExpr1 m_e1;
			TExpr2 m_e2;
		};
		*/

		/*
		template < typename Expr >
		class TExprSqrt
		{
		public:
			EXPR_RESULT_TYPE(double);
			TExprSqrt(const Expr &e) : m_e(e) {}
			EXPRFUNC_1ARG(operator(), return std::sqrt(m_e(i1)); );
			EXPRFUNC_2ARG(operator(), return std::sqrt(m_e(i1, i2)); );
			EXPRFUNC_3ARG(operator(), return std::sqrt(m_e(i1, i2, i3)); );

		private:
			Expr m_e;
		};
		*/

		
		/// Evaluate one expression and then the other
		template < typename Expr1, typename Expr2 >
		class TExprCouple
		{
		public:
			// TODO: okay, maybe you want to return the value of the second expression...
			EXPR_RESULT_TYPE( void );
			TExprCouple(const Expr1 &e1, const Expr2 &e2) : m_e1(e1), m_e2(e2) {}
			EXPRFUNC_1ARG(operator(), ( m_e1(i1),			m_e2(i1)); );
			EXPRFUNC_2ARG(operator(), ( m_e1(i1, i2),		m_e2(i1, i2)); );
			EXPRFUNC_3ARG(operator(), ( m_e1(i1, i2, i3),	m_e2(i1, i2, i3)); );
      
		private:
			Expr1 m_e1;
			Expr2 m_e2;
		};

	} // namespace expr
} // namespace til

#endif

