#ifndef TIL_TEXPR_OPERATORS_H
#define TIL_TEXPR_OPERATORS_H

/// \file Part of the TExpr package -- do not include this file
// TODO: put package files (i.e. files that should not be included directly by library
// users, but that are only included inside the library by the package it belongs to) in
// a separate physical folder, so that the user is not overwhelmed with include files
// when looking at the directory...

namespace til { namespace expr {

    //-------------------//
   //  Unary operators  //
  //-------------------//

#define TIL_DEFINITION_UNARY_OPERATOR_EXPRS(opName, functor)    \
  template < typename Expr >                                    \
  inline                                                        \
  TExpr< TExprUnaryOperator< Expr, functor > >                  \
  opName ( const TExpr<Expr> &e)                                \
  {                                                             \
    typedef TExprUnaryOperator< Expr, functor > TExprRet;       \
    return TExpr<TExprRet>(TExprRet(e.getExpr(), functor ()));  \
  }                                                             \

	TIL_DEFINITION_UNARY_OPERATOR_EXPRS(operator-, functor::Negate);
	TIL_DEFINITION_UNARY_OPERATOR_EXPRS(abs, functor::Abs);
	TIL_DEFINITION_UNARY_OPERATOR_EXPRS(operator*, functor::Deref);

#undef TIL_DEFINITION_UNARY_OPERATOR_EXPRS


    //-------------------------------//
   //  Binary operators with values //
  //-------------------------------//

	
#define TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE_LEFT(opName, functor)                 \
template < typename Expr, typename T >                                                  \
inline                                                                                  \
TExpr< TExprBinaryOperator< Expr, TExprConstant<T>, functor > >                         \
opName ( const TExpr<Expr> & e, const T & value )                                       \
{                                                                                       \
typedef TExprBinaryOperator<Expr, TExprConstant<T>, functor > TExprRet;	                \
return TExpr<TExprRet>(TExprRet(e.getExpr(), TExprConstant<T>(value), functor()));			\
}

#define TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE_RIGHT(opName, functor)                \
template < typename Expr, typename T >                                                  \
inline                                                                                  \
TExpr< TExprBinaryOperator< TExprConstant<T>, Expr, functor > >                         \
opName ( const T & value, const TExpr<Expr> & e )                                       \
{                                                                                       \
typedef TExprBinaryOperator<TExprConstant<T>, Expr, functor > TExprRet;                 \
return TExpr<TExprRet>(TExprRet(TExprConstant<T>(value), e.getExpr(), functor()));      \
}

#define TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE(opName, functor)                \
TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE_LEFT(opName, functor)                   \
TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE_RIGHT(opName, functor)

  TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE(operator+, functor::Plus);
  TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE(operator-, functor::Minus);
  TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE(operator*, functor::Multiplies);
  TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE(operator/, functor::Divides);
  TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE(operator<, functor::Less);
  TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE(operator<=, functor::Less_Equal);
  TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE(operator==, functor::Equal_To);
  TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE(operator!=, functor::Not_Equal_To);

/*
	TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE_LEFT(operator+, functor::Plus);
	TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE_LEFT(operator-, functor::Minus);
	TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE_LEFT(operator*, functor::Multiplies);
	TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE_LEFT(operator/, functor::Divides);
	TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE_LEFT(operator<, functor::Less);
	TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE_LEFT(operator<=, functor::Less_Equal);
	TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE_LEFT(operator==, functor::Equal_To);
	TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE_LEFT(operator!=, functor::Not_Equal_To);
*/

#undef TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE_LEFT
#undef TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE_RIGHT
#undef TIL_DEFINITION_BINARY_OPERATOR_EXPR_VALUE

  //----------------------------------//
 //  Binary operators between TExpr  //
//----------------------------------//

#define TIL_DEFINITION_BINARY_OPERATOR_EXPRS(opName, functor)             \
template < typename Expr1, typename Expr2 >                               \
inline                                                                    \
TExpr< TExprBinaryOperator< Expr1, Expr2, functor > >                     \
opName ( const TExpr<Expr1> &e1, const TExpr<Expr2> &e2)                  \
{                                                                         \
typedef TExprBinaryOperator<Expr1, Expr2, functor > TExprRet;             \
return TExpr<TExprRet>(TExprRet(e1.getExpr(), e2.getExpr(), functor()));  \
}

	TIL_DEFINITION_BINARY_OPERATOR_EXPRS(operator+, functor::Plus);
	TIL_DEFINITION_BINARY_OPERATOR_EXPRS(operator-, functor::Minus);
	TIL_DEFINITION_BINARY_OPERATOR_EXPRS(operator*, functor::Multiplies);
	TIL_DEFINITION_BINARY_OPERATOR_EXPRS(operator/, functor::Divides);
	TIL_DEFINITION_BINARY_OPERATOR_EXPRS(operator<, functor::Less);
	TIL_DEFINITION_BINARY_OPERATOR_EXPRS(operator<=, functor::Less_Equal);
	TIL_DEFINITION_BINARY_OPERATOR_EXPRS(operator==, functor::Equal_To);
	TIL_DEFINITION_BINARY_OPERATOR_EXPRS(operator!=, functor::Not_Equal_To);

	// NB: ?= assign-operators are defined inside TExpr class
	//TIL_DEFINITION_BINARY_OPERATOR_EXPRS(operator+=, functor::AddTo);

#undef TIL_DEFINITION_BINARY_OPERATOR_EXPRS


	/// Sequential evaluation of two template expressions
	template < typename Expr1, typename Expr2 >
	TExpr<TExprCouple<Expr1, Expr2> >
	operator,(
	const TExpr<Expr1> &e1,	///< This expression is evaluated first
	const TExpr<Expr2> &e2	///< This expression is evaluated last
	)
	{
		typedef TExprCouple<Expr1, Expr2> TExprRet;
		return TExpr<TExprRet>(TExprRet(e1.getExpr(), e2.getExpr()));
	}

	
}} // namespace til::expr

#endif

