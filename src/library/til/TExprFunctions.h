#ifndef TIL_TEXPR_FUNCTIONS_H
#define TIL_TEXPR_FUNCTIONS_H


/// \file Belongs to the TExpr package
/// Do not include directly, include til/TExpr.h instead

namespace til { namespace expr {


	/// if/then block for expressions
	// TODO: put some checking in there. One could check that
	// the return type of the first expression is a bool. One could
	template < typename ExprIf, typename ExprThen >
// The following don't work, because void::value_type is not defined, but if we
// want TEXpr to work even for non-iteratos, we would just have to use
// traits for iterator types (probably endemic to TExpr only, so we can put
// that in traits), and never use Iterator::value_type but value_type<Iterator>
// that would be defined for numeric types as well (or rather pointers on
// numeric types). To do that we need to have an IsIteratorType I guess
//			enable_if<are_same_types<bool, typename ExprIf::TypeStruct<void,void,void>::Type>,
		TExpr<TExprIfThen<ExprIf, ExprThen> >
//			>
		if_then(
		const TExpr<ExprIf> &eIf,		///< This expression has to return a boolean
		const TExpr<ExprThen> &eThen	///< Executed if eIf returns true
		)
	{
		typedef TExprIfThen<ExprIf, ExprThen> TExprRet;
		return TExpr<TExprRet>(TExprRet(eIf.getExpr(), eThen.getExpr()));
	}

	/// If/then/else structure for template expressions
	template < typename ExprIf, typename ExprThen, typename ExprElse >
	TExpr<TExprIfThenElse<ExprIf, ExprThen, ExprElse> >
	if_then_else(
	const TExpr<ExprIf>	  &eIf,		///< This expression has to return a boolean
	const TExpr<ExprThen> &eThen,	///< Executed if eIf returns true
	const TExpr<ExprElse> &eElse	///< Executed if eIf returns false
	)
	{
		typedef TExprIfThenElse<ExprIf, ExprThen, ExprElse> TExprRet;
		return TExpr<TExprRet>(TExprRet(eIf.getExpr(), eThen.getExpr(), eElse.getExpr()));
	}


  /// Square root.
	template < typename Expr >
	//TExpr<TExprSqrt<Expr> >
	TExpr< TExprUnaryOperator< Expr, functor::Sqrt> >
	sqrt(const TExpr<Expr> &e)
	{
		typedef TExprUnaryOperator< Expr, functor::Sqrt> TExprRet;
		return TExpr<TExprRet>(TExprRet(e.getExpr(), functor::Sqrt()));
	}

  // dunno why, but this one cannot be in this namespace...
  /*
  /// static type casting.
  template < typename TTo, typename Expr >
	inline
  TExpr<TExprUnaryOperator< Expr, functor::Cast<TTo> > >
	cast (const TExpr<Expr> &e)
	{
		typedef TExprUnaryOperator<Expr, functor::Cast<TTo> > TExprRet;
		return TExpr<TExprRet>(TExprRet(e.getExpr(), functor::Cast<TTo>()));
	}
  */

  /// Insert a left value in a template expression.
  template < typename T >
	TExpr<TExprLValue<T> >
	var(T & lvalue)
	{
		return TExpr<TExprLValue<T> >(TExprLValue<T>(lvalue));
	}

  template < typename TFunctor >
  typename boost::enable_if<
    is_detemplated_functor<TFunctor>,
    expr::TExprFunctorHelper<TFunctor>
  >::type
  func(TFunctor & f)
  {
    return expr::TExprFunctorHelper<TFunctor>(f);
  }
  
  template < typename TFunctor >
  typename boost::disable_if<
    is_detemplated_functor<TFunctor>,
    expr::TExprFunctorHelper<expr::functor::Wrap<TFunctor> >
  >::type
  func(TFunctor & f)
  {
    return expr::TExprFunctorHelper<expr::functor::Wrap<TFunctor> >(f);
  }

  /*

	/// Bind a detemplated functor.
	template < typename TFunctor, typename Expr>
	typename boost::enable_if<
		is_detemplated_functor<TFunctor>
		,
		expr::TExpr<expr::TExprUnaryOperator<Expr, TFunctor> >
	>::type
	bind(TFunctor & functor, const expr::TExpr<Expr> & e)
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
	bind(TFunctor & functor, const expr::TExpr<Expr> & e)
	{
		typedef expr::TExprUnaryOperator<Expr, expr::functor::Wrap<TFunctor> > TExprRet;
		return expr::TExpr<TExprRet>(TExprRet(e.getExpr(), expr::functor::Wrap<TFunctor>(functor)));
	}
  */

  /// static type casting.
  template < typename TTo, typename Expr >
	inline
  TExpr<TExprUnaryOperator< Expr, functor::Cast<TTo> > >
	cast (const TExpr<Expr> &e)
	{
		typedef TExprUnaryOperator<Expr, functor::Cast<TTo> > TExprRet;
		return TExpr<TExprRet>(TExprRet(e.getExpr(), functor::Cast<TTo>()));
	}
  
  
  template < typename Expr1, typename Expr2 >
  inline
  TExpr<TExprBinaryOperator< Expr1, Expr2, functor::CastTo > >
  castTo(TExpr<Expr1> e1, TExpr<Expr2> e2)
  {
    typedef TExprBinaryOperator< Expr1, Expr2, functor::CastTo > TExprRet;
    return TExpr<TExprRet>(TExprRet(e1.getExpr(), e2.getExpr(), functor::CastTo()));
  }
  
  }// namespace expr

} // namespace til

#endif

