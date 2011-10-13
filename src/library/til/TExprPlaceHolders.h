#ifndef TIL_TEXPR_PLACE_HOLDERS_H
#define TIL_TEXPR_PLACE_HOLDERS_H

/// \file Belongs to TExpr package
/// Do not include directly, include til/TExpr.h instead

// includes from TIL library
#include "til/TExprMacros.h"

namespace til { namespace expr {

  /// Placeholder for the first argument.
	class FirstArgument
	{
	public:

		EXPR_RESULT_TYPE( Iterator1 );
		EXPR_ITERATORTYPE( Iterator1 );

    template < class Iterator1 >
    typename TypeStruct<Iterator1>::Type
    operator()(Iterator1 &i1)
    { return i1; }
    
    template < class Iterator1, class Iterator2>
    typename TypeStruct<Iterator1, Iterator2>::Type
    operator()(Iterator1 & i1, Iterator2 &)
    { return i1; }
    
    template < class Iterator1, class Iterator2, class Iterator3 >
    typename TypeStruct<Iterator1, Iterator2, Iterator3>::Type
    operator()(Iterator1 & i1, Iterator2 &, Iterator3 &)
    { return i1; }

    /*
		EXPRFUNC_1ARG( operator(), return i1; );
		EXPRFUNC_2ARG( operator(), return i1; );
		EXPRFUNC_3ARG( operator(), return i1; );		
    */
	};


	/// Placeholder for the second argument.
	class SecondArgument
	{
	public:

		EXPR_RESULT_TYPE( Iterator2 );
		EXPR_ITERATORTYPE( Iterator2 );

    template < class Iterator1, class Iterator2>
    typename TypeStruct<Iterator1, Iterator2>::Type
    operator()(Iterator1 &, Iterator2 & i2)
    { return i2; }
    
    template < class Iterator1, class Iterator2, class Iterator3 >
    typename TypeStruct<Iterator1, Iterator2, Iterator3>::Type
    operator()(Iterator1 &, Iterator2 & i2, Iterator3 &)
    { return i2; }

    /*
		EXPRFUNC_2ARG( operator(), return i2; )
		EXPRFUNC_3ARG( operator(), return i2; )
    */
    
		EXPRFUNC_2ARG_RET( getLValue, return i2; , const Iterator2 & );
		EXPRFUNC_3ARG_RET( getLValue, return i2; , const Iterator2 & );
	};

	/// Placeholder for the third argument.
	class ThirdArgument
	{
	public:

		EXPR_RESULT_TYPE( Iterator3 );
		EXPR_ITERATORTYPE( Iterator3 );

    template < class Iterator1, class Iterator2, class Iterator3 >
    typename TypeStruct<Iterator1, Iterator2, Iterator3>::Type
    operator()(Iterator1 &, Iterator2 &, Iterator3 & i3)
    { return i3; }

		//EXPRFUNC_3ARG( operator(), return i3; )

		EXPRFUNC_3ARG_RET( getLValue, return i3; , const Iterator3 & );
	};

}} // namespace til::expr


#endif

