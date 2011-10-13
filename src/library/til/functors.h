#ifndef TIL_FUNCTORS_H
#define TIL_FUNCTORS_H

/// \file This file contains small basic functors.
/// Some of them wrap standard functions (e.g. std::sqrt) into functors, which is usefull for template
/// expression.
/// Some others wrap arithmetic operators. The standard already define classes for them (e.g. std::plus), but
/// taking arguments all having the same type. Always having consistent types when doing numerical operations
/// is probably an excellent idea; However, it is necessary here to allow different types for these operators
/// to allow higher level operations in template expressions, e.g. a multiplication between a matrix and a vector.

// includes from STL
#include <cmath>		    // sqrt, exp
#include <functional>	  // unary_function, binary_function

// includes from BOOST
#include <boost/call_traits.hpp>

// includes from TIL
#include "til/traits.h"

namespace til { namespace functor
{


  /// Absolute value.
 	template < typename T >
	struct Abs : public std::unary_function<T,T>
	{
    typedef std::unary_function<T,T>        Base;
    typedef typename Base::argument_type    argument_type;
    typedef typename Base::result_type      result_type;
        
		inline result_type operator()(argument_type x) const
		{
			return std::abs(x);//( x>0? x : -x );
		}
	};


  /// Square
  /*
  template < typename T >
  struct Square : public std::unary_function<T,T>
  {
    typedef std::unary_function<T,T>        Base;
    typedef typename Base::argument_type    argument_type;
    typedef typename Base::result_type      result_type;

    inline result_type operator()(argument_type x) const
    {
      return x*x;
    }
  };
  */
  
  template < typename T, unsigned int D >
  struct Pow
    : public std::unary_function<T,T>
  {
    typedef std::unary_function<T,T>        Base;
    typedef typename Base::argument_type    argument_type;
    typedef typename Base::result_type      result_type;

    inline result_type operator()(argument_type x) const
    {
      if (D%2)
        return square(Pow<T,D/2>()(x)) * x;
      else
        return square(Pow<T,D/2>()(x));
    }
  };

  template < typename T >
  struct Pow<T,1>
    : public std::unary_function<T,T>
  {
    typedef std::unary_function<T,T>        Base;
    typedef typename Base::argument_type    argument_type;
    typedef typename Base::result_type      result_type;

    inline result_type operator()(argument_type x) const
    {
      return x;
    }
  };
  
  template < typename T >
  struct Pow<T,0>
    : public std::unary_function<T,T>
  {
    typedef std::unary_function<T,T>        Base;
    typedef typename Base::argument_type    argument_type;
    typedef typename Base::result_type      result_type;

    class InvalidArgument : public std::exception {};
    
    inline result_type operator()(argument_type x) const
    {
      if (x == 0) throw InvalidArgument();
      else return 1;
    }
  };

  /// Square root.
	template < typename T >
	struct Sqrt
    : public std::unary_function<T,T>
	{
    typedef std::unary_function<T,T>        Base;
    typedef typename Base::argument_type    argument_type;
    typedef typename Base::result_type      result_type;

		inline result_type operator()(argument_type x) const
		{
			return std::sqrt(x);
		}
	};

  /// Exponential.
	template < typename T >
	struct Exp
    : public std::unary_function<T,T>
	{
    typedef std::unary_function<T,T>        Base;
    typedef typename Base::argument_type    argument_type;
    typedef typename Base::result_type      result_type;
		result_type operator()(argument_type x) const
		{
			return std::exp(x);
		}
	};

  /// Static cast functor.
  template < typename TTo, typename TFrom  >
  struct CastTo
    : public std::binary_function<TTo &, typename boost::call_traits<TFrom>::param_type, void>
  {
    typedef std::binary_function<TTo &, typename boost::call_traits<TFrom>::param_type, void> Base;
    typedef typename Base::first_argument_type    first_argument_type;
    typedef typename Base::second_argument_type   second_argument_type;
    typedef typename Base::result_type            result_type;
    void operator()(first_argument_type x, second_argument_type y) const { x = y; }
  };
  
  /// Two-step static class functor.
  template < typename TFrom >
  struct Convertor
  {
    Convertor(typename boost::call_traits<TFrom>::param_type from) : m_from(from) {}
    
    template < typename TTo >
    void into(TTo & to)
    {
      CastTo<TTo, TFrom>()(to, m_from);
    }
    
    typename boost::call_traits<TFrom>::param_type m_from;
  };
  
  /// Static cast functor.
  template < typename TTo, typename TFrom >
  struct Cast
    : public std::unary_function<typename boost::call_traits<TFrom>::param_type, TTo>
	{
    typedef std::unary_function<typename boost::call_traits<TFrom>::param_type,TTo>        Base;
    typedef typename Base::argument_type    argument_type;
    typedef typename Base::result_type      result_type;
		result_type operator()(argument_type x) const
		{
			// Actually, what is the real difference between static casting and
			// constructing.
			// return static_cast<TTo>(x);
			// return TTo(x);
      // Do not call castTo because they should rely on different mechanisms:
      // Cast rely on constructors
      // castTo implements operator=
      TTo y;
      CastTo<TTo,TFrom>()(y,x);
      return y;
		}
	};

  /// Dereferenciation functor.
	template < typename T >
	struct Deref
    : public std::unary_function<T, typename deref<T>::type>
	{
    typedef std::unary_function<T, typename deref<T>::type> Base;
    typedef typename Base::argument_type    argument_type;
    typedef typename Base::result_type      result_type;
		result_type operator()(argument_type x) const
		{
			return *x;
		}
	};

  /// In-place addition functor.
	template < typename T1, typename T2 >
	struct AddTo
    : public std::binary_function<typename boost::add_reference<T1>::type,typename boost::call_traits<T2>::param_type,void>
	{
    typedef std::binary_function<typename boost::add_reference<T1>::type,typename boost::call_traits<T2>::param_type,void> Base;
    typedef typename Base::first_argument_type    first_argument_type;
    typedef typename Base::second_argument_type   second_argument_type;
    typedef typename Base::result_type            result_type;
		result_type operator()(first_argument_type x, second_argument_type y) const
		{
			x += y;
		}
	};

  /// In-place subtraction functor.
	template < typename T1, typename T2 >
	struct SubTo
    : public std::binary_function<typename boost::add_reference<T1>::type,typename boost::call_traits<T2>::param_type,void>
	{
    typedef std::binary_function<typename boost::add_reference<T1>::type,typename boost::call_traits<T2>::param_type,void> Base;
    typedef typename Base::first_argument_type    first_argument_type;
    typedef typename Base::second_argument_type   second_argument_type;
    typedef typename Base::result_type            result_type;
		result_type operator()(first_argument_type x, second_argument_type y) const
		{
			x -= y;
		}
	};

  /// In-place multiplication functor.
	template < typename T1, typename T2 >
	struct MulTo
    : public std::binary_function<typename boost::add_reference<T1>::type,typename boost::call_traits<T2>::param_type,void>
	{
    typedef std::binary_function<typename boost::add_reference<T1>::type,typename boost::call_traits<T2>::param_type,void> Base;
    typedef typename Base::first_argument_type    first_argument_type;
    typedef typename Base::second_argument_type   second_argument_type;
    typedef typename Base::result_type            result_type;
		result_type operator()(first_argument_type x, second_argument_type y) const
		{
			x *= y;
		}
	};

  /// In-place division functor.
	template < typename T1, typename T2 >
	struct DivTo
    : public std::binary_function<typename boost::add_reference<T1>::type,typename boost::call_traits<T2>::param_type,void>
	{
    typedef std::binary_function<typename boost::add_reference<T1>::type,typename boost::call_traits<T2>::param_type,void> Base;
    typedef typename Base::first_argument_type    first_argument_type;
    typedef typename Base::second_argument_type   second_argument_type;
    typedef typename Base::result_type            result_type;
		result_type operator()(first_argument_type x, second_argument_type y) const
		{
			x /= y;
		}
	};


  /// Addition functor.
  // TODO: this stupid combine doesn't mean anything. Or what I mean is that it is relevant to numerical types
  // only. Probably we should have a default_add_trait struct somewhere, that itself could default to combine.
	template < typename T1, typename T2, typename TRes = typename combine<T1,T2>::type >
	struct Add
    : public std::binary_function<typename boost::call_traits<T1>::param_type, typename boost::call_traits<T2>::param_type, TRes>
	{
    typedef std::binary_function<typename boost::call_traits<T1>::param_type, typename boost::call_traits<T2>::param_type, TRes> Base;
    typedef typename Base::first_argument_type    first_argument_type;
    typedef typename Base::second_argument_type   second_argument_type;
    typedef typename Base::result_type            result_type;
		result_type operator()(first_argument_type x, second_argument_type y) const
		{
			return x + y;
		}
	};

  /// Subtraction functor.
	template < typename T1, typename T2, typename TRes = typename combine<T1,T2>::type >
	struct Sub
    : public std::binary_function<typename boost::call_traits<T1>::param_type, typename boost::call_traits<T2>::param_type, TRes>
	{
    typedef std::binary_function<typename boost::call_traits<T1>::param_type, typename boost::call_traits<T2>::param_type, TRes> Base;
    typedef typename Base::first_argument_type    first_argument_type;
    typedef typename Base::second_argument_type   second_argument_type;
    typedef typename Base::result_type            result_type;
		result_type operator()(first_argument_type x, second_argument_type y) const
		{
			return x - y;
		}
	};

  /// Multiplication functor.
	template < typename T1, typename T2, typename TRes = typename combine<T1,T2>::type >
	struct Mul
    : public std::binary_function<typename boost::call_traits<T1>::param_type, typename boost::call_traits<T2>::param_type, TRes>
	{
    typedef std::binary_function<typename boost::call_traits<T1>::param_type, typename boost::call_traits<T2>::param_type, TRes> Base;
    typedef typename Base::first_argument_type    first_argument_type;
    typedef typename Base::second_argument_type   second_argument_type;
    typedef typename Base::result_type            result_type;
		result_type operator()(first_argument_type x, second_argument_type y) const
		{
			return x * y;
		}
	};

  /// Division functor.
	template < typename T1, typename T2, typename TRes = typename combine<T1,T2>::type >
	struct Div
    : public std::binary_function<typename boost::call_traits<T1>::param_type, typename boost::call_traits<T2>::param_type, TRes>
	{
    typedef std::binary_function<typename boost::call_traits<T1>::param_type, typename boost::call_traits<T2>::param_type, TRes> Base;
    typedef typename Base::first_argument_type    first_argument_type;
    typedef typename Base::second_argument_type   second_argument_type;
    typedef typename Base::result_type            result_type;
		result_type operator()(first_argument_type x, second_argument_type y) const
		{
      // TODO: should we use til::Fraction instead of raw '/' here?
			return x / y;
		}
	};

  /// Function call functor.
  template < typename TFunctor, typename T >
  struct Call
    : public std::binary_function<TFunctor, T, typename TFunctor::result_type>
  {
    typedef std::binary_function<TFunctor, T, typename TFunctor::result_type> Base;
    typedef typename Base::first_argument_type    first_argument_type;
    typedef typename Base::second_argument_type   second_argument_type;
    typedef typename Base::result_type            result_type;
    result_type operator()(first_argument_type x, second_argument_type y) const
    {
      return x(y);
    }
  };

  /// Assignment functor.
  template < typename T >
	struct Assign
    : public std::binary_function<typename boost::add_reference<T>::type,typename boost::call_traits<T>::param_type,void>
	{
    typedef std::binary_function<typename boost::add_reference<T>::type,typename boost::call_traits<T>::param_type,void> Base;
    typedef typename Base::first_argument_type    first_argument_type;
    typedef typename Base::second_argument_type   second_argument_type;
    typedef typename Base::result_type            result_type;    
		result_type operator()(first_argument_type x, second_argument_type y)
		{
			x = y;
		}
	};

  /*
	/// Accumulation functor.
	/// Very similar indeed to std::accumulate function.
	// NB: no default is given for TAccumulation, because this is an important choice
	// that depends on the context and should not be ignored.
	template < typename T, typename TAccumulation, typename TBinaryOperator = Add<TAccumulation,T> >
	class Accumulate : public std::unary_function<T,void>
	{
  public: // typedefs
    typedef std::unary_function<T,void>     Base;
    typedef typename Base::argument_type    argument_type;
    typedef typename Base::result_type      result_type;

	public: // constructors & destructor
		Accumulate() : Base(), m_accumulation() {}
    
	public: // functions
		TAccumulation get() const { return m_accumulation; }
    
	public: // operators
		result_type operator()(argument_type x)
		{
			m_accumulation = TBinaryOperator()(m_accumulation, x);
		}
    
	private: // data
		TAccumulation m_accumulation;
	};
  */

}} // namespace til::functor


  //-------------------------------------------------------------------------------------------------
  
    //--------------------//
   //  helper functions  //
  //--------------------//

namespace til
{
  template < typename TRes, typename T1, typename T2 >
  inline TRes add(T1 x, T2 y)
  { return functor::Add<T1,T2,TRes>()(x,y); }

  template < typename TRes, typename T1, typename T2 >
  inline TRes sub(T1 x, T2 y)
  { return functor::Sub<T1,T2,TRes>()(x,y); }

  template < typename TRes, typename T1, typename T2 >
  inline TRes mul(T1 x, T2 y)
  { return functor::Mul<T1,T2,TRes>()(x,y); }

  template < typename TRes, typename T1, typename T2 >
  inline TRes div(T1 x, T2 y)
  { return functor::Div<T1,T2,TRes>()(x,y); }

  template < typename TTo, typename TFrom >
  inline void convert(TTo & x, const TFrom & y)
  { functor::CastTo<TTo,TFrom>()(x,y); }
  
  template < typename TTo, typename TFrom >
  inline TTo convert(const TFrom & y)
  { return functor::Cast<TTo, TFrom>()(y); }
  
  template < typename T >
  inline functor::Convertor<T>
  convert2(const T & from)
  { return functor::Convertor<T>(from); }
}

#endif


