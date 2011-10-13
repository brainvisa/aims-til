#ifndef TIL_BASIC_FUNCTORS_H
#define TIL_BASIC_FUNCTORS_H

// includes from STL
#include <cassert>
#include <limits>

namespace til {
namespace functor {

	/// The following functions are necessary for the coming functors. They should probably
	/// be the only cases when those ternary functors (that are really doing two things:
	/// operate THEN assign) are simply a call to z = x OP y.

template < typename T1, typename T2, typename T3 >
inline
typename boost::enable_if_c<
  std::numeric_limits<T1>::is_specialized &&
  std::numeric_limits<T2>::is_specialized
>::type
add(T1 x, T2 y, T3& z)
{
  z = x + y;
}

template < typename T1, typename T2, typename T3 >
inline
typename boost::enable_if_c<
  std::numeric_limits<T1>::is_specialized &&
  std::numeric_limits<T2>::is_specialized
>::type
sub(T1 x, T2 y, T3& z)
{
  z = x - y;
}

template < typename T1, typename T2, typename T3 >
inline
typename boost::enable_if_c<
  std::numeric_limits<T1>::is_specialized &&
  std::numeric_limits<T2>::is_specialized
>::type
mul(T1 x, T2 y, T3& z)
{
  z = x * y;
}

template < typename T1, typename T2, typename T3 >
inline
typename boost::enable_if_c<
  std::numeric_limits<T1>::is_specialized &&
  std::numeric_limits<T2>::is_specialized
>::type
div(T1 x, T2 y, T3& z)
{
  assert(y != 0);
  z = x / y;
}
/*

  template < typename T1 = None, typename T2 = None, typename T3 = typename combine<T1, T2>::type >
  struct Add
  {
  private: // classes
	struct FlagNone {};
	struct FlagOther {};

	template < typename T > struct Flag { typedef FlagOther type; };
	template <> struct Flag<None> { typedef FlagNone type; };

  public: // operators

	template < typename X1, typename X2 >
		void operator() (X1 &x, const X2 & y)
	{
		dispatch(x, y, typename Flag<T1>::type(), typename Flag<T2>::type());
	}

  private: // functions
	template < typename X1, typename X2 >
	void dispatch(X1 x, X2 y, FlagNone, FlagNone)
	{
		x += y;
	}

	void dispatch(T1 x1, T2 x2, FlagOther, FlagOther)
	{
		x += y;
	}

    //t T1 T2 add(T1 T2)

    template <typename _T1, typename T2>
    inline void operator()(T1 &x, const T2 &y) { x += y; }

    template <typename T1, typename T2, typename T3>
    inline void operator()(const T1 &x, const T2 &y, T3 & z) { add(x,y,z); }
  };
*/

  /*
  struct Add
  {
    template <typename T1, typename T2>
    inline void operator()(T1 &x, const T2 &y) { x += y; }

    template <typename T1, typename T2, typename T3>
    inline void operator()(const T1 &x, const T2 &y, T3 & z) { add(x,y,z); }
  };

  struct Sub
  {
    template <typename T1, typename T2>
    inline void operator()(T1 &x, const T2 &y) { x -= y; }

    template <typename T1, typename T2, typename T3>
    inline void operator()(const T1 &x, const T2 &y, T3 & z) { sub(x,y,z); }
  };

  struct Div
  {
    template < typename T1, typename T2 >
    inline void operator()(T1 & x, const T2 & y) { x /= y); }
    
    template < typename T1, typename T2, typename T3>
    inline void operator()(const T1 &x, const T2 &y, T3 & z) { div(x,y,z); }
  };
  
  // Actually it might be a better idea to use *=/ * operators in there rather
  // than mul. Because these functors should be used at the numeric point,
  // meaning these operators should be defined on these types.
  // No. the problem stems from z = x*y which is a bad idea (constructor and copy
  // called). So mul has to be called. And thus, for consistency, also for the
  // *= operation.
  // TODO: It just might be that we never need a binary form. Since ternary forms
  // are inlined, it may see that two of the elements are the same and thus
  // simplify the whole thing.
  struct Mul
  {
    template < typename T1, typename T2 >
    inline void operator()(T1 & x, const T2 & y) { x *= y; }

    template < typename T1, typename T2, typename T3 >
    inline void operator()(const T1 & x, const T2 & y, T3 & z) { mul(x,y,z); }
  };
  
  / *
  template < typename T1 = void, typename T2 = T1 >
  struct Sqrt
  {

	//TODO: use a numeric_cast
    template < typename LT1, typename LT2 >
    inline void operator()(const LT1 & x, LT2 & y) { y = static_cast<LT2>(std::sqrt<double>(x)); }
    
    inline void operator()(const T1 & x, T2 & y) { y = static_cast<T2>(std::sqrt<double>(x)); }
  };
  */

} // namespace functor
} // namespace til

#endif

