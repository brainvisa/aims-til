#ifndef TIL_NUMERIC_ARRAY_OPERATORS_H
#define TIL_NUMERIC_ARRAY_OPERATORS_H

/// \file Belongs to numeric array package.
/// Do not include directly, include til/numeric_array instead

// includes from TIL
#include "til/functors.h"
// #include "til/sandbox.h"

namespace til
{

  //-----------------------------------------------------------------------------------------------
  
    //-----------------------//
   //  Arithmetic functors  //
  //-----------------------//

  namespace functor
  {
    // NB: I don't know if there is a case when we might have a better choice of iterator
    // for res than the random access through operator[] ??
#define TIL_DEFINE_TEMPLATE_OP(name, op)                                              \
    template < typename T1, typename T2, typename TRes, std::size_t D >               \
    struct name <numeric_array<T1,D>, numeric_array<T2,D>, TRes >                     \
      : public std::binary_function<const numeric_array<T1,D> &, const numeric_array<T2,D> &, TRes >      \
    {                                                                                 \
      typedef std::binary_function<const numeric_array<T1,D> &, const numeric_array<T2,D> &, TRes > Base; \
      typedef typename Base::first_argument_type    first_argument_type;              \
      typedef typename Base::second_argument_type   second_argument_type;             \
      typedef typename Base::result_type            result_type;                      \
      TRes operator()(first_argument_type x, second_argument_type y)                  \
      {                                                                               \
        typedef typename value_type_of<TRes>::type value_type;                        \
        TRes res;                                                                     \
        for (std::size_t i = 0; i < D; ++i)                                           \
          res[i] = static_cast<value_type>(x[i]) op static_cast<value_type>(y[i]);    \
        return res;                                                                   \
      }                                                                               \
    };                                                                                \
    template < typename T, typename TRes, std::size_t D >                             \
    struct name <T, numeric_array<T,D>, TRes >                                        \
      : public std::binary_function<T, const numeric_array<T,D> &, TRes >             \
    {                                                                                 \
      typedef std::binary_function<T, const numeric_array<T,D> &, TRes > Base;        \
      typedef typename Base::first_argument_type    first_argument_type;              \
      typedef typename Base::second_argument_type   second_argument_type;             \
      typedef typename Base::result_type            result_type;                      \
      TRes operator()(first_argument_type x, second_argument_type y)                  \
      {                                                                               \
        TRes res;                                                                     \
        for (std::size_t i = 0; i < D; ++i) res[i] = x op y[i];                       \
        return res;                                                                   \
      }                                                                               \
    };                                                                                \
    template < typename T, typename TRes, std::size_t D >                             \
    struct name <numeric_array<T,D>, T, TRes >                                        \
      : public std::binary_function<const numeric_array<T,D> &, T, TRes >             \
    {                                                                                 \
      typedef std::binary_function<const numeric_array<T,D> &, T, TRes > Base;        \
      typedef typename Base::first_argument_type    first_argument_type;              \
      typedef typename Base::second_argument_type   second_argument_type;             \
      typedef typename Base::result_type            result_type;                      \
      TRes operator()(first_argument_type x, second_argument_type y)                  \
      {                                                                               \
        TRes res;                                                                     \
        for (std::size_t i = 0; i < D; ++i) res[i] = x[i] op y;                       \
        return res;                                                                   \
      }                                                                               \
    };                                                                                \

    TIL_DEFINE_TEMPLATE_OP(Add, +)
    TIL_DEFINE_TEMPLATE_OP(Sub, -)
    TIL_DEFINE_TEMPLATE_OP(Mul, *)
    // TODO: don't we want some safeguard here for division? Say using til::fraction?
    TIL_DEFINE_TEMPLATE_OP(Div, /)

#undef TIL_DEFINE_TEMPLATE_OP


  //-----------------------------------------------------------------------------------------------

    //----------//
   //  CastTo  //
  //----------//

    template < typename T1, typename T2, std::size_t D >
    struct CastTo<numeric_array<T1,D>, numeric_array<T2,D> >
      : public std::binary_function<numeric_array<T1,D> &, const numeric_array<T2,D> &, void>
    {
      typedef std::binary_function<numeric_array<T1,D> &, const numeric_array<T2,D> &, void> Base;
      typedef typename Base::first_argument_type    first_argument_type;
      typedef typename Base::second_argument_type   second_argument_type;
      void operator()(first_argument_type x, second_argument_type y)
      {
        for (std::size_t i = 0; i < D; ++i) x[i] = y[i];
      }
    };

  } // namespace functor


  //-----------------------------------------------------------------------------------------------

    //-------------------------//
   //  artithmetic operators  //
  //-------------------------//


  // Definition of array-on-array operations
  // NB: I am not sure whether we should allow different value types here, because
  // then the problem occurs of the return type. I think a functor should probably
  // take care of these cases.

/*
#define TIL_OPERATOR(op)
  template < typename T1, typename T2, std::size_t D >                        \
  numeric_array<T,D>                                                          \
  operator op (const numeric_array<T1,D> & x, const numeric_array<T2,D> & y)  \
  {                                                                           \
    numeric_array<T,D> res;                                                   \
    for (std::size_t i = 0; i < D; ++i) res[i] = x[i] op y[i];                \
  }                                                                           \
*/
#define TIL_OPERATOR(op, func)                                                        \
  template < typename T, std::size_t D >                                              \
  inline                                                                              \
  numeric_array<T,D>                                                                  \
  operator op (const numeric_array<T,D> & x, const numeric_array<T,D> & y)            \
  {                                                                                   \
    return func <numeric_array<T,D>,numeric_array<T,D>,numeric_array<T,D> >() (x,y);  \
  }                                                                                   \

  TIL_OPERATOR(+, functor::Add);
  TIL_OPERATOR(-, functor::Sub);
  TIL_OPERATOR(*, functor::Mul);
  TIL_OPERATOR(/, functor::Div);

#undef TIL_OPERATOR


  // Operations with constants.
  // NB: I did not used boost::call_traits::param_type here. This is not because it
  // is always better to pass by values: the numerical type might actually not be numerical
  // at all, say std::complex, or even bigger. But I think the compiler can do more
  // intelligent stuff when a value is passed and this value is known at compile time.
  // So there's definitely a trade off here.
  // TODO: check how call_traits::param_type really works, and use it at least if it passes by
  // value all numerical types.

  // NB: for the second operator, we dont call x op y, because the operator (esp.
  // multiplication) might not be commutative.
#define TIL_OPERATOR(op, func)                                                  \
  template < typename T , std::size_t D >                                       \
  inline                                                                        \
  numeric_array<T,D>                                                            \
  operator op (const numeric_array<T,D> & x, T y)                               \
  {                                                                             \
    return func <numeric_array<T,D>,T,numeric_array<T,D> >()(x,y);              \
  }                                                                             \
  template < typename T, std::size_t D >                                        \
  inline                                                                        \
  numeric_array<T,D>                                                            \
  operator op (T y, const numeric_array<T,D> & x)                               \
  {                                                                             \
    return func <T,numeric_array<T,D>,numeric_array<T,D> >()(y,x);              \
  }                                                                             \
  
  TIL_OPERATOR(+, functor::Add);
  TIL_OPERATOR(-, functor::Sub);
  TIL_OPERATOR(*, functor::Mul);
  // Purposedly not implemented
  TIL_OPERATOR(/, functor::Div)
#undef TIL_OPERATOR



  //-----------------------------------------------------------------------------------------------


    //------------------------//
   //  Comparison operators  //
  //------------------------//


#define TIL_BOOL_VALUE_OPERATOR(opname, op)                                 \
  template < typename T1, typename T2, std::size_t D >                      \
  inline bool opname (const numeric_array<T1,D> & x,                        \
                      const numeric_array<T2,D> & y)                        \
  {                                                                         \
    for (std::size_t i = 0; i < D; ++i)                                     \
    {                                                                       \
      if (!(x[i] op y[i])) return false;                                    \
    }                                                                       \
    return true;                                                            \
  }                                                                         \
  template < typename T, std::size_t D >                                    \
  inline bool opname (const numeric_array<T,D> & x,                         \
                      typename boost::call_traits<T>::param_type y)         \
  {                                                                         \
    for (std::size_t i = 0; i < D; ++i)                                     \
    {                                                                       \
      if (!(x[i] op y)) return false;                                       \
    }                                                                       \
    return true;                                                            \
  }                                                                         \

  TIL_BOOL_VALUE_OPERATOR(operator==, ==);
  TIL_BOOL_VALUE_OPERATOR(operator!=, !=);  

  // The problem redefining operator< to do a all< is that there is no more
  // ordering on our objects. The problem is that these operators are taken
  // by default for maps and so on. So it's better to point out the potential
  // problem and give an explicit name to these.

  TIL_BOOL_VALUE_OPERATOR(all_less, <);
  TIL_BOOL_VALUE_OPERATOR(all_less_equal, <=);
  TIL_BOOL_VALUE_OPERATOR(all_greater, >);
  TIL_BOOL_VALUE_OPERATOR(all_greater_equal, >=);

#undef TIL_BOOL_VALUE_OPERATOR


  //-----------------------------------------------------------------------------------------------

    //----------------------------//
   //  Bad comparison operators  //
  //----------------------------//


#define BAD(op)                                                                       \
  template < typename T, std::size_t D >                                              \
  bool operator op (const numeric_array<T,D> & x, const numeric_array<T,D> & y)       \
  {                                                                                   \
    return x.do_not_use_these_operators();                                            \
  }                                                                                   \

  /// Do not use this, because it is ambiguous -- use all_XXX or any_XXX instead.
  BAD(<)
  /// Do not use this, because it is ambiguous -- use all_XXX or any_XXX instead.
  BAD(<=)
  /// Do not use this, because it is ambiguous -- use all_XXX or any_XXX instead.
  BAD(>)
  /// Do not use this, because it is ambiguous -- use all_XXX or any_XXX instead.
  BAD(>=)

#undef BAD

} // namespace til

#endif


