#ifndef TIL_POINT_OPERATORS_H_
#define TIL_POINT_OPERATORS_H_

namespace til
{
  namespace functor
  {
    /// By default, bounce any cast to the underlying container.
    template < typename X, typename T, std::size_t D, typename TStorage >
    struct CastTo<Point<T,D,TStorage>,X> : std::binary_function<Point<T,D,TStorage> &, const X &, void>
    {
      typedef std::binary_function<Point<T,D,TStorage> &, const X &, void> Base;
      typedef typename Base::first_argument_type    first_argument_type;
      typedef typename Base::second_argument_type   second_argument_type;
      void operator()(first_argument_type p, second_argument_type x)
      {
        CastTo<TStorage,X>()(p.data(),x);
      }
    };
  }


#define TIL_BOUNCE_BOOL_OP_TO_DATA(op)                                                                \
  template < typename T1, typename T2, std::size_t D, typename TStorage1, typename TStorage2 >        \
  inline bool operator op (const Point<T1,D,TStorage1> & x, const Point<T2,D,TStorage2> & y)          \
  { return x.data() op y.data(); }                                                                    \

  TIL_BOUNCE_BOOL_OP_TO_DATA(==)
  TIL_BOUNCE_BOOL_OP_TO_DATA(!=)

#undef TIL_BOUNCE_BOOL_OP_TO_DATA

#define TIL_BOUNCE_VEC_OP_TO_DATA(op, arg1, arg2, res, functor)                                             \
  template < typename T, std::size_t D, typename TStorage >                                                 \
  inline res <T,D,TStorage> operator op (const arg1 <T,D,TStorage> & x, const  arg2 <T,D,TStorage> & y)    \
  { return functor<TStorage,TStorage,res <T,D,TStorage> >()(x.data(), y.data()); }

  TIL_BOUNCE_VEC_OP_TO_DATA(+, Point, Vector, Point, functor::Add);
  TIL_BOUNCE_VEC_OP_TO_DATA(+, Vector, Point, Point, functor::Add);
  TIL_BOUNCE_VEC_OP_TO_DATA(-, Point, Vector, Point, functor::Sub);
  TIL_BOUNCE_VEC_OP_TO_DATA(-, Point, Point, Vector, functor::Sub);

#undef TIL_BOUNCE_VEC_OP_TO_DATA


#define TIL_BOUNCE_SCALAR_OP_TO_DATA(op, functor)                                            \
  template < typename T, std::size_t D, typename TStorage >                         \
  inline Point<T,D,TStorage> operator op (const Point<T,D,TStorage> & x, T y)       \
  { return functor <TStorage,T,Point<T,D,TStorage> >()(x.data(), y); }              \
  template < typename T, std::size_t D, typename TStorage >                         \
  inline Point<T,D,TStorage> operator op (T y, const Point<T,D,TStorage> & x)       \
  { return functor <T,TStorage,Point<T,D,TStorage> >()(y, x.data()); }             \


  TIL_BOUNCE_SCALAR_OP_TO_DATA(+, functor::Add)
  TIL_BOUNCE_SCALAR_OP_TO_DATA(-, functor::Sub)

#undef TIL_BOUNCE_SCALAR_OP_TO_DATA


  template < typename T, std::size_t D, typename TStorage >
  std::ostream & operator<<(std::ostream & os, const Point<T,D,TStorage> & p)
  {
    return os << p.data();
  }

}

#endif /*POINT_OPERATORS_H_*/
