#ifndef TIL_CUBIC_SPLINE_INTERPOLATION_TPP_
#define TIL_CUBIC_SPLINE_INTERPOLATION_TPP_

namespace til
{
  
  // small helper functions
  namespace
  {
    template < typename T >
    inline T cubic_1(T x)
    { return -0.25*cube(x-1); }
    
    template < typename T >
    inline T cubic_2(T x)
    { return 0.75*(x-3)*square(x-1)+1; }

    template < typename T >
    inline T cubic_3(T x)
    { return 0.75*(x-2)*square(x)+1; }

    template < typename T >
    inline T cubic_4(T x)
    { return 0.25*cube(x); }
  }
  
  template < typename T >
  inline T CubicSplineInterpolation<T>::compute(T f1, T f2, T f3, T f4, T x)
  { return f1*cubic_1(x) + f2*cubic_2(x) + f3*cubic_3(x) + f4*cubic_4(x); }
  
} // namespace til

#endif /*CUBIC_SPLINE_INTERPOLATION_TPP_*/
