#ifndef TIL_KEYS_INTERPOLATION_TPP_
#define TIL_KEYS_INTERPOLATION_TPP_

namespace til
{
  namespace
  {
    template < typename T >
    static inline T lapl(T x1, T x2, T x3)
    { return x1 + x3 - 2*x2; }

    template < typename T >
    static inline T bilapl(T x1, T x2, T x3, T x4, T x5)
    { return x1 + x5 + 6*x3 - 4*(x2+x4); }

    template < typename T >
    static inline T keysPo1(T x)
    { return x*(x*x-1)/T(6.0); }

    template < typename T >
    static inline T keysPo2(T x)
    { return -x*x*(x-1)/T(12.0);}
  }
  
  
  template < typename T >
  T KeysInterpolation<T>::compute(T f1, T f2, T f3, T f4, T f5, T f6, T x) 
  {
    return
      x     * f4 + keysPo1(x)   * lapl(f3, f4, f5) + keysPo2(x)   * bilapl(f2, f3, f4, f5, f6) +
      (1-x) * f3 + keysPo1(1-x) * lapl(f2, f3, f4) + keysPo2(1-x) * bilapl(f1, f2, f3, f4, f5);
  }
   
} // namespace til


#endif /*KEYS_INTERPOLATION_TPP_*/
