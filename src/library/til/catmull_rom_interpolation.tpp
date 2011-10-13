#ifndef TIL_CATMULL_ROM_INTERPOLATION_TPP_
#define TIL_CATMULL_ROM_INTERPOLATION_TPP_

namespace til
{
  namespace
  {
    template < typename T >
    static inline T CRPo(T x)
    { return x*x*(x-1)/T(2.0); }
  
    template < typename T >
    static inline T lapl_(T x1, T x2, T x3)
    { return x1 + x3 - 2*x2; }
  }
  
  
  template < typename T > 
  T CatmullRomInterpolation<T>::compute(T f1, T f2, T f3, T f4, T x)
  {
    return 
      CRPo(x)   * lapl_(f2, f3, f4) + x     * f3 +
      CRPo(1-x) * lapl_(f1, f2, f3) + (1-x) * f2;
  }
  
} // namespace til

#endif /*CATMULL_ROM_INTERPOLATION_TPP_*/
