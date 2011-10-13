#ifndef TIL_MISC_SCALAR_FUNCTIONS_H_
#define TIL_MISC_SCALAR_FUNCTIONS_H_

/// \file Simple, standard, small numerical functions

namespace til
{
  ///< Square of a number
  template < typename T > INLINE const T square(const T &a) { return a*a;}  
  ///< Cube of a number
  template < typename T > INLINE const T cube(const T &a) { return a*a*a;}
  ///< Max of two numbers
  template < typename T > INLINE const T max(const T &a, const T &b) { return ((a>b)?a:b);}
  ///< Min of two numbers
  template < typename T > INLINE const T min(const T &a, const T &b) { return ((a<b)?a:b);}
  ///< Max of three numbers
  template < typename T > INLINE const T max(const T &a, const T &b, const T &c) { return max(max(a,b),c);}
  ///< Min of three numbers
  template < typename T > INLINE const T min(const T &a, const T &b, const T &c) { return min(min(a,b),c);}
  ///< Min of four numbers
  template < typename T > INLINE const T min(const T &a, const T &b, const T &c, const T &d) { return min(min(min(a,b),c),d);}
  ///< Force x into [a, b]
  template < typename T > 
  inline const T crop(const T & x, const T & a, const T & b) { return min(max(x, a), b); }
  ///< Force x into [0, dim-1]
  template < typename T, typename TInt >
  inline const T rangecrop(const T & x, TInt dim) { return crop<T>(x, 0, dim-1); }
  ///< Absolute value of a number
  // NB: I commented it out because std::abs should be used instead, it may be faster because working at the bit level.
  //template < typename T > INLINE const T abs(const T &a) { return max(a, -a); }
  ///< Euclidean norm of 2D vector (a, b)
  template < typename T > INLINE double norm(const T &a, const T &b) { return std::sqrt(double(a*a+b*b)); }
  ///< Euclidean norm of 3D vector (a, b, c)
  template < typename T > INLINE double norm(const T &a, const T &b, const T &c) { return std::sqrt(double(a*a+b*b+c*c)); }


  template < typename T >
  inline
  T integer_pow(T x, unsigned int n)
  {
    T res = (n%2 ? x : 1);
    while (n >>= 1)
    {
      x *= x;
      if (n%2)
      {
        res *= x;
      }
    }
    return res;
  }

} // namespace til

#endif /*MISC_SCALAR_FUNCTIONS_H_*/

