

namespace til
{
  
  //---------------------------------------------------------------------------

  template < typename T >
  T IntegerRoot<T>::operator()(T n, std::size_t N)
  {
    // this split is just to avoid annoying warning that the n<0 test always fails for unsigned types.
    return this->compute(n, N, cat2type<is_signed, T>());      
  }
  
  //---------------------------------------------------------------------------

  template < typename T >
  T IntegerRoot<T>::compute(T n, std::size_t N, label::Passed<is_signed>)
  {
    if (n < 0)
    {
      if (N%2)
      {
        throw NoRoot();
      }
      else
      {
        return -this->solve(-n, N);
      }
    }
    else
    {
      return this->solve(n, N);
    }
  }

  //---------------------------------------------------------------------------

  template < typename T >
  T IntegerRoot<T>::compute(T n, std::size_t N, label::Failed<is_signed>)
  {
    return this->solve(n, N);
  }

  //---------------------------------------------------------------------------

  template < typename T >
  T IntegerRoot<T>::solve(T n, std::size_t N)
  {
    for (std::size_t i = 0;;++i)
    {
      T p = integer_pow(i, N);
      if (p == n) return i;
      else if (p > n) throw NoRoot();
    }
    // should never go here
    throw NoRoot();
  }
  
  //---------------------------------------------------------------------------

  /// Sort a, b, c in decreasing order (a>=b>=c).
  template < typename T >
  void sort(T & a, T & b, T & c)
  {
    if (a>=b) {     // a >= b
      if (a>=c) {   // a >= b c
        if (b>=c) { // a >= b >= c
          return;
        } else {  // a >= c > b
          std::swap(b,c);
        }
      } else {    // c > a >= b
        shift_right(a,b,c,a);
      }
    } else {      // b > a
      if (a<c) {    // b c > a
        if (b<=c) { // c > b > a
          std::swap(a,c);
        } else {    // b > c > a
          shift_right(c,b,a,c);
        }
      } else {    // b > a > c
        std::swap(b,a);
      }
    }
  }

} // namespace til

