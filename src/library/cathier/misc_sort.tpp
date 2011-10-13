

namespace til
{
  //---------------------------------------------------------------------------
  
  template < typename TVector2D, typename T >
  inline
  TVector2D
  sortedVector(T a, T b)
  {
    if (a > b) std::swap(a,b);
    TVector2D res;
    res[0] = a;
    res[1] = b;
    return res;
  }

  //---------------------------------------------------------------------------  

  template < typename TVector3D, typename T >
  TVector3D
  sorted_vector(T a, T b, T c)
  {
    til::sort(a, b, c);
    TVector3D res;
    res[0] = a;
    res[1] = b;
    res[2] = c;
    return res;
  }

  //---------------------------------------------------------------------------  

  /// NB: the argument is passed by value *on purpose*.
  template < typename TVector3D >
  TVector3D
  sorted_vector(TVector3D v)
  {
    til::sort(v[0], v[1], v[2]);
    return v;
  }

  //---------------------------------------------------------------------------

  template < typename T >
  inline
  std::pair<T,T>
  make_sorted_pair(T a, T b)
  {
    if (a >= b)
      return std::make_pair(a,b);
    else
      return std::make_pair(b,a);
  }

  //---------------------------------------------------------------------------  

} // namespace til

