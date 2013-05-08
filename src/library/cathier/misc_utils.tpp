
#include <cmath>

namespace til
{
  //---------------------------------------------------------------------------
  
  /// Computes the centroid of a list of vertices.
  template < typename TRes, typename TIterator >
  void mean(TIterator begin, TIterator end, TRes & res)
  {
    res = TRes();
    // TODO: OK, counting is not necessary for random access iterators. Maybe some partial specialization could help.
    std::size_t count = 0;
    for (; begin != end; ++begin)
    {
      res += *begin;
      ++count;
    }  
    // Can't use accumulate because I can't force my functors to cast to othr type
    // there is another reason: a = a + x is inefficient for classes.
    //res = std::accumulate(c.begin(), c.end(), TPoint3D());
  
    //mul(res, static_cast<typename value_type<TRes>::type>(1.0/size(c)));
    res *= static_cast<typename value_type_of<TRes>::type>(1.0/count);
  }

  //---------------------------------------------------------------------------

  template < typename TPoint, typename TPointCollection >
  void stdev(const TPointCollection & c, TPoint & res)
  {
    res = TPoint();
    {
      typename TPointCollection::const_iterator iPoint = c.begin();
      for (; iPoint != c.end(); ++iPoint)
      {
        typename TPoint::iterator i = res.begin();
        typename TPointCollection::value_type::const_iterator j = iPoint->begin();
        for (; i != res.end(); ++i, ++j)
        {
          *i += square(*j);
        }
      }
    }  
    {
      typename TPoint::iterator i = res.begin();
      for (; i != res.end(); ++i)
      {
        *i *= 1.0/size(c);
      }
    }
    //mul(res, 1.0/size(c));
    //res.data() *= 1.0/size(c);
    TPoint center;
    centroid(c, center);
    res[0] = std::sqrt(res[0] - square(center[0]));
    res[1] = std::sqrt(res[1] - square(center[1]));
    res[2] = std::sqrt(res[2] - square(center[2]));
  }

  //---------------------------------------------------------------------------

  // TODO: add a test so that it works only for signed integers... or is >> working only on signed integers anyway?
  template < typename T >
  bool is_dyadic(T n)
  {
    if (n < 2) return false;
    while (n >>= 1)
    {
      if (n%2 == 1 && n != 1) return false;
    }
    return true;
  };

  //---------------------------------------------------------------------------

  /// Returns the largest m so that n >= 2^m.
  /// Only for signed integer... TODO: add a check for that.
  // NB: reason for using int: well, because int is largely sufficient to cover all possible
  // exponent, so I don't want to put extra computing e.g log2(i)-2.
  template < typename T >
  inline int log2(T n)
  {
    int res = 0;
    for (; n > 0; n >>= 1) ++res;
    return res;
  }

  //---------------------------------------------------------------------------

  // TODO: problem is, if an int is passed, it is converted silently into an unsigned int... so the security
  // of using an unsigned argument is kind of lost... Better accept signed input and check for sign explicitely
  // inside the function it seems :/
  template < typename T >
  inline T exp2(unsigned int n)
  {
    T res = 1;
    for (; n > 0; --n) res *= 2;
    return res;
  }
  /*
  {
    int res = 1;
    while (--n) res <<= 1;
    return res;
  }
  */

  //---------------------------------------------------------------------------

  template < typename TIndexCollection >
  //typename boost::enable_if<boost::is_same<typename TIndexCollection::value_type::value_type, std::size_t>,
    std::vector<std::list<std::size_t> >
  //>::type
  invertIndices(const std::vector<TIndexCollection> & c)
  {
    // Do a first pass to get maximum index
    std::size_t mx = 0;
    for (std::size_t i = 0; i < size(c); ++i)
      for (std::size_t j = 0; j < size(c[i]); ++j)
        mx = max(c[i][j], mx);

    // allocate result
    std::vector<std::list<std::size_t> > res(mx+1);
    
    // Fill result
    for (std::size_t i = 0; i < size(c); ++i)
      for (std::size_t j = 0; j < size(c[i]); ++j)
        res[c[i][j]].push_back(i);
    
    return res;
  }

  //---------------------------------------------------------------------------

  template < typename T >
  inline std::size_t max_index(T x0, T x1, T x2)
  {
    if (x0>=x1) {
      if (x0>=x2)
        return 0;
      else
        return 2;
    } else {
      if (x1>=x2)
        return 1;
      else
        return 2;
    }
  }

  //---------------------------------------------------------------------------

  template < typename T >
  inline std::size_t min_index(T x0, T x1, T x2)
  {
    if (x0<=x1) {
      if (x0<=x2)
        return 0;
      else
        return 2;
    } else {
      if (x1<=x2)
        return 1;
      else
        return 2;
    }
  }

  //---------------------------------------------------------------------------

  template < typename T >
  inline bool is_nan(T x)
  {
    // TODO: COMPILER_MAINTENANCE : this has to be regularly checked for compilers :(
#ifdef __GNUC__
    // GNU has its own implementation
    return std::isnan(x);
#else
    // The IEEE-compliant way -- but I guess it is compiler-dependant.
    return x != x;
#endif
  }

  //---------------------------------------------------------------------------

  template < typename T >
  inline T lower_dyadic(T n)
  {
    T res = 1;
    while (res <= n) res *= 2;
    return res /= 2;
  };

    
    
} // namespace til

