#ifndef TIL_SANDBOX_H
#define TIL_SANDBOX_H

namespace til
{

  /*
  template < typename TContainer >
  struct All
  {
    All(const TContainer & c) : m_c(c) {}
    const TContainer & m_c;
  };

  template < typename TContainer >
  inline All<TContainer> all(const TContainer & c) { return All<TContainer>(c); }

  template < typename TContainer1, typename TContainer2 >
  inline bool operator<(const All<TContainer1> & c1, const All<TContainer2> & c2)
  {
  }
  */


  namespace 
  {
    // Small helper for the following
    template < typename T >
    inline void minto(T & x, T y)
    {
      if (x>y) x=y;
    }
  }


  template < typename TRange >
  inline
  typename boost::enable_if<is_range<TRange>, typename TRange::value_type>::type
  min(TRange r)
  {
    typedef typename TRange::value_type value_type;
    // We don't want to use numeric_limits because this function might be useful for
    // non numerical types.
    if (!r.ok()) return value_type();
    value_type res = *r;
    for (; r.ok(); ++r)
    {
      // The function call is simply to avoid calling *r twice, as well
      // as to explicitely creating a temporary variable if not needed.
      // TODO: helper func is not working, maybe there is a pb in the typedefs
      // in range.
      if (res > *r) res = *r;
    }
    return res;
  }
}


#endif

