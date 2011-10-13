#ifndef TIL_LOOP_H
#define TIL_LOOP_H

// includes from BOOST
#include <boost/call_traits.hpp>

// includes from TIL
#include "til/is_traits.h"
#include "til/stditerator.h"

namespace til {

namespace detail
{
  template < class TContainer, typename T, class TBinaryFunctor >
  void
  loop_cv(TContainer & c, T v, TBinaryFunctor f)
  {
    for (typename stditerator<TContainer>::type iC = c.begin(); iC != c.end(); ++iC)
    {
      loop(*iC, v, f);
    }
  }
}

/// Apply a binary functor to a collection and a constant
template < class TContainer, typename T, class TBinaryFunctor >
typename boost::enable_if_c<
  is_container<TContainer>::value &&
  !is_container<T>::value
>::type
loop(TContainer & c, T v, TBinaryFunctor f)
{
  detail::loop_cv(c, v, f);
}

namespace detail
{
  template < typename T1, typename T2, class TBinaryFunctor >
  inline void
  loop_vv(T1 & x, T2 & y, TBinaryFunctor f)
  {
    f(x,y);
  }
}

template < typename T1, typename T2, class TBinaryFunctor >
inline
typename boost::enable_if_c<
  !is_container<T1>::value &&
  !is_container<T2>::value
>::type
loop(T1 & x, T2 & y, TBinaryFunctor f)
{
  detail::loop_vv(x, y, f);
}

namespace detail
{
  template < class TContainer1, class TContainer2, typename TBinaryFunctor >
  void
  loop_cc(TContainer1 & c1, TContainer2 &res, TBinaryFunctor f)
  {
    assert(size(c1) == size(res));
    typename stditerator<TContainer1>::type iC1 = c1.begin();
    typename stditerator<TContainer2>::type iRes = res.begin();
    for (; iC1 != c1.end(); ++iC1, ++iRes)
    {
      loop(*iC1, *iRes, f);
    } 
  }
}

/// Apply a binary functor to two collections
template < class TContainer1, class TContainer2, typename TBinaryFunctor >
typename boost::enable_if_c<
  is_container<TContainer1>::value &&
  is_container<TContainer2>::value
>::type
loop(const TContainer1 & c1, TContainer2 &res, TBinaryFunctor f)
{
  detail::loop_cc(c1, res, f);
}


namespace detail
{
  template < typename T1, typename T2, typename T3, class TBinaryFunctor >
  inline void
  loop_vvv(T1 x, T2 y, T3 & z, TBinaryFunctor f)
  {
    f(x,y,z);
  }
}

template < typename T1, typename T2, typename T3, class TBinaryFunctor >
inline
typename boost::enable_if_c<
  !is_container<T1>::value &&
  !is_container<T2>::value &&
  !is_container<T3>::value
>::type
loop(T1 x, T2 y, T3 & z, TBinaryFunctor f)
{
  detail::loop_vvv(x, y, z, f);
}

namespace detail
{
  template < class TContainer1, class TContainer2, class TContainerRes, class TBinaryFunctor >
  void
  loop_ccc(const TContainer1 & c1, const TContainer2 &c2, TContainerRes &res, TBinaryFunctor f)
  {
    assert(size(c1) == size(c2));
    assert(size(c1) == size(res));
    typename TContainer1::const_iterator iC1 = c1.begin();
    typename TContainer1::const_iterator cEnd = c1.end();
    typename TContainer2::const_iterator iC2 = c2.begin();
    typename TContainerRes::iterator iRes = res.begin();
    for (; iC1 != cEnd; ++iC1, ++iC2, ++iRes)
    {
      loop(*iC1, *iC2, *iRes, f);
    }
  }
}

template < class TContainer1, class TContainer2, class TContainerRes, class TBinaryFunctor >
typename boost::enable_if_c<
  is_container<TContainer1>::value &&
  is_container<TContainer2>::value &&
  is_container<TContainerRes>::value
>::type
loop(const TContainer1 & c1, const TContainer2 &c2, const TContainerRes &res, TBinaryFunctor f)
{
  detail::loop_ccc(c1, c2, res, f);
}

namespace detail
{
  template < class TContainer, class TBinaryFunctor >
  void
  loop_cvc(const TContainer & c1, typename boost::call_traits<typename TContainer::value_type>::param_type v, TContainer & res, TBinaryFunctor f)
  {
    assert(size(c1) == size(res));
    typename TContainer::const_iterator iC1 = c1.begin();
    typename TContainer::iterator iRes = res.begin();
    for (; iC1 != c1.end(); ++iC1, ++iRes)
    {
      loop(*iC1, v, *iRes, f);
    }
    //std::transform(c1.begin(), c1.end(), res.begin(), std::bind2nd(std::minus<typename TContainer::value_type>(), v));
  }
}

template < class TContainer, typename T, class TBinaryFunctor >
typename boost::enable_if_c<
  is_container<TContainer>::value &&
  !is_container<T>::value
>::type
loop(const TContainer & c1, T v, TContainer &res, TBinaryFunctor f)
{
  detail::loop_cvc(c1, v, res, f);
}







/*



  namespace label
  {
    template < std::size_t D >
    struct ecd_fixed  {};
    struct ecd_dense  {};
    struct ecd_sparse {};
    struct ecd_rle    {};
  };


  template < typename T >
  struct loop_label {};
  
  template < typename T, typename D >
  struct loop_label<std::pair<numeric_array<T,D>::iterator,numeric_array<T,D>::iterator> > : public ecd_fixed<D> {};
  
  template < typename T, typename D >
  struct loop_label<numeric_array<T,D>::const_iterator> : public ecd_fixed<D> {};

  namespace detail
  {

    template < typename TFunctor, typename TIterator >
    void loop_standard(TFunctor f, std::pair<TIterator, TIterator> range)
    {
      for (; range.first != range.second; ++range.first)
      {
        f(range.first);
      }
    }
    
    template < typename TFunctor, typename TIterator >
    void newloop(TFunctor f, std::pair<TIterator,TIterator> range, ::til::label::ecd_dense)
    {
      detail::loop_standard(f, range);
    }


    template < typename TFunctor, typename TIterator1, typename TIterator2, std::size_t D >
    void loop(TFunctor f, TIterator i, ecd_fixed<D>, TIterator2 j, ecd_fixed<D>)
    {
      for (std::size_t i = 0; i < D; ++i)
      {
        f(i,j);
      }
    }
  }

*/


} // namespace til

#endif

