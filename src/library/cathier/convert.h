#ifndef _CONVERT_H_
#define _CONVERT_H_

// includes from BOOST
//#include <boost/numeric/conversion/converter.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

// includes from TIL
//#include "til/is_traits.h"
#include "til/templateTools.h"
#include "til/traits.h"

// includes from TIL
#include "convert.h"
#include "globalTraits.h"

namespace til
{
  //////////////////// Conversion functions ////////////////////
  
  /*
  /// Convert from type TFrom to type TTo
  template < typename TTo, typename TFrom >
  inline
  void
  convert(const TFrom &x, TTo &y)
  {
    /// This function juste merely call the adequate conversion function
    typename ConvertTraits<TFrom, TTo>::Tag Tag;
    convert(x, y, Tag());
  }
  
  
  /// Default conversion: static cast
  template < typename TTo, typename TFrom >
  inline
  void
  convert(const TFrom & x, TTo & y, convertTag::Default)
  {
    y = static_cast<TTo>(x);  
  }
  
  
  /// For numerical values: numerical cast
  template < typename TTo, typename TFrom >
  inline
  void
  convert(const TFrom & x, TTo & y, convertTag::Numeric)
  {
    y = boost::numeric::converter<TTo,TFrom>::convert(x);
  }
  
  
  
  template < typename TPoint3DTo, typename TPoint3DFrom >
  inline
  typename boost::enable_if_c<Is_3DPoint<TPoint3DFrom>::value && Is_3DPoint<TPoint3DTo>::value>::type
  convert(const TPoint3DFrom & x, TPoint3DTo & y)
  {
    y[0] = x[0];
    y[1] = x[1];
    y[2] = x[2];
  }
  */
  
  /// default conversion between types: static conversion
  template < typename TTo, typename TFrom >
  inline
  typename boost::enable_if_c<
    // static_cast must work: even the default conversion excludes some cases!
    boost::is_convertible<TFrom,TTo>::value
    // numerical types
    && !(std::numeric_limits<TFrom>::is_specialized || std::numeric_limits<TTo>::is_specialized)
    && !is_same_naked<TFrom, TTo>::value
  >::type
  convert(const TFrom & x, TTo &y)
  {
    y = static_cast<TTo>(x);
  }
  
  /// conversion between numerical types
  struct Convert_numerical
  {
    template < typename TTo, typename TFrom >
    inline void operator()(const TFrom & x, TTo & y) const
    {
      // TODO: this is too slow because it also does range checking.
      //y = boost::numeric::converter<TTo, TFrom>::convert(x);
      y = x;
    }
  };
  
  template < typename TTo, typename TFrom >
  inline
  typename boost::enable_if_c<
       std::numeric_limits<TFrom>::is_specialized
    && std::numeric_limits<TTo>::is_specialized
    // test for gcc3.3
    && !is_same_naked<TFrom, TTo>::value
  >::type
  convert(const TFrom & x, TTo &y)
  {
    Convert_numerical()(x,y);
  }
  
  namespace detail
  {
    template < typename T >
    inline void
    convert_operatorEqual(const T & x, T & y)
    {
      y = x;
    }
  }
  
  /// By default, if input and output types are the same, use operator=.
  /// NB: this could be dangerous
  template < typename T1, typename T2 >
  inline
  typename boost::enable_if<is_same_naked<T1,T2> >::type
  convert(const T1 & x, T2 & y)
  {
    detail::convert_operatorEqual(x,y);
  }
  
    
  namespace detail
  {
    template < int D, typename T1, typename T2 >
    inline
    void convert_fixedLoop(const T1 & x, T2 & y)
    {
      for (int i=0; i<D; ++i) convert(x[i], y[i]);
    }  
  }
  
  
  inline void convert(const Point3df & x, Point3df & y)
  {
    detail::convert_fixedLoop<3>(x,y);
  }
  
  /// TODO: The idea here is to have no convert with a loop; use for_each(..., Convert()) instead
  struct Convert
  {
    template < typename T1, typename T2 >
    void
    operator()(const T1 & x, T2 & y) const
    {
      convert(x,y);
    }
  };
  
  template < typename T1, typename T2 >
  struct Convert2
  {
    void
    operator()(const T1 & x, T2 & y) const
    {
      convert(x,y);
    }
  };
  

  template < typename TTo, typename TFrom >
  TTo convertTo(const TFrom & x)
  {
    TTo y;
    convert(x,y);
    return y;
  }

  //---------------------------------------------------------------------------

  
  // This is just plain crap, it should at least be a recursive stuff...
  // obviously a quick fix...
  template < typename T, typename B >
  void convert_allocate
  (
    std::vector<til::sparse_vector<T,B> > const & m
  , std::vector<std::vector<std::pair<std::size_t, T> > > & v)
  {
    typedef std::vector<til::sparse_vector<T,B> > T1;
    typedef std::vector<std::vector<std::pair<std::size_t, T> > > T2;

    v.resize(m.size());
    typename T1::const_iterator im = m.begin();
    typename T2::iterator iv = v.begin();
    for (; iv != v.end(); ++im, ++iv)
    {
      iv->resize(im->getMap().size());
      typename T1::value_type::Map::const_iterator iim = im->getMap().begin();
      typename T2::value_type::iterator iiv = iv->begin();
      for (; iiv != iv->end(); ++iim, ++iiv)
      {
        *iiv = *iim;
      }
    }
  }
  
  
  //---------------------------------------------------------------------------
  
  /*
  /// Converting a map<key,value> into a vector<pair<key, value> >.
  /// This is very usefull when the part filling the map is over and we are
  /// left with a sparse array that is now constant. Transforming it into a
  /// fixed vector speeds up access a lot.
  // TODO: shouldn't this go into a Convert/Cast class?
  template < typename TKey, typename TMap, typename TComp, typename TAlloc1, typename TAlloc2 >
  void convert(const std::map<TKey, TMap, TComp, TAlloc1> & m, std::vector<std::pair<TKey, TMap>, TAlloc2> & v)
  {
    detail::loop_cc(m, v, Convert());
  }
  */
  
  //---------------------------------------------------------------------------
  
} // namespace til

#endif //_CONVERT_H_

