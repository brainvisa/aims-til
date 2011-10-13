#ifndef _GLOBALTRAITS_H_
#define _GLOBALTRAITS_H_

// includes from STL
#include <algorithm>  // transform
#include <cassert>
#include <cmath>      // sqrt
#include <functional> // mul, bind2nd
#include <limits>     // numeric_limits
#include <vector>

// includes from BOOST
#include <boost/array.hpp>
#include <boost/call_traits.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

// includes from TIL
#include "til/til_declarations.h"
#include "til/is_traits.h"
#include "til/templateTools.h"
#include "til/TExpr.h"
#include "til/traits.h"

// includes from TIL
#include "stditerator.h"


/// A macro to apply a macro for all numerical types
/// (and one thought templates would obsolete this kind of crap ;)
#define TIL_FOR_ALL_NUMERICAL_TYPES(macro)  \
macro (double);                             \
macro (float);                              \
macro (long);                               \
macro (int);                                \
macro (short);                              \
macro (unsigned short);                     \
macro (char);                               \
macro (unsigned char);


namespace til
{
  
  namespace functor
  {
    struct Size
    {
      typedef std::size_t result_type;
  
      template < typename TCollection >
      std::size_t operator()(const TCollection & c)
      {
        return size(c);
      }
    };
  }
  
  namespace detail
  {
    template < class TContainer1, class TContainer2, typename TBinaryFunctor >
    void
    loop_cc(const TContainer1 & c1, TContainer2 &res, TBinaryFunctor f)
    {
      assert(size(c1) == size(res));
      typename TContainer1::const_iterator iC1 = c1.begin();
      typename TContainer2::iterator iRes = res.begin();
      for (; iC1 != c1.end(); ++iC1, ++iRes)
      {
        loop(*iC1, *iRes, f);
      } 
    }
  }
  
  
  
  
  /// Apply a functor to two sparse vectors only where they both have data.
  template < typename TMap1, typename TMap2, typename TFunctor >
  class Loop_mapPairs
  {
  public: // operators
  
    void operator()(const TMap1 & map1, const TMap2 & map2, TFunctor f)
    {
      typename TMap1::const_iterator i1 = map1.begin();
      typename TMap2::const_iterator i2 = map2.begin();
  
      for(;;)
      {
        if (i1->first < i2->first)
        {
          if (++i1 == map1.end()) return;        
        } 
        else if (i1->first > i2->first)
        {
          if (++i2 == map2.end()) return;        
        }
        else
        {
          //res.insert(res.end(), std::make_pair(i1->first, i1->second * i2->second));
          //f(i1->first, i1->second, i2->second);
          f(*i1, *i2);
          if (++i1 == map1.end()) return;
          if (++i2 == map2.end()) return;
        }
      }
    }
    
  private: // checks
  
    typedef typename boost::enable_if<is_map<typename naked_type<TMap1>::type > >::type CheckTMap1isMap;
    typedef typename boost::enable_if<is_map<typename naked_type<TMap2>::type > >::type CheckTMap2isMap;
  };
  
  
  
  template < typename TMap1, typename TMap2, typename TFunctor >
  //typename boost::enable_if<boost::is_stateless<TFunctor> >::type
  typename boost::enable_if<boost::is_empty<TFunctor> >::type
  loop_mapPairs(TMap1 & map1, TMap2 & map2, TFunctor f)
  {
    Loop_mapPairs<TMap1, TMap2, TFunctor>()(map1, map2, f);
  }
  
  
  template < typename TMap1, typename TMap2, typename TFunctor >
  //typename boost::disable_if<boost::is_stateless<TFunctor> >::type
  typename boost::disable_if<boost::is_empty<TFunctor> >::type
  loop_mapPairs(TMap1 & map1, TMap2 & map2, TFunctor & f)
  {
    Loop_mapPairs<TMap1, TMap2, TFunctor&>()(map1, map2, f);
  }
  
  
  /// Apply a functor to a pair of sparse vectors where either one has data.
  template < typename TMap1, typename TMap2, typename TFunctor >
  class Loop_mapEach
  {
  public: // operators
    
    //TODO: I think it should really use iterator arguments instead. And those iterators probably
    // should be private data so that no argument passing between functions
    void operator()(TMap1 & map1, TMap2 & map2, TFunctor f)
    {
      typename TMap1::const_iterator i1 = map1.begin();
      typename TMap2::const_iterator i2 = map2.begin();
  
      if (i1 == map1.end()) goto map1over_checkmap2;
      if (i2 == map2.end()) goto map2over;
  
      // NB: not using goto's would me to have additional end checks at each if's. This function is
      // very low-level thus time critical, so let's stick with goto's.      
      for (;;)
      {
        if (i1->first < i2->first)
        {
          //res.insert(res.end(), std::make_pair(i1->first, f(i1->second, T())));
          //f(i1->first, i1->second, T());
          // TODO: This is probably no good. Imagine that mapped_type is something heavy to allocate.
          // Probably better to instanciate a default somewhere. Policy here?
          f(*i1, typename TMap2::mapped_type());
          if (++i1 == map1.end())
            //return this->map1over(map2, f);
            goto map1over;
        }
        else if (i1->first > i2->first)
        {
          //res.insert(res.end(), std::make_pair(i2->first, f(i2->second, T())));
          //f(i2->first, i2->second, T());
          f(typename TMap1::mapped_type(), *i2);
          if (++i2 == map2.end())
            //return this->map2over(map1, f);
            goto map2over;
        }
        else
        {
          //res.insert(res.end(), std::make_pair(i1->first, f(i1->second, i2->second)));
          //f(i1->first, i1->second, i2->second);
          f(*i1, *i2);
          if (++i1 == map1.end())
          {
            ++i2; 
            //return this->map1over_checkmap2(map2, f);
            goto map1over_checkmap2;
          }
          if (++i2 == map2.end()) 
            //return this->map2over(map1, f);
            goto map2over;
        }
      }
  
      return;
      
  map1over_checkmap2:
      if (i2 == map2.end()) return;
  
  map1over:
      do
      {
        f(typename TMap1::mapped_type(), *i2);
      }
      while (++i2 != map2.end());
      return;
  
  map2over:
      do
      {
        f(*i1, typename TMap2::mapped_type());
      }
      while (++i1 != map1.end());
      return;
    }

  /*
  private: // functions
  
    void map1over_checkmap2(TMap2 & map2, TFunctor f)
    {
      if (++i2 == map2.end()) return;
      this->map1over(map2, f);
    }
    
    void map1over(TMap2 & map2, TFunctor f)
    {
      do
      {
        f(typename TMap1::mapped_type(), *i2);
      }
      while (++i2 != map2.end());
    }

    void map2over(TMap1 & map1, TFunctor f)
    {
      do
      {
        f(*i1, typename TMap2::mapped_type());
      }
      while (++i1 != map1.end());
    }
    */
  };
  
  
  template < typename TMap1, typename TMap2, typename TFunctor >
  //typename boost::enable_if<boost::is_stateless<TFunctor> >::type
  typename boost::enable_if<boost::is_empty<TFunctor> >::type
  loop_mapEach(TMap1 & map1, TMap2 & map2, TFunctor f)
  {
    Loop_mapEach<TMap1, TMap2, TFunctor>()(map1, map2, f);
  }
  
  
  template < typename TMap1, typename TMap2, typename TFunctor >
  //typename boost::disable_if<boost::is_stateless<TFunctor> >::type
  typename boost::disable_if<boost::is_empty<TFunctor> >::type
  loop_mapEach(TMap1 & map1, TMap2 & map2, TFunctor & f)
  {
    Loop_mapEach< TMap1, TMap2, TFunctor & >()(map1, map2, f);
  }
  
  
  
  
  template < typename TMap1, typename TMap2, typename TFunctor >
  /*typename boost::enable_if_c<
    is_map<TMap1>::value && is_map<TMap2>::value
  >::type*/
  class Find_mapEach
  {
  public:
    void operator()(TMap1 &map1, TMap2 &map2, TFunctor f)
    {
      typename TMap1::const_iterator i1 = map1.begin();
      typename TMap2::const_iterator i2 = map2.begin();
    
      // NB: not using goto's would me to have additional end checks at each if's. This function is
      // very low-level thus time critical, so let's stick with goto's.      
      for (;;)
      {
        if (i1->first < i2->first)
        {
          //res.insert(res.end(), std::make_pair(i1->first, f(i1->second, T())));
          //f(i1->first, i1->second, T());
          // TODO: This is probably no good. Imagine that mapped_type is something heavy to allocate.
          // Probably better to instanciate a default somewhere. Policy here?
          if (f(*i1, typename TMap2::mapped_type())) return;
          if (++i1 == map1.end()) goto map1over;
        }
        else if (i1->first > i2->first)
        {
          //res.insert(res.end(), std::make_pair(i2->first, f(i2->second, T())));
          //f(i2->first, i2->second, T());
          f(typename TMap1::mapped_type(), *i2);
          if (++i2 == map2.end()) goto map2over;
        }
        else
        {
          //res.insert(res.end(), std::make_pair(i1->first, f(i1->second, i2->second)));
          //f(i1->first, i1->second, i2->second);
          f(*i1, *i2);
          if (++i1 == map1.end()) goto map1over_checkmap2;
          if (++i2 == map2.end()) goto map2over;
        }
      }
      return;
    
    map1over_checkmap2:
      
      if (++i2 == map2.end()) return;
    
    map1over:
      do
      {
        f(typename TMap1::mapped_type(), *i2);
      }
      while (++i2 != map2.end());
      return;
    
    map2over:
      do
      {
        f(*i1, typename TMap2::mapped_type());
      }
      while (++i1 != map1.end());
      return;
    }
  };
  
  
  // Actually this should not be made available as soone as type T has operator*
  // This is precisely because operator* sucks (return value) that this mul
  // functions are introduced.
  template < typename T1, typename T2, typename T3 >
  inline
  typename boost::enable_if_c<
    std::numeric_limits<T1>::is_specialized &&
    std::numeric_limits<T2>::is_specialized
  >::type
  sub(T1 x, T2 y, T3& z)
  {
    z = x - y;
  }
  
  template < typename T >
  inline
  typename boost::enable_if_c<std::numeric_limits<T>::is_specialized>::type
  mul(T &x, T y)
  {
    x *= y;
  }
  
  // Actually this should not be made available as soone as type T has operator*
  // This is precisely because operator* sucks (return value) that this mul
  // functions are introduced.
  template < typename T1, typename T2, typename T3 >
  inline
  typename boost::enable_if_c<
    std::numeric_limits<T1>::is_specialized &&
    std::numeric_limits<T2>::is_specialized
  >::type
  mul(T1 x, T2 y, T3& z)
  {
    z = x * y;
  }
  
  template < typename T1, typename T2 >
  typename boost::enable_if_c<is_container<T1>::value>::type
  mul(T1 & x, const T2 & y)
  {
    loop(x,y,functor::Mul<T1,T2>());
  }
  
  template < typename T1, typename T2, typename T3 >
  typename boost::enable_if_c<is_container<T1>::value>::type
  mul(const T1 & x, const T2 & y, T3 & z)
  {
    loop(x,y,z,functor::Mul<T1,T2>());
  }
      
  template < typename TContainer >
  inline void resize(TContainer & c, std::size_t n)
  {
    c.resize(n);
  }
  
  template < typename T, std::size_t D >
  inline void resize(boost::array<T, D>, std::size_t n)
  {
    if (D != n)
    {
      throw std::length_error("boost::array can't be resized to a different size");
    }
  }

} // namespace til

#endif //_GLOBALTRAITS_H_
