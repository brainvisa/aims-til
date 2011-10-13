#ifndef SPARSEVECTOR_H_
#define SPARSEVECTOR_H_

// includes from STL
#include <functional>
#include <limits>
#include <map>

// includes from BOOST
#include "boost/type_traits.hpp"

// includes from TIL
#include "til/sparse_vector.h"
#include "til/templateTools.h"

//includes from TIL
#include "globalTraits.h"

namespace til
{

  /// A mathematical vector, using sparse storage policy.
  template < typename T >
  class SparseVector : public sparse_vector<T>
  {
  public: // typedefs
    
    typedef sparse_vector<T> Base;
    
  public: // constructors & destructor
  
    /// Create a null vector of size 0
    SparseVector() : Base() {}
    /// Create a null vector of length d.
    SparseVector(std::size_t d) : Base(d) {}
    
    // TODO: I am not sure whether we want to enfore the default value to be zero?
  };
  

  namespace policy
  {
      
    /// Label class for StoreFunctor storage policy classes
    class StoreFunctorStoragePolicy_label {};
    
    /// Default storage policy for StoreFunctor.
    /// The default here is assignment. Works with numerical types
    /*
    template < typename T >
    struct DefaultStoragePolicy : public StoreFunctoreStoragePolicy_label
    {
      public: // operators
      void operator()(T & storage, T value, T) const
      {
        storage = value;
      }
    };
    */
    
    /// StoreFunctor storage policy : accumulation.
    template < typename T, typename TAccumulation >
    class Accumulate : public StoreFunctorStoragePolicy_label
    {
    public: // typedef
      typedef Accumulate<T, TAccumulation> Self;
      
    public: // operators
    
      template < typename X >
      void operator()(TAccumulation & storage, T value, X)        const { Self::store(storage, value); }
      template < typename X >
      void operator()(TAccumulation & storage, T value, X, X)     const { Self::store(storage, value); }
      template < typename X >
      void operator()(TAccumulation & storage, T value, X, X, X)  const { Self::store(storage, value); }
  
    private: // functions
      static void store(TAccumulation & storage, T value)
      {
        storage += TAccumulation(value);
      }
    };
    
    /// StoreFunctor storage policy for SparseVectors.    
    template < typename T >
    class AppendLast : public StoreFunctorStoragePolicy_label
    {
    public: // operators
      void operator()(SparseVector<T> & sv, T value, const std::pair<std::size_t, T> & pair) const
      {
        sv.getMap().insert(sv.sparse_end(), std::make_pair(pair.first, value));
      }
      void operator()(SparseVector<T> & sv, T value, const std::pair<std::size_t, T> & pair, T) const
      {
        sv.getMap().insert(sv.sparse_end(), std::make_pair(pair.first, value));
      }
      void operator()(SparseVector<T> & sv, T value, T, const std::pair<std::size_t, T> & pair) const
      {
        sv.getMap().insert(sv.sparse_end(), std::make_pair(pair.first, value));
      }
      void operator()(SparseVector<T> & sv, T value, const std::pair<std::size_t, T> & pair, std::pair<std::size_t, T>) const
      {
        sv.getMap().insert(sv.sparse_end(), std::make_pair(pair.first, value));
      }
    };
  }
  
  
  template < typename TFunctor >
  class IndexWrap
  {
  };
  
  /// Wrap functors so that it can accept pairs as input -- it then apply the functor on the second argument
  template < typename TFunctor, typename TOutput = typename TFunctor::result_type >
  class WrapFunctorForPairs
  {
  public: // typedefs
  
    typedef TOutput return_type;
  
  private: // Tag classes
  
    struct Pair_Tag {};
    struct Default_Tag {};
    
    template < typename _T >
    struct Tag { typedef Default_Tag type; };
    template < typename _TKey, typename _TValue >
    struct Tag<std::pair<_TKey, _TValue> > { typedef Pair_Tag type; };

    template < typename _T > typename Tag<typename naked_type<_T>::type>::type tagOf(_T) { return typename Tag<typename naked_type<_T>::type>::type (); }

  public: // constuctors & destructor

    WrapFunctorForPairs() : m_functor() {}


  public: // operators

    /*
    template < typename T1 >
    TOutput operator()(T1 x1)
    {
      return m_functor(x1);
    }
    
    
    template < typename TK1, typename T1 >
    TOutput operator()(std::pair<TK1, T1> & pair1)
    {
      return m_functor(pair1->second());
    }
   */
   
    template < typename T1, typename T2 >
    TOutput operator()(T1 x1, T2 x2)
    {
      return this->wrap(x1, tagOf(x1), x2, tagOf(x2));
    }

    /*
    template < typename T1, typename T2 >
    TOutput operator()(T1 x1, T2 x2)
    {
      return m_functor(x1, x2);
    }
    */
    
    /*
    template < typename TK1, typename T1, typename T2 >
    TOutput operator()(std::pair<TK1, T1> & pair1, T2 x2)
    {
      return m_functor(pair1->second, x2);
    }

    template < typename T1, typename TK2, typename T2 >
    TOutput operator()(T1 x1, std::pair<TK2, T2> & pair2)
    {
      return m_functor(x1, pair2->second);
    }

    

    template < typename TK1, typename T1, typename TK2, typename T2 >
    TOutput operator()(std::pair<TK1, T1> & pair1, std::pair<TK2, T2> & pair2) const 
    {
    }
    */

  private: // functions
  
    template < typename T1, typename T2 >
    TOutput wrap(T1 x1, Default_Tag, T2 x2, Default_Tag)
    {
      return m_functor(x1, x2);
    }

    template < typename T1, typename T2 >
    TOutput wrap(T1 x1, Default_Tag, T2 & pair2, Pair_Tag)
    {
      return m_functor(x1, pair2.second);
    }

    template < typename T1, typename T2 >
    TOutput wrap(T1 & pair1, Pair_Tag, T2 x2, Default_Tag)
    {
      return m_functor(pair1.second, x2);
    }

    template < typename T1, typename T2 >
    TOutput wrap(T1 & pair1, Pair_Tag, T2 & pair2, Pair_Tag)
    {
      return m_functor(pair1.second, pair2.second);
    }

  private: // data
    TFunctor m_functor;
  };

  /// Label for loop functors
  class LoopFunctor_label {};
  
  /*
  template < typename T >
  struct is_loop_functor : public boost::is_base_of<LoopFunctor_label, T> {};
  */
  
  template < typename TFunctor, typename TOutput = typename naked_type<TFunctor>::type::return_type >
  class LoopWraper : public LoopFunctor_label
  {
  public: // typedefs
  
    typedef TOutput return_type;
    
  public: // constructors & destructors

    LoopWraper() : m_functor() {}
    LoopWraper(TFunctor & f) : m_functor(f) {}

  public: // operators
  
    template < typename T1 >
    TOutput operator()(T1 x1)
    { return m_functor(x1); }

    template < typename T1, typename T2 >
    TOutput operator()(T1 x1, T2 x2)
    { return m_functor(x1, x2); }

    template < typename T1, typename T2, typename T3 >
    TOutput operator()(T1 x1, T2 x2, T3 x3)
    { return m_functor(x1, x2, x3); }

  public: // functions
    bool proceed() { return true; }
  private:
    TFunctor m_functor;
  };

  template < typename TFunctor >
  typename boost::enable_if<boost::is_stateless<TFunctor>, LoopWraper<TFunctor> >::type
  loopWraper(TFunctor f) { return LoopWraper<TFunctor>(); }

  template < typename TFunctor >
  typename boost::disable_if<boost::is_stateless<TFunctor>, LoopWraper<TFunctor&> >::type
  loopWraper(TFunctor &f) { return LoopWraper<TFunctor&>(f); }

  
  /// A functor wrapper to transform a functor with output in an inplace functor.
  /// Note that the idea here is to have a silent storage, i.e. StoreFunctor uses the same
  /// number of arguments than the functor it wraps. However, each time it is called, the
  /// functor store the value away, in a manner that depends on the storage type and on the
  /// storage policy on that type.
  template < typename TFunctor, typename TStorage, typename TStoragePolicy >
  class StoreFunctor
  {
  public: // typedefs
  
    typedef void return_type;
    
  public: // constructors & destructor
    StoreFunctor() : m_functor(), m_data(), m_storagePolicy() {}
    StoreFunctor(const TStorage & init) : m_functor(), m_data(init), m_storagePolicy() {}

  public: // set & get
    const TStorage & get() const { return m_data; }

  public: // operators
  
    /*
    template < typename T1 >
    void operator()(std::size_t pos, T1 x)
    {
      m_data.getMap().insert(m_data.end(), std::make_pair(pos, m_functor(x)));
    }
    
    template < typename T1, typename T2 >
    void operator()(std::size_t pos, T1 x, T2 y)
    {
      m_data.getMap().insert(m_data.end(), std::make_pair(pos, m_functor(x,y)));
    }
    */

    template < typename T1 >
    void operator()(T1 & x1)
    { this->store<T1&>(x1); }

    template < typename T1 >
    void operator()(const T1 & x1)
    { this->store<const T1&>(x1); }

    template < typename T1, typename T2 >
    void operator()(T1 & x1, T2 & x2)
    { this->store<T1&,T2&>(x1,x2); }

    template < typename T1, typename T2 >
    void operator()(T1 & x1, const T2 & x2)
    { this->store<T1&,const T2&>(x1,x2); }

    template < typename T1, typename T2 >
    void operator()(const T1 & x1, T2 & x2)
    { this->store<const T1&,T2&>(x1,x2); }

    template < typename T1, typename T2 >
    void operator()(const T1 & x1, const T2 & x2)
    { this->store<const T1&,const T2&>(x1,x2); }

    template < typename T1, typename T2, typename T3 >
    void operator()(T1 & x1, T2 & x2, T3 & x3)
    {
      m_storagePolicy(m_data, m_functor(x1, x2, x3), x1, x2, x3);
    }

  private: // functions
  
    template < typename T1 >
    void store(T1 x1)
    {
      m_storagePolicy(m_data, m_functor(x1), x1);
    }

    template < typename T1, typename T2 >
    void store(T1 x1, T2 x2)
    {
      m_storagePolicy(m_data, m_functor(x1, x2), x1, x2);
    }
  
  
  private: // data
    TFunctor m_functor;
    TStorage m_data;
    TStoragePolicy m_storagePolicy;
  };
  
  

  
  
  /// Pointwise multiplication.
  template < typename T >
  SparseVector<T> operator*(const SparseVector<T> &sv1, const SparseVector<T> &sv2)
  {
    //StoreFunctor<std::multiplies<T>, SparseVector<T> > f(SparseVector<T>(sv1.size()));
    StoreFunctor<WrapFunctorForPairs<std::multiplies<T> >, SparseVector<T>, policy::AppendLast<T> > f(SparseVector<T>(sv1.size()));
    loop_mapPairs(sv1.getMap(), sv2.getMap(), f);
    return f.get();
  }
  
  /// Pointwise addition
  template < typename T >
  SparseVector<T> operator+(const SparseVector<T> &sv1, const SparseVector<T> &sv2)
  {
    //StoreFunctor<std::plus<T>, SparseVector<T> > f(SparseVector<T>(sv1.size()));
    StoreFunctor<WrapFunctorForPairs<std::plus<T> >, SparseVector<T>, policy::AppendLast<T> > f(SparseVector<T>(sv1.size()));
    loop_mapEach(sv1.getMap(), sv2.getMap(), f);
    return f.get();
  }

  /// Pointwise subtraction
  template < typename T >
  SparseVector<T> operator-(const SparseVector<T> &sv1, const SparseVector<T> &sv2)
  {
    StoreFunctor<WrapFunctorForPairs<std::minus<T> >, SparseVector<T>, policy::AppendLast<T> > f(SparseVector<T>(sv1.size()));
    //StoreFunctor<std::minus<T>, SparseVector<T> > f(SparseVector<T>(sv1.size()));
    loop_mapEach(sv1.getMap(), sv2.getMap(), f);
    return f.get();
  }
  
  /// Pointwise division
  template < typename T >
  SparseVector<T> operator/(const SparseVector<T> &sv1, const SparseVector<T> &sv2)
  {
    //StoreFunctor<std::divides<T>, SparseVector<T> > f(SparseVector<T>(sv1.size()));
    StoreFunctor<WrapFunctorForPairs<std::divides<T> >, SparseVector<T>, policy::AppendLast<T> > f(SparseVector<T>(sv1.size()));
    loop_mapEach(sv1.getMap(), sv2.getMap(), f);
    return f.get();
  }
  
  /// Dot product between two sparse vectors.
  /// The template argument TPrecision is the type on which we wish to make the computations.
  template < typename TPrecision, typename T >
  TPrecision dot(const SparseVector<T> &sv1, const SparseVector<T> &sv2)
  {
    //StoreFunctor<std::multiplies<T>, T > f(0.0);
    StoreFunctor<WrapFunctorForPairs<std::multiplies<TPrecision> >, TPrecision, policy::Accumulate<TPrecision,TPrecision> > f(0.0);
    loop_mapPairs(sv1.getMap(), sv2.getMap(), f);
    return f.get();
  }
   
  /* 
  template < typename T >
  bool operator==(const SparseVector<T> &sv1, const SparseVector<T> &sv2)
  {
  }
  */

  template < typename T >
  std::ostream& operator<<(std::ostream& os, const SparseVector<T> v)
  {
    typename SparseVector<T>::sparse_const_iterator i = v.sparse_begin();
    for (; i != v.sparse_end(); ++i)
    {
      os << "( " << i->first << " : " << i->second << " ) ";
    }
    return os;
  }
  
} // namespace til


#endif /*SPARSEVECTOR_H_*/
