#ifndef TIL_SPARSE_VECTOR_H_
#define TIL_SPARSE_VECTOR_H_

// includes from TIL
//#include "til/miscTools.h"

// includes from TIL
#include "sparse_vector_policies.h"
#include "std_wrap.h"
#include "value_proxy.h"

namespace til
{

  //----------------------------------------------------------------------------------

  
  //namespace {
    
    /// Apply a functor to a pair of sparse vectors where either one has data.
    // NB: the design of Loop_mapEach is a bit crappy because it requires the functor to be highly non standard,
    // accepting combinations of pairs/T, plus the assignment is a problem, because when the first argument
    // passed is a number, of course, the container cannot be updated. I wonder if there is a general way
    // to do that.
    // Here, the functor is standard, it is a stl-like functor returning a value, and this value
    // is assigned to the first container. So, restricted use, but cleaner.
    // I also restricted it to sparse_vector class because I want to use the set() functionality, which
    // is useful to check for non-zero entries. TODO: Maybe this should be a policy instead
    template < typename TMap1, typename TMap2, typename TFunctor >
    class Loop_mapEachAssign
    {
    public: // operators
      
      //TODO: I think it should really use iterator arguments instead. And those iterators probably
      // should be private data so that no argument passing between functions
      void operator()(TMap1 & map1, const TMap2 & map2, TFunctor f)
      {
        typename TMap1::iterator i1 = map1.begin();
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
            i1->second = f(i1->second, typename TMap2::mapped_type());
            if (++i1 == map1.end())
              //return this->map1over(map2, f);
              goto map1over;
          }
          else if (i1->first > i2->first)
          {
            //res.insert(res.end(), std::make_pair(i2->first, f(i2->second, T())));
            //f(i2->first, i2->second, T());
            // create a new entry in map1
            map1[i2->first] = f(typename TMap1::mapped_type(), i2->second);
            if (++i2 == map2.end())
              //return this->map2over(map1, f);
              goto map2over;
          }
          else
          {
            //res.insert(res.end(), std::make_pair(i1->first, f(i1->second, i2->second)));
            //f(i1->first, i1->second, i2->second);
            i1->second = f(i1->second, i2->second);
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
          map1[i2->first] = f(typename TMap1::mapped_type(), i2->second);
        }
        while (++i2 != map2.end());
        return;
    
    map2over:
        do
        {
          i1->second = f(i1->second, typename TMap2::mapped_type());
        }
        while (++i1 != map1.end());
        return;
      }
    };
    
    template < typename TMap1, typename TMap2, typename TFunctor >
    //typename boost::enable_if<boost::is_stateless<TFunctor> >::type
    void
    loop_mapEachAssign(TMap1 & map1, TMap2 & map2, TFunctor f)
    {
      Loop_mapEachAssign<TMap1, TMap2, TFunctor>()(map1, map2, f);
    }
  //}
  

  //----------------------------------------------------------------------------------

    //-----------------//
   //  sparse_vector  //
  //-----------------//


  /// A class that mimic the behavior of std::vector but with a storage policy focused on sparse data.
  /// Basically, a sparse vector has a default value, and stores only values that are not equal to this
  /// default value. Conceptually just a std::map with a default value, if you like. But then, everything
  /// starts being more complicated because of the apparent similarity with std::vector...
  /// Note that this is a container -- if you are looking for a more mathematical object, have a look
  /// at SparseVector, the sparse counterpart of Vector.

  template < typename T, typename BaselinePolicy = policy::SVBaseline_Value<T> >
  class sparse_vector
  {  
  private: // classes
  
    typedef value_proxy<sparse_vector<T> >                ValueProxy;
    typedef const_value_proxy<const sparse_vector<T> >    ConstValueProxy;
    //class ValueProxy;
    //class ConstValueProxy;
    
  public: // classes

    /// An iterator on non-default values only
    class sparse_iterator;    
    /// A const iterator on non-default values only
    class sparse_const_iterator;
    /// A const iterator on all values, STL-like
    class const_iterator;
    class iterator;

  public: // typedefs

    typedef sparse_vector<T>                Self;
    typedef std::map<std::size_t, T>        Map;
    //typedef typename Map::key_type          key_type;
    typedef typename Map::mapped_type       value_type;
    typedef T                               const_reference;
    typedef ValueProxy                      reference;
    typedef const T *                       const_pointer;

  public: // constructors & destructor
  
    /// Create a null vector of size 0
    sparse_vector() : 
        m_data()
      , m_size()
      , m_baselinePolicy()
      , m_proxyIndex(std::numeric_limits<std::size_t>::max())
    {}

    /// Create a null vector of length d.
    sparse_vector(std::size_t d) :
        m_data()
      , m_size(d)
      , m_baselinePolicy()
      , m_proxyIndex(std::numeric_limits<std::size_t>::max())
    {}
    
  public: // set & get
  
    /// Read-write access to i-th element.
    /// NB: never use this operator if you are using map iterators at the same time. This is because
    /// operator[] may remove elements from the map and invalidate an existing iterator pointing to the
    /// removed element.
    ValueProxy operator[](std::size_t i)
    {
      assert( i < m_size );
      return ValueProxy(i, *this);
    }    

    // I don't really understand why we should mess with a ConstValueProxy here?
    // We can understand that writing is costly, but reading surely has to be as
    // fast as possible.
    /*
    ConstValueProxy operator[](std::size_t i) const
    {
      assert( i < m_size );
      return ConstValueProxy(i, this);
    }
    */

    /// Get value of i-th element.
    T operator[](std::size_t i) const
    {
      return this->get(i);
    }
  
  public: //functions

    // NB: these functions bring nothing more to the user than the more standard 
    // operator[], so it's advisable to use the later instead. The only reason these
    // are public is that (1) proxies need them, (2) I would rather not mess around
    // with friends if I can, and (3) these functions having the same functionality
    // as operator[], they are as harmless (or harmful). The only reason to set them
    // to private would be to reduce noise in the API (which is a very good reason btw).

    /// Get i-th value.
    // NB: we have to return a T and not a const T & since return value might be created
    // inside the function
    // NB: no non-const counterpart is given, because we don't want to initialize a value
    // to default value.
    T get(std::size_t i) const
    {
      assert( i < m_size);
      typename Map::const_iterator pos ( m_data.find(i) );
      if (pos != m_data.end())
      {
        return pos->second;
      }
      else
      {
        return m_baselinePolicy(i);
      }
    }

    /*
    /// Get a reference on the n-th value
    T & get(std::size_t i)
    {
      assert( i < m_size);
        //return m_data[i];
        
        // NB: The following code is doing the same thing as the commented single line above.
        // There is no performance gain or loss in doing that -- actually the running time
        // seems to be exactly the same. However, the following code runs with multimap, while
        // the previous line doesn't. So I keep it for later comparisons.
        
      typename Map::iterator pos ( m_data.find(i) );
      if (pos != m_data.end())
      {
        return pos->second;
      }
      else
      {
        m_data.insert(std::make_pair(i,value));
      }
    }
    */
        
    /// Set n-th value
    void set(std::size_t i, const T & value)
    {
      assert(i < m_size);
  
      if (value != m_baselinePolicy(i))
      {
        // Set element to non-default value

        //m_data[i] = value;
        
        // NB: The following code is doing the same thing as the commented single line above.
        // There is no performance gain or loss in doing that -- actually the running time
        // seems to be exactly the same. However, the following code runs with multimap, while
        // the previous line doesn't. So I keep it for later comparisons.
        sparse_iterator pos = m_data.find(i);

        if (pos != m_data.end())
        {
          pos->second = value;
        }
        else
        {
          m_data.insert(std::make_pair(i,value));
        }        
      }
      else
      {
        // Remove element set to default value.
        
        sparse_iterator pos = m_data.find(i);
        if (pos != m_data.end())
        {
          m_data.erase(pos);
        }
      }
    }
    

  /*  
  public: // friend functions

    template < typename X, typename TFunctor >
    void applyFunctor_all(const SparseVector<X> &sv1, const SparseVector<X> &sv2, TFunctor & f);

    template < typename X, typename TFunctor >
    void applyFunctor_pairs(const SparseVector<X> &sv1, const SparseVector<X> &sv2, TFunctor & f);
  
  private: // set & get
  */  

    /// Returns the internal data.
    // TODO: put this in private and add friends
    Map       & getMap()       { return m_data;}
    Map const & getMap() const { return m_data;}

    //T getDefaultValue() const { return m_defaultValue; }    
    //void setDefaultValue(typename boost::call_traits<T>::param_type value) { m_defaultValue = value; }

    BaselinePolicy baselinePolicy() const { return m_baselinePolicy; }
    BaselinePolicy & baselinePolicy () { return m_baselinePolicy; }

  public: // functions
  
    /// Alike STL container's "size".
    /// Return the size of the vector (this is *not* the number of non-default elements, but the
    /// real size of the underlying vector).
    std::size_t size() const { return m_size; }

    /// Alike STL container's "clear".
    void clear()
    {
      m_data.clear();
      m_size = 0;
    }

    /// Alike STL container's "empty".
    bool empty() const { return m_size == 0; }

    /// Alike STL container's "erase"
    void erase(const sparse_iterator & i) { m_data.erase(i); }

    ///
    bool is_null() const
    {
      for(typename Map::const_iterator i = m_data.begin(); i != m_data.end(); ++i)
      {
        if ((i->second)!=0) return false;
      }
      return true;
    }
  /*
  public: // operators used for proxy operations -- should be private, have to be public :(
    // All of these operations are highly risky, and should be used only by generic functions that
    // do not know they deal with sparse_vectors
    
    // TODO: check that we really gain in speed by *not* using a proxy class. If not, there is not
    // reason to get stuck with this very unsafe way of doing things...
    
    /// Purposedly not publicly documented.
    // Assign a value to proxy-ed element
    void operator=(typename boost::call_traits<T>::param_type value)
    {
      this->set(m_proxyIndex, value);
      //m_proxyIndex = std::numeric_limits<std::size_t>::max();
    }

    /// Purposedly not publicly documented.
    // Get the value of the proxy-ed element
    operator T ()
    {
      return this->get(m_proxyIndex);
    }
    
  */
  
  public: // iterators
  
    // Iterators on the data. Note that these iterators iterate only on the non-zero
    // values of the vector; plus, they return a pair (index, value), not a simple value.
    // Thus, they are not compatible and alike other vector iterators -- hence the non-standard names.
  
    sparse_iterator       sparse_begin()       { return m_data.begin(); }
    sparse_iterator       sparse_end  ()       { return m_data.end  (); }
    sparse_const_iterator sparse_begin() const { return m_data.begin(); }
    sparse_const_iterator sparse_end  () const { return m_data.end  (); }
             
    const_iterator begin() const { return const_iterator(0, *this); }
    const_iterator end()   const { return const_iterator(this->size(), *this); }
    iterator begin() { return iterator(0, *this); }
    iterator end() { return iterator(this->size(), *this); }

  public: // exceptions
  
    class OutOfRange : public std::out_of_range
    {
    public: // constructors & destructor
      OutOfRange(const std::string & msg) : std::out_of_range(msg) {}
    };
  
  public: // operators
  
    void operator+=(const sparse_vector<T> & a)
    { loop_mapEachAssign(this->getMap(), a.getMap(), std::plus<T>()); }
    void operator-=(const sparse_vector<T> & a)
    { loop_mapEachAssign(this->getMap(), a.getMap(), std::minus<T>()); }
    void operator*=(const sparse_vector<T> & a)
    { loop_mapEachAssign(this->getMap(), a.getMap(), std::multiplies<T>()); }
    void operator/=(const sparse_vector<T> & a)
    { loop_mapEachAssign(this->getMap(), a.getMap(), std::divides<T>()); }

    void operator+=(const T & value)
    { for (typename Map::iterator i = m_data.begin(); i != m_data.end(); ++i) i->second += value; }
    void operator-=(const T & value)
    { for (typename Map::iterator i = m_data.begin(); i != m_data.end(); ++i) i->second -= value; }
    void operator*=(const T & value)
    { for (typename Map::iterator i = m_data.begin(); i != m_data.end(); ++i) i->second *= value; }
        
  private: // data

    // Contains the actual data
    Map m_data;
    // (virtual) size of our vector
    std::size_t m_size;
    // our default values
    BaselinePolicy m_baselinePolicy;
    // index used during proxy access
    std::size_t m_proxyIndex;
  };


  //----------------------------------------------------------------------------------


    //-----------------------------//
   //  sparse_vector::ValueProxy  //
  //-----------------------------//


  /*
  template < typename T, typename BaselinePolicy >
  class sparse_vector<T>::ValueProxy : public value_proxy<sparse_vector<T> >
  {
  public: // typedefs
    typedef value_proxy<sparse_vector<T> > Base;
  
  public: // constructors & destructor    
      ValueProxy(std::size_t i, sparse_vector<T> & sv) : Base(i, sv) {}

  public: // operators
      // unfortunately we have to redefine this function by hand :(
      void operator=(typename boost::call_traits<T>::param_type value) { this->Base::operator=(value); }
  };
  

  template < typename T, typename BaselinePolicy >
  class sparse_vector<T>::ConstValueProxy : public const_value_proxy<const sparse_vector<T> >
  {
  public: // typedefs
    typedef const_value_proxy<const sparse_vector<T> > Base;
  
  public: // constructors & destructor    
      ConstValueProxy(std::size_t i, const sparse_vector<T> * pSparseVector) : Base(i, pSparseVector) {}
  public: // operators
      // unfortunately we have to redefine this function by hand :(
      void operator=(typename boost::call_traits<T>::param_type value) { this->Base::operator=(value); }
  };
  */


  // NB: these have no interest at all; they are just wrappers around Map::iterators
  // so that they are unique and distinguisable from Map::iterators, so that we can
  // specialize stl functions for sparse_vector iterators only.
  // TODO: change this so that it has the look'n'feel of a standard container.
  // Actually, it probably need to be an indexed_iterator, namely an iterator yielding
  // its index, so make it a sparse_indexed_iterator

  // TODO: actually, sparse_iterators should be much smarter than that. In particular, they should
  // know whether to destroy current entry or not if it is set to zero, based on the container
  // policy.
  template < typename T, typename BaselinePolicy >
  class sparse_vector<T,BaselinePolicy>::sparse_iterator : public Map::iterator
  {
  public: // typedefs
    typedef typename Map::iterator Base;
  public: // constructors & destructor
    sparse_iterator() : Base() {}
    sparse_iterator(const Base & b) : Base(b) {}
  };
  
  
  template < typename T, typename BaselinePolicy >
  class sparse_vector<T,BaselinePolicy>::sparse_const_iterator : public Map::const_iterator
  {
  public: // typedefs
    typedef typename Map::const_iterator Base;
  public: // constructors & destructor
    sparse_const_iterator() : Base() {}
    sparse_const_iterator(const Base & b) : Base(b) {}
    sparse_const_iterator(const sparse_iterator & i) : Base(i) {}
  };

  namespace detail
  {
    template < typename TSparseVector, typename Proxy >
    class sparse_vector_iterator_base
    {
    public: // typedefs
      typedef sparse_vector_iterator_base         Self;
      typedef typename TSparseVector::value_type  value_type;
    public: // constructors & destructor
      //const_iterator(Self * pSparseVector) : m_valueProxy(0, pSparseVector) {}
      sparse_vector_iterator_base(std::size_t i, TSparseVector & sv) : m_valueProxy(i, sv) {}
      template < typename XSparseVector, typename XProxy >
      sparse_vector_iterator_base(sparse_vector_iterator_base<XSparseVector,XProxy> it)
        : m_valueProxy(it.proxy()) {}
    public: // set & get
      const Proxy & proxy() const { return m_valueProxy; }
      Proxy & proxy() { return m_valueProxy; }
    public: // operators
      //const ValueProxy & operator*() { return m_valueProxy; }
      //reference operator*() { return *m_valueProxy; }
      void operator++() { ++(m_valueProxy.index()); }
    public: // friends
      //friend bool operator==(const typename sparse_vector<T>::const_iterator & i1, const typename sparse_vector<T>::const_iterator & i2) { return i1.m_valueProxy.index() == i2.m_valueProxy.index(); }
      //friend bool operator!=(const typename sparse_vector<T>::const_iterator & i1, const typename sparse_vector<T>::const_iterator & i2) { return i1.m_valueProxy.index() != i2.m_valueProxy.index(); }
      /*
      template < typename T1, typename P1, typename T2, typename P2 >
      friend bool operator == (const sparse_vector_iterator_base<T1,P1> & i1,const sparse_vector_iterator_base<T2,P2> & i2);
      template < typename T1, typename P1, typename T2, typename P2 >
      friend bool operator != (const sparse_vector_iterator_base<T1,P1> & i1,const sparse_vector_iterator_base<T2,P2> & i2);
      */


      //void operator=(const Self & other) { m_valueProxy = other.m_valueProxy; }


    protected: // data
      Proxy m_valueProxy;
    };

    template < typename T1, typename P1, typename T2, typename P2 >
    inline bool operator ==
    (
     const sparse_vector_iterator_base<T1,P1> & i1,
     const sparse_vector_iterator_base<T2,P2> & i2
    )
    { return i1.proxy().index() == i2.proxy().index(); }

    template < typename T1, typename P1, typename T2, typename P2 >
    inline bool operator !=
    (
     const sparse_vector_iterator_base<T1,P1> & i1,
     const sparse_vector_iterator_base<T2,P2> & i2
    )
    { return i1.proxy().index() != i2.proxy().index(); }

  }

  template < typename T, typename BaselinePolicy >
  class sparse_vector<T,BaselinePolicy>::const_iterator
    : public detail::sparse_vector_iterator_base<const sparse_vector<T, BaselinePolicy>, ConstValueProxy >
  {
  public: // typedefs
    typedef detail::sparse_vector_iterator_base<const sparse_vector<T, BaselinePolicy>, ConstValueProxy > Base;
    typedef const_reference reference;
    typedef const T * pointer;
  public: // constructors
    const_iterator(std::size_t i, const sparse_vector<T,BaselinePolicy> & p) : Base(i,p) {}
    const_iterator(const iterator & it) : Base(it) {}
    const_reference operator*() { return const_reference(this->proxy()); }
    pointer operator->() { return Base::m_valueProxy.operator->(); }
  };

  template < typename T, typename BaselinePolicy >
  class sparse_vector<T,BaselinePolicy>::iterator
    : public detail::sparse_vector_iterator_base<sparse_vector<T,BaselinePolicy>, ValueProxy>
  {
  public: // typedefs
    typedef detail::sparse_vector_iterator_base<sparse_vector<T,BaselinePolicy>, ValueProxy> Base;
    typedef ValueProxy & reference;
    typedef T * pointer;
  public: // constructors
    iterator(std::size_t i, sparse_vector<T,BaselinePolicy> & p) : Base(i,p) {}
    ValueProxy & operator*() { return this->proxy(); }
    pointer operator->() { return Base::m_valueProxy.operator->(); }


    //void operator=(const iterator & i) { return this->Base::operator=(i); }
  };

  /*
  template < typename T >
  class sparse_vector<T>::const_iterator
  {
  public: // typedefs
    typedef const_iterator Self;
  public: // constructors & destructor
    //const_iterator(Self * pSparseVector) : m_valueProxy(0, pSparseVector) {}
    const_iterator(std::size_t i, const Self * pSparseVector) : m_valueProxy(i, pSparseVector) {}
  public: // operators
    //const ValueProxy & operator*() { return m_valueProxy; }
    T operator*() { return T(m_valueProxy); }
    void operator++() { ++(m_valueProxy.index()); }
  public: // friends
    //friend bool operator==(const typename sparse_vector<T>::const_iterator & i1, const typename sparse_vector<T>::const_iterator & i2) { return i1.m_valueProxy.index() == i2.m_valueProxy.index(); }
    //friend bool operator!=(const typename sparse_vector<T>::const_iterator & i1, const typename sparse_vector<T>::const_iterator & i2) { return i1.m_valueProxy.index() != i2.m_valueProxy.index(); }
    friend inline bool operator==(const Self & i1, const Self & i2) { return i1.m_valueProxy.index() == i2.m_valueProxy.index(); }
    friend inline bool operator!=(const Self & i1, const Self & i2) { return i1.m_valueProxy.index() != i2.m_valueProxy.index(); }
  private: // data
    ConstValueProxy m_valueProxy;
  };
  */

  /// Specialized fill for sparse_vector.
  /// Note that the success of the fill on sparse vector is not guarenteed and
  /// depends on the actual baseline policy. Some baseline policies won't compile
  /// this function. Others might throw an exception if the value is invalid.
  template < typename T, typename BaselinePolicy >
  inline void fill(sparse_vector<T,BaselinePolicy> & v, typename boost::call_traits<T>::param_type value)
  {
    v.getMap().clear();
    v.baselinePolicy().setValue(value);
  }  
}


// Package includes
#include "til/sparse_vector_operators.h"
#include "til/sparse_vector_tools.h"

#endif /*SPARSE_VECTOR_H_*/
