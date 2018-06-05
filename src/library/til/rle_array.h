#ifndef TIL_RLE_ARRAY_H
#define TIL_RLE_ARRAY_H

// includes from STL
#include <list>

// includes from TIL
#include "til/traits.h"
#include "til/value_proxy.h"

namespace til 
{
  /// A run-length encoded array.
  /// Run-length encoding is a simple way to compress data. It stores data using two
  /// parameters: a value, and the number of time this value is repeated.
  /// The range of usefulness of this scheme is therefore limited to arrays where
  /// successive values have a very high probability to be strictly equal. This is
  /// generally true for smooth segmentation/label images, where this encoding is most 
  /// useful.
  /// Due to the simplicity of its compression scheme, reading from an rle_array is
  /// fast, and extremely so when iterating; writing is a more complex operation, 
  /// but still tractable.
  template < typename T, typename TCount = unsigned int >
  class rle_array
  {
  public: // classes
    
    class sparse_iterator;
    class const_sparse_iterator;
    class iterator;
    class const_iterator;
    class Segment;

 private: // typedefs

    typedef std::list<Segment> Data;
    typedef value_proxy<rle_array<T,TCount>, std::pair<typename Data::iterator, TCount> > ValueProxy;

  public: // typedefs
  
    typedef rle_array<T,TCount>   Self;
    typedef T                     value_type;
    typedef const T &             const_reference;
    typedef const T *             const_pointer;
    typedef ValueProxy            reference;
    typedef TCount                count_type;

  public: // constructors & destructor
  
    /// Default constructor, yield a zero-size container.
    rle_array()
      : m_data()
      , m_size(0)
    {}

    /// Standard constructor providing a size and an optional fill value.
    /// If not value is provided, the default value of type T is used.
    rle_array(std::size_t i, T value = T())
      : m_data(1, Segment(TCount(i), value))
      , m_size(i)
    {}
  
  public: // iterators

    const_iterator begin() const
    { return const_iterator(m_data.begin(), 0, 0); }
    const_iterator end() const
    { return const_iterator(m_data.end(), 0, m_size); }
    iterator begin()
    { return iterator(m_data.begin(), 0, 0, *this); }
    iterator end()
    { return iterator(m_data.end(), 0, m_size, *this); }

  public: // operators

    const T & operator[](std::size_t i) const;
    ValueProxy operator[](std::size_t i);

  public: // set & get for proxies
  
    void set(std::pair<typename Data::iterator, TCount> & index, const T & value);
    template < typename XIterator >
    const T & get(std::pair<XIterator, TCount> & index) const
    { return index.first->value(); }

    std::size_t size() const { return m_size; }
  public: // functions

  void print() const
	{
		for (typename Data::const_iterator iRepValue = m_data.begin(); iRepValue != m_data.end(); ++iRepValue)
		{
			std::cout << "Value : " <<  iRepValue->value() << "  ";
			std::cout << "Repeat : " << iRepValue->length() << ",  ";
		}
		std::cout << std::endl;
	}

  private: // data

    Data          m_data;
    std::size_t   m_size;
  };

  //---------------------------------------------------------------------------

  /// A Segment is composed of a value and a length.
  /// The RLE array is simply a collection of those segments.
  template < typename T, typename TCount >
  class rle_array<T,TCount>::Segment
  {
  public: // typedefs
    typedef TCount count_type;
  public: // constructors & destructor
    
    Segment() : m_value(), m_length(0) {};
    // follow std container convention for container constructors: first length, then value
    Segment(TCount r, const T &v) : m_value(v), m_length(r) {};

  public: // set & get

    // I think a good reason not to return a TCount for const is that if we
    // do ++(x.length()) we can't be sure if the const has not actually been called!
    const TCount & length() const { return m_length; }
    TCount & length() { return m_length; }
    const T & value() const { return m_value; }
    T & value() { return m_value; }
    //void setValue(T v) { m_value = v; }

  private: // data

    /// A value
    T       m_value;
    /// Number of contiguous occurences of the value
    TCount  m_length;
  };

  //---------------------------------------------------------------------------

  template < typename T, typename TCount >
  class rle_array<T,TCount>::const_sparse_iterator : public Data::const_iterator
  {
  public: // typedefs
    typedef typename Data::const_iterator Base;
  public: // constructors & destructor
    //sparse_iterator() : Base() {}
    const_sparse_iterator(const Base & b) : Base(b) {}
  };

  template < typename T, typename TCount >
  class rle_array<T,TCount>::sparse_iterator : public Data::iterator
  {
  public: // typedefs
    typedef typename Data::iterator Base;
  public: // constructors & destructor
    sparse_iterator() : Base() {}
    sparse_iterator(const Base & b) : Base(b) {}
  };

  namespace detail
  {

    //-------------------------------------------------------------------------
    
    // TODO: I think the whole crap can be redesigned so that a value proxy is used
    // for const and non-const without penalty, so that index is always stored
    // inside value_proxy and there is not value_proxy of value& kind of crap.
    template < typename TRLEArray, typename TIterator >
    class rle_array_iterator_base
    {
    public: // typedefs
      
      typedef rle_array_iterator_base<TRLEArray, TIterator>     Self;
      typedef TRLEArray                                         container;
      typedef typename TRLEArray::count_type                    count_type;
      typedef std::pair<TIterator, count_type>                  index_type;
      // NB: the access policy has to pass a reference on the index here, since the index
      // may be modified by the set function.
      typedef value_proxy<TRLEArray, index_type, policy::VPAccess_Default<container, index_type&> >                
                                                                ValueProxy;
      typedef typename container::value_type                    value_type;
      //typedef const T & reference;
      //typedef typename TIterator::reference                     reference;
      //typedef typename TIterator::pointer                       pointer;
      typedef typename ValueProxy &                             reference;
      typedef typename ValueProxy::const_pointer                pointer;

    public: // constructors & destructor

      rle_array_iterator_base(const TIterator & i, count_type segindex, std::size_t pos, container & c)
        : m_valueProxy(index_type(i, segindex), c)
        , m_i(pos)
      {
      }

      template < typename XRLEArray, typename XIterator >
      rle_array_iterator_base(const rle_array_iterator_base<XRLEArray, XIterator> & b)
        : m_valueProxy(b.proxy())
        , m_i(b.pos())
      {
      }

    public: // set & get

      /// Return the current position of the iterator in the array.
      std::size_t pos() const { return m_i; }
      //index_type & index() { return m_index; }
      //const index_type & index() const { return m_index; }
      const ValueProxy & proxy() const { return m_valueProxy; }

    public: // operators

      // NB: It is essential that a reference on ValueProxy is passed, because
      // an assignment might affect the iterator hold by the value proxy, due to
      // list insertions.
      //const T & operator*() { return m_index.first->value(); }
      reference operator*() { return m_valueProxy; }

      //pointer operator->() { return m_index.first->operator->(); }
      pointer operator->() { return m_valueProxy.operator->(); }

      void operator++()
      {
        ++m_i;
        // Are we at the end of a segment?
        //if (++(m_index.second) == m_index.first->length())
        if (++(m_valueProxy.index().second) == m_valueProxy.index().first->length())
        {
          // go to next segment
          //++m_index.first;
          //m_index.second = 0;
          ++(m_valueProxy.index().first);
          m_valueProxy.index().second = 0;
        }
      }

      void print()
      {
        std::cout << "Pointing on " << m_valueProxy.index().first->value() <<" "<< m_valueProxy.index().first->length() << std::endl;
        std::cout << "Local pos: " << m_valueProxy.index().second << std::endl;
        std::cout << "Absolute pos: " << m_i << std::endl;
      }

      /*
      void operator=(const Self & other)
      {
        m_valueProxy = other.m_valueProxy;
        m_i = other.m_i;
      }
      */

    private: // data
      
      ValueProxy m_valueProxy;
      //index_type m_index;
      /// Current position in the array
      std::size_t m_i;
    };
    
    //-------------------------------------------------------------------------

    template < typename TRLEArray1, typename TRLEArray2, typename TIterator1, typename TIterator2 >
    inline
    bool operator==
    (
     const rle_array_iterator_base<TRLEArray1,TIterator1> & i1,
     const rle_array_iterator_base<TRLEArray2,TIterator2> & i2
    )
    { return (i1.proxy().index().first == i2.proxy().index().first) && (i1.proxy().index().second == i2.proxy().index().second); }
    
    //-------------------------------------------------------------------------

    template < typename TRLEArray1, typename TRLEArray2, typename TIterator1, typename TIterator2 >
    inline
    bool operator!=
    (
     const rle_array_iterator_base<TRLEArray1,TIterator1> & i1,
     const rle_array_iterator_base<TRLEArray2,TIterator2> & i2
    )
    { return !(i1 == i2); }

    //-------------------------------------------------------------------------

  } // namespace

  //-------------------------------------------------------------------------

  template < typename T, typename TCount >
  class rle_array<T,TCount>::const_iterator : public detail::rle_array_iterator_base<const rle_array<T,TCount>, typename Data::const_iterator>
  {
  public: // typedefs
    typedef detail::rle_array_iterator_base<const rle_array<T,TCount>, typename Data::const_iterator> Base;
    typedef typename Base::value_type   value_type;
    typedef typename Base::reference    reference;
    typedef typename Base::pointer      pointer;
  public: // constructors & destructor
    const_iterator(const typename Data::const_iterator & i, TCount segindex, std::size_t pos, const rle_array<T,TCount> & c)
      : Base(i, segindex, pos, c)
    {}
    /// This is only to convert from an rle_array<T,TCount>::iterator, so that one can
    /// write const_iterator i = a.begin() even when a is non-const. Convertions from
    /// other types are unsafe and probably won't compile.
    template < typename XRLEArray, typename XIterator >
    const_iterator(const detail::rle_array_iterator_base<XRLEArray,XIterator> & b) : Base(b) {}
  };

  template < typename T, typename TCount >
  class rle_array<T,TCount>::iterator : public detail::rle_array_iterator_base<rle_array<T,TCount>,typename Data::iterator>
  {
  public: // typedefs
    typedef detail::rle_array_iterator_base<rle_array<T,TCount>,typename Data::iterator> Base;
    typedef typename Base::value_type   value_type;
    typedef typename Base::reference    reference;
    typedef typename Base::pointer      pointer;
  public: // constructors & destructor
    iterator(const typename Data::iterator & i, TCount segindex, std::size_t pos, rle_array<T,TCount> & c)
      : Base(i, segindex, pos, c)
    {}


    //void operator=(const iterator & other) { this->Base::operator=(other); }

  };

} // namespace til

// package include
#include "rle_array.tpp"

#endif

