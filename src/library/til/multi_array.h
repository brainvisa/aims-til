#ifndef TIL_MULTI_ARRAY_H
#define TIL_MULTI_ARRAY_H

// includes from STL
#include <cassert>
#include <numeric>

// includes from BOOST
#include "boost/array.hpp"

// includes from TIL
#include "til/basic_range.h"
#include "til/labels.h"
#include "til/meta.h"
#include "til/numeric_array.h"
#include "til/traits.h"

namespace til
{

  // Idea: we first define a multi array that can combine containers, i.e.
  // have a container of container (or not), but in anycase that is 2D.
  // then, in higher dim, simply iterate?
  //template < typename TContainer, std::size_t D >
  //class multi_conc;
/*
  namespace detail
  {
    template < typename TMultiArray, std::size_t D >
    class multi_conc_two
    {
    public: // typedef

      typedef TMultiArray                                 outside_container;
      typedef typename value_type_of<TMultiArray>::type   inside_container;
    
    public: // classes

      class range;
      class const_range;

    public: // constructors & destructor

      multi_conc_two() : m_data() {}
      multi_conc(const boost::array<std::size_t, 2> & dims)
        : m_data(dims[0], inside_container(dims[1]))
      {
      }
    
    public: // range

      range whole_range() { return range(til::whole_range(m_data), til::whole_range(*til::whole_range(m_data))); };

      //range begin() { return range(m_data.begin(), m_data.front().begin()); }
      //range end()  { return iterator(m_data.end(), m_data.back().end()); }

    public: // functions

      boost::array<std::size_t, 2> dims() const
      { return m_dims; }

    private: // data
      TContainer m_data;
      boost::array<std::size_t, 2> m_dims;
    };
  } // namespace detail

  template < typename TContainer, std::size_t D >
  class multi_conc : public detail::multi_conc_two<TContainer, D >
  {
  };

  template < typename TContainer, std::size_t D >
  const int multi_conc<TContainer, D>::d = D;
*/
/*
  template < typename TContainer >
  class multi_conc<TContainer, 1>
  {
  public: // typedefs
    
    typedef typename TContainer::value_type       value_type;
    typedef typename TContainer::reference        reference;
    typedef typename TContainer::const_reference  const_reference;
  
  public: // constants

    static const int d;

  public: // classes

    class range : public basic_range<typename TContainer::iterator>
    {
    private: // typedefs
      typedef typename TContainer::iterator iterator;
    public: // typedefs
      typedef basic_range<iterator> Base;
    public: // constructors
      range(iterator begin, iterator end) : Base(begin,end) {}
    };
    
    class const_range : public basic_range<typename TContainer::const_iterator>
    {
    private: // typedefs
      typedef typename TContainer::const_iterator iterator;
    public: // typedefs
      typedef basic_range<iterator> Base;
    public: // constructors
      const_range(iterator begin, iterator end) : Base(begin,end) {}
    };
  
  public: // constructors

    multi_conc() : m_data() {}
    explicit multi_conc(std::size_t d) : m_data(d) {}
    multi_conc(std::size_t d, value_type v) : m_data(d,v) {}
    multi_conc(const TContainer & data) : m_data(data) {}

  public: // ranges
    range whole_range() { return range(m_data.begin(), m_data.end()); }
    const_range whole_range() const { return range(m_data.begin(), m_data.end()); }

  public: // operators

    reference operator[](std::size_t d) { return m_data[d]; }
    const_reference operator[](std::size_t d) const { return m_data[d]; }

    reference operator()(std::size_t d) { return m_data[d]; }
    const_reference operator()(std::size_t d) const { return m_data[d]; }

  private: // data

    TContainer m_data;
  };

  */
/*
  template < typename TContainer >
  class multi_conc<TContainer, 2>
  {
  public: // typedef

    typedef TContainer outside_container;
    // TODO: check that this is actually a container?
    typedef typename value_type_of<TContainer>::type inside_container;
    typedef typename TContainer::reference reference;
    typedef typename TContainer::const_reference const_reference;
  public: // classes

    class range;
    class const_range;

  public: // constructors & destructor

    multi_conc() : m_data() {}
    multi_conc(const numeric_array<std::size_t, 2> & dims)
      : m_data(dims[0], inside_container(dims[1]))
      , m_dims(dims)
    {
      assert(size(m_data) == dims[0]);
    }
  
  public: // range

    range whole_range() { return range(til::whole_range(m_data), til::whole_range(*til::whole_range(m_data))); };

    //range begin() { return range(m_data.begin(), m_data.front().begin()); }
    //range end()  { return iterator(m_data.end(), m_data.back().end()); }

  public: // operators

    typename TContainer::reference operator[](std::size_t i) { return m_data[i]; }
    typename TContainer::const_reference operator[](std::size_t i) const { return m_data[i]; }

    reference operator()(std::size_t i, std::size_t j) { return m_data[i][j]; }
    const_reference operator()(std::size_t i, std::size_t j) const { return m_data[i][j]; }

  public: // set & get

    boost::array<std::size_t, 2> dims() const
    { return m_dims; }

    std::size_t size() const { return m_dims[0]*m_dims[1]; }

  private: // data
    TContainer m_data;
    numeric_array<std::size_t, 2> m_dims;
  };

*/
  namespace detail
  {
  /*
    /// Composition of two ranges.
    template < typename TOuterRange, typename TInnerRange >
    class range_compo : public range_label
    {
    public: // typedefs

      typedef TOuterRange outer_range;
      typedef TInnerRange inner_range;
      typedef typename inner_range::reference reference;
      typedef typename inner_range::value_type value_type;

    public: // constructors

      range_compo(const TOuterRange & outer, const TInnerRange & inner) : 
            m_outer(outer), m_inner(inner)
      {}

    public: // operators
      
      reference operator*() { return *m_inner; }

      void operator++() { ++m_inner;};

      bool ok()
      {
        if (m_inner.ok()) return true;
        ++m_outer;
        if (m_outer.ok()) 
        {
          m_inner = whole_range(*m_outer);
          return true;
        }
        return false;
      }

    public: //friends

      template < typename T1, typename T2, typename U1, typename U2 >
      friend bool operator!=(const range_compo<T1,T2> & i1, const range_compo<U1,U2> & i2);

    private: // data
      outer_range m_outer;
      inner_range m_inner;
    };

    template < typename T1, typename T2, typename U1, typename U2 >
    bool operator!=(const range_compo<T1,T2> & i1, const range_compo<U1,U2> & i2)
    {
      return (i1.m_inner != i2.m_inner);
    }
    */

    //----------------------------------------------------------------------------

    /*
    template < typename TIterator >
    class volumetric_iterator_base
    {
    public: // operators
      reference operator*() { return *m_inner; }
      void operator++() { ++m_inner; }
    protected: // data
      TIterator m_i;
    };
    */

    //----------------------------------------------------------------------------

      //---------------------------//
     //  inner_voliterator_compo  //
    //---------------------------//

    template < typename TOuterIterator, std::size_t D, typename TInnerVolIterator >
    class inner_voliterator_compo
    {
    public: // typedefs
      typedef TOuterIterator outer_iterator;
      typedef TInnerVolIterator inner_iterator;
      typedef typename inner_iterator::reference reference;
      typedef typename inner_iterator::value_type value_type;
      typedef meta::int_type<D> dim;
    protected: // typedefs
      typedef numeric_array<std::size_t,D> Vector;
    public: // constructors
      inner_voliterator_compo(TOuterIterator i, const Vector & size)
        : m_i(i)
        , m_size(size)
        , m_inner(i->begin())
      {}
    public: // operators
      reference operator*() { return *m_inner; }
      void operator++() { ++m_inner; }
    public: // functions
      template < typename TPosIterator >
      void from_pos(TOuterIterator begin, TPosIterator pbegin, TPosIterator pend)
      {
        m_i = begin;
        TPosIterator pmid = pbegin + D;
        std::advance(m_i, pos2offset(pbegin, pmid, m_size.begin()));
        m_inner.from_pos(m_i->begin(), pmid, pend);
      }
    private: // data
      TOuterIterator m_i;
      numeric_array<std::size_t,D> m_size;
      TInnerVolIterator m_inner;
    };
  }
  
  //----------------------------------------------------------------------------

    //---------------------//
   //  voliterator_compo  //
  //---------------------//

  template < typename TOuterIterator, std::size_t D, typename TInnerVolIterator >
  class voliterator_compo
    : public detail::inner_voliterator_compo<TOuterIterator,D,TInnerVolIterator>
  {
  public: // typedefs
    typedef detail::inner_voliterator_compo<TOuterIterator,D,TInnerVolIterator> Base;
    typedef typename Base::outer_iterator  outer_iterator;
    typedef typename Base::inner_iterator  inner_iterator;
    typedef typename Base::reference       reference;
    typedef typename Base::value_type      value_type;
  protected: // typedefs
    typedef typename Base::Vector          Vector;
  public: // constructors
    voliterator_compo(TOuterIterator begin, const Vector & size)
      : Base(begin, size)
      , m_begin(begin)
    {}
  public: // functions
    template < std::size_t D2 >
    void from_pos(const numeric_array<std::size_t,D2> & pos)
    { this->from_pos(pos.begin(), pos.end()); }
    template < typename TPosIterator >
    void from_pos(TPosIterator pbegin, TPosIterator pend)
    { this->from_pos(m_begin, pbegin, pend); }
  private: // data
    TOuterIterator m_begin;
    TOuterIterator m_i;
    numeric_array<std::size_t,D> m_size;
    TInnerVolIterator m_inner;
  };

  namespace detail
  {
    //----------------------------------------------------------------------------

    /// Composition of an iterator and a range
    template < typename TOuterIterator, typename TInnerRange >
    class iterator_range_compo : public range_label
    {
    public: // typedefs

      typedef TOuterIterator outer_iterator;
      typedef TInnerRange inner_range;
      typedef typename inner_range::reference reference;
      typedef typename inner_range::value_type value_type;

    public: // constructors

      iterator_range_compo(
        TOuterIterator  begin,
        TOuterIterator  end
        )
        : m_i(begin)
        , m_end(end)
        , m_inner(begin->whole_range())
      {}

    public: // operators
      
      reference operator*() { return *m_inner; }
      void operator++() { ++m_inner; }
      bool ok()
      {
        if (m_inner.ok()) return true;
        ++m_i;
        if (m_i != m_end)
        {
          m_inner = m_i->whole_range();
          return true;
        }
        return false;
      }

    public: //friends

      template < typename T1, typename T2, typename U1, typename U2 >
      friend bool operator!=(const iterator_range_compo<T1,T2> & i1, const iterator_range_compo<U1,U2> & i2);

    private: // data

      outer_iterator m_i;
      outer_iterator m_end;
      inner_range m_inner;
    };


    template < typename T1, typename T2, typename U1, typename U2 >
    bool operator!=(const iterator_range_compo<T1,T2> & i1, const iterator_range_compo<U1,U2> & i2)
    {
      return (i1.m_inner != i2.m_inner);
    }

  } // namespace detail


/*
  template < typename TContainer >
  class multi_conc<TContainer,2>::range
    : public detail::range_compo<typename range_of<TContainer>::type, typename range_of<typename value_type_of<TContainer>::type>::type>
  {
  public: // typedefs
    typedef detail::range_compo<typename range_of<TContainer>::type, typename range_of<typename value_type_of<TContainer>::type>::type> Base;
  public: // constructors
    range(const outer_range & outer, const inner_range & inner) : Base(outer,inner) {}
  };
*/
  /*
  template < typename TIterator >
  inline
  typename TIterator::value_type
  prod(TIterator begin, const TIterator & end)
  {
    typename TIterator::value_type res = 1;
    for (; begin != end; ++begin) res *= *begin;
    return res;
  }
  */

  namespace
  {
    template < typename TIterator1, typename TIterator2 >
    inline std::size_t
    offset_prod(TIterator1 pos, TIterator1 posend, TIterator2 dims)
    {
      std::size_t res = *pos;
      while (++pos != posend) res = *(++dims) * res + *pos;
      return res;
      /*
      ++pos;
      ++dims;
      for(; pos!= posend; ++dims, ++pos) res = *dims * res + *pos;
      return res;
      */
      /*
      do {
        if (++pos == posend) return res;
        ++dims;
        res = *dims * res+ *pos;
      }
      */
    }

    template < typename T, std::size_t D, typename TIterator >
    numeric_array<T,D>
    slice(TIterator i)
    {
      numeric_array<T,D> res;
      std::copy<TIterator>(i, i+D, res.begin());
      return res;
    }
  }

  namespace detail
  {

    //------------------------------------------------------------------------------

      //--------------------//
     //  multi_array_base  //
    //--------------------//


    /// Wraps a container containing multi_array inside a multi_array, bringing
    /// D extra dimensions to the dimensions alreay spanned by the contained multi-array.
    template < typename TContainer, std::size_t D, bool b >
    class multi_array_base : public multi_array_label
    {
    public: // typedefs
      typedef TContainer                                  outer_container;
      typedef typename value_type_of<TContainer>::type    inner_container;
      typedef typename inner_container::value_type        value_type;
      typedef typename inner_container::reference         reference;
      typedef typename inner_container::const_reference   const_reference;
      typedef meta::int_type<D> my_d_type;
      typedef typename meta::add_type<my_d_type, typename inner_container::d_type>::type d_type;
      typedef numeric_array<std::size_t, d_type::value> coord_type;
    private: // helper
      static inline std::size_t size_helper(const coord_type & dims)
      {
        return std::accumulate(dims.rbegin(), dims.rbegin()+D, std::size_t(1), std::multiplies<std::size_t>());
      }
    public: // classes

      class range
        : public iterator_range_compo<typename outer_container::iterator, typename inner_container::range>
      {
      public:
        typedef iterator_range_compo<typename outer_container::iterator, typename inner_container::range> Base;
        typedef typename Base::outer_iterator outer_iterator;
        range(outer_iterator begin, outer_iterator end) : Base(begin, end) {}
      };
      class const_range
        : public iterator_range_compo<typename outer_container::const_iterator, typename inner_container::const_range>
      {
      public:
        typedef iterator_range_compo<typename outer_container::const_iterator, typename inner_container::const_range> Base;
        typedef typename Base::outer_iterator outer_iterator;
        const_range(outer_iterator begin, outer_iterator end) : Base(begin, end) {}
      };
    public: // constructors
      multi_array_base() : m_data() {}
      //explicit multi_array_base(const numeric_array<std::size_t, DTotal> & dims)
      explicit multi_array_base(const coord_type & dims)
        : m_data(size_helper(dims), 
                 inner_container(slice<std::size_t, d_type::value-D>(dims.begin())))
        , m_dims(dims)
      {
        assert(std::accumulate(dims.begin(), dims.end(), std::size_t(1), std::multiplies<std::size_t>())
          == til::size(m_data) * til::size(m_data[0]));
      }
      multi_array_base(const coord_type & dims, value_type v)
        // TODO: look for an alternative to this stupid slice function.
        : m_data(size_helper(dims), 
          inner_container(slice<std::size_t, d_type::value-D>(dims.begin()), v))
        , m_dims(dims)
      {
        assert(std::accumulate(dims.begin(), dims.end(),std::size_t(1),std::multiplies<std::size_t>())
          == size(m_data) * size(m_data[0]));
      }
      multi_array_base(const TContainer & data) : m_data(data) {}
    public: // range
      range whole_range() { return range(m_data.begin(), m_data.end()); }
      const_range whole_range() const { return const_range(m_data.begin(), m_data.end()); }
    public: //  set & get
      coord_type dims() const { return m_dims; }
      std::size_t size() const { return til::size(m_data) * til::size(m_data.front()); } 
      //std::size_t D() const { return d_type::value; }
    public: // operators
      typename TContainer::reference operator[](std::size_t i) { return m_data[i]; }
      typename TContainer::const_reference operator[](std::size_t i) const{ return m_data[i]; }
      reference operator()(const coord_type & p)
      { return this->operator()(p.rbegin(), p.rend()); }
      const_reference operator()(const coord_type & p) const
      { return this->operator()(p.rbegin(), p.rend()); }
      template < typename TIterator >
      reference operator()(TIterator begin, TIterator end)
      {
        return m_data[offset_prod(begin, begin+D, m_dims.rbegin())](begin+D,end);
          //(slice<std::size_t, d_type::value-D>(p.begin()));
          //[offset_prod(p.rbegin()+D, p.rend(), m_dims.rbegin()+D)];
      }
      template < typename TIterator >
      const_reference operator()(TIterator begin, TIterator end) const
      {
        return m_data[offset_prod(begin, begin+D, m_dims.rbegin())](begin+D,end);
          //(slice<std::size_t, d_type::value-D>(p.begin()));
          //[offset_prod(p.rbegin()+D, p.rend(), m_dims.rbegin()+D)];
      }
    private: // data
      TContainer m_data;
      coord_type m_dims;
    };


    //------------------------------------------------------------------------------

      //--------------------//
     //  multi_array_base  //
    //--------------------//

    /// Wraps a simple, linear container into a multidimensional array of dimension D.
    template < typename TContainer, std::size_t D >
    class multi_array_base<TContainer, D, false> : public multi_array_label
    {
    public: // typedefs
      typedef multi_array_base<TContainer,D,false>    Self;
      typedef typename TContainer::value_type         value_type;
      typedef typename TContainer::reference          reference;
      typedef typename TContainer::const_reference    const_reference;
      typedef meta::int_type<D> d_type;
      typedef numeric_array<std::size_t, D> coord_type;
    private: // helper
      static inline std::size_t size_helper(const coord_type & dims)
      {
        return std::accumulate(dims.begin(), dims.end(), std::size_t(1), std::multiplies<std::size_t>());
      }
    private: // typedefs
      typedef typename TContainer::iterator iterator;
      typedef typename TContainer::const_iterator const_iterator;
    public: // classes
      /*
      class range : public basic_range<typename TContainer::iterator>
      {
      public: 
        typedef basic_range<iterator> Base;
        class range(iterator begin, iterator end) : Base(begin, end) {}
      };
      class const_range : public basic_range<typename TContainer::const_iterator>
      {
      public:
        typedef basic_range<typename TContainer::const_iterator> Base;
        class const_range(iterator begin, iterator end) : Base(begin, end) {}
      };
      */
      typedef basic_range<iterator>                         range;
      typedef basic_range<const_iterator>                   const_range;
      typedef basic_volumetric_iterator<iterator,D>         volumetric_iterator;
      typedef basic_volumetric_iterator<const_iterator,D>   const_volumetric_iterator;
    public: // constructors
      multi_array_base() : m_data() {}
      explicit multi_array_base(const coord_type & dims)
        : m_data(Self::size_helper(dims))
        //prod(dims.begin(), dims.end()))
        , m_dims(dims)
      {
        assert(til::size(m_data) == std::accumulate(m_dims.begin(), m_dims.end(),std::size_t(1),std::multiplies<std::size_t>()));
      }
      multi_array_base(const coord_type & dims, value_type v)
        : m_data(size_helper(dims))
        //prod(dims.begin(), dims.end()), v)
        , m_dims(dims)
      {
        assert(size(m_data) == prod(m_dims.begin(), m_dims.end()));
      }
      multi_array_base(const TContainer & data) : m_data(data)
      {
        assert(size(m_data) == size(data));
      }
    public: // range
      range whole_range() { return range(m_data.begin(), m_data.end()); }
      const_range whole_range() const { return const_range(m_data.begin(), m_data.end()); }
    public: // set & get
      std::size_t size() const { return til::size(m_data); }
      coord_type dims() const { return m_dims; }
    public: // operators
      reference operator[](std::size_t i) { return m_data[i]; }
      const_reference operator[](std::size_t i) const { return m_data[i]; }
      reference operator()(const coord_type & p)
      {
        return this->operator()(p.rbegin(), p.rend());
      }
      const_reference operator()(const coord_type & p) const
      {
        return this->operator()(p.rbegin(), p.rend());
      }
      template < typename TIterator >
      reference operator()(TIterator begin, TIterator end)
      {
        return m_data[offset_prod(begin, end, m_dims.rbegin())];
      }
      template < typename TIterator >
      const_reference operator()(TIterator begin, TIterator end) const
      {
        return m_data[offset_prod(begin, end, m_dims.rbegin())];
      }
    private: // data
      TContainer m_data;
      coord_type m_dims;
    };

  } // namespace detail


  //------------------------------------------------------------------------------


    //---------------//
   //  multi_array  //
  //---------------//


  /// A multi-dimensional container.
  /// This class has two purposes. The first is to create a multidimensional
  /// array from a simple container. The second is to concatenate these containers
  /// into yet higher dimensional containers.
  /// Caveat: the container multi_array is based on must have random access, i.e.
  /// operator[]. This makes sense for most applications I guess.
  template < typename TContainer, std::size_t D >
  class multi_array
    : public detail::multi_array_base<TContainer,D,is_multi_array<typename value_type_of<TContainer>::type>::value>
  {
  public: // typedefs
    typedef detail::multi_array_base<TContainer, D, is_multi_array<typename value_type_of<TContainer>::type>::value > Base;
    typedef typename Base::coord_type   coord_type;
    typedef typename Base::value_type   value_type;
  public: // constructors
    multi_array() : Base() {}
    explicit multi_array(const coord_type & dims) : Base(dims) {}
    multi_array(const coord_type & dims, value_type v) : Base(dims, v) {}
    //multi_array(const TContainer & data) : Base(data) {}
  };

  //-----------------------------------------------------------------------------

  /*
  /// Concatenations of multi-arrays
  template < typename TContainer >
  class multi_conc_2<TContainer>
  {
  public: // typedef

    typedef TContainer outside_container;
    typedef typename value_type_of<TContainer>::type inside_container;
    typedef typename TContainer::reference reference;
    typedef typename TContainer::const_reference const_reference;
  public: // classes

    class range;
    class const_range;

  public: // constructors & destructor

    multi_conc() : m_data() {}
    multi_conc(const numeric_array<std::size_t, 2> & dims)
      : m_data(dims[0], inside_container(dims[1]))
      , m_dims(dims)
    {
      assert(size(m_data) == dims[0]);
    }
  
  public: // range

    range whole_range() { return range(til::whole_range(m_data), til::whole_range(*til::whole_range(m_data))); };

    //range begin() { return range(m_data.begin(), m_data.front().begin()); }
    //range end()  { return iterator(m_data.end(), m_data.back().end()); }

  public: // operators

    typename TContainer::reference operator[](std::size_t i) { return m_data[i]; }
    typename TContainer::const_reference operator[](std::size_t i) const { return m_data[i]; }

    reference operator()(std::size_t i, std::size_t j) { return m_data[i][j]; }
    const_reference operator()(std::size_t i, std::size_t j) const { return m_data[i][j]; }

  public: // functions

    boost::array<std::size_t, 2> dims() const
    { return m_dims; }

  private: // data
    TContainer m_data;
    numeric_array<std::size_t, 2> m_dims;
  };
  */

} // namespace til

#endif


