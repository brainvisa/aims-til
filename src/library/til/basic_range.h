#ifndef TIL_BASIC_RANGE_H
#define TIL_BASIC_RANGE_H

// includes from STL
#include <vector>

// includes from TIL
#include "til/labels.h"
#include "til/meta.h"
#include "til/miscTools.h"
#include "til/traits.h"

namespace til
{

  /// A simple range made out of a couple of iterators.
  template < typename TIterator >
  class basic_range
    : public range_label
  {
  public: // typedefs
    typedef TIterator iterator;
    typedef typename TIterator::value_type  value_type;
    typedef typename TIterator::reference   reference;
    typedef typename TIterator::pointer     pointer;
  public: // constructors
    basic_range(TIterator begin, TIterator end) : m_i(begin), m_end(end) {}
  public: // operators
    // NB: there is no operator--. operator++ should be really interpreted as "next".
    // If you want to scan backwards, use reverse iterators.
    void operator++() { ++m_i; }
    pointer operator->() { return m_i->operator->(); }
    reference operator*() { return *m_i; }
    bool ok() const { return m_i != m_end; }

    /*
    void operator=(const basic_range<TIterator> & other)
    {
      m_i = other.m_i;
      m_end = other.m_end;
    }
    */

  private: // data
    TIterator m_i;
    TIterator m_end;
  };


  namespace detail
  {
    template < typename TIterator, std::size_t D >
    class inner_basic_volumetric_iterator
      : public TIterator
    {
    public: // typedefs
      typedef TIterator Base;
      typedef meta::int_type<D> dim;
    protected: // typedefs
      // TODO: en fait, ca n'a aucun sens d'avoir des typedef private, autant tous les mettre en protected, non?
      typedef numeric_array<std::size_t, D> Vector;
    public: // constructors
      inner_basic_volumetric_iterator(TIterator i, const Vector & size)
        : Base(i)
        , m_size(size)
      {}
    public: // set & get
      Vector & size() { return m_size; }
    public: // functions
      template < typename TPosIterator >
      void from_pos(TIterator begin, TPosIterator pbegin, TPosIterator pend)
      {
        *this = begin;
        std::advance(*this, pos2offset(pbegin, pend, m_size.begin()));
      }
    private: // data
      Vector m_size;
    };
  }

  /*
  template < typename TIterator, std::size_t D, typename TPosIterator >
  void advance(detail::inner_basic_volumetric_iterator<TIterator,D> & volI, TPosIterator pbegin, TPosIterator pend)
  {
    std::advance(volI, pos2offset(pbegin, pend, volI.size().begin()));
  }
  */

  template < typename TIterator, std::size_t D >
  class basic_volumetric_iterator
    : public detail::inner_basic_volumetric_iterator<TIterator,D>
  {
  public: // typedefs
    typedef detail::inner_basic_volumetric_iterator<TIterator,D> Base;
  protected: // typedefs
    typedef typename Base::Vector Vector;
  public: // constructors
    basic_volumetric_iterator(TIterator begin, const Vector & size)
      : Base(begin, size)
      , m_begin(begin)
    {}
  public: // functions

    void from_pos(const Vector & pos)
    { this->from_pos(pos.begin(), pos.end()); }

    template < typename TPosIterator >
    void from_pos(TPosIterator pbegin, TPosIterator pend)
    { this->from_pos(m_begin, pbegin, pend); }

  private: // functions
    /*
    void initOffset(const Vector & size)
    {
      m_offset[0] = 1;
      for (std::size_t i = 1; i < D; ++i)
      {
        m_offset[i] = m_offset[i-1] * size[i-1];
      }
    }
    */
  private: // data
    TIterator m_begin;
  };

/*
  template < typename TIterator >
  class volumetric_range : public basic_range<TIterator>
  {
  public: // typedefs
    typedef basic_range<TIterator> Base;
    typedef typename Base::reference reference;
  public: // operators
    reference getValue(
  };
*/

  /*
  template < typename TContainer >
  typename range_of<TContainer>::type
  whole_range(TContainer & c)
  { return c.whole_range(); }

  template < typename TContainer >
  typename const_range_of<TContainer>::type
  whole_range(const TContainer & c)
  { return c.whole_range(); }
  */

} // namespace til

#endif

