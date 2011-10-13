#ifndef TIL_ORDERED_ITERATOR_H_
#define TIL_ORDERED_ITERATOR_H_

// includes from STL
#include <iterator>

// includes from TIL
#include "til/traits.h"

namespace til
{
  
  /// Iterate through a random access container according to some predefined order.
  template < typename TIterator, typename TOrderIterator >
  class Ordered_iterator
    : public std::iterator<
      typename std::iterator_traits<TOrderIterator>::iterator_category,
      typename std::iterator_traits<TIterator>::value_type,
      typename std::iterator_traits<TIterator>::difference_type,
      typename std::iterator_traits<TIterator>::pointer,
      typename std::iterator_traits<TIterator>::reference>
  {
  public: // typedefs
    typedef Ordered_iterator<TIterator, TOrderIterator> Self;
    typedef std::iterator<
      typename std::iterator_traits<TOrderIterator>::iterator_category,
      typename std::iterator_traits<TIterator>::value_type,
      typename std::iterator_traits<TIterator>::difference_type,
      typename std::iterator_traits<TIterator>::pointer,
      typename std::iterator_traits<TIterator>::reference> Base;
    typedef typename Base::reference reference;
    
  public: // constructors
  
    /// Construct an Ordered_iterator from a random-access iterator on the beginning of the data and
    /// an iterator on the beginning of the indices.
    Ordered_iterator(TIterator databegin, TOrderIterator orderbegin)
      : Base()
      , m_databegin(databegin)
      , m_orderit(orderbegin)
    {}
    
  public: // operators
  
    void operator++() { ++m_orderit; }
    reference operator*() { return *this->current_iter(); }
    reference operator->() { return *this->current_iter(); }
    bool operator==(const Ordered_iterator & it) const { return this->current_iter() == it.current_iter(); }
    bool operator!=(const Ordered_iterator & it) const { return !(*this == it); }

  private: // functions
  
    // TODO: adapt this for maps...
    TIterator current_iter() { return m_databegin + *m_orderit; }

  private: // checks
  
    // check that TIterator is a random access iterator
    typedef enable_if<is_same<typename std::iterator_traits<TIterator>::iterator_category, std::random_access_iterator_tag> >
      Check_is_random_access_iterator;
  
  private: // data, input
  
    TIterator m_databegin;
    TOrderIterator m_orderit;
  };
  
  template < typename TCollection, typename TOrderIterator >
  Ordered_iterator<typename TCollection::iterator, TOrderIterator>
  ordered_iterator(TOrderIterator obegin, TCollection & coll)
  { return Ordered_iterator<typename TCollection::iterator, TOrderIterator>(coll.begin(), obegin); }
  
  template < typename TCollection, typename TOrderIterator >
  Ordered_iterator<typename TCollection::const_iterator, TOrderIterator>
  ordered_iterator(TOrderIterator obegin, const TCollection & coll)
  { return Ordered_iterator<typename TCollection::const_iterator, TOrderIterator>(coll.begin(), obegin); }
  
} // namespace til


#endif /*PERMUATION_ITERATOR_H_*/
