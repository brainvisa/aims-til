#ifndef TIL_BASIC_ITERATOR_H
#define TIL_BASIC_ITERATOR_H

// includes from STL
#include <iterator>

// includes from BOOST
#include "boost/type_traits.hpp"

// includes from TIL
#include "til/til_common.h"

namespace til
{
  /// A simple wrapper around T*.
  /// This class is just there to encapsulate a T* iterator in a class... Nothing to
  /// be proud of :/
  /// Not fully compliant because other features not needed yet.
  // TODO: it seems that std::iterator does not define const_reference...
  // is it standard or just MSVC?
  template < typename T >
  class basic_iterator : public std::iterator<std::random_access_iterator_tag, T>
  {
  public://
    typedef std::iterator<std::random_access_iterator_tag, T> Base;
    //TODO: why is this crap necessary?? Should peep into some of the stl algo
    // and look, if indeed value_type of const_iterators are const, how they
    // deal with that.
    typedef typename boost::remove_const<typename Base::value_type>::type value_type;
  public: // constructors & destructor
    basic_iterator() : Base(), m_pointer() {}
    basic_iterator(T * p) : Base(), m_pointer(p) {}
    template < typename X >
    basic_iterator(basic_iterator<X> i) : Base(), m_pointer((T*)i) {}
  public: // operators
    void operator++() { ++m_pointer; }
    // NB: it probably doesn't make sense to have a const operator*. That's probably
    // why iterators don't define const_reference...
    //const T & operator*() const { return *m_pointer; }
    T & operator*() { return *m_pointer; }
    operator T* () { return m_pointer; }
    T * operator->() { return m_pointer; }
  public: // friends
    template < typename T1, typename T2 >
    friend bool operator==(const basic_iterator<T1> & i1, const basic_iterator<T2> & i2);
    template < typename T1, typename T2 >
    friend bool operator!=(const basic_iterator<T1> & i1, const basic_iterator<T2> & i2);
  private: // data
    T * m_pointer;
  };

  // NB: actually, the fact that we have two template parameters T1 and T2 won't
  // enable us to compare pointers on different object. This is just a mean to compare
  // const and non-const pointers.
  template < typename T1, typename T2 >
  inline
  bool operator==(const basic_iterator<T1> & i1, const basic_iterator<T2> & i2)
  {
    //return static_cast<T1*>(i1) == static_cast<T2*>(i2);
    return i1.m_pointer == i2.m_pointer;
  }

  template < typename T1, typename T2 >
  inline
  bool operator!=(const basic_iterator<T1> & i1, const basic_iterator<T2> & i2)
  {
    //return static_cast<T1*>(i1) != static_cast<T2*>(i2);
    return i1.m_pointer != i2.m_pointer;
  }

} // namespace til

#endif

