#ifndef TIL_CAT2TYPE_H_
#define TIL_CAT2TYPE_H_

// includes from BOOST
#include <boost/utility/enable_if.hpp>

// includes from TIL
#include "til/til_common.h"

///\file Label for categories. Kind of the same thing as type2type, except that here
/// we want a type to label a whole 'category', not a single type. A category can be 
/// any group of types. Generally, a category is defines via a is_a<...> template.
// Right now, only binary tests -- I dunno if it's wise to go further... :/
// Actually, maybe the only safe way to generalize to more tests is to use a combination
// of binary tests. And that can be used straight away, by using several labels as
// function parameters.

namespace til
{
  namespace label
  {
    /// A label corresponding to a positive answer of a specific is_a test
    template < template <typename> class TIS_A >
    struct Passed
    {};

    /// A label corresponding to a negative answer of a specific is_a test
    template < template <typename> class TIS_A >
    struct Failed
    {};
  }
  
  namespace detail
  {
    template < template <typename> class TIS_A, typename T, bool B >
    struct test_whether_impl
    { typedef label::Failed<TIS_A> type; };

    template < template <typename> class TIS_A, typename T >
    struct test_whether_impl<TIS_A, T, true>
    { typedef label::Passed<TIS_A> type; };
  }    
  
  /// A structure used for is_a tests.
  /// Result of the test is coded in the 'type' typedef.
  template < template <typename> class TIS_A, typename T >
  struct test_whether
    : public detail::test_whether_impl<TIS_A, T, TIS_A<T>::value>
  {};
  
  /// Do an is_a test.
  template < template <typename> class TIS_A, typename T >
  inline 
  typename test_whether<TIS_A,T>::type
  cat2type()
  {
    typedef typename test_whether<TIS_A,T>::type TReturn;
    return TReturn();
  }
  
  /*
  template < template <typename> class TIS_A, typename T >
  inline
  typename boost::enable_if<TIS_A<T>, label::Passed<TIS_A> >::type
  cat2type()
  {
    return label::Passed<TIS_A>();
  }

  template < template <typename> class TIS_A, typename T >
  inline
  typename boost::disable_if<TIS_A<T>, label::Failed<TIS_A> >::type
  cat2type()
  {
    return label::Failed<TIS_A>();
  }
  */
  
} // namespace til


#endif /*CAT2TYPE_H_*/
