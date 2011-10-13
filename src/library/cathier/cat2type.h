#ifndef CAT2TYPE_H_
#define CAT2TYPE_H_

// includes from BOOST
#include <boost/utility/enable_if.hpp>

/// \file Label for categories. Kind of the same thing as type2type, except that here
/// we want a type to label a whole 'category', not a single type. A category can be 
/// any group of types. Generally, a category is defines via a is_a<...> template.
// Right now, only binary tests -- I dunno if it's wise to go further... :/
// Actually, maybe the only safe way to generalize to more tests is to use a combination
// of binary tests. And that can be used straight away, by using several labels as
// function parameters.

namespace til
{
  //---------------------------------------------------------------------------

  /// Namespace for labels, i.e. empty structures that are use only for their meaning.
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
  
  //---------------------------------------------------------------------------

  template < template <typename> class TIS_A, typename T>
  inline
  typename boost::enable_if<TIS_A<T>, label::Passed<TIS_A> >::type
  cat2type()
  {
    return label::Passed<TIS_A>();
  }

  //---------------------------------------------------------------------------

  template < template <typename> class TIS_A, typename T>
  inline
  typename boost::disable_if<TIS_A<T>, label::Failed<TIS_A> >::type
  cat2type()
  {
    return label::Failed<TIS_A>();
  }

  //---------------------------------------------------------------------------

} // namespace


#endif /*CAT2TYPE_H_*/
