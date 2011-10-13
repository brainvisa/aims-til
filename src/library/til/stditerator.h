#ifndef TIL_STDITERATOR_H
#define TIL_STDITERATOR_H

// includes from BOOST
#include <boost/type_traits.hpp>

namespace til{
namespace detail
{
  /// Technical implementation detail for stditerator
  template < typename T, bool b >
  struct stdit { typedef typename T::iterator type; };
  
  template < typename T >
  struct stdit<T, true> { typedef typename T::const_iterator type; };
}

/// Select T2::const_iterator T is const, T2::iterator otherwise.
template < typename T, typename T2 = T >
struct stditerator : public detail::stdit<T, boost::is_const<T2>::value > {};

} // namespace til

#endif

