#ifndef TIL_VALUE_PROXY_POLICY_H
#define TIL_VALUE_PROXY_POLICY_H

// includes from TIL
#include "til/traits.h"

namespace til { namespace policy
{

  /// Policy for value_proxy that defines how it should access the container elements.
  template < typename TContainer, typename TIndex >
  class VPAccess_Default
  {
  public: // typedefs

    typedef typename TContainer::value_type           value_type;
    //typedef typename reference_of<TContainer>::type     reference;
    typedef typename TContainer::const_reference      const_reference;
    typedef typename TContainer::const_pointer        const_pointer;

  public: // functions

    /// Get an element value
    //static value_type get(TContainer * pContainer, TIndex i)
    //{ return pContainer->get(i); }
    static const_reference get(TContainer & container, TIndex i)
    { return container.get(i); }
          
    static const_pointer getPointer(TContainer & container, TIndex i)
    { return container.getPointer(i); }

    /// Set a value to an element
    //static void set(TContainer * pContainer, TIndex i, value_type value)
    //{ return pContainer->set(i, value); }
    static void set(TContainer & container, TIndex i, value_type value)
    { container.set(i, value); }
  };

}} // namespace til::policy

#endif

