#ifndef TIL_VALUE_PROXY_H_
#define TIL_VALUE_PROXY_H_

// includes from STL
#include <cstdlib>

// includes from BOOST
#include "boost/call_traits.hpp"
#include "boost/type_traits.hpp"

// includes from TIL
#include "til/traits.h"
#include "til/value_proxy_policies.h"

namespace til
{
  
  //------------------------------------------------------------------------------------------------

    //---------------------//
   //  const_value_proxy  //
  //---------------------//

  /// The value_proxy is designed to provide a simple access to container values
  /// that are in fact complicated (e.g. when data is compressed).
  /// The value proxy do not rely on standard access interface (e.g. operator[]) because
  /// they are designed precisely to help building them.
  /// Therefore they necessarily rely on a non-standard interface, namely get(index) and
  /// set(index, value) functions. This is just the default interface though -- this can 
  /// be changed to suit your needs through the AccessPolicy template parameter.
  
  template < typename TContainer, typename TIndex = std::size_t, typename AccessPolicy = policy::VPAccess_Default<TContainer,TIndex> > //, typename T = typename value_type<TContainer>::type >
  class const_value_proxy
  {
  public: // typedefs

    typedef const_value_proxy<TContainer, TIndex, AccessPolicy>   Self;
    //typedef TContainer                                 container;
    typedef typename value_type_of<TContainer>::type    value_type;
    typedef TIndex                                      index_type;
    typedef typename AccessPolicy::const_reference      const_reference;
    typedef typename AccessPolicy::const_pointer        const_pointer;

  public: // constructors & destructor
    
    //const_value_proxy(TIndex i, TContainer * pContainer) : m_i(i), m_pContainer(pContainer) {}
    const_value_proxy(TIndex i, TContainer & container) : m_i(i), m_pContainer(&container) {}
    template < typename XContainer, typename XIndex, typename XAccessPolicy>
    const_value_proxy(const const_value_proxy<XContainer, XIndex, XAccessPolicy> & other)
      : m_i(other.index())
      , m_pContainer(&other.container())
    {}

  public: // operators

    const_pointer operator->() const
    {
      return AccessPolicy::getPointer(*m_pContainer, m_i);
    }
  
    //operator const value_type () const
    operator const_reference ()
    {
      //return m_pContainer->get(m_i);
      return AccessPolicy::get(*m_pContainer, m_i);
    }
    
  public: // set & get
    
    // TIndex might be a reference. Therefore, manipulation of
    // TIndex has to go through such objects.
    typename boost::call_traits<TIndex>::const_reference index() const { return m_i; }
    typename boost::call_traits<TIndex>::reference       index()       { return m_i; }  
    //TContainer *  pContainer() { return m_pContainer; }
    TContainer & container() { return *m_pContainer; }
    const TContainer & container() const { return *m_pContainer; }

    /*
    void operator=(const Self & other)
    {
      m_i = other.m_i;
      m_container = other.m_container;
    }
    */

  private: // data

    TIndex        m_i;
    TContainer *  m_pContainer;
    //TContainer &  m_container;
  };
  

  //------------------------------------------------------------------------------------------------

    //---------------//
   //  value_proxy  //
  //---------------//

  /// Adds some left-value operations to const_value_proxy.
  template < typename TContainer, typename TIndex = std::size_t, typename AccessPolicy = policy::VPAccess_Default<TContainer,TIndex> > //, typename T = typename value_type<TContainer>::type >
  class value_proxy : public const_value_proxy<TContainer, TIndex, AccessPolicy>
  {
  public: // typedefs
  
    typedef value_proxy<TContainer, TIndex, AccessPolicy>         Self;
    typedef const_value_proxy<TContainer, TIndex, AccessPolicy>   Base;
    typedef typename Base::value_type                             value_type;

  public: // constructors & destructor
    
    //value_proxy(TIndex i, TContainer * pContainer) : Base(i, pContainer) {}
    value_proxy(TIndex i, TContainer & container) : Base(i, container) {}
    template < typename XContainer, typename XIndex, typename XAccessPolicy >
    value_proxy(const value_proxy<XContainer,XIndex,XAccessPolicy> & other)
      : Base(other)
    {}

  public: // operators
  
    // NB: actually, I guess all these operations could be done using a value_type& conversion
    // operator. The problem doing so is that when getting a reference on a value that was not
    // initialized before, the container has to initialize it to the default value. Which might
    // be bad, because we might be breaking the rule that no value should be initialized to
    // the default value.
  
    void operator=(typename boost::call_traits<value_type>::param_type value)
    {
      // NB: we have to use set, and not getref, because we may be inserting a new value.
      //m_pContainer->set(m_i, value);
      //AccessPolicy::set(this->pContainer(), this->index(), value);
      AccessPolicy::set(this->container(), this->index(), value);
    }
    
    void operator+=(typename boost::call_traits<value_type>::param_type value)
    {
      this->assignOp(value, std::plus<value_type>());
    }

    void operator-=(typename boost::call_traits<value_type>::param_type value)
    {
      this->assignOp(value, std::minus<value_type>());
    }

    void operator*=(typename boost::call_traits<value_type>::param_type value)
    {
      this->assignOp(value, std::multiplies<value_type>());
    }

    void operator/=(typename boost::call_traits<value_type>::param_type value)
    {
      this->assignOp(value, std::divides<value_type>());
    }
    
    /*
    operator value_type & () 
    {
      //return m_pContainer->getref(m_i);
      return AccessPolicy::getref(this->pContainer(), this->index());
    }
    */
  

    /*
    void operator=(const Self & other)
    {
      this->Base::operator=(other);
    }
    */

  private: // functions
  
    template < typename TFunctor >
    void assignOp(typename boost::call_traits<value_type>::param_type value, TFunctor f)
    {
      //AccessPolicy::set(this->pContainer(), this->index(),
      //  f(AccessPolicy::get(this->pContainer(), this->index()), value));
      AccessPolicy::set(this->container(), this->index(),
        f(AccessPolicy::get(this->container(), this->index()), value));
    }
  };
  
} // namespace til


#endif /*VALUE_PROXY_H_*/
