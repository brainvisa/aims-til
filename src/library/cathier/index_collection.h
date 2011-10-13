#ifndef INDEX_COLLECTION_H_
#define INDEX_COLLECTION_H_

// includes from STL
#include <cassert>

// includes from BOOST
#include <boost/call_traits.hpp>
#include <boost/type_traits.hpp>

// includes from TIL
#include "til/templateTools.h"

// includes from TIL
#include "declarations_external.h"
#include "globalTraits.h"

namespace til
{

  //////////////KEEP THIS/////////////
  // I highly suspect it fails because of stupid GCC 3.3
  
  // NB: we template on TPointer since it can be an iterator
  template < typename TIndex, typename T, typename TAlloc, typename TPointer >
  inline
  typename boost::enable_if_c<
    boost::is_same<TIndex, std::size_t>::value &&
    !boost::is_same<TPointer, std::size_t>::value
    ,
    TIndex
  >::type
  getIndex(const std::vector<T, TAlloc> & c, TPointer pElem)
  {
    return std::distance(&*(c.begin()), &*pElem);
  }
    
  template < typename TIndex, typename TCollection >
  inline
  typename boost::enable_if<
    boost::is_same<TIndex, std::size_t>
    ,
    TIndex
  >::type
  getIndex(const TCollection &, std::size_t i)
  {
    return i;
  }

  template < typename TIndex, typename T, typename TAlloc, typename TPointer >
  inline
  typename boost::enable_if_c<
    boost::is_same<TIndex, T*>::value &&
    !boost::is_same<TPointer, std::size_t>::value
    ,
    TIndex
  >::type
  getIndex(const std::vector<T, TAlloc> & c, TPointer pElem)
  {
    return &*pElem;
  }
  
  template < typename TIndex, typename TContainer >
  inline
  typename boost::enable_if<
    boost::is_same<TIndex, typename value_type_of<TContainer>::type*>
    ,
    TIndex
  >::type
  getIndex(const TContainer & c, std::size_t i)
  {
    return &(c[i]);
  }
  

  /*
  namespace detail
  {
    template < typename TContainer, typename TPointer, typename TReturn >
    TReturn
    indexee_pointer(TContainer &, TPointer p)
    {
      return *p;
    }

    template < typename TContainer, typename TScalarIndex, typename TReturn >
    TReturn
    indexee_scalar(const TContainer & c, TScalarIndex i)
    {
      return c[i];
    }      
  }
  */
  
  /*
  template < typename TContainer >
  inline
  typename type_if< boost::is_const<TContainer>::value,
    typename boost::call_traits<typename value_type<TContainer>::type>::const_reference,
    typename boost::call_traits<typename value_type<TContainer>::type>::reference
  >::type
  indexee(TContainer & c, typename value_type<TContainer>::type * p)
  {
    return *p;
  }
  

  template < typename TContainer >
  inline
  typename type_if< boost::is_const<TContainer>::value,
    typename boost::call_traits<typename value_type<TContainer>::type>::const_reference,
    typename boost::call_traits<typename value_type<TContainer>::type>::reference
  >::type
  indexee(TContainer & c, std::size_t i)
  {
    return c[i];
  }
  */
  
  namespace detail
  {
    
    template < typename TContainer, typename TPointer >
    class index_collection_pointer
    {
    public: // typedefs

      typedef TPointer                         index_type;
      typedef typename TContainer::value_type  indexed_type;

    public: // constructors & destructor

      index_collection_pointer() {}
      index_collection_pointer(TContainer &) {}

    public: // functions

      typename boost::call_traits<indexed_type>::const_reference
      get(typename boost::call_traits<TPointer>::param_type p) const
      {
        return *p;
      }

      typename boost::call_traits<indexed_type>::reference
      get(typename boost::call_traits<TPointer>::param_type p)
      {
        return *p;
      }
    };
    
    template < typename TContainer, typename TScalarIndex >
    class index_collection_scalar
    {
    public: // typedefs

      typedef TScalarIndex                              index_type;
      typedef typename TContainer::value_type           indexed_type;
      typedef typename TContainer::reference            reference;
      typedef typename TContainer::const_reference      const_reference;

    public: // constructors & destructor

      index_collection_scalar() : m_pContainer() {}
      index_collection_scalar(TContainer & c) : m_pContainer(&c) {}

    public: // set & get

      void setContainer(TContainer & c) { m_pContainer = &c; }
      TContainer & getContainer() { assert(m_pContainer != 0); return *m_pContainer; }
      const TContainer & getContainer() const { assert(m_pContainer != 0); return *m_pContainer; }

    public: // functions

      const_reference get(TScalarIndex i) const
      {
        assert(m_pContainer != 0);
        return (*m_pContainer)[i];
      }

      reference get(TScalarIndex i)
      {
        assert(m_pContainer != 0);
        return (*m_pContainer)[i];
      }
      
    private: // data
    
      // TODO: why using a shared_ptr? It's a pain because it requires the user to use it as well.
      // OK, but we need to use a pointer, because we need to be able to default-initialize the class
      TContainer * m_pContainer;
      //shared_ptr<TContainer> m_pContainer;
    };
  } // namespace detail
  
  
  template < typename TContainer, typename TIndex >
  class index_collection;
  
  
  template < typename TContainer, typename T >
  class index_collection< TContainer, T* > 
    : public detail::index_collection_pointer<TContainer, T*>
  {
  public: // typedefs

    typedef detail::index_collection_pointer<TContainer, T*> Base;
    
  public: // constructors & destructor
  
    index_collection() : Base() {}
    index_collection(TContainer & c) : Base(c) {}
  };
  

  template < typename TContainer >
  class index_collection< TContainer, std::size_t >
    : public detail::index_collection_scalar<TContainer, std::size_t >
  {
  public: // typedefs
  
    typedef detail::index_collection_scalar<TContainer, std::size_t > Base;

  public: // constructors & destructor
  
    index_collection() : Base() {}
    index_collection(TContainer & c) : Base(c) {}
  };



  // NB: TIndexCollection could actually be a reference
  template < typename TIndexCollection >
  class virtual_collection
  {
  public: // typedefs
    typedef typename TIndexCollection::reference reference;
    
  public: //
    // NB: we have to use call_traits as TIndexCollection might be a reference
    virtual_collection(typename boost::call_traits<TIndexCollection>::param indices)
      : m_indices(indices)
    {}
  public: //
    reference operator[](std::size_t i) { return this->get(i); }
  private: // data
    TIndexCollection m_indices;
  };


} // namespace til


#endif /*INDEX_COLLECTION_H_*/
