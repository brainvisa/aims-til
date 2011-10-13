#ifndef TIL_ACCESSORS_H_
#define TIL_ACCESSORS_H_

/// \file Definition of basic accessors.
/// Accessors are a different concept than iterator, because their aim is to be able to relate
/// different information; that is, the very same index yields different information when passed
/// to different accessors.
/// For example, for meshes, the index could be a vertex index. When given this index, one accessor could
/// return the vertex position, another one could return the neighbors of this vertex, etc. So while the index
/// clearly represent a given vertex, there is no limit to the information that can be attached to it, nor do
/// we say how this information is stored in memory.


namespace til { 
/// namespace for accessors.
namespace xsr
{
  
  //-------------------------------------------------------------------------------------
  
    //------------------//
   //  Integer_access  //
  //------------------//

  /// Integer access over a random-access container.
  // TODO: add requirement that TCollection is random-access.
  template < typename TCollection, typename TIndex = std::size_t >
  class Integer_access
  {
  public: // typedefs
    typedef TIndex                                    index_type;
    typedef typename TCollection::value_type          value_type;
    typedef typename TCollection::reference           reference;
    typedef typename TCollection::const_reference     const_reference;

  public: // constructors
  
    // NB: Actually, it can be quite comfortable to have it non explicit. We can pass collections rather
    // than constructing explicitely an accessor.
    Integer_access(TCollection & data) : m_data(data) {}
    
  public: // operators
    reference operator()(index_type i) { return m_data[i]; }
    const_reference operator()(index_type i) const { return m_data[i]; }
  private: // data
    TCollection & m_data;
  };
  
  //-------------------------------------------------------------------------------------
    
  /*
  template < typename TReference, typename TValueType = naked_type<TReference>::type >
  struct GetVertexField
  {
    typedef TValueType    value_type;
    typedef TReference    reference;
    template < typename TIterator >
    static TResult get(TIterator i) { return i->vertex(); }
  };
  */
 
 
  //-------------------------------------------------------------------------------------

    //-------------------//
   //  Iterator_access  //
  //-------------------//
 
  /// Return value pointed at by an iterator.
  template < typename TIterator > //, typename Getter >
  struct Iterator_access
  {
  public: // typedefs
    typedef typename TIterator::value_type  value_type;
    typedef TIterator                       index_type;
  public: // static accessors
    value_type & operator()(index_type i) const
    { return *i; }
  };

  
}} // namespace til::xsr

#endif /*ACCESSORS_H_*/
