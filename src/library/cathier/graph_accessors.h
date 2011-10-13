#ifndef TIL_GRAPH_ACCESSORS_H_
#define TIL_GRAPH_ACCESSORS_H_

/// \File Belongs to graph package.
/// Do not include directly, include "graph.h" instead.

namespace til { namespace xsr {
  
  
  /// namespace for graph accessors
  namespace graph
  {
    template < typename TIterator >
    struct Position
    {
    public: // typedefs
      typedef typename TIterator::value_type::Vertex        value_type;
      typedef typename TIterator::value_type::Vertex &      reference;
      typedef typename TIterator::value_type::VertexIndex   index_type;
    public: // static accessors
      // NB: it is not stupid to fix the reference at value_type, because the index_type, hence the iterator,
      // is also fixed inside the meshvertexnode class. 
      reference
      operator()(index_type i) const
      { return i->pos(); }
    };
    
    template < typename TIterator >
    struct Neighbors
    {
    public: // typedefs
      typedef typename TIterator::value_type::VertexIndexCollection     value_type;
      typedef typename TIterator::value_type::VertexIndexCollection &   reference;
      typedef typename TIterator::value_type::VertexIndex               index_type;
    public: // static accessors
      reference
      operator()(index_type i) const
      { return i->neighbors(); }
    };
  }

  
  
}} // namespace til::xsr


#endif /*GRAPH_ACCESSORS_H_*/
