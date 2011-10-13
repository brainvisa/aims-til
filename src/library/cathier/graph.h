#ifndef TIL_GRAPH_H_
#define TIL_GRAPH_H_

namespace til
{
  /*
  template < typename T >
  struct OrientedGraphNode;
  
  template < typename T >
  struct OrientedGraphEdge
  {
    std::list<OrientedGraphNode<T> >::iterator from;
    std::list<OrientedGraphNode<T> >::iterator to;
  };
  */
  
  //---------------------------------------------------------------------------
  
  //NB: SIMPLE because edges are encoded in nodes. This is more efficient, on the other hand less general, since
  // edges cannot have properties of their own (e.g. such as a cost) as there have no label to name them.
  template < typename T >
  struct SimpleOrientedGraphNode
  {
    explicit SimpleOrientedGraphNode(const T & v) : value(v) {}
    std::list<typename std::list<SimpleOrientedGraphNode<T> >::iterator> from;
    std::list<typename std::list<SimpleOrientedGraphNode<T> >::iterator> to;
    T value;
  };
  
  //---------------------------------------------------------------------------
  
  template < typename TMeshVertexNode >
  struct MeshFaceNodeX
  {
    boost::array<typename std::list<TMeshVertexNode>::iterator,3> face;
  };

  //---------------------------------------------------------------------------
  
  /*
  template < typename T >
  struct DefaultPointer
  {
    typedef typename std::list<T>::iterator type;
  };
  template < typename TAttribute, template < typename > class TPointerPolicy = DefaultPointer >
  */

  //---------------------------------------------------------------------------
  
  // TODO: What we would like here is to be able to derive from MeshVertexNode. The problem comes from the fact
  // that right now, MeshFaceNode must have the knowledge of MeshVertexNode and vice versa. It means that
  // these classes cannot be templated over the other class. So the stuff to think about is how to get rid
  // of this dependancy.
  template < typename TAttribute >
  class MeshVertexNodeX
  {
  public: // typedefs
  
    typedef MeshVertexNodeX<TAttribute>                                     Self;
    typedef til::numeric_array<float,3>                                     Vertex;
    typedef std::list<Self>                                                 VertexCollection;
    typedef typename VertexCollection::iterator                             VertexIndex;
    typedef std::list<VertexIndex>                                          VertexIndexCollection;
  
    typedef MeshFaceNodeX<Self>                   Face;
    typedef std::list<Face>                       FaceCollection;
    typedef typename FaceCollection::iterator     FaceIndex;
    typedef std::list<FaceIndex>                  FaceIndexCollection;
  
    //typedef std::list<typename TPointerPolicy<Self>::type>                  NeighborCollection;
    //typedef std::list<typename TPointerPolicy<MeshFaceNodeX<Self> >::type>  FaceCollection;
  
  public: // set & get
   
    Vertex       & pos()       { return m_pos; }
    Vertex const & pos() const { return m_pos; }
  
    VertexIndexCollection       & neighbors()       { return m_neighbors; }
    VertexIndexCollection const & neighbors() const { return m_neighbors; }
  
    FaceIndexCollection       & faces()       { return m_faces; }
    FaceIndexCollection const & faces() const { return m_faces; }
  
    TAttribute       & attribute()       { return m_attribute; }
    TAttribute const & attribute() const { return m_attribute; }
  
  private: // data
  
    Vertex                  m_pos;
    VertexIndexCollection   m_neighbors;
    FaceIndexCollection     m_faces;
    TAttribute              m_attribute;
  };
  
  //---------------------------------------------------------------------------
  
  template < typename TAttribute, typename TEdgeAttribute >
  class DirectedNode
  {
  public: // typedefs
    typedef DirectedNode<TAttribute, TEdgeAttribute> Self;
  private: // data
    TAttribute m_attribute;
    std::list<std::pair<Self, TEdgeAttribute> > m_neighbors;
  };
  
  template < typename TAttribute >
  class UndirectedNode
  {
  private: // data
    TAttribute m_attribute;
    std::list<Edge> m_edges;
  };

  //---------------------------------------------------------------------------

} // namespace til


// Package include
#include "graph_accessors.h"

#endif /*GRAPH_H_*/
