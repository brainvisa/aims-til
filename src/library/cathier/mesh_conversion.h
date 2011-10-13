#ifndef TIL_MESH_CONVERSION_H_
#define TIL_MESH_CONVERSION_H_

namespace til
{
  
  //---------------------------------------------------------------------------
  
    //----------------------------//
   //  Graph2ListMeshConvertor2  //
  //----------------------------//
  
  template < typename TVertexCollection, typename TFaceCollection, typename TGraphVertices, typename TGraphFaces >
  class Graph2ListMeshConvertor2
  {
  public: // typedefs
  
    //typedef std::map<typename std::list<TVertexNode>::const_iterator, std::size_t, ItComp<typename std::list<TVertexNode>::const_iterator> > VertexIndexMap;
    typedef std::map<const void*, std::size_t> VertexIndexMap;
  
  public: // set & get
  
    VertexIndexMap & index_vertex() { return m_index_vertex; }
    //shared_ptr<TVertexCollection> vertices() { return m_vertices; }
    //shared_ptr<TFaceCollection> faces() { return m_faces; }
    
  public: // operators
  
    //template < typename TGraphVertices, typename TGraphFaces >
    void operator()
    (
      const TGraphVertices & graph_vertices,
      const TGraphFaces & graph_faces,
      TVertexCollection & vertices,
      TFaceCollection & faces
    );
    
  private: // data, output
  
    //shared_ptr<TVertexCollection>   m_vertices;
    //shared_ptr<TFaceCollection>     m_faces;
    VertexIndexMap                  m_index_vertex;
  };
  

  //---------------------------------------------------------------------------

    //---------------------------//
   //  List2GraphMeshConvertor  //
  //---------------------------//

  /// Converts a (vertexCollection, faceCollection) mesh into a
  /// (vertexList, faceList) mesh graph.
  template < typename TVertexNode, typename TFaceNode, typename TVertexCollection, typename TFaceCollection, typename TNeighborCollection >
  class List2GraphMeshConvertor
  {
  public: // operators
  
    //template < typename TVertexCollection, typename TFaceCollection, typename TNeighborCollection >
    void operator()
    (
      TVertexCollection const &                     vertices,             ///< [input] mesh vertices
      TFaceCollection const &                       faceIndices,          ///< [input] mesh faces
      std::vector<std::list<std::size_t> > const &  invertedFaceIndices,  ///< [input] inverted face indices
      TNeighborCollection const &                   neighbors,            ///< [input] vertex neighborhoods
      std::list<TVertexNode> &                      graph_vertices,
      std::list<TFaceNode> &                        graph_faces
    );
    
  public: // set & get
  
    /*
    std::list<TVertexNode>       & graph_vertices()       { return m_graph_vertices; }
    std::list<TVertexNode> const & graph_vertices() const { return m_graph_vertices; }
  
    std::list<TFaceNode>       & graph_faces()       { return m_graph_faces; }
    std::list<TFaceNode> const & graph_faces() const { return m_graph_faces; }
    */
    /*
    std::vector<typename std::list<TVertexNode>::iterator>       & index_vertex()       { return m_index_vertex; }
    std::vector<typename std::list<TVertexNode>::iterator> const & index_vertex() const { return m_index_vertex; }
  
    std::vector<typename std::list<TFaceNode>::iterator>       & index_face()       { return index_face; }
    std::vector<typename std::list<TFaceNode>::iterator> const & index_face() const { return index_face; }
    */
    
  private: // data, internal
    //std::list<TVertexNode>    m_graph_vertices;
    //std::list<TFaceNode>      m_graph_faces;
    std::vector<typename std::list<TVertexNode>::iterator>  m_index_vertex;
    std::vector<typename std::list<TFaceNode>::iterator>    m_index_face;
  };
  
  //--  shortcuts  --//
  
  // TODO: this is wrong, because l2g is destroyed. EIther use smart pointers, OR pass along the
  // reference to vertexOut/facesOut to l2g.
  
  template < typename TVerticesIn, typename TFacesIn, typename TNeighbors, typename TVertexNode, typename TFaceNode >
  void list2graph_mesh_conversion
  (
    const TVerticesIn & verticesIn,
    const TFacesIn & facesIn,
    const TNeighbors & neighc,
    std::list<TVertexNode> & verticesOut,
    std::list<TFaceNode> & facesOut
  )
  {
    std::vector<std::list<std::size_t> > faceIndices = til::invertIndices(facesIn);
    List2GraphMeshConvertor<TVertexNode, TFaceNode, TVerticesIn, TFacesIn, TNeighbors> l2g;
    l2g(verticesIn, facesIn, faceIndices, neighc, verticesOut, facesOut);
  }

  template < typename TVerticesIn, typename TFacesIn, typename TVertexNode, typename TFaceNode >
  void list2graph_mesh_conversion
  (
    const TVerticesIn & verticesIn,
    const TFacesIn & facesIn,
    std::list<TVertexNode> & verticesOut,
    std::list<TFaceNode> & facesOut
  )
  {
    shared_ptr<std::vector<std::vector<std::size_t> > > neighc = circular_neighborhoods(verticesIn, facesIn);
    list2graph_mesh_conversion(verticesIn, facesIn, *neighc, verticesOut, facesOut);
  }
  
  //---------------------------------------------------------------------------
  
  // TODO: remove this crap
  
  template < typename TMesh >
  class Graph2ListMeshConvertor
  {
  public: // typedefs
  
    //typedef std::map<typename std::list<TVertexNode>::const_iterator, std::size_t, ItComp<typename std::list<TVertexNode>::const_iterator> > VertexIndexMap;
    typedef std::map<const void*, std::size_t> VertexIndexMap;
  
  public: // set & get
  
    VertexIndexMap & index_vertex() { return m_index_vertex; }
    TMesh & mesh() { return m_mesh; }
    
  public: // operators
  
    template < typename TVertexNode, typename TFaceNode >
    void operator()
    (
      const std::list<TVertexNode> & graph_vertices,
      const std::list<TFaceNode> & graph_faces
    )
    {
      {
        {
          getVertices(m_mesh).resize(graph_vertices.size());
          std::size_t c = 0;
          for (typename std::list<TVertexNode>::const_iterator i = graph_vertices.begin(); i != graph_vertices.end(); ++i)
          {
            getVertices(m_mesh)[c] = i->pos();
            m_index_vertex.insert(std::make_pair(static_cast<const void*>(&*i),c));
            ++c;
          }
        }
        {
          int c = 0;
          getFaceIndices(m_mesh).resize(graph_faces.size());
          for (typename std::list<TFaceNode>::const_iterator i = graph_faces.begin(); i != graph_faces.end(); ++i)
          {
            for (int j = 0; j < 3; ++j)
              getFaceIndices(m_mesh)[c][j] = m_index_vertex[static_cast<const void*>(&*(i->face[j]))];
            ++c;
          }
        }
      }
    }
    
  private: // data
    TMesh m_mesh;
    VertexIndexMap m_index_vertex;
  };
    
} // namespace til

#include "mesh_conversion.tpp"


#endif /*MESH_CONVERSION_H_*/
