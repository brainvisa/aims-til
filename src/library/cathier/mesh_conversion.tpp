#ifndef TIL_MESH_CONVERSION_TPP_
#define TIL_MESH_CONVERSION_TPP_

namespace til
{
  
  //---------------------------------------------------------------------------
  
  template < typename TVertexCollection, typename TFaceCollection, typename TGraphVertices, typename TGraphFaces >
  //template < typename TVertexCollection, typename TFaceCollection >
  //template < typename TGraphVertices, typename TGraphFaces >
  //void Graph2ListMeshConvertor2<TVertexCollection, TFaceCollection>
  void Graph2ListMeshConvertor2 < TVertexCollection, TFaceCollection, TGraphVertices, TGraphFaces >
  ::operator()
  (
    const TGraphVertices & graph_vertices,
    const TGraphFaces & graph_faces,
    TVertexCollection & vertices,
    TFaceCollection & faces
  )
  {
    // TODO: use accessors instead
    {
      vertices.resize(graph_vertices.size());
      std::size_t c = 0;
      for (typename TGraphVertices::const_iterator i = graph_vertices.begin(); i != graph_vertices.end(); ++i)
      {
        vertices[c] = i->pos();
        m_index_vertex.insert(std::make_pair(static_cast<const void*>(&*i),c));
        ++c;
      }
    }
    {
      int c = 0;
      faces.resize(graph_faces.size());
      for (typename TGraphFaces::const_iterator i = graph_faces.begin(); i != graph_faces.end(); ++i)
      {
        for (int j = 0; j < 3; ++j)
          faces[c][j] = m_index_vertex[static_cast<const void*>(&*(i->face[j]))];
        ++c;
      }
    }
  }
  
  //---------------------------------------------------------------------------

  template < typename TVertexNode, typename TFaceNode, typename TVertexCollection, typename TFaceCollection, typename TNeighborCollection >
  //template < typename TVertexNode, typename TFaceNode >
  //template < typename TVertexCollection, typename TFaceCollection, typename TNeighborCollection >
  //void List2GraphMeshConvertor<TVertexNode, TFaceNode>::
  void List2GraphMeshConvertor<TVertexNode, TFaceNode, TVertexCollection, TFaceCollection, TNeighborCollection>::
  operator()
  (
      TVertexCollection                     const & vertices,             ///< [input] mesh vertices
      TFaceCollection                       const & faceIndices,          ///< [input] mesh faces
      std::vector<std::list<std::size_t> >  const & invertedFaceIndices,  ///< [input] inverted face indices
      TNeighborCollection                   const & neighbors,            ///< [input] vertex neighborhoods
      std::list<TVertexNode>                      & graph_vertices,
      std::list<TFaceNode>                        & graph_faces
  )
  {
    assert(vertices.size() == invertedFaceIndices.size());
    assert(vertices.size() == neighbors.size());

    m_index_vertex.resize(vertices.size());
    m_index_face.resize(faceIndices.size());
    {
      // push all vertices into the vertex list and get an index2iterator translation
      for (std::size_t i = 0; i < vertices.size(); ++i)
      {
        // create a node with current position
        TVertexNode m;
        m.pos() = vertices[i];
        // push this node at the end of the vertex list
        m_index_vertex[i] = graph_vertices.insert(graph_vertices.end(), m);
      }
      
      // now push all faces into the face list
      for (std::size_t i = 0; i < til::size(faceIndices); ++i)
      {
        // create a face with current face indices
        TFaceNode f;
        for (int j = 0; j < 3; ++j)
        {
          assert(vertices.size() > faceIndices[i][j]);
          f.face[j] = m_index_vertex[faceIndices[i][j]];
        }
        // push this face at the end of the face list
        m_index_face[i] = graph_faces.insert(graph_faces.end(), f);
      }
      
      // now we need a second pass on the vertices to fill in the neighbor
      // and the faces structure.
      for (std::size_t i = 0; i < vertices.size(); ++i)
      {
        for (std::list<std::size_t>::const_iterator j = invertedFaceIndices[i].begin(); j != invertedFaceIndices[i].end(); ++j)
        {
          assert(faceIndices.size() > *j);
          m_index_vertex[i]->faces().push_back(m_index_face[*j]);
        }
        for (typename TNeighborCollection::value_type::const_iterator iN = neighbors[i].begin();
             iN != neighbors[i].end(); ++iN)
        {
          assert(vertices.size() > *iN);
          m_index_vertex[i]->neighbors().push_back(m_index_vertex[*iN]);
        }
      }
    }
  }
}

#endif /*MESH_CONVERSION_TPP_*/
