#ifndef _MESHTRAITS_H_
#define _MESHTRAITS_H_

// includes from STL
#include <vector>

// includes from BOOST
#include <boost/array.hpp>
#include <boost/type_traits.hpp>

// include from TIL
#include "cyclic_iterator.h"
#include "Mesh.h"
#include "misc_sort.h"
#include "miscUtils.h"

// declarations
#include "til/til_declarations.h"


namespace til
{
  /// Default values for MeshTraits
  struct MeshTraits_default
  {
    enum { has_edge_indices       = 0 };
    enum { has_neighbor_indices   = 0 };
    enum { has_faces_indices      = 0 };
  };
  
  /// Traits for meshes.
  /// Mesh algorithms rely on these traits only. This allows to adapt the
  /// algorithms to any mesh type by simply specializing traits for our favorite
  /// mesh class
  template < class TMesh >
  struct MeshTraits {};
  
  
  
  
  /*
  template < typename TMesh >
  const typename MeshTraits<TMesh>::Vertex & 
  getVertex(std::size_t i, const TMesh & mesh)
  {
    return getVertices(mesh)[i];
  }
  
  template < typename TMesh >
  const typename MeshTraits<TMesh>::Vertex &
  getVertex(const typename MeshTraits<TMesh>::Vertex *v, const TMesh &)
  {
    return *v;
  }
  */
  
  
  /////////////////////////// MESH //////////////////////////////////////////////
  
  /// Mesh traits for meshes of type Mesh.
  template < typename TParam >
  // NB: Deriving from TParam simply enable us to get all the typedefs.
  struct MeshTraits<Mesh<TParam> > : public MeshTraits_default, public TParam
  {
    // constants
    enum { has_faces_indices = 1 };
    
    // typedefs
    //typedef MeshFaceCollection<TParam> FaceCollection;
  };
  
  template < typename TParam >
  inline
  const typename MeshTraits<Mesh<TParam> >::VertexCollection &
  getVertices(const Mesh<TParam> & mesh) { return mesh.getVertices(); }
  
  template < typename TParam >
  inline
  typename MeshTraits<Mesh<TParam> >::VertexCollection &
  getVertices(Mesh<TParam> & mesh) { return mesh.getVertices(); }
  
  /*
  template < typename TParam >
  inline
  const typename MeshTraits<Mesh<TParam> >::FaceIndexCollection &
  getFaceIndices(const Mesh<TParam> & mesh) { return mesh.getFaceIndices(); }
  
  template < typename TParam >
  inline
  typename MeshTraits<Mesh<TParam> >::FaceIndexCollection &
  getFaceIndices(Mesh<TParam> & mesh) { return mesh.getFaceIndices(); }
  */
  
  inline
  const MeshTraits<Mesh1>::FaceIndexCollection &
  getFaceIndices(const Mesh1 & mesh) { return mesh.getFaceIndices(); }
  
  inline
  const MeshTraits<Mesh2>::FaceIndexCollection &
  getFaceIndices(const Mesh2 & mesh) { return mesh.getFaceIndices(); }
  
  inline
  MeshTraits<Mesh1>::FaceIndexCollection &
  getFaceIndices(Mesh1 & mesh) { return mesh.getFaceIndices(); }
  
  inline
  MeshTraits<Mesh2>::FaceIndexCollection &
  getFaceIndices(Mesh2 & mesh) { return mesh.getFaceIndices(); }
  
  /*
  template < typename TParam >
  class MeshFaceCollection
  {
  public: // constructors & destructor
    
    MeshFaceCollection(Mesh<TParam> & mesh) :
      m_mesh(mesh)
    {    
    }
  
  public: // functions
  
  
    class ConstIterator
    {
    / *    
    public: // Face
  
      class Face
      {
      public: // constructors & destructor
        Face(const AimsSurface<D,T> & mesh, typename MeshTraits<AimsSurface<D,T> >::FaceIndex &f)
          : m_mesh(mesh), m_f(f) {}
    
      public: // functions
        typename MeshTraits<AimsSurface<D,T> >::Vertex &
          operator[](int i)
        {
          return getVertices(m_mesh)[m_f[i]];
        }
    
      private: // data
        const AimsSurface<D,T> & m_mesh;
        typename MeshTraits<AimsSurface<D,T> >::FaceIndex & m_f;
      };
    * /
    
    public: // constructors & destructor
      ConstIterator(const Mesh<TParam> & mesh, typename MeshTraits<Mesh<TParam> >::FaceIndexCollection::iterator i)
        : m_vertices(getVertices(mesh)), m_i(i) {}
    
    public: // functions
      inline void operator++() { ++m_i; }
      inline const ConstIterator & operator*() const { return (*this); }
      inline bool operator!=(const ConstIterator &i) const { return m_i != i.m_i; }
  
      inline
      const typename MeshTraits<Mesh<TParam> >::Vertex &
      operator[](int i) const
      {
        return m_vertices[(*m_i)[i]];
      }
      
  
    private: // data
  
      const typename MeshTraits<Mesh<TParam> >::VertexCollection & m_vertices;
      //const AimsSurface<D,T> & m_mesh;
      typename MeshTraits<Mesh<TParam> >::FaceIndexCollection::iterator m_i;
      //Face m_face;
    };
    
    
    typedef ConstIterator const_iterator;
  
    ConstIterator begin()  { return ConstIterator(m_mesh, getFaceIndices(m_mesh).begin()); }
    ConstIterator end()    { return ConstIterator(m_mesh, getFaceIndices(m_mesh).end()); }
  
  private: // data
    Mesh<TParam> & m_mesh;
  };
  */
  
  /*
  template < >
  class MeshFaceCollection<MeshParam2>
  {
  public: // constructors & destructor
    
    MeshFaceCollection(Mesh<MeshParam2> & mesh) :
      m_mesh(mesh)
    {    
    }
  
  public: // functions
  
  
    class ConstIterator
    {
    / *    
    public: // Face
  
      class Face
      {
      public: // constructors & destructor
        Face(const AimsSurface<D,T> & mesh, typename MeshTraits<AimsSurface<D,T> >::FaceIndex &f)
          : m_mesh(mesh), m_f(f) {}
    
      public: // functions
        typename MeshTraits<AimsSurface<D,T> >::Vertex &
          operator[](int i)
        {
          return getVertices(m_mesh)[m_f[i]];
        }
    
      private: // data
        const AimsSurface<D,T> & m_mesh;
        typename MeshTraits<AimsSurface<D,T> >::FaceIndex & m_f;
      };
    * /
    
    public: // constructors & destructor
      ConstIterator(const Mesh<MeshParam2> & mesh,  MeshTraits<Mesh<MeshParam2> >::FaceIndexCollection::iterator i)
        : m_vertices(getVertices(mesh)), m_i(i) {}
    
    public: // functions
      inline void operator++() { ++m_i; }
      inline const ConstIterator & operator*() const { return (*this); }
      inline bool operator!=(const ConstIterator &i) const { return m_i != i.m_i; }
  
  
      inline
      const  MeshTraits<Mesh<MeshParam2> >::Vertex &
      operator[](int i) const
      {
        return *((*m_i)[i]);
      }
  
      / *
      inline
      typename boost::enable_if<boost::is_same<MeshParam2, MeshParam1>, 
        const typename MeshTraits<Mesh<MeshParam2> >::Vertex &>::type
      get(int i) const
      {
        return m_vertices[(*m_i)[i]];
      }
  
      inline
      typename boost::enable_if<boost::is_same<MeshParam2, MeshParam2>, 
        const typename MeshTraits<Mesh<MeshParam2> >::Vertex &>::type
      get(int i) const
      {
        return *((*m_i)[i]);
      }
      * /
  
    private: // data
  
      const  MeshTraits<Mesh<MeshParam2> >::VertexCollection & m_vertices;
      //const AimsSurface<D,T> & m_mesh;
       MeshTraits<Mesh<MeshParam2> >::FaceIndexCollection::iterator m_i;
      //Face m_face;
    };
    
    typedef ConstIterator const_iterator;
  
    ConstIterator begin()  { return ConstIterator(m_mesh, getFaceIndices(m_mesh).begin()); }
    ConstIterator end()    { return ConstIterator(m_mesh, getFaceIndices(m_mesh).end()); }
  
  private: // data
    Mesh<MeshParam2> & m_mesh;
  };
  */
  
  /*
  template < typename TParam >
  inline
  MeshFaceCollection<TParam>
  getFaces(Mesh<TParam> & mesh)
  {
    return MeshFaceCollection<TParam>(mesh);
  }
  */
  
  /*
  template < typename TMesh, typename TEdgeCollection >
  class MeshEdgeExtractor
  {
  private: // typedefs
    typedef typename MeshTraits<TMesh>::FaceIndex::value_type   VertexIndex;
    typedef numeric_array<VertexIndex,2>                        Edge;
    typedef std::set<Edge, Lexicographical_compare<Edge> >       UniqueEdges;

  public: // functions
    void operator()(const TMesh & mesh, TEdgeCollection & edges)
    {
      UniqueEdges uniqueEdges;
      this->_getEdges(mesh, uniqueEdges);
      edges.resize(uniqueEdges.size());
      std::copy(uniqueEdges.begin(), uniqueEdges.end(), edges.begin());
    }
    
  private: // functions

     void _getEdges(const TMesh & mesh, UniqueEdges & res)
    {
      typename MeshTraits<TMesh>::FaceIndexCollection::const_iterator iFic = getFaceIndices(mesh).begin();
      for (; iFic != getFaceIndices(mesh).end(); ++iFic)
      {
        typename MeshTraits<TMesh>::FaceIndex::const_iterator iFi1 = iFic->begin();
        const_cyclic_iterator<typename MeshTraits<TMesh>::FaceIndex> iFi2(*iFic, iFi1+1);
        //const_cyclic_iterator<typename MeshTraits<Mesh<TParam> >::FaceIndex> iFi2(iFic, iFi1+1);
        //const_cyclic_iterator<typename TParam::FaceIndexCollection::value_type> iFi2(iFic, iFi1+1);
        //const_cyclic_iterator<boost::array<Point<float,3>*, 3> > iFi2(iFic, iFi1+1);
        for (;iFi1 != iFic->end(); ++iFi1, ++iFi2)
        {
          // We sort vertex by index to input only edges (a,b) where
          // a < b. This is to ensure that the same edge is not
          // input twice as (a,b) and (b,a).
          //res.insert((*iFi1 < *iFi2 ? Edge(*iFi1, *iFi2) : Edge(*iFi2, *iFi1)));
          res.insert(til::sortedVector<Edge>(*iFi1, *iFi2));
        }
      }
    }
  };
  
  template < typename TParam >
  std::set<numeric_array<typename TParam::FaceIndex::value_type,2>,
    Lexicographical_compare<numeric_array<typename TParam::FaceIndex::value_type,2> > >
  getEdges(const Mesh<TParam> & mesh)
  {
    typedef typename TParam::FaceIndex::value_type VertexIndex;
    typedef numeric_array<VertexIndex,2> Edge;
    std::set<Edge, Lexicographical_compare<Edge> > res;
    typename TParam::FaceIndexCollection::const_iterator iFic = getFaceIndices(mesh).begin();
    for (; iFic != getFaceIndices(mesh).end(); ++iFic)
    {
      typename TParam::FaceIndex::const_iterator iFi1 = iFic->begin();
      const_cyclic_iterator<typename TParam::FaceIndex> iFi2(*iFic, iFi1+1);
      //const_cyclic_iterator<typename MeshTraits<Mesh<TParam> >::FaceIndex> iFi2(iFic, iFi1+1);
      //const_cyclic_iterator<typename TParam::FaceIndexCollection::value_type> iFi2(iFic, iFi1+1);
      //const_cyclic_iterator<boost::array<Point<float,3>*, 3> > iFi2(iFic, iFi1+1);
      for (;iFi1 != iFic->end(); ++iFi1, ++iFi2)
      {
        // We sort vertex by index to input only edges (a,b) where
        // a < b. This is to ensure that the same edge is not
        // input twice as (a,b) and (b,a).
        //res.insert((*iFi1 < *iFi2 ? Edge(*iFi1, *iFi2) : Edge(*iFi2, *iFi1)));
        res.insert(til::sortedVector<Edge>(*iFi1, *iFi2));
      }
    }
    return res;
  }

  template < typename TResult, typename TParam >
  void
  getEdges(const Mesh<TParam> & mesh, TResult & res)
  {
    MeshEdgeExtractor<Mesh<TParam>, TResult> edgeExtractor;
    edgeExtractor(mesh, res);
  }
  */
  
  //////////////////////////////////////////////////////////////////////////////
  
  template < >
  struct MeshTraits<Mesh_N> : public MeshTraits<Mesh1>
  {
    // constants
    enum { has_neighbor_indices = 1 };
    
    // typedefs
    typedef std::vector<std::size_t>        NeighborIndex;
    typedef std::vector<NeighborIndex>      NeighborIndexCollection;
  };
  
  
  inline
  MeshTraits<Mesh_N>::NeighborIndexCollection &
  getNeighborIndices(Mesh_N & mesh)
  {
    return mesh.getNeighborIndices();
  }
  
  inline
  const MeshTraits<Mesh_N>::NeighborIndexCollection &
  getNeighborIndices(const Mesh_N & mesh)
  {
    return mesh.getNeighborIndices();
  }
  
  template < >
  struct MeshTraits<Mesh2_N> : public MeshTraits<Mesh2>
  {
    // constants
    enum { has_neighbor_indices = 1 };
    
    // typedefs
    typedef Mesh2_N::NeighborIndex                 NeighborIndex;
    typedef Mesh2_N::NeighborIndexCollection       NeighborIndexCollection;
  };
  
  
  inline
  MeshTraits<Mesh2_N>::NeighborIndexCollection &
  getNeighborIndices(Mesh2_N & mesh)
  {
    return mesh.getNeighborIndices();
  }
  
  inline
  const MeshTraits<Mesh2_N>::NeighborIndexCollection &
  getNeighborIndices(const Mesh2_N & mesh)
  {
    return mesh.getNeighborIndices();
  }
  
  
  //////////////////////////////////////////////////////////////////////////////
  // Adding and removing attributes to a mesh
  //////////////////////////////////////////////////////////////////////////////
  
  
  /// A trait for knowing the return type of a Mesh when adding
  /// redundant attributes such as normals, neighbors, etc.
  // TODO: right now this kind of stuff is really a pain in the asset, and
  // it will only get more painful as the number of possible attribute grows, but
  // I can't really find an alternative where all of this could be done
  // automatically. Namely, the problem is to have only one possible class for a 
  // given set of added attributes, i.e. we don't want to have both AddNormal<AddNeighbors>
  // and AddNeighbors<AddNormal>. Even more important, we don't want to add twice
  // the same attributes, like AddNormal<AddNeighbors<AddNormal>>.
  template < class TMesh >
  struct MeshAttributes
  {
  };
  
  /// Mesh attributes traits for Mesh
  /*
  template <typename TMeshParam>
  struct MeshAttributes<Mesh<TMeshParam> >
  {
    typedef typename detail::AddNeighborIndiceAttribute<Mesh<TMeshParam> >     AddNeighborIndices;
    typedef typename Mesh<TMeshParam>                                          RemoveNeighborIndices;
    typedef typename detail::AddNormalAttribute<Mesh<TMeshParam> >             AddNormals;
    typedef typename Mesh<TMeshParam>                                          RemoveNormals;
  };
  */
  
  /// Mesh attributes traits for Mesh_N
  template <>
  struct MeshAttributes<Mesh_N>
  {
    typedef Mesh_N    AddNeighborIndices;
    typedef Mesh1      RemoveNeighborIndices;
    typedef Mesh_NNo  AddNormals;
    typedef Mesh_N    RemoveNormals;
  };
  
  template <>
  struct MeshAttributes<Mesh_No>
  {
    typedef Mesh_NNo  AddNeighborIndices;
    typedef Mesh_No   RemoveNeighborIndices;
    typedef Mesh_NNo  AddNormals;
    typedef Mesh_No   RemoveNormals;
  };
  
  template <>
  struct MeshAttributes<Mesh_NNo>
  {
    typedef Mesh_NNo  AddNeighborIndices;
    typedef Mesh_No   RemoveNeighborIndices;
    typedef Mesh_NNo  AddNormals;
    typedef Mesh_N    RemoveNormals;
  };
  
  
  //////////////// Get methods
  
  
  /// Returns the i-th vertex of a face of a mesh of type Mesh1.
  inline
  const MeshTraits<Mesh1>::Vertex &
  getFaceVertex
  (
   const Mesh1 & mesh,                        ///< The mesh
   const std::vector<numeric_array<std::size_t, 3> >::const_iterator & iFC,    ///< An iterator pointing on the face
   int i                                      ///< The number of the face point
  )
  {
    return mesh.getVertices()[(*iFC)[i]];
  }
  
  /// Returns the i-th vertex of a face of a mesh of type Mesh2.
  inline
  const MeshTraits<Mesh2>::Vertex &
  getFaceVertex
  (
   const Mesh2 &,                           ///< The mesh (actually unused here)
   const MeshTraits<Mesh2>::FaceIndexCollection::const_iterator & iFC, ///< An iterator pointing on the face
   int i                                    ///< The number of the face
  )
  {
    return *((*iFC)[i]);
  }
  
  
  /// Returns the vertex neighbors of a mesh of type Mesh_N given by its index.
  inline
  const MeshTraits<Mesh_N>::Vertex &
  getVertexNeighbor
  (
   const Mesh_N & mesh,                                       ///< The mesh
   const MeshTraits<Mesh_N>::NeighborIndex::value_type & iNi  ///< The neighbor index
  )
  {
    return mesh.getVertices()[iNi];
  }
  
  /// Returns the vertex neighbors of a mesh of type Mesh_N given by its index.
  inline
  const numeric_array<float,3> &
  getVertexNeighbor
  (
   const Mesh2_N &,                                           ///< The mesh (actually unused here)
   const MeshTraits<Mesh2_N>::NeighborIndex::value_type & iNi ///< The neighbor index
  )
  {
    return *iNi;
  }
  
  
  /// Return the face indices of a mesh of type AimsTimeSurface.
  /// Conveniance macro for looping through all faces of a mesh
  #define TIL_FOR_ALL_CONST_FACES(mesh)                                                       \
  typename MeshTraits<TMesh>::FaceIndexCollection::const_iterator iFace;                      \
  for (iFace = getFaceIndices(mesh).begin(); iFace != getFaceIndices(mesh).end(); ++iFace)    \
  
  /// Conveniance macro for looping through all vertices of a mesh
  #define TIL_FOR_ALL_CONST_VERTICES(mesh)                                                    \
  typename MeshTraits<TMesh>::VertexCollection::const_iterator iVertex;                       \
  for (iVertex = getVertices(mesh).begin(); iVertex != getVertices(mesh).end(); ++iVertices)  \
  
} // namespace til

#endif //_MESHTRAITS_H_
