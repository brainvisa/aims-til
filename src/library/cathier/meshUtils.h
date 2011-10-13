#ifndef _MESHUTILS_H_
#define _MESHUTILS_H_

// includes from STL
#include <cmath>
#include <numeric>    // accumulate


// includes from BOOST
#include <boost/shared_ptr.hpp>
using boost::shared_ptr;

// includes from AIMS
#include "aims/mesh/surface.h"
#include "aims/vector/vector.h"

// includes from TIL
#include "til/functors.h"
#include "til/is_traits.h"
#include "til/loop.h"
#include "til/numeric_array.h"
#include "til/TExpr.h"
#include "til/traits.h"
//#include "til/Vector3.h"

// includes from TIL
#include "cyclic_iterator.h"
//#include "convert.h"
#include "geometrics.h"
#include "globalTraits.h"
#include "Mesh.h"
#include "MeshTraits.h"
#include "miscUtils.h"


// declarations
template < typename T > class MeshTraits;




AimsSurfaceTriangle * 
makeSphere
(
 const Point3df &center, ///< the center of the sphere
 float radius,           ///< the radius of the spere
 int iter                ///< the number of subdivision iterations
 );

namespace til
{
  
  /// Return the index of vertex v[i] -- which is, obviously, i.
  template < typename TMesh >
  inline
  std::size_t getVertexNumber(const TMesh &, std::size_t i) { return i; }
  
  /// Return the index of the vertex pointed at by *p.
  template < typename TMesh, typename TVertex >
  // make sure that Vertex Collection is contiguous in memory so that fast
  // arithmetic can be used
  //typename boost::enable_if<boost::is_same<typename MeshTraits<TMesh>::VertexCollection::difference_type, std::ptrdiff_t>, std::size_t >::type
  inline
  typename boost::enable_if<is_same_naked<typename MeshTraits<TMesh>::Vertex, TVertex>, std::size_t>::type
  getVertexNumber(const TMesh & mesh, TVertex *p)
  {
    // NB: it seems that for some reason, we don't have to divide by sizeof.
    // I would like to be sure this is standard...
    //return (p - &(getVertices(mesh)[0])) / sizeof(typename MeshTraits<TMesh>::Vertex);
    return p - &(getVertices(mesh)[0]);
  
  }
  
  
  namespace detail
  {
/*    template < class TMeshFrom, class TMeshTo >
    void convert_mesh_1(const TMeshFrom & meshFrom, TMeshTo & meshTo)
    {
      // allocate space for vertices
      getVertices(meshTo).resize(getVertices(meshFrom).size());
      // convert vertices
      using namespace til::expr;
      std::transform(getVertices(meshFrom).begin(), getVertices(meshFrom).end(), getVertices(meshTo).begin(), 
        ::til::functor::Cast<typename MeshTraits<TMeshTo>::Vertex, typename MeshTraits<TMeshFrom>::Vertex>());
      //loop(getVertices(meshTo), getVertices(meshFrom), Convert2<typename MeshTraits<TMeshFrom>::Vertex, typename MeshTraits<TMeshTo>::Vertex>());
      //convert(getVertices(meshFrom), getVertices(meshTo));
      // allocate space for face indices
      getFaceIndices(meshTo).resize(getFaceIndices(meshFrom).size());
      // convert face indices
      detail::loop_xx(castTo(*_2, *_1), getFaceIndices(meshFrom), getFaceIndices(meshTo));
      //std::transform(getFaceIndices(meshFrom).begin(), getFaceIndices(meshFrom).end(), getFaceIndices(meshTo).begin(),
      //loop(getFaceIndices(meshTo), getFaceIndices(meshFrom), Convert());
      //convert(getFaceIndices(meshFrom), getFaceIndices(meshTo));
    }
*/

    inline void convert_aimsmeshTomesh1(const AimsTimeSurface<3, Void> & meshFrom, til::Mesh1 & meshTo)
    {
      // allocate space for vertices
      int size = meshFrom.vertex().size();
      getVertices(meshTo).resize(size);

      // convert vertices
      using namespace til::expr;
      std::transform(meshFrom.vertex().begin(), meshFrom.vertex().end(), getVertices(meshTo).begin(), 
        ::til::functor::Cast< til::MeshTraits< til::Mesh1 >::Vertex, Point3df >());
      size = meshFrom.polygon().size();
      getFaceIndices(meshTo).resize(size);

      // convert face indices
      detail::loop_xx(castTo(*_2, *_1), meshFrom.polygon(), getFaceIndices(meshTo));
    }

    inline void convert_mesh1Toaimsmesh(const til::Mesh1 & meshFrom, AimsTimeSurface<3, Void> & meshTo)
    {
      // allocate space for vertices
      int size = getVertices(meshFrom).size();
      meshTo.vertex().resize(size);

      // convert vertices
      using namespace til::expr;
      std::transform(getVertices(meshFrom).begin(), getVertices(meshFrom).end(), meshTo.vertex().begin(), 
        ::til::functor::Cast< Point3df, til::MeshTraits< til::Mesh1 >::Vertex >());
      size = getFaceIndices(meshFrom).size();
      meshTo.polygon().resize(size);

      // convert face indices
      detail::loop_xx(castTo(*_2, *_1), getFaceIndices(meshFrom), meshTo.polygon());
    }
  
    ///  When indices of from are numbers and two are pointers
    template < class TMeshFrom, class TMeshTo >
    inline void convert_mesh_2(const TMeshFrom & meshFrom, TMeshTo & meshTo)
    {
      // allocate space for vertices
      getVertices(meshTo).resize(getVertices(meshFrom).size());
      // convert vertices
      //convert(getVertices(meshFrom), getVertices(meshTo));
      //std::transform(getVertices(meshFrom).begin(), getVertices(meshFrom).end(), getVertices(meshTo).begin,
      //  functor::Cast<typename TMeshTo::Vertex, typename TMeshFrom::Vertex>());
      {
        using namespace til::expr;
        detail::loop_xx(castTo(*_1, *_2), getVertices(meshTo), getVertices(meshFrom));
      }
      //loop(getVertices(meshTo), getVertices(meshFrom), Convert());
      {
        //allocate space for face indices
        getFaceIndices(meshTo).resize(getFaceIndices(meshFrom).size());
        // copy face indices
        typename MeshTraits<TMeshFrom>::FaceIndexCollection::const_iterator iff = getFaceIndices(meshFrom).begin();
        typename MeshTraits<TMeshTo>::FaceIndexCollection::iterator ift = getFaceIndices(meshTo).begin();
        for (; iff != getFaceIndices(meshFrom).end(); ++iff, ++ift)
        {                
          (*ift)[0] =  &(getVertices(meshTo)[(*iff)[0]]);
          (*ift)[1] =  &(getVertices(meshTo)[(*iff)[1]]);
          (*ift)[2] =  &(getVertices(meshTo)[(*iff)[2]]);
        } 
      }
    }
    template < class TMeshFrom, class TMeshTo >
    inline void convert_mesh_3(const TMeshFrom & meshFrom, TMeshTo & meshTo)
    {
      // allocate space for vertices
      getVertices(meshTo).resize(getVertices(meshFrom).size());
      // copy vertices
      //convert(getVertices(meshFrom), getVertices(meshTo));
      //loop(getVertices(meshTo), getVertices(meshFrom), Convert());
      using namespace til::expr;
      detail::loop_xx(castTo(*_1, *_2), getVertices(meshTo), getVertices(meshFrom));
      //std::transform(getVertices(meshFrom).begin(), getVertices(meshFrom).end(), getVertices(meshTo).begin(),
      //  functor::Cast<typename TMeshTo::Vertex, typename TMeshFrom::Vertex>());
      {
        //allocate space for face indices
        getFaceIndices(meshTo).resize(getFaceIndices(meshFrom).size());
        // copy face indices
        typename MeshTraits<TMeshFrom>::FaceIndexCollection::const_iterator iff = getFaceIndices(meshFrom).begin();
        typename MeshTraits<TMeshTo>::FaceIndexCollection::iterator ift = getFaceIndices(meshTo).begin();
        for (; iff != getFaceIndices(meshFrom).end(); ++iff, ++ift)
        {
          (*ift)[0] =  getVertexNumber(meshFrom, (*iff)[0]);
          (*ift)[1] =  getVertexNumber(meshFrom, (*iff)[1]);
          (*ift)[2] =  getVertexNumber(meshFrom, (*iff)[2]);
        } 
      }
    }
  }
  
  
  /*
  /// Computes the centroid of a flat polygon in space
  template < typename TVertexCollection >
  Point3dd
  getCentroid(const TVertexCollection &vertices)
  {
    typename TVertexCollection::const_iterator iVertex;
    typename TVertexCollection::const_iterator iNextVertex;
    Point3dd res(0.0, 0.0, 0.0);
    double area, totalArea = 0.0;
    
    // Center all computations around average for numerical precision
    Point3dd m = getAverage(vertices);
    
    // Loop through all pair of vertices
    iVertex = vertices.begin();
    iNextVertex = vertices.begin();
    ++iNextVertex;
    for (; iNextVertex != vertices.end(); ++iVertex, ++iNextVertex)
    {
      // compute the (signed) area of the triangle (Center, Vertex, NextVertex)
      area = dot(*iNextVertex - *iVertex, *iVertex - m);
      // add it to the total area
      totalArea += area;
      // add weighted centroid of this triangle to the result
      res += area * (*iVertex + *iNextVertex + m) / 3.0;
    }
    // Do the same thing for the last pair
    area = dot(*(vertices.begin()) - *iNextVertex, *iNextVertex - m);
    totalArea += area;
    res += area * (*(vertices.begin()) + *iNextVertex + m) / 3.0;
    // Divide by the total area
    res /= totalArea;
    // Add center
    res += m;
    return res;
  }
  */

  //---------------------------------------------------------------------------

  template < typename TFaceCollection >
  class Get_edges_and_faces
  {
  public: // typedefs
    typedef std::size_t index_type;
    typedef std::pair<index_type,index_type> Edge;
    typedef std::map<Edge, std::vector<index_type> > Edges;
  public: // functions
    void operator()
    (
      TFaceCollection const & faces,
      std::vector<Edge> & res,
      std::vector<std::vector<index_type> > & newfaces
    );
  private: // functions
    void insert_edge
    (
      Edges & edges,
      Edge edge,
      std::size_t i
    );
  };
  
  //---------------------------------------------------------------------------

  template < typename TFaceCollection >
  inline void Get_edges_and_faces<TFaceCollection>::insert_edge
  (
    Edges & edges,
    Edge edge,
    std::size_t i
  )
  {
    typename Edges::iterator pos = edges.find(edge);
    if (pos != edges.end())
    {
      pos->second.push_back(i);
    }
    else
    {
      edges.insert(typename Edges::value_type(edge, std::vector<std::size_t>(1,i)));
    }
  }

  //---------------------------------------------------------------------------

  template < typename TFaceCollection >
  inline void Get_edges_and_faces<TFaceCollection>::operator()
  (
    TFaceCollection const & faces,
    std::vector<Edge> & res,
    std::vector<std::vector<index_type> > & newfaces
  )
  {
    Edges edges;
    typename TFaceCollection::const_iterator iFace;
    std::size_t i;
    for (i = 0, iFace = faces.begin(); iFace != faces.end(); ++iFace, ++i)
    { 
      // We sort vertex by index to input only edges (a,b) where
      // a < b. This is to ensure that the same edge is not
      // input twice as (a,b) and (b,a).
      insert_edge(edges, make_sorted_pair((*iFace)[0], (*iFace)[1]), i);
      insert_edge(edges, make_sorted_pair((*iFace)[1], (*iFace)[2]), i);
      insert_edge(edges, make_sorted_pair((*iFace)[2], (*iFace)[0]), i);      
    }
    // convert result into a std::vector
    res.resize(edges.size());
    //std::copy(edges.begin(), edges.end(), res.begin());      
    //std::transform(edges.begin(), edges.end(), res.begin(), Return_pair1<Edge, std::vector<index_type> >());
    
    newfaces.resize(faces.size());
    typename Edges::const_iterator ie;
    for (i = 0, ie = edges.begin(); ie != edges.end(); ++ie, ++i)
    {
      res[i] = ie->first;
      assert(ie->second.size() == 2);
      for (std::size_t j = 0; j < ie->second.size(); ++j)
      {
        newfaces[ie->second[j]].push_back(i);
      }
    }
  }

  //---------------------------------------------------------------------------

  template < typename TFaceCollection >
  inline void
  //shared_ptr<std::vector<std::pair<std::size_t, std::size_t> > >
  get_edges_and_faces
  (
    TFaceCollection const & faces,
    std::vector<std::pair<std::size_t, std::size_t> > & res,
    std::vector<std::vector<std::size_t> > & newfaces
  )
  {
    Get_edges_and_faces<TFaceCollection>()(faces, res, newfaces);
  }

  //---------------------------------------------------------------------------
  
  /// Collects edges from a set of faces
  template < typename TFaceCollection >
  inline shared_ptr<std::vector<std::pair<std::size_t, std::size_t> > >
  faces2edges(TFaceCollection const & faces)
  {
    typedef std::size_t index_type;
    typedef std::pair<index_type,index_type> Edge;
    
    std::set<Edge> edges;
    typename TFaceCollection::const_iterator iFace = faces.begin();
    for (; iFace != faces.end(); ++iFace)
    { 
      // We sort vertex by index to input only edges (a,b) where
      // a < b. This is to ensure that the same edge is not
      // input twice as (a,b) and (b,a).
      edges.insert(make_sorted_pair((*iFace)[0], (*iFace)[1]));
      edges.insert(make_sorted_pair((*iFace)[1], (*iFace)[2]));
      edges.insert(make_sorted_pair((*iFace)[2], (*iFace)[0]));      
    }
    // convert result into a std::vector
    shared_ptr<std::vector<Edge> > res(new std::vector<Edge>(edges.size()));
    std::copy(edges.begin(), edges.end(), res->begin());
    return res;
  }

  //---------------------------------------------------------------------------

  /// Collects edges from a set of neighbors.
  /// NB: This assumes that the set of neighbors of a point does *not* contain the point.
  template < typename TNeighborhoodCollection >
  inline void
  neighborhoods2edges
  (
    TNeighborhoodCollection const & neigh,                      ///< [input] The collection of neighborhoods
    std::vector<std::pair<std::size_t, std::size_t> > & edges   ///< [output] The computed edges
  )
  {
    typedef std::pair<std::size_t, std::size_t> Edge;
    std::set<Edge> my_edges;
    for (std::size_t i = 0; i < neigh.size(); ++i)
    {
      for (std::size_t j = 0; j < neigh[i].size(); ++j)
      {
        my_edges.insert(make_sorted_pair(i, neigh[i][j]));
      }
    }
    edges.resize(edges.size());
    std::copy(my_edges.begin(), my_edges.end(), edges.begin());
    // TODO: remove this test because it does not *have* to be the case, because neighborhoods may
    // not be symmetrical, even though I wish they are.
    assert(edges.size() * 2 == neigh.size());
  }
  
  //---------------------------------------------------------------------------

  /// Get point neighbors so that they form a loop around points.
  /// This is truely useful only for meshes that have a fixed number of vertices per faces, because
  /// the circular neighbor then also contains the neighboring face information in a natural and compact fashion.
  /// Note that if the triangles are always turning in the same direction (i.e. their normals are consistent
  /// surface-wise), then so will the circular neighbors, and in the same direction as the triangles.
  // TODO: add templation on return type and use "if (MeshTraits::has_face_indices)" etc.
  template < typename TFaceCollection >
  shared_ptr<std::vector<std::vector<typename TFaceCollection::value_type::value_type> > >
  circular_neighborhoods(TFaceCollection const & faces, std::size_t nVertices);

  //---------------------------------------------------------------------------

  template < typename TVertexCollection, typename TFaceCollection >
  inline
  shared_ptr<std::vector<std::vector<typename TFaceCollection::value_type::value_type> > >
  circular_neighborhoods(const TVertexCollection & vertices, const TFaceCollection & faces)
  { return circular_neighborhoods(faces, vertices.size()); }

  //---------------------------------------------------------------------------
  
  // TODO: add templation on return type and use "if (MeshTraits::has_face_indices)" etc.
  template < typename TMesh >
  inline shared_ptr<std::vector<std::vector<typename MeshTraits<TMesh>::FaceIndex::value_type> > >
  getNeighborIndices(const TMesh & mesh)
  {
    typedef typename MeshTraits<TMesh>::FaceIndex::value_type VertexIndex;
  
    // Allocate the result -- it should have as many elements as
    // the number of vertices of the mesh.
    std::vector<std::set<VertexIndex> > res(size(getVertices(mesh)));
    // for all faces
    for (typename MeshTraits<TMesh>::FaceIndexCollection::const_iterator iFic = getFaceIndices(mesh).begin();
         iFic != getFaceIndices(mesh).end(); ++iFic)
    {
      // for all couple of points on the face
      typename MeshTraits<TMesh>::FaceIndex::const_iterator iFi1 = iFic->begin();
      typename MeshTraits<TMesh>::FaceIndex::const_iterator iFi2 = iFi1+1;
      for (;iFi2 != iFic->end(); ++iFi1, ++iFi2)
      {
        res[getVertexNumber(mesh, *iFi2)].insert(*iFi1);
        res[getVertexNumber(mesh, *iFi1)].insert(*iFi2);
      }
      // we use iFi1 that should now points to the last element.
      // do not use container.end() instead! This points nowhere.
      res[getVertexNumber(mesh, *(iFic->begin()))].insert(*iFi1);
      res[getVertexNumber(mesh, *iFi1)].insert(*(iFic->begin()));
    }  
    //convert the result into a vector
    shared_ptr<std::vector<std::vector<VertexIndex> > > res2(new std::vector<std::vector<VertexIndex> >);
    allocate_sameSize(res, *res2);
    //detail::loop_cc(*res2, res, Convert());
    {
      using namespace til::expr;
      for (std::size_t i = 0; i < size(res); ++i)
      {
        detail::loop_xx(castTo(*_1, *_2), (*res2)[i], res[i]);      
      }
      //detail::loop_xx(castTo(*_1, *_2), *res2, res);
    }
    //loop(res, *res2, Convert());
    //convert(res, *res2);
    //convert<std::vector<std::set<std::size_t> >, std::vector<std::vector<std::size_t> > >(res, *res2);
    return res2;
  }

  //---------------------------------------------------------------------------

  template < typename TVertexCollection, typename TFaceIndexCollection >
  inline void getNeighborIndices
  (
    const TVertexCollection & vertices,               ///< Input mesh vertices
    const TFaceIndexCollection & faces,               ///< Input mesh faces
    std::vector<std::vector<std::size_t> > & neigh    ///< Output neighborhoods
  )
  {
    typedef typename TFaceIndexCollection::value_type FaceIndices;
    typedef std::size_t VertexIndex;
    // Allocate the result -- it should have as many elements as
    // the number of vertices of the mesh.
    std::vector<std::set<VertexIndex> > res(vertices.size());
    // for all faces
    for (typename TFaceIndexCollection::const_iterator iFic = faces.begin(); iFic != faces.end(); ++iFic)
    {
      // for all couple of points on the face
      typename FaceIndices::const_iterator iFi1 = iFic->begin();
      typename FaceIndices::const_iterator iFi2 = iFi1+1;
      for (;iFi2 != iFic->end(); ++iFi1, ++iFi2)
      {
        res[*iFi2].insert(*iFi1);
        res[*iFi1].insert(*iFi2);
      }
      // we use iFi1 that should now points to the last element.
      // do not use container.end() instead! This points nowhere.
      res[*(iFic->begin())].insert(*iFi1);
      res[*iFi1].insert(*(iFic->begin()));
    }
    //convert the result into a vector
    neigh.resize(res.size());
    {
      using namespace til::expr;
      for (std::size_t i = 0; i < size(res); ++i)
      {
        neigh[i].resize(res[i].size());
        detail::loop_xx(castTo(*_1, *_2), neigh[i], res[i]);
      }
    }
  }

  //---------------------------------------------------------------------------
  
  template < typename TParam >
  inline detail::AddNeighborIndexAttribute<Mesh<TParam> >
  addNeighborsToMesh(const Mesh<TParam> & mesh)
  {
    shared_ptr<typename detail::AddNeighborIndexAttribute<Mesh<TParam> >::NeighborIndexCollection > nbi = getNeighborIndices(mesh);
    return detail::AddNeighborIndexAttribute<Mesh<TParam> >(mesh, nbi);
  }
  
  
  /*
  MeshAttributes<Mesh1>::AddNeighborIndices
  addNeighborsToMesh(const Mesh1 & mesh);
  */
  
  namespace detail
  {
    /// Apply a functor for each pair (vertex, neighbor_of_vertex), given a mesh
    /// and another container.
    /// The functor should expect calls of the form f(Vertex, NeighboringVertex, ExtraVariable).
    template < typename TMesh, typename TContainer, typename TFunctor >
    inline void 
    for_each_neighbors_N(const TMesh & mesh, TContainer & c, TFunctor f)
    {
      // Allocate output with the same structure of neighbor indices of mesh
      // if the main container has the wrong size.
      /*
      if (size(c) != size(getNeighborIndices(mesh)))
      {
        allocate_sameSize(getNeighborIndices(mesh), c);
      }
      */
      assert(size(c) == size(getNeighborIndices(mesh)));
    
      // for all vertices
      typename MeshTraits<TMesh>::VertexCollection::const_iterator          iV        = getVertices(mesh).begin();
      typename MeshTraits<TMesh>::NeighborIndexCollection::const_iterator   iNic      = getNeighborIndices(mesh).begin();
      typename TContainer::iterator                                         iC        = c.begin();
      for (; iNic != getNeighborIndices(mesh).end(); ++iNic, ++iC, ++iV)
      {
        // for all neighbors of current vertex
        typename MeshTraits<TMesh>::NeighborIndex::const_iterator iNi = iNic->begin();
        typename TContainer::value_type::iterator iL = iC->begin();
        for (; iNi != iNic->end(); ++iNi, ++iL)
        {
          f(*iV, getVertexNeighbor(mesh, *iNi), *iL);
        }
      }  
    }
  
    /// Apply a functor for each pair (vertex, neighbor_of_vertex), given a mesh
    /// and another container.
    /// The functor should expect calls of the form f(Vertex, NeighboringVertex, ExtraVariable).
    template < typename TMesh, typename TContainer, typename TFunctor >
    inline void 
    for_each_neighbors_V(const TMesh & mesh, TContainer & c, TFunctor f)
    {
      /*
      // Allocate output with the same structure of neighbor indices of mesh
      // if the main container has the wrong size.
      if (size(c) != size(getNeighborIndices(mesh)))
      {
        allocate_sameSize(getNeighborIndices(mesh), c);
      }
      */
      assert(size(c) == size(getNeighborIndices(mesh)));
    
      // for all vertices
      typename MeshTraits<TMesh>::VertexCollection::const_iterator          iV        = getVertices(mesh).begin();
      typename MeshTraits<TMesh>::NeighborIndexCollection::const_iterator   iNic      = getNeighborIndices(mesh).begin();
      typename TContainer::iterator                                         iC        = c.begin();
      for (; iNic != getNeighborIndices(mesh).end(); ++iNic, ++iC, ++iV)
      {
        // for all neighbors of current vertex
        typename MeshTraits<TMesh>::NeighborIndex::const_iterator iNi = iNic->begin();
        for (; iNi != iNic->end(); ++iNi)
        {
          f(*iV, getVertexNeighbor(mesh, *iNi), *iC);
        }
      }  
    }
  }
  
  //---------------------------------------------------------------------------

  template < typename TMesh, typename TContainer, typename TFunctor >
  inline typename boost::enable_if_c<
    is_container<typename TContainer::value_type>::value
  >::type
  for_each_neighbors(const TMesh & mesh, TContainer & c1, TFunctor f)
  {
    detail::for_each_neighbors_N(mesh, c1, f);
  }
  
  //---------------------------------------------------------------------------

  template < typename TMesh, typename TContainer, typename TFunctor >
  inline typename boost::enable_if_c<
    !is_container<typename TContainer::value_type>::value
  >::type
  for_each_neighbors(const TMesh & mesh, TContainer & c1, TFunctor f)
  {
    detail::for_each_neighbors_V(mesh, c1, f);
  }
  
  //---------------------------------------------------------------------------
  
  
  namespace detail
  {
    /// Apply a functor for each pair (vertex, neighbor_of_vertex), given a mesh
    /// and two extra containers.
    /// The functor should expect calls of the form f(Vertex, NeighboringVertex, ExtraVariable1, ExtraVariable2).
    /// Here, both containers follow the structure of NeighborIndices (i.e. a container of
    /// container).
    template < typename TMesh, typename TContainer1, typename TContainer2, typename TFunctor >
    inline void
    for_each_neighbors_NN(const TMesh & mesh, TContainer1 & c1, TContainer2 & c2, TFunctor f)
    {
      /*
      // Allocate output with the same structure of neighbor indices of mesh
      // if the main container has the wrong size.
      if (size(c1) != size(getNeighborIndices(mesh)))
      {
        allocate_sameSize(getNeighborIndices(mesh), c1);        
      }
      if (size(c2) != size(getNeighborIndices(mesh)))
      {
        allocate_sameSize(getNeighborIndices(mesh), c2);        
      }
      */
      assert(size(c1) == size(c2));
      assert(size(c1) == size(getNeighborIndices(mesh)));
      
      
      // for all vertices
      typename MeshTraits<TMesh>::VertexCollection::const_iterator          iV        = getVertices(mesh).begin();
      typename MeshTraits<TMesh>::NeighborIndexCollection::const_iterator   iNic      = getNeighborIndices(mesh).begin();
      typename TContainer1::iterator                                        iC1       = c1.begin();
      typename TContainer2::iterator                                        iC2       = c2.begin();
      for (; iNic != getNeighborIndices(mesh).end(); ++iNic, ++iC1, ++iC2, ++iV)
      {
        // for all neighbors of current vertex
        typename MeshTraits<TMesh>::NeighborIndex::const_iterator iNi = iNic->begin();
        typename TContainer1::value_type::iterator iL1 = iC1->begin();
        typename TContainer2::value_type::iterator iL2 = iC2->begin();
        for (; iNi != iNic->end(); ++iNi, ++iL1, ++iL2)
        {
          f(*iV, getVertexNeighbor(mesh, *iNi), *iL1, *iL2);
        }
      }
    }
    
    /// Apply a functor for each pair (vertex, neighbor_of_vertex), given a mesh
    /// and two extra containers.
    /// Here, the first container should follow NeighborIndex structure (container of container)
    /// and the second the vertex structure (container)
    template < typename TMesh, typename TContainer1, typename TContainer2, typename TFunctor >
    inline void
    for_each_neighbors_NV(const TMesh & mesh, TContainer1 & c1, TContainer2 & c2, TFunctor f)
    {
      // Allocate output with the same structure of neighbor indices of mesh
      // if the main container has the wrong size.
      /*
      if (size(c1) != size(getNeighborIndices(mesh)))
      {
        allocate_sameSize(getNeighborIndices(mesh), c1);        
      }
      if (size(c2) != size(getVertices(mesh)))
      {
        allocate_sameSize(getVertices(mesh), c2);
      }
      */
      assert(size(c1) == size(c2));
      assert(size(c1) == size(getNeighborIndices(mesh)));
          
      // for all vertices
      typename MeshTraits<TMesh>::VertexCollection::const_iterator          iV        = getVertices(mesh).begin();
      typename MeshTraits<TMesh>::NeighborIndexCollection::const_iterator   iNic      = getNeighborIndices(mesh).begin();
      typename stditerator<TContainer1>::type                               iC1       = c1.begin();
      typename stditerator<TContainer2>::type                               iC2       = c2.begin();
      for (; iNic != getNeighborIndices(mesh).end(); ++iNic, ++iC1, ++iC2, ++iV)
      {
        // for all neighbors of current vertex
        typename MeshTraits<TMesh>::NeighborIndex::const_iterator iNi = iNic->begin();
        typename stditerator<typename TContainer1::value_type, TContainer1>::type iL1 = iC1->begin();
        //typename stditerator<typename TContainer1::value_type>::type iL1 = iC1->begin();
        for (; iNi != iNic->end(); ++iNi, ++iL1)
        {
          f(*iV, getVertexNeighbor(mesh, *iNi), *iL1, *iC2);
        }
        /*
        if (TFunctor::PerVertex) {
          f.perVertex(v[*iNi], *iC2, size(*iNic));
        }*/
      }  
    }
  }
  
  //---------------------------------------------------------------------------
  
  /// Apply a functor for each pair (vertex, neighbor_of_vertex), given a mesh
  /// and two extra containers.
  /// The functor should expect calls of the form f(Vertex, NeighboringVertex, ExtraVariable1, ExtraVariable2).
  /// Here, both containers follow the structure of NeighborIndices (i.e. a container of
  /// container).
  template < typename TMesh, typename TContainer1, typename TContainer2, typename TFunctor >
  inline typename boost::enable_if_c<
    is_container<typename TContainer1::value_type>::value &&
    is_container<typename TContainer2::value_type>::value
  >::type
  for_each_neighbors(const TMesh & mesh, TContainer1 & c1, TContainer2 & c2, TFunctor f)
  {
    detail::for_each_neighbors_NN(mesh, c1, c2, f);
  }
  
  /// Apply a functor for each pair (vertex, neighbor_of_vertex), given a mesh
  /// and two extra containers.
  /// Here, the first container should follow NeighborIndex structure (container of container)
  /// and the second the vertex structure (container)
  template < typename TMesh, typename TContainer1, typename TContainer2, typename TFunctor >
  inline typename boost::enable_if_c<
    is_container<typename TContainer1::value_type>::value &&
    !is_container<typename TContainer2::value_type>::value
  >::type
  for_each_neighbors(const TMesh & mesh, TContainer1 & c1, TContainer2 & c2, TFunctor f)
  {
    detail::for_each_neighbors_NV(mesh, c1, c2, f);
  }
  
  //---------------------------------------------------------------------------
  
  /// Computes edge lengths of a mesh.
  /// Actually, the lengths are computed between a point and its neighbors, not
  /// between edges. This means that the same length is computed twice. This also
  /// explains why the second argument is an array of arrays, not a simple array.
  /// NB: their is only partial checking for the allocation of the output. For
  /// efficiency reasons, we allocate the output only if the main container has
  /// the wrong size. If not, we do not check for indiviual sizes of all the other
  /// containers.
  template < typename TPrecision, typename TMesh >
  inline typename boost::enable_if_c<
    // Mesh must have neighbor indices
    MeshTraits<TMesh>::has_neighbor_indices
  >::type
  getEdgeLengths
  (
   const TMesh &mesh,                               ///< input mesh
   std::vector<std::vector<TPrecision> > & lengths  ///< output containing edge lenghts
  )
  {
    if (size(lengths) != size(getNeighborIndices(mesh)))
    {
      allocate_sameSize(getNeighborIndices(mesh), lengths);
    }
    for_each_neighbors(mesh, lengths, functor::Diff());
  }

  //---------------------------------------------------------------------------
  
  template < typename TVertexXsr, typename TCircularNeighborXsr, typename prec_type >
  class Mesh_curvature2
  {
  public: // typedefs
  
    //typedef typename precision<typename TVertexCollection::value_type>::type prec_type;
    //typedef typename TCircularNeighborIndices::value_type   Neighborhood;
    //typedef typename TVertexCollection::value_type          Vertex;

    typedef typename TCircularNeighborXsr::value_type  Neighborhood;
    typedef typename TCircularNeighborXsr::reference   NeighborhoodRef;
    //typedef typename TVertexXsr::value_type            Vertex;
    typedef typename TVertexXsr::reference             VertexRef;
    typedef typename TVertexXsr::index_type            index_type;

  public: // constructors

    Mesh_curvature2(TVertexXsr vertexXsr, TCircularNeighborXsr neighXsr)
      : m_vertexXsr(vertexXsr)
      , m_neighXsr(neighXsr)
    {}

  public: // set & get
  
    /// Return computed (signed) Gaussian curvature at vertex.
    prec_type gaussianCurvature() const { return m_orientedGaussianCurvature; }
    /// Return computed (signed) mean curvature at vertex.
    prec_type meanCurvature() const { return m_orientedMeanCurvature; }
    /// Return computed (signed) principal curvatures at vertex.
    /// the first in the pair is always the one with highest norm.
    std::pair<prec_type, prec_type> principalCurvatures() const { return m_principalCurvatures; }
    /// Return computed voronoi area at vertex.
    prec_type voronoiArea() const { return m_voronoiArea; }
    /// Return computed normal.
    const numeric_array<prec_type, 3> & normal() const { return m_normal; }
    
    /// Return unoriented Gaussian curvature at vertex.
    prec_type unorientedGaussianCurvature() const { return m_gaussianCurvature; }
    /// Return unoriented mean curvature at vertex.
    prec_type unorientedMeanCurvature() const { return m_meanCurvature; }
    

  public: // functions
  
    /// Computes all the good stuff at the i-th vertex.
    /// NB: nothing is returned here, use the appropriate get functions after calling compute.
    void process(index_type i)
    {
      NeighborhoodRef nh = m_neighXsr(i);
      // Curvature cannot be computed if vertex has not at least three neighbors
      assert(size(nh) >= 2);
      VertexRef p = m_vertexXsr(i);

      // initializations
      m_voronoiArea = 0;
      // mean curvature vector
      //Vector<prec_type, 3> m;
      std::fill(m_normal.begin(), m_normal.end(), prec_type(0));
      // Gaussian curvature
      m_gaussianCurvature = 0;
      
      numeric_array<prec_type, 3> e1, e2, e3;
      //numeric_array<prec_type, 3> e1, e2, e3;
      // Loop through all neighboring triangles
      typename Neighborhood::const_iterator iNh1 = nh.begin();
      const_cyclic_iterator<Neighborhood> iNh2 (nh, iNh1+1);      
      for (; iNh1 != nh.end(); ++iNh1, ++iNh2)
      {
        // compute the edges, their squared length and their length
        e1 = p - m_vertexXsr(*iNh1);
        prec_type n21 = norm2<prec_type>(e1);
        prec_type n1 = std::sqrt(n21);
        
        e2 = p - m_vertexXsr(*iNh2);
        prec_type n22 = norm2<prec_type>(e2);
        prec_type n2 = std::sqrt(n22);

        e3 = m_vertexXsr(*iNh2) - m_vertexXsr(*iNh1);
        prec_type n23 = norm2<prec_type>(e3);
        //prec_type n3 = std::sqrt(n23);
        //e3 *= 1/norm<prec_type>(e3);

        prec_type d1 = dot(e1, e3);
        prec_type d2 = dot(e2, e3);
        
        prec_type biarea = std::sqrt(n21*n23 - d1*d1);

        // Subtracting the angle between the two side edges to the gaussian curvature
        m_gaussianCurvature -= std::acos( max(prec_type(-1), min(prec_type(1), dot(e1, e2) / (n1*n2))));

        //cot1 =  d1 / norm<prec_type>(e1 - d1 * e3);
        //cot2 = -d2 / norm<prec_type>(e2 - d2 * e3);
        
        //cot1 =  d1 / dist<prec_type>(e1, d1 * e3);
        //cot2 = -d2 / dist<prec_type>(e2, d2 * e3);

        // NB: aren't the denominator all proportional to area? In that case
        // needed only once?
        /*
        prec_type cot1 = d1 / biarea;
        prec_type cot2 = -d2 / biarea;
        m += cot1 * e2 + cot2 * e1;
        */
        m_normal += (d1 * e2 - d2 * e1) * (prec_type(1)/biarea);
        
        if (dot(e1, e2) <= 0)
        {
          //n3 = dist<prec_type>(m_vertices[*iNh1], m_vertices[*iNh2]);
          //area += heronsFormula(n1, n2, n3) / 2;
          m_voronoiArea += biarea / 4;
        }
        else
        {
          if (d1 <= 0 || d2 >= 0)
          {
            //n3 = dist<prec_type>(m_vertices[*iNh1], m_vertices[*iNh2]);
            //area += heronsFormula(n1, n2, n3) / 4;
            m_voronoiArea += biarea / 8;
          }
          else
          {
            m_voronoiArea += 1.0/8.0 * (n22 * d1 -  n21 * d2 ) / biarea;
          }
        }
      }
      
      // finishing computation of gaussian curvature
      m_gaussianCurvature += 2*M_PI;
      m_gaussianCurvature /= m_voronoiArea;
      // finishing computation of mean curvature
      // NB: [Meyer et al.] the mean curvature is half the norm of K, hence the factor 4.
      m_meanCurvature = norm<prec_type>(m_normal) / (4*m_voronoiArea);
      // computes the principal curvatures from the unsigned mean and Gaussian curvatures
      prec_type delta = std::sqrt(max(prec_type(0), square(m_meanCurvature) - m_gaussianCurvature));
      m_principalCurvatures.first = m_meanCurvature + delta;
      m_principalCurvatures.second = m_meanCurvature - delta;

      // takes into account the sign
      // TODO: this test is actually not always accurate -- try to find a more robust criterion.
      if (dot(m_normal, cross(e1, e2)) < 0)
      {
        m_orientedMeanCurvature =  -m_meanCurvature;
        m_orientedGaussianCurvature = -m_gaussianCurvature;
        m_principalCurvatures.first = -m_principalCurvatures.first;
        m_principalCurvatures.second = -m_principalCurvatures.second;
        m_normal *= -1/norm(m_normal);
      }
      else
      {
        m_orientedMeanCurvature =  m_meanCurvature;
        m_orientedGaussianCurvature = m_gaussianCurvature;
        m_normal *= 1/norm(m_normal);
      }
    }
    
  private: // data, input
  
    TVertexXsr             m_vertexXsr;
    TCircularNeighborXsr   m_neighXsr;
    
  private: // data, output

    prec_type m_gaussianCurvature;
    prec_type m_meanCurvature;
    prec_type m_orientedGaussianCurvature;
    prec_type m_orientedMeanCurvature;
    prec_type m_voronoiArea;
    std::pair<prec_type, prec_type> m_principalCurvatures;
    numeric_array<prec_type, 3> m_normal;
  };

  //---------------------------------------------------------------------------

  template < typename TVertexCollection, typename TCircularNeighborIndices, typename prec_type >
  class Mesh_curvature
  {
  public: // typedefs
  
    //typedef typename precision<typename TVertexCollection::value_type>::type prec_type;
    typedef typename TCircularNeighborIndices::value_type   Neighborhood;
    typedef typename TVertexCollection::value_type          Vertex;

  public: // constructors

    Mesh_curvature(const TVertexCollection & vertices, const TCircularNeighborIndices & neighs)
    : m_vertices(vertices), m_neighs(neighs) {}

  public: // set & get
  
    /// Return computed (signed) Gaussian curvature at vertex.
    prec_type gaussianCurvature() const { return m_orientedGaussianCurvature; }
    /// Return computed (signed) mean curvature at vertex.
    prec_type meanCurvature() const { return m_orientedMeanCurvature; }
    /// Return computed (signed) principal curvatures at vertex.
    /// the first in the pair is always the one with highest norm.
    std::pair<prec_type, prec_type> principalCurvatures() const { return m_principalCurvatures; }
    /// Return computed voronoi area at vertex.
    prec_type voronoiArea() const { return m_voronoiArea; }
    /// Return computed normal.
    const numeric_array<prec_type, 3> & normal() const { return m_normal; }
    
    /// Return unoriented Gaussian curvature at vertex.
    prec_type unorientedGaussianCurvature() const { return m_gaussianCurvature; }
    /// Return unoriented mean curvature at vertex.
    prec_type unorientedMeanCurvature() const { return m_meanCurvature; }
    

  public: // functions
  
    /// Computes all the good stuff at the i-th vertex.
    /// NB: nothing is returned here, use the appropriate get functions after calling compute.
    void process(std::size_t i)
    {
      Neighborhood nh = m_neighs[i];
      // Curvature cannot be computed if vertex has not at least three neighbors
      assert(size(nh) >= 2);  
      const Vertex & p = m_vertices[i];

      // initializations
      m_voronoiArea = 0;
      // mean curvature vector
      //Vector<prec_type, 3> m;
      std::fill(m_normal.begin(), m_normal.end(), prec_type(0));
      // Gaussian curvature
      m_gaussianCurvature = 0;
      
      numeric_array<prec_type, 3> e1, e2, e3;
      //numeric_array<prec_type, 3> e1, e2, e3;
      // Loop through all neighboring triangles
      typename Neighborhood::const_iterator iNh1 = nh.begin();
      const_cyclic_iterator<Neighborhood> iNh2 (nh, iNh1+1);      
      for (; iNh1 != nh.end(); ++iNh1, ++iNh2)
      {
        // compute the edges, their squared length and their length
        e1 = p - m_vertices[*iNh1];
        prec_type n21 = norm2<prec_type>(e1);
        prec_type n1 = std::sqrt(n21);
        
        e2 = p - m_vertices[*iNh2];
        prec_type n22 = norm2<prec_type>(e2);
        prec_type n2 = std::sqrt(n22);

        e3 = m_vertices[*iNh2] - m_vertices[*iNh1];
        prec_type n23 = norm2<prec_type>(e3);
        //prec_type n3 = std::sqrt(n23);
        //e3 *= 1/norm<prec_type>(e3);

        prec_type d1 = dot(e1, e3);
        prec_type d2 = dot(e2, e3);
        
        prec_type biarea = std::sqrt(n21*n23 - d1*d1);

        // Subtracting the angle between the two side edges to the gaussian curvature
        m_gaussianCurvature -= std::acos( max(prec_type(-1), min(prec_type(1), dot(e1, e2) / (n1*n2))));

        //cot1 =  d1 / norm<prec_type>(e1 - d1 * e3);
        //cot2 = -d2 / norm<prec_type>(e2 - d2 * e3);
        
        //cot1 =  d1 / dist<prec_type>(e1, d1 * e3);
        //cot2 = -d2 / dist<prec_type>(e2, d2 * e3);

        // NB: aren't the denominator all proportional to area? In that case
        // needed only once?
        /*
        prec_type cot1 = d1 / biarea;
        prec_type cot2 = -d2 / biarea;
        m += cot1 * e2 + cot2 * e1;
        */
        m_normal += (d1 * e2 - d2 * e1) * (prec_type(1)/biarea);
        
        if (dot(e1, e2) <= 0)
        {
          //n3 = dist<prec_type>(m_vertices[*iNh1], m_vertices[*iNh2]);
          //area += heronsFormula(n1, n2, n3) / 2;
          m_voronoiArea += biarea / 4;
        }
        else
        {
          if (d1 <= 0 || d2 >= 0)
          {
            //n3 = dist<prec_type>(m_vertices[*iNh1], m_vertices[*iNh2]);
            //area += heronsFormula(n1, n2, n3) / 4;
            m_voronoiArea += biarea / 8;
          }
          else
          {
            m_voronoiArea += 1.0/8.0 * (n22 * d1 -  n21 * d2 ) / biarea;
          }
        }
      }
      
      // finishing computation of gaussian curvature
      m_gaussianCurvature += 2*M_PI;
      m_gaussianCurvature /= m_voronoiArea;
      // finishing computation of mean curvature
      // NB: [Meyer et al.] the mean curvature is half the norm of K, hence the factor 4.
      m_meanCurvature = norm<prec_type>(m_normal) / (4*m_voronoiArea);
      // computes the principal curvatures from the unsigned mean and Gaussian curvatures
      prec_type delta = std::sqrt(max(prec_type(0), square(m_meanCurvature) - m_gaussianCurvature));
      m_principalCurvatures.first = m_meanCurvature + delta;
      m_principalCurvatures.second = m_meanCurvature - delta;

      // takes into account the sign
      // TODO: this test is actually not always accurate -- try to find a more robust criterion.
      if (dot(m_normal, cross(e1, e2)) < 0)
      {
        m_orientedMeanCurvature =  -m_meanCurvature;
        m_orientedGaussianCurvature = -m_gaussianCurvature;
        m_principalCurvatures.first = -m_principalCurvatures.first;
        m_principalCurvatures.second = -m_principalCurvatures.second;
        m_normal *= -1/norm(m_normal);
      }
      else
      {
        m_orientedMeanCurvature =  m_meanCurvature;
        m_orientedGaussianCurvature = m_gaussianCurvature;
        m_normal *= 1/norm(m_normal);
      }
      
    }
    
  private: // data, input
  
    const TVertexCollection         &  m_vertices;
    const TCircularNeighborIndices  &  m_neighs;
    
  private: // data, output

    prec_type m_gaussianCurvature;
    prec_type m_meanCurvature;
    prec_type m_orientedGaussianCurvature;
    prec_type m_orientedMeanCurvature;
    prec_type m_voronoiArea;
    std::pair<prec_type, prec_type> m_principalCurvatures;
    numeric_array<prec_type, 3> m_normal;
  };

  //---------------------------------------------------------------------------

  namespace policy
  {

    //-------------------------------------------------------------------------

    /// Mesh policy for vertices indexed with an integer. This therefore assumes that the
    /// vertex collection is a vector of such.
    /// NB: TVertexCollection should offer random access.
    template < typename TCollection >
    class IntegerIndexing
    {
    public: // typedefs
      typedef std::size_t                       index_type;
      typedef typename TCollection::value_type  value_type;

    public: // constructors
      // NB: Actually, it can be quite comfortable to have it non explicit. We can pass collections rather
      // than constructing explicitely an accessor.
      IntegerIndexing(TCollection & data) : m_data(data) {}
      
    public: // static accessors
      typename TCollection::reference operator()(index_type i)
      { return m_data[i]; }
      typename TCollection::const_reference operator()(index_type i) const 
      { return m_data[i]; }
      
    private: // data
      TCollection & m_data;
    };
    
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
    
    //-------------------------------------------------------------------------
   
    /// Return value pointed at by an iterator.
    template < typename TIterator > //, typename Getter >
    struct IteratorIndexing
    {
    public: // typedefs
      typedef typename TIterator::value_type  value_type;
      typedef TIterator                       index_type;
    public: // static accessors
      value_type & operator()(index_type i) const
      { return *i; }
    };
  
    //-------------------------------------------------------------------------

    template < typename TIterator > //, typename Getter >
    //struct IteratorIndexing
    struct GraphIteratorIndexing_Position
    {
    public: // typedefs
      typedef typename TIterator::value_type::Vertex        value_type;
      typedef typename TIterator::value_type::VertexIndex   index_type;
    public: // static accessors
      // NB: it is not stupid to fix the reference at value_type, because the index_type, hence the iterator,
      // is also fixed inside the meshvertexnode class. 
      value_type & operator()(index_type i) const
      { return i->pos(); }
      //{ return Getter::get(i); }
    };
    
    //-------------------------------------------------------------------------
    
    template < typename TIterator >
    struct GraphIteratorIndexing_Neighbors
    {
    public: // typedefs
      typedef typename TIterator::value_type::VertexIndexCollection   value_type;
      typedef typename TIterator::value_type::VertexIndex             index_type;
    public: // static accessors
      value_type & operator()(index_type i) const
      { return i->neighbors(); }
    };
    
    //-------------------------------------------------------------------------
    
  } // namespace policy

  //---------------------------------------------------------------------------

  template < typename TVertexAccessPolicy, typename TCircularNeighborhoodAccessPolicy, typename prec_type >
  class MeshCurvature2
  {
  public: // typedefs
  
    //typedef typename precision<typename TVertexCollection::value_type>::type prec_type;
    typedef typename TCircularNeighborhoodAccessPolicy::value_type   Neighborhood;
    //typedef typename TVertexCollection::value_type          Vertex;
    typedef typename TVertexAccessPolicy::value_type        Vertex;
    typedef typename TVertexAccessPolicy::index_type        VertexIndex;
    //typedef typename TVertexAccessPolicy::collection_type   VertexCollection;

  public: // constructors

    MeshCurvature2(TVertexAccessPolicy vertexAccess, TCircularNeighborhoodAccessPolicy neighborAccess)
    : m_vertexAccess(vertexAccess), m_neighborAccess(neighborAccess) {}

  public: // set & get
  
    /// Return computed (signed) Gaussian curvature at vertex.
    prec_type gaussianCurvature() const { return m_orientedGaussianCurvature; }
    /// Return computed (signed) mean curvature at vertex.
    prec_type meanCurvature() const { return m_orientedMeanCurvature; }
    /// Return computed (signed) principal curvatures at vertex.
    /// the first in the pair is always the one with highest norm.
    std::pair<prec_type, prec_type> principalCurvatures() const { return m_principalCurvatures; }
    /// Return computed voronoi area at vertex.
    prec_type voronoiArea() const { return m_voronoiArea; }
    /// Return computed normal.
    const numeric_array<prec_type, 3> & normal() const { return m_normal; }
    
    /// Return unoriented Gaussian curvature at vertex.
    prec_type unorientedGaussianCurvature() const { return m_gaussianCurvature; }
    /// Return unoriented mean curvature at vertex.
    prec_type unorientedMeanCurvature() const { return m_meanCurvature; }
    

  public: // functions
  
    /// Computes all the good stuff at the i-th vertex.
    /// NB: nothing is returned here, use the appropriate get functions after calling compute.
    void process(VertexIndex i)
    {
      const Neighborhood & nh = m_neighborAccess(i);
      // Curvature cannot be computed if vertex has not at least three neighbors
      assert(size(nh) >= 2);  
      const Vertex & p = m_vertexAccess(i);

      // initializations
      m_voronoiArea = 0;
      // mean curvature vector
      numeric_array<prec_type, 3> m;
      // Gaussian curvature
      m_gaussianCurvature = 0;
      
      numeric_array<prec_type, 3> e1, e2, e3;
      // Loop through all neighboring triangles
      typename Neighborhood::const_iterator iNh1 = nh.begin();
      const_cyclic_iterator<Neighborhood> iNh2 (nh, nh.begin()); ++iNh2;
      for (; iNh1 != nh.end(); ++iNh1, ++iNh2)
      {
        // compute the edges, their squared length and their length
        e1 = p - m_vertexAccess(*iNh1);
        prec_type n21 = norm2<prec_type>(e1);
        prec_type n1 = std::sqrt(n21);
        
        e2 = p - m_vertexAccess(*iNh2);
        prec_type n22 = norm2<prec_type>(e2);
        prec_type n2 = std::sqrt(n22);

        e3 = m_vertexAccess(*iNh2) - m_vertexAccess(*iNh1);
        prec_type n23 = norm2<prec_type>(e3);
        //prec_type n3 = std::sqrt(n23);
        //e3 *= 1/norm<prec_type>(e3);

        prec_type d1 = dot(e1, e3);
        prec_type d2 = dot(e2, e3);
        
        prec_type biarea = std::sqrt(n21*n23 - d1*d1);

        // Subtracting the angle between the two side edges to the gaussian curvature
        m_gaussianCurvature -= std::acos( max(prec_type(-1), min(prec_type(1), dot(e1, e2) / (n1*n2))));

        //cot1 =  d1 / norm<prec_type>(e1 - d1 * e3);
        //cot2 = -d2 / norm<prec_type>(e2 - d2 * e3);
        
        //cot1 =  d1 / dist<prec_type>(e1, d1 * e3);
        //cot2 = -d2 / dist<prec_type>(e2, d2 * e3);

        // NB: aren't the denominator all proportional to area? In that case
        // needed only once?
        /*
        prec_type cot1 = d1 / biarea;
        prec_type cot2 = -d2 / biarea;
        m += cot1 * e2 + cot2 * e1;
        */
        m_normal += (d1 * e2 - d2 * e1) * (prec_type(1)/biarea);
        
        if (dot(e1, e2) <= 0)
        {
          //n3 = dist<prec_type>(m_vertices[*iNh1], m_vertices[*iNh2]);
          //area += heronsFormula(n1, n2, n3) / 2;
          m_voronoiArea += biarea / 4;
        }
        else
        {
          if (d1 <= 0 || d2 >= 0)
          {
            //n3 = dist<prec_type>(m_vertices[*iNh1], m_vertices[*iNh2]);
            //area += heronsFormula(n1, n2, n3) / 4;
            m_voronoiArea += biarea / 8;
          }
          else
          {
            m_voronoiArea += 1.0/8.0 * (n22 * d1 -  n21 * d2 ) / biarea;
          }
        }
      }
      
      // finishing computation of gaussian curvature
      m_gaussianCurvature += 2*M_PI;
      m_gaussianCurvature /= m_voronoiArea;
      // finishing computation of mean curvature
      // NB: [Meyer et al.] the mean curvature is half the norm of K, hence the factor 4.
      m_meanCurvature = norm<prec_type>(m_normal) / (4*m_voronoiArea);
      // computes the principal curvatures from the unsigned mean and Gaussian curvatures
      prec_type delta = std::sqrt(max(prec_type(0), square(m_meanCurvature) - m_gaussianCurvature));
      m_principalCurvatures.first = m_meanCurvature + delta;
      m_principalCurvatures.second = m_meanCurvature - delta;

      // takes into account the sign
      if (dot(m_normal, cross(e1, e2)) < 0)
      {
        m_orientedMeanCurvature =  -m_meanCurvature;
        m_orientedGaussianCurvature = -m_gaussianCurvature;
        m_principalCurvatures.first = -m_principalCurvatures.first;
        m_principalCurvatures.second = -m_principalCurvatures.second;
        m_normal *= -1/norm(m_normal);
      }
      else
      {
        m_orientedMeanCurvature =  m_meanCurvature;
        m_orientedGaussianCurvature = m_gaussianCurvature;
        m_normal *= 1/norm(m_normal);
      }
    }
    
  private: // data, input

    TVertexAccessPolicy                 m_vertexAccess;
    TCircularNeighborhoodAccessPolicy   m_neighborAccess;
    
  private: // data, output

    prec_type m_gaussianCurvature;
    prec_type m_meanCurvature;
    prec_type m_orientedGaussianCurvature;
    prec_type m_orientedMeanCurvature;
    prec_type m_voronoiArea;
    std::pair<prec_type, prec_type> m_principalCurvatures;
    numeric_array<prec_type, 3> m_normal;
  };

//-----------------------------------------------------------------------------

  //----------------------//
 //  GonzalezClustering  //
//----------------------//

/// Gonzalez clustering
template < typename TData, typename TPrec, typename TDist = SquaredEuclideanDist<TPrec> >
class GonzalezClustering
{
private: // typedefs
  typedef typename TData::const_iterator iterator;

public: // typedefs
  typedef TPrec                       prec_type;
  typedef typename TData::value_type  value_type;
  
public: // constructors
  
  GonzalezClustering(const TData & data)
    : m_begin(data.begin())
    , m_end(data.end())
    , m_dist()
    , m_quant(new std::vector<iterator>(1, data.begin()))
    , m_labels(new std::vector<std::size_t>(data.size(), 0))
  {}

public: // set & get

  /// Returns data quantization, i.e. a subset of the original data of representative of its clusters.
  shared_ptr<std::vector<iterator> > quantization() { return m_quant; }

  shared_ptr<std::vector<std::size_t> > labels() { return m_labels; }

public: // functions, main

  /// Clustering, until max cluster diameter is less than a threshold
  void clusterize_maxDiam(prec_type maxDiam)
  {
    unsigned int niter = 0;
    for (;;)
    {
      ++niter;
      //std::cout << "iter " << niter << std::endl;
      
      this->computeClusters();
      this->moveClusterSeeds();

      // Stop criterion
      std::size_t count = 0;
      bool flag = true;
      {
        for (iterator i = m_begin; i != m_end; ++i, ++count)
        {
          std::size_t label = (*m_labels)[count];
          if (m_dist(*i, *(*m_quant)[label]) > maxDiam)
          {
            flag = false;
            break;
          }
        }
      }
      if (flag) break;
      
      this->addNewClusterSeed();

      if (m_quant->size() != (niter+1))
      {
        std::cerr << "Warning: early termination of clustering" << std::endl;
        break;
      }
    }
  }
    
public: // functions, detail

  /// update labels to reflect current clustering.
  void computeClusters()
  {
    // for all points
    std::size_t count = 0;
    for (iterator i = m_begin; i != m_end; ++i, ++count)
    {
      // look for the quantized vector closest to point
      prec_type min_dist = std::numeric_limits<prec_type>::max();
      for (std::size_t j = 0; j < m_quant->size(); ++j)
      {
        prec_type tmp = m_dist(*i, *(*m_quant)[j]);
        if (tmp < min_dist)
        {
          min_dist = tmp;
          (*m_labels)[count] = j;
        }
      }
    }
  }
  
  /// Move cluster quantized vector as close to cluster center as possible
  void moveClusterSeeds()
  {
    std::vector<value_type> centers(m_quant->size(), value_type());
    {
      // for all points
      std::vector<unsigned int> counts(m_quant->size(), 0);
      std::size_t count = 0;
      for (iterator i = m_begin; i != m_end; ++i, ++count)
      {
        std::size_t label = (*m_labels)[count];
        // accumulate center points
        centers[label] += *i;
        ++counts[label];
      }
      // finish computation of centers
      for (std::size_t i = 0; i < centers.size(); ++i)
      {
        assert(counts[i] > 0);
        centers[i] *= prec_type(1)/counts[i];
        //(*m_quant)[i] = centers[i] / counts[i];
      }
    }
    std::vector<prec_type> dists(m_quant->size(), std::numeric_limits<prec_type>::max());
    // for all points
    std::size_t count = 0;
    for (iterator i = m_begin; i != m_end; ++i, ++count)
    {
      std::size_t label = (*m_labels)[count];
      // set point closest to its cluster centroid as cluster quantization
      prec_type tmp = m_dist(*i, centers[label]);
      if (tmp < dists[label])
      {
        dists[label] = tmp;
        (*m_quant)[label] = i;
      }
    }
  }
  
  /// Add a new cluster to the list
  void addNewClusterSeed()
  {
    prec_type max_dist = 0.0;
    iterator newSeed = m_begin;
    // for all points
    std::size_t count = 0;
    for (iterator i = m_begin; i != m_end; ++i, ++count)
    {
      std::size_t label = (*m_labels)[count];
      // look for the point with largest distance to its cluster quantized vector
      prec_type tmp = m_dist(*i, *(*m_quant)[label]);
      if (tmp > max_dist)
      {
        max_dist = tmp;
        newSeed = i;
      }
    }
    // Check that point can be added
    if (max_dist == 0.0)
    {
      std::cerr << "Seed cannot be added" << std::endl;
      return;
    }
#ifndef NDEBUG
    // Check that added vector quant is not already there
    {
      typename std::vector<iterator>::iterator i = std::find(m_quant->begin(), m_quant->end(), newSeed);
      assert( i == m_quant->end() );
    }
#endif
    // add point
    m_quant->push_back(newSeed);
  }
  
  
private: //data, input

  iterator m_begin;
  iterator m_end;
  TDist m_dist;
  
private: // data, output

  shared_ptr<std::vector<iterator> > m_quant;
  
private: // data, internals

  //std::vector<prec_type> m_dist;
  shared_ptr<std::vector<std::size_t> > m_labels;
};



/// Gonzalez clustering.
template < typename TCollection, typename TPrec >
inline shared_ptr<std::vector<typename TCollection::const_iterator> >
gonzalez_clustering(const TCollection & data, TPrec maxDiam)
{
  GonzalezClustering<TCollection, TPrec> gc(data);
  gc.clusterize_maxDiam(maxDiam);
  return gc.quantization();
}


/*
template < typename TVertexCollection, typename TFaceCollection >
void remove_vertex(typename TVertexCollection::iterator i, const TVertexCollection & v, const TFaceCollection & f)
{
  for (std::list<std::list<MeshFaceNode>::iterator>::iterator j = i->faces.begin(); j != i->faces.end(); ++j)
  {
    // Remove index to this face
    for (std::size_t k = 0; k < 3; ++k)
    {
      if ((*j)->face[k] == i) continue;
      std::size_t tmp = til::size((*j)->face[k]->faces);
      (*j)->face[k]->faces.remove(*j);
    }
    // remove face itself
    graph_faces.erase(*j);
  }
  // remove point in neighbor's list
  for (std::list<std::list<MeshVertexNode>::iterator>::iterator j = i->neighbors.begin(); j != i->neighbors.end(); ++j)
  {
    (*j)->neighbors.remove(i);
  }
  // remove point
  i = graph_vertices.erase(i);
  //++i;
}
*/

/// Simple flattening of a point and its neighborhood.
/// The neighborhood is assumed to be given in (counter) clockwise order.
/// This flattening tries to be as exact as possible from the point of view of the center point: distances to neighbors are preserved,
/// and proportions between angles between two neighbors are preserved.
  /*
  class SimpleNeighborFlattening
  {
  public: // typedef
    typedef std::vector<Point<float,2> > Res;
  public: // constructors
    SimpleNeighborFlattening() : m_res() {}
  public: // set & get
    const Res & res() const { return m_res; }
    Res & res() { return m_res; }
  public: // operators
    /// Computes the flattened neighbor starting from a point and its (counter-)clockwise neighborhood.
    void operator()(
      const Point<float,3>                                & point,
      const std::vector<std::list<MeshVertexNode>::iterator>   & neighbors
    )
    {
      m_angles.clear();
      m_norms.clear();
      m_angles.reserve(til::size(neighbors));
      m_norms.reserve(til::size(neighbors));
      
      std::vector<std::list<MeshVertexNode>::iterator>::const_iterator n = neighbors.begin();
      const_cyclic_iterator<std::vector<std::list<MeshVertexNode>::iterator> > n2(neighbors, ++neighbors.begin());
      for (; n != neighbors.end(); ++n, ++n2)
      {
        m_norms.push_back(norm<double>((*n)->pos-point));
        m_angles.push_back(std::acos(dot((*n2)->pos-point, (*n)->pos-point) / (norm<double>((*n2)->pos-point) * norm<double>((*n)->pos-point))));
      }
      std::partial_sum(m_angles.begin(), m_angles.end(), m_angles.begin());
      std::transform(m_angles.begin(), m_angles.end(), m_angles.begin(), std::bind2nd(std::multiplies<double>(), 2*M_PI/angles.back()));
      
      m_res.clear();
      //m_res.resize(til::size(neighbors));
      m_res.reserve(til::size(neighbors));
      for (std::size_t i = 0; i < size(res); ++i)
      {
        m_res.push_back(Point<float,2>(std::cos(angles[i]) * norms[i], std::sin(angles[i]) * norms[i]));
      }      
    }
  private: // data
    Res m_res;
  private: // computation variables
    std::vector<double> m_angles;
    std::vector<double> m_norms;
  };

  */

  //---------------------------------------------------------------------------
  
  template < typename TPrec, typename TVertexIndex, typename TVertexCollection, typename TNeighborhood, typename TVertex >
  inline TPrec
  dist2_surface(const TVertex & p, TVertexIndex i, const TVertexCollection & vertices, const TNeighborhood & neighc)
  {
    typedef typename TVertex::value_type prec_type;
    prec_type dmin = std::numeric_limits<prec_type>::max();
    numeric_array<TPrec, 3> tmp;
    typename TNeighborhood::const_iterator ineigh = neighc.begin();
    const_cyclic_iterator<TNeighborhood> ineigh2(neighc, neighc.begin());
    ++ineigh2;
    for (; ineigh != neighc.end(); ++ineigh, ++ineigh2)
    {
      geo::Project::PointOnTriangle3D<typename TVertexCollection::value_type, numeric_array<TPrec, 3> >()(p, vertices[i], vertices[*ineigh], vertices[*ineigh2], tmp);
      prec_type d = dist2(p, tmp, prec<prec_type>());
      if (d < dmin)
      {
        dmin = d;
      }
    }
    return dmin;
  }

  //---------------------------------------------------------------------------

  /// Returns the index of the neighborhood that is the closest.
  template < typename TPrec, typename TVertexIndex, typename TVertexCollection, typename TNeighborhood, typename TVertex >
  inline numeric_array<double,3>
  closest_normal(const TVertex & p, TVertexIndex i, const TVertexCollection & vertices, const TNeighborhood & neighc)
  {
    typedef typename TVertex::value_type prec_type;
    prec_type dmin = std::numeric_limits<prec_type>::max();
    numeric_array<TPrec, 3> tmp;
    typename TNeighborhood::const_iterator ineigh = neighc.begin();
    const_cyclic_iterator<TNeighborhood> ineigh2(neighc, neighc.begin());
    numeric_array<double, 3> res;
    ++ineigh2;
    for (; ineigh != neighc.end(); ++ineigh, ++ineigh2)
    {
      geo::Project::PointOnTriangle3D<typename TVertexCollection::value_type, numeric_array<TPrec, 3> >()(p, vertices[i], vertices[*ineigh], vertices[*ineigh2], tmp);
      prec_type d = dist2(p, tmp, prec<prec_type>());
      if (d < dmin)
      {
        geo::triangle_normal(vertices[i], vertices[*ineigh], vertices[*ineigh2], res);
        dmin = d;
      }
      else if (d == dmin)
      {
        numeric_array<double, 3> tmp;    
        geo::triangle_normal(vertices[i], vertices[*ineigh], vertices[*ineigh2], tmp);
        res += tmp;
      }
    }
    return res;
  }
  
  // --------------------------------------------------------------------------
  
  template < typename TVertexAccessPolicy, typename TNeighborhoodAccessPolicy, typename TVertexAccessPolicy2 >
  class ScaleIndependantLaplacianSmoothing
  {
  public: // typedefs
  
    typedef typename TNeighborhoodAccessPolicy::value_type  Neighborhood;
    typedef typename TVertexAccessPolicy::value_type        Vertex;
    typedef typename TVertexAccessPolicy::index_type        VertexIndex;
    typedef typename TVertexAccessPolicy2::index_type       VertexIndex2;
    typedef typename Vertex::value_type  prec_type;

  public: // constructors

    ScaleIndependantLaplacianSmoothing
    (
      TVertexAccessPolicy vertexAccess, 
      TNeighborhoodAccessPolicy neighborAccess,
      TVertexAccessPolicy2 vertexAccess2,
      prec_type lambda
    )
      : m_vertexAccess(vertexAccess)
      , m_neighborAccess(neighborAccess)
      , m_vertexAccess2(vertexAccess2)
      , m_lambda(lambda)
    {}

  public: // functions

    void operator()(VertexIndex begin, VertexIndex end, VertexIndex2 begin2)
    {
      for (; begin != end; ++begin, ++begin2)
      {
        m_buffer[0] = m_buffer[1] = m_buffer[2] = 0;
        prec_type E = 0;
        m_vertexAccess2(begin2) = m_vertexAccess(begin);
        //prec_type coeff = m_lambda / m_neighborAccess(begin).size();
        for (typename Neighborhood::const_iterator iNeigh = m_neighborAccess(begin).begin(); iNeigh != m_neighborAccess(begin).end(); ++iNeigh)
        {
          m_tmp = m_vertexAccess(*iNeigh) - m_vertexAccess(begin);
          prec_type d = norm(m_tmp); //dist(m_vertexAccess(*iNeigh), m_vertexAccess(begin));
          E += d;
          m_tmp *= 1/d;
          m_buffer += m_tmp;
          //m_vertexAccess2(begin2) += (m_vertexAccess(*iNeigh) - m_vertexAccess2(begin2));
        }
        m_vertexAccess2(begin2) += m_lambda * (2/E) * m_buffer;
      }
    }
      
  private: // data, input
  
    TVertexAccessPolicy         m_vertexAccess;
    TNeighborhoodAccessPolicy   m_neighborAccess;
    TVertexAccessPolicy2        m_vertexAccess2;
    prec_type                   m_lambda;
    
  private: // data, internals

    numeric_array<prec_type, 3> m_buffer;
    numeric_array<prec_type, 3> m_tmp;
  };

  //---------------------------------------------------------------------------
  
    //----------------------//
   //  LaplacianSmoothing  //
  //----------------------//
  
  /// Discrete Laplacian smoothing.
  /// The neat thing about this class is that it can be used to smooth both a mesh and data attached to it.
  template < typename TInputAccessPolicy, typename TOutputAccessPolicy, typename TNeighborhoodAccessPolicy >
  class LaplacianSmoothing
  {
  public: // typedefs
  
    typedef typename TNeighborhoodAccessPolicy::value_type  Neighborhood;
    typedef typename TInputAccessPolicy::value_type         Data;
    typedef typename TInputAccessPolicy::index_type         InputIndex;
    typedef typename TOutputAccessPolicy::index_type        OutputIndex;
    typedef typename precision<Data>::type                  prec_type;

  public: // constructors

    LaplacianSmoothing
    (
      TInputAccessPolicy  inputAccess, 
      TOutputAccessPolicy outputAccess,
      TNeighborhoodAccessPolicy neighborAccess,
      prec_type lambda
    )
      : m_inputAccess(inputAccess)
      , m_outputAccess(outputAccess)
      , m_neighborAccess(neighborAccess)
      , m_lambda(lambda)
    {}

  public: // functions

    void operator()(InputIndex begin, InputIndex end, OutputIndex begin2)
    {
      for (; begin != end; ++begin, ++begin2)
      {
        // Check that now to avoid division by zero later 
        assert(m_neighborAccess(begin).size() > 0);

        m_outputAccess(begin2) = m_inputAccess(begin) * (1 - m_lambda);
        prec_type coeff = m_lambda / m_neighborAccess(begin).size();
        typename Neighborhood::const_iterator iNeigh = m_neighborAccess(begin).begin();
        for (; iNeigh != m_neighborAccess(begin).end(); ++iNeigh)
        {
          m_outputAccess(begin2) += coeff * m_inputAccess(*iNeigh);

          // NB: this could be in a policy
          if (is_nan(m_outputAccess(begin2)))
          {
            std::cout << "(LS:is-nan)";
          }
        }
      }
    }
    
  private: // data, input
  
    TInputAccessPolicy          m_inputAccess;
    TOutputAccessPolicy         m_outputAccess;
    TNeighborhoodAccessPolicy   m_neighborAccess;
    prec_type m_lambda;
  };


  /// Laplacian smoothing -- helper function
  template < typename TVertex >
  inline void laplacian_smoothing
  (
    std::vector<TVertex> & vertices,
    std::vector<std::vector<std::size_t> > & neighbors,
    unsigned int nstep,
    double coeff
  )
  {
    typedef std::vector<TVertex>                                                        TVertexCollection;
    typedef policy::IntegerIndexing<TVertexCollection>                             VertexIndexing;
    typedef policy::IntegerIndexing<std::vector<std::vector<std::size_t> > >       NeighborIndexing;
    typedef LaplacianSmoothing<VertexIndexing, VertexIndexing, NeighborIndexing>   Smoother;

    std::size_t n = vertices.size();

    std::vector<TVertex> verticesBuff(n);
    VertexIndexing indVertex(vertices);
    VertexIndexing indVertexBuff(verticesBuff);
    NeighborIndexing indNeigh(neighbors);
    Smoother smoothToBuffer(indVertex, indVertexBuff, indNeigh, coeff);
    Smoother smoothFromBuffer(indVertexBuff, indVertex, indNeigh, coeff);
    for (unsigned int i = 0; i < nstep; ++i)
    {
      smoothToBuffer(0, n, 0);
      smoothFromBuffer(0, n, 0);
    }
  }

  //---------------------------------------------------------------------------

  
  template < typename TVertexAccessPolicy, typename TNeighborhoodAccessPolicy > //, typename TPrec >
  class LambdaMuSmoothing
  {
  public: // typedefs ---------------------------------------------------------
    
    //typedef TPrec prec_type;
    typedef typename TNeighborhoodAccessPolicy::value_type  Neighborhood;
    typedef typename TVertexAccessPolicy::value_type        Vertex;
    typedef typename TVertexAccessPolicy::index_type        VertexIndex;
    typedef typename Vertex::value_type  prec_type;
    
  public: // constructors -----------------------------------------------------
    
    LambdaMuSmoothing
    (
      TVertexAccessPolicy vertexAccess,
      TNeighborhoodAccessPolicy neighborAccess,
      const prec_type lambda,
      const prec_type mu,
      int bufferSize
    )
      : m_vertexAccess(vertexAccess)
      , m_neighborAccess(neighborAccess)
      , m_lambda(lambda)
      , m_mu(mu)
    {
      this->setBufferSize(bufferSize);
    }
    
    
  public: // set & get --------------------------------------------------------
  
    prec_type & lambda() { return m_lambda; }
    prec_type & mu() { return m_mu; }
    void setBufferSize(int bufferSize) { m_vertices.resize(bufferSize); }
        
  public: // functions --------------------------------------------------------

    void operator()(VertexIndex begin, VertexIndex end)
    {
      this->lambda_step(begin, end);
      this->mu_step(begin, end);
    }
    
  private: // functions -------------------------------------------------------
  
    void lambda_step(VertexIndex begin, VertexIndex end)
    {
      for (std::size_t i = 0; begin != end; ++begin, ++i)
      {
        m_vertices[i] = m_vertexAccess(begin) * (1 - m_lambda);
        prec_type coeff = m_lambda / m_neighborAccess(begin).size();
        typename Neighborhood::const_iterator iNeigh = m_neighborAccess(begin).begin();
        for (; iNeigh != m_neighborAccess(begin).end(); ++iNeigh)
        {
          m_vertices[i] += coeff * m_vertexAccess(*iNeigh);
        }
      }
    }
    
    void mu_step(VertexIndex begin, VertexIndex end)
    {
      for (std::size_t i = 0; begin != end; ++begin, ++i)
      {
        m_vertexAccess(begin) = m_vertices[i] * (1 + m_mu);
        prec_type coeff = m_mu / m_neighborAccess(begin).size();
        typename Neighborhood::const_iterator iNeigh = m_neighborAccess(begin).begin();
        for (; iNeigh != m_neighborAccess(begin).end(); ++iNeigh)
        {
          m_vertexAccess(begin) -= coeff * m_vertexAccess(*iNeigh);
        }
      }
    }
    
  private: // data, input -----------------------------------------------------

    TVertexAccessPolicy         m_vertexAccess;
    TNeighborhoodAccessPolicy   m_neighborAccess;
    prec_type                   m_lambda;
    prec_type                   m_mu;
    
  private: // data, internals -------------------------------------------------
  
    std::vector<Vertex> m_vertices;    
  };

  //----------------------------------------------------------------------------

  /// Computes the n-neighborhood, i.e. collection of neighbors that are at most n jumps away.
  template < typename NeighborCollection, typename NeighborCollectionOut >
  void
  get_n_neighborhood
  (
    NeighborCollection const & neigh      ///< [input] The 1-neighborhoods
  , NeighborCollectionOut & neigh_n_out   ///< [output] The computed n-neighborhood 
  , unsigned int n                        ///< [input] the desired neighborhood width
  );

  //---------------------------------------------------------------------------

  /// Convert a set of neighborhoods into a set of edges.
  /// If the neighborhood is symmetric, the size of the edge set should be half of those
  /// of the neighbors.
  template < typename T >
  void
  neighbors2edges
  (
    std::vector<std::vector<T> >  const & neighbors
  , std::vector<std::pair<T, T> >       & edges
  );

  //---------------------------------------------------------------------------

  template < typename T >
  class Neighboring_faces
  {
  public: // exceptions
    struct InconsistentArguments : public std::exception {};
  public: // operator
    /// Computes the index of the neighboring faces.
    /// The neighboring faces follow the same order than the neighbors, i.e.
    /// the first face corresponds to the first two neighbors, etc.
    void operator()
    (
      std::vector<std::vector<T> > const & cneighs              ///< [input] Circular neighborhoods
    , std::vector<til::numeric_array<T, 3> > const & faces      ///< [input] Faces
    , std::vector<std::vector<T> > & neighborfaces              ///< [output] neighbor faces
    );
  };

  //---------------------------------------------------------------------------

  /// Computes the index of the neighboring faces.
  /// The neighboring faces follow the same order than the neighbors, i.e.
  /// the first face corresponds to the first two neighbors, etc.
  template < typename T >
  inline void
  neighboring_faces
  (
    std::vector<std::vector<T> > const & cneighs            ///< [input] Circular neighborhoods
  , std::vector<til::numeric_array<T, 3> > const & faces    ///< [input] Faces
  , std::vector<std::vector<T> > & neighborfaces            ///< [output] neighbor faces
  )
  {
    Neighboring_faces<T>()(cneighs, faces, neighborfaces);
  }

  //---------------------------------------------------------------------------

} // namespace til


// package include
#include "mesh_utils.tpp"

#endif //_MESHUTILS_H_
