#ifndef TIL_MESH_DECIMATION
#define TIL_MESH_DECIMATION

// includes from STL

// includes from TIL
#include "graph.h"
#include "mesh_conversion.h"
#include "meshUtils.h"

namespace til
{

  template < typename TVertexXsr, typename TNeighborhoodXsr, typename TPrec >
  struct MeshWaveletEnergy2
  {
  public: // typedefs
  
    typedef typename TVertexXsr::index_type        index_type;
    typedef typename TVertexXsr::value_type        Vertex;
    typedef typename TNeighborhoodXsr::value_type  Neighborhood;
    typedef typename TNeighborhoodXsr::reference   NeighborhoodRef;

  public: // constructors
    MeshWaveletEnergy2(const TVertexXsr & vertexXsr, const TNeighborhoodXsr & neighXsr)
      : m_neighXsr(neighXsr)
      , m_vertexXsr(vertexXsr)
      , m_mc(vertexXsr, neighXsr)
    {}

  public: // functions
  
    TPrec operator()(index_type i)
    {
      m_mc.process(i);
      Vertex c;
      NeighborhoodRef neigh = m_neighXsr(i);
      for (typename Neighborhood::const_iterator j = neigh.begin(); j != neigh.end(); ++j)
      //for (std::size_t j = 0; j < m_neighc[i].size(); ++j)
      {
        //std::cout << m_vertices[m_neighc[i][j]] << " ";
        c += m_vertexXsr(*j);
        //c += m_vertices[m_neighc[i][j]];
      }
      //std::cout << std::endl;
      c *= TPrec(1.0) / neigh.size();
      //std::cout << "mean: " << c << std::endl;
      //centroid(m_neighc[i], c);
      //return dist<TPrec>(m_vertices[i], centroid<Vertex>(m_neighc[i])) / m_mc.voronoiArea();
      //m_error = m_vertices[i] - c
      return dist(m_vertexXsr(i), c, prec<TPrec>()) / m_mc.voronoiArea();
    }
  
  private: // data, input
  
    TNeighborhoodXsr                                       m_neighXsr;
    TVertexXsr                                             m_vertexXsr;

  private: // data, internals

    Mesh_curvature2<TVertexXsr, TNeighborhoodXsr, TPrec>    m_mc;
  };

  //---------------------------------------------------------------------------


  /// Distance of a vertex to the mean of its neighbors divided by the area of
  /// its Voronoi cell.
  template < typename TVertexCollection, typename TNeighborhoods, typename TPrec >
  struct MeshWaveletEnergy
  {
  public: // typedefs
    typedef typename TVertexCollection::value_type Vertex;
  public: // constructors
    MeshWaveletEnergy
    (
      TVertexCollection const & vertices
    , TNeighborhoods const & neighc
    )
      : m_vertices(vertices)
      , m_neighc(neighc)
      , m_mc(vertices, neighc)
    {}
      
  public: // functions
  
    TPrec operator()(std::size_t i)
    {
      m_mc.process(i);
      Vertex c;
      for (std::size_t j = 0; j < m_neighc[i].size(); ++j)
      {
        c += m_vertices[m_neighc[i][j]];
      }
      c *= TPrec(1.0) / m_neighc[i].size();
      //centroid(m_neighc[i], c);
      //return dist<TPrec>(m_vertices[i], centroid<Vertex>(m_neighc[i])) / m_mc.voronoiArea();
      //m_error = m_vertices[i] - c
      return dist(m_vertices[i], c, prec<TPrec>()) / m_mc.voronoiArea();
    }

  private: // data, input
    const TVertexCollection & m_vertices;
    const TNeighborhoods & m_neighc;  
  private: // data, internal
    Mesh_curvature<TVertexCollection, TNeighborhoods, TPrec > m_mc;
  };
  
  //---------------------------------------------------------------------------

  namespace
  {
    
    
    /// Compare the address of the pointee of two iterators.
    /*
    template < typename TIterator >
    struct ItComp
    {
      bool operator()
      (
       //const TIterator & i,
       //const TIterator & j
       TIterator i,
       TIterator j
      )
      {
        return &*i < &*j;
      }
    };
    */
    
    template < typename T >
    struct GroComp
    {
      //typedef std::list<MeshVertexNode>::iterator T;
    
      GroComp(T a, T b, T c)
        : m_a(a)
        , m_b(b)
        , m_c(c)
      {}
      
      boost::array<T, 3>
      helper(T a, T b, T c)
      {
        boost::array<T,3> tmp;
        tmp[0] = a;
        tmp[1] = b;
        tmp[2] = c;
        return tmp;
      }
      
      template < typename TFace >
      bool operator()(const TFace & f)
      {
        return
          f->face == helper(m_a, m_b, m_c) ||
          f->face == helper(m_c, m_a, m_b) ||
          f->face == helper(m_b, m_c, m_a) ||
          f->face == helper(m_b, m_a, m_c) ||
          f->face == helper(m_c, m_b, m_a) ||
          f->face == helper(m_a, m_c, m_b)
          ;
      }
      
      T m_a;
      T m_b;
      T m_c;
    };

  }

  //---------------------------------------------------------------------------

  namespace
  {
    struct Temp
    {
      double angle;
      std::size_t index;
      friend bool operator<(const Temp & x, const Temp & y) { return x.angle < y.angle; }
    };
    struct TempComp
    {
      bool operator()(const Temp & x, const Temp & y) { return x.angle < y.angle; }
    };
  }

  //---------------------------------------------------------------------------

  /// Given a list of 2D points, compute a bounding equilateral triangle.
  /// The three resulting points are writen starting from the out iterator.
  // TODO: move to geo
  template < typename TPointIterator, typename TOutIterator >
  TOutIterator
  bounding_triangle(TPointIterator pbegin, TPointIterator pend, TOutIterator out)
  {
    float miny, k1, k2;
    k1 = k2 = -std::numeric_limits<float>::max();
    miny = std::numeric_limits<float>::max();
    for (; pbegin != pend; ++pbegin)
    {
      miny = min(miny, (*pbegin)[1]);
      max_helper(k1, (*pbegin)[1] - std::sqrt(3.0f) * (*pbegin)[0]);
      max_helper(k2, (*pbegin)[1] + std::sqrt(3.0f) * (*pbegin)[0]);
    }
    
    til::numeric_array<float,2> M;
    M[0] = (k2 - k1) / (2*std::sqrt(3.0f));
    M[1] = 2*miny/3 + (k1+k2)/6;
      
    til::numeric_array<float,2> A;
    A[0] = (miny - k1) * (1/std::sqrt(3.0f));
    A[1] = miny;
    til::numeric_array<float,2> B;
    B[0] = (k2 - miny) * (1/std::sqrt(3.0f));
    B[1] = miny;
    til::numeric_array<float,2> C;
    C[0] = (k2 - k1) / (2*std::sqrt(3.0f));
    C[1] = (k1+k2)/2;

    A += 2.0f*(A-M);
    B += 2.0f*(B-M);
    C += 2.0f*(C-M);
    
    *(out++) = A;
    *(out++) = B;
    *(out++) = C;
    return out;
  }

  //---------------------------------------------------------------------------


    //-------------------------------//
   //  SimpleDelaunayTriangulation  //
  //-------------------------------//

  /*
  template < typename TFace >
  class SimpleDelaunayTriangulation
  {
  public: // operators
    std::list<boost::array<std::size_t,3> >
    operator()(std::vector<numeric_array<float, 2> > points);
  };
  */
  
  inline
  std::list<boost::array<std::size_t,3> >
  simple_delaunay_triangulation(std::vector<numeric_array<float, 2> > points)
  {
    typedef boost::array<std::size_t,3> Face;
    typedef std::list<Face> FaceList;
    
    // Get size now before points are added
    std::size_t n = points.size();
    
    // Don't waste time calling this for less than 4 points.
    assert(points.size() >= 4);
    
    FaceList faces;
  
    // compute bounding triangle
    //bounding_triangle(points.begin(), points.end(), std::back_inserter(points));
    
    {
      float miny, k1, k2;
      k1 = k2 = -std::numeric_limits<float>::max();
      miny = std::numeric_limits<float>::max();
      for (std::size_t i = 0; i < points.size(); ++i)
      {
        miny = min(miny, points[i][1]);
        k1 = max(k1, points[i][1] - std::sqrt(3.0f)*points[i][0]);
        k2 = max(k2, points[i][1] + std::sqrt(3.0f)*points[i][0]);
      }
      
      numeric_array<float,2> M;
      M[0] = (k2 - k1) / (2*std::sqrt(3.0f));
      M[1] = 2*miny/3 + (k1+k2)/6;
        
      numeric_array<float,2> A;
      A[0] = (miny - k1) * (1/std::sqrt(3.0f));
      A[1] = miny;
      numeric_array<float,2> B;
      B[0] = (k2 - miny) * (1/std::sqrt(3.0f));
      B[1] = miny;
      numeric_array<float,2> C;
      C[0] = (k2 - k1) / (2*std::sqrt(3.0f));
      C[1] = (k1+k2)/2;
  
      A += 2.0f*(A-M);
      B += 2.0f*(B-M);
      C += 2.0f*(C-M);
      
      points.push_back(A);
      points.push_back(B);
      points.push_back(C);
      
      boost::array<std::size_t,3> face = { { points.size()-1, points.size()-2, points.size()-3 } };
      faces.push_back(face);
    }
    
  
    //Face face = {0,1,2};
    //faces.push_back(face);
    //for (std::size_t i = 3; i < size(points); ++i)
    for (std::size_t i = 0; i < n; ++i)
    {
      //std::cout << "another point :" << points[i] << std::endl;
  
      // Remove triangles if we lie in their circumcircle, and collect their vertex.
      std::set<std::size_t> neighbors;
      //e.insert(0);
      //e.insert(i-1);
      {
        FaceList::iterator f = faces.begin();
        while (f != faces.end())
        {
          //std::cout << "face " << points[(*f)[0]] << " " <<  points[(*f)[1]] << " " <<  points[(*f)[2]] << ":";
          if (geo::is_in_circumcircle<double>(points[i], points[(*f)[0]], points[(*f)[1]], points[(*f)[2]]))
          {
            //std::cout << "inside" << std::endl;
            neighbors.insert((*f)[0]);
            neighbors.insert((*f)[1]);
            neighbors.insert((*f)[2]);
            f = faces.erase(f);
          }
          else
          {
            //std::cout << "outside circle" << std::endl;
            ++f;
          }
        }
      }
      
      assert(neighbors.size() >= 3);
      
      // order neighbors according to the angle they make with the current point
      
      std::vector<Temp> tmp;
      tmp.reserve(neighbors.size());
      for (std::set<std::size_t>::iterator iNeighbor = neighbors.begin(); iNeighbor != neighbors.end(); ++iNeighbor)
      {
        Temp t;
        t.index = *iNeighbor;
        numeric_array<float,2> v = points[*iNeighbor]-points[i];
        t.angle = geo::angle(v[0], v[1]);
        tmp.push_back(t);
      }
      //std::sort(tmp.begin(), tmp.end(), TempComp());
      std::sort(tmp.begin(), tmp.end());
      
      // Form new faces
      for (std::size_t it = 0; it < tmp.size(); ++it)
      {
        //std::cout << "Adding face " << i << " " << tmp[it].index << " " << tmp[ (it+1) % size(tmp) ].index << std::endl;
        Face face = { {i, tmp[it].index, tmp[ (it+1) % tmp.size() ].index } };
        faces.push_back(face);
      }
    }
    /*
      {
        std::cout << "faces" << std::endl;
        for (FaceList::iterator f = faces.begin(); f != faces.end(); ++f)
        {
          std::cout << (*f)[0] << " " << (*f)[1] << " " << (*f)[2] << std::endl;
        }
        std::cout << "endfaces" << std::endl;
      }
  */
  
    // Remove extra faces
    {
      FaceList::iterator f = faces.begin();
      while (f != faces.end())
      {
        // Remove faces with supertriangle vertex
        if ((*f)[0] >= n ||
            (*f)[1] >= n ||
            (*f)[2] >= n)
        {
          f = faces.erase(f);
        }
        else
        {
          ++f;
        }
      }
    }
  /*
      {
        FaceList::iterator f = faces.begin();
        std::cout << "cross2D: " << cross(points[(*f)[1]] - points[(*f)[0]], points[(*f)[2]] - points[(*f)[0]]) << std::endl;
      }
  */
  
    // Check that normals are ok
    {
      for (FaceList::iterator f = faces.begin(); f != faces.end(); ++f)
      {
        if (cross(points[(*f)[1]] - points[(*f)[0]], points[(*f)[2]] - points[(*f)[0]]) < 0)
        {
          std::cout << "wS!";
        }
      }
    }
  
  
    // Check if we have to resort to normal testing to remove faces
    if (faces.size() != n-2)
    {
      /*
      std::cout << "NI!" << faces.size() << "-" << n << std::endl;
      {
        std::cout << "faces2" << std::endl;
        for (FaceList::iterator f = faces.begin(); f != faces.end(); ++f)
        {
          std::cout << (*f)[0] << " " << (*f)[1] << " " << (*f)[2] << std::endl;
        }
        std::cout << "endfaces" << std::endl;
      }
      */
          
      // Remove faces with bad normal
      /*
      {
        FaceList::iterator f = faces.begin();
        while (f != faces.end())
        {
          if (cross(points[(*f)[1]] - points[(*f)[0]], points[(*f)[2]] - points[(*f)[0]]) < 0)
          {
            f = faces.erase(f);
          }
          else
          {
            ++f;
          }
        }
      }
      */
  
      {
        FaceList::iterator f = faces.begin();
        while (f != faces.end())
        {
          bool flagdel = false;
          int count = 0;
          for (int i = 0; i < 3; ++i)
          {
            if ((*f)[i] > (*f)[(i+1)%3])
            {
              if (++count > 1)
              {
                flagdel = true;
                //f = faces.erase(f);
                break;
              }
            }
          }
          if (flagdel) f = faces.erase(f);
          else ++f;
        }
      }
      
      // Check again that we have the expected number of faces
      if (faces.size() != n-2)
      {
        std::cout << "NI2!" << faces.size() << "-" << n << std::endl;
        {
          std::cout << "faces2" << std::endl;
          for (FaceList::iterator f = faces.begin(); f != faces.end(); ++f)
          {
            std::cout << (*f)[0] << " " << (*f)[1] << " " << (*f)[2] << std::endl;
          }
          std::cout << "endfaces" << std::endl;
          std::cout << "points" << std::endl;
          for (std::size_t i = 0; i < n; ++i) std::cout << points[i] << std::endl;
          std::cout << "end points" << std::endl;
        }
      }
    }
    
    // check that all outer edges is missing
    {
      std::set<boost::array<std::size_t,2> > outer_edges;
      for (FaceList::iterator f = faces.begin(); f != faces.end(); ++f)
      {
        for (std::size_t i = 0; i < 3; ++i)
        {
          std::size_t delta = std::abs(int((*f)[i]) - int((*f)[(i+1)%3]));
          if (delta == 1 || delta == (n-1))
          {
            boost::array<std::size_t,2> tmp = { {min((*f)[i], (*f)[(i+1)%3]), max((*f)[i], (*f)[(i+1)%3])} };
            outer_edges.insert(tmp);
          }
        }
      }
      if (outer_edges.size() != n)
      {
        std::cout << "wU!";
      }
    }
    
    //exit(0);
    
    /*
    std::cout << "nfaces: " << faces.size() << std::endl;
    for (FaceList::const_iterator i = faces.begin(); i != faces.end(); ++i)
    {
      std::cout << (*i)[0] << " " << (*i)[1] << " " << (*i)[2] << std::endl;
    }
    */  
    return faces;
  }

  //---------------------------------------------------------------------------


  /// Remove a vertex from a mesh.
  /// Does all the bookkeeping crap: removing the vertex means removing associated faces, and warning neghbors
  /// about the missing of their fellow.
  /// The function does *not* remesh the surface, i.e. if nothing else is done, the surface will present a whole
  /// where the vertex were lying.
  template < typename TVertexNode, typename TFaceNode >
  typename std::list<TVertexNode>::iterator
  remove_vertex
  (
   typename std::list<TVertexNode>::iterator    i,
   std::list<TVertexNode> &                     graph_vertices,
   std::list<TFaceNode> &                       graph_faces
  )
  {
    // remove faces point belongs to
    for (typename TVertexNode::FaceIndexCollection::iterator j = i->faces().begin(); j != i->faces().end(); ++j)
    {
      // Remove index to this face
      for (std::size_t k = 0; k < 3; ++k)
      {
        // don't remove face from facelist of point i itself -- because j is iterating on this very list,
        // plus, it will be deleted with i anyway.
        if ((*j)->face[k] == i) continue;
  #ifndef NDEBUG
        std::size_t tmp = (*j)->face[k]->faces().size();
  #endif
        (*j)->face[k]->faces().remove(*j);
        assert((*j)->face[k]->faces().size() == tmp-1);
      }
      // remove face itself
  #ifndef NDEBUG
      std::size_t tmp = graph_faces.size();
  #endif
      graph_faces.erase(*j);
      assert(graph_faces.size() == tmp-1);
    }
  
    // remove point in neighbor's list
    for (typename TVertexNode::VertexIndexCollection::iterator j = i->neighbors().begin(); j != i->neighbors().end(); ++j)
    {
      (*j)->neighbors().remove(i);
    }
    // remove point
    return graph_vertices.erase(i);
  }


  //---------------------------------------------------------------------------

  // NB: neighbors should be a std::vector so far
  template < typename TNeighbors >
  std::vector<numeric_array<float,2> >
  simple_neighborhood_flattening
  (
   const numeric_array<float,3>                                & point,
  // Point<float,3>                                point,
   const TNeighbors   & neighbors
  )
  {
    //std::cout << "nneighbors " << neighbors.size() << std::endl;
    //std::cout << "point " << point << std::endl;
    //std::cout << "neighbors " << std::endl;
    //for (std::size_t i = 0; i < neighbors.size(); ++i) std::cout << neighbors[i]->pos() << std::endl;
    
  
    std::vector<double> angles;
    std::vector<double> norms;
    angles.reserve(neighbors.size());
    norms.reserve(neighbors.size());
    
    typename TNeighbors::const_iterator n = neighbors.begin();
    const_cyclic_iterator<TNeighbors> n2(neighbors, ++neighbors.begin());
    //std::cout << "point: " << point << std::endl;
    for (; n != neighbors.end(); ++n, ++n2)
    {
      //std::cout << (*n2)->pos() << " " << (*n)->pos() << std::endl;
      //std::cout << (*n2)->pos()-point<< " " << (*n)->pos()-point << std::endl;
      norms.push_back(norm((*n)->pos()-point, prec<double>()));
      /*
      double tmp4 = norm<double>((*n2)->pos()-point);
      double tmp3 = norm<double>((*n)->pos()-point);
      double tmp1 = tmp4 * tmp3;
      double tmp2 = dot((*n2)->pos() - point, (*n)->pos() - point);
      double tmp5 = std::acos( tmp2 / tmp1);
      angles.push_back(tmp5);
      */
      angles.push_back(std::acos(dot((*n2)->pos()-point, (*n)->pos()-point) / (norm((*n2)->pos()-point, prec<double>()) * norm((*n)->pos()-point, prec<double>()))));
  
      //std::cout << angles.size() << " " << angles.capacity() << std::endl;
    }
  
    /*
    std::cout << "norms: ";
    std::copy(norms.begin(), norms.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
  
    std::cout << "angles: ";
    std::copy(angles.begin(), angles.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    */
    
    std::partial_sum(angles.begin(), angles.end(), angles.begin());
  
    /*
    std::cout << "angles summed: ";
    std::copy(angles.begin(), angles.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    */
    
    //double totalAngle = std::accumulate(angles.begin(), angles.end(), 0.0);
    std::transform(angles.begin(), angles.end(), angles.begin(), std::bind2nd(std::multiplies<double>(), 2*M_PI/angles.back()));
  
    /*
    std::cout << "angles normalized: ";
    std::copy(angles.begin(), angles.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    */
    
    std::vector<numeric_array<float, 2> > res(neighbors.size());
    for (std::size_t i = 0; i < size(res); ++i)
    {
      res[i][0] = std::cos(angles[i]) * norms[i];
      res[i][1] = std::sin(angles[i]) * norms[i];
      //res[i][0] = std::cos(angles[i]);
      //res[i][1] = std::sin(angles[i]);
    }
    /*
    std::cout << "newpoints" << std::endl;
    for (std::size_t i = 0; i < res.size(); ++i)
    {
      std::cout << res[i] << std::endl;
    }
    std::cout << "end newpoints" << std::endl;
    */
    
    return res;
  }


  //---------------------------------------------------------------------------


  namespace
  {
    template < typename T >
    struct ValueFinder
    {
      ValueFinder(T value) : m_value(value) {}
      bool operator()(const SimpleOrientedGraphNode<T> & i) { return i.value == m_value; }
      T m_value;
    };
    
    template < typename T >
    ValueFinder<T> valueFinder(T value) { return ValueFinder<T>(value); }
  }
  
  /// A class to remove a vertex from a mesh, and remeshing the hole.
  /// NB: remeshing is not always possible, therefore removing a given point in a mesh is not always possible.
  /// Therefore, it may be important to check the return value of operator() to check if the point has actually been
  /// removed.
  template < typename TVertexNode, typename TFaceNode >
  class Vertex_remover
  {
  public: // enum
  
    enum Error {
      none = 0,
      invalid_neighborhood,
      triangulation_failure,
    };
    
  private: // typedefs
  
    typedef SimpleOrientedGraphNode<typename std::list<TVertexNode>::iterator> NeighborGraphNode;
    typedef std::list<NeighborGraphNode> NeighborGraph;
    typedef typename TVertexNode::VertexIndexCollection::value_type VertexPointer;
    
  public: // constructors
  
    Vertex_remover(std::list<TVertexNode> & graph_vertices, std::list<TFaceNode> & graph_faces)
      : m_graph_vertices(graph_vertices)
      , m_graph_faces(graph_faces)
      , m_error(none)
    {}
  
  public: // set & get
  
    std::vector<VertexPointer> & neighbors() { return m_neighbors; }
    Error error() const { return m_error; }
    
  
  private: // functions
  
    /// Triangularize the neighbors of i.
    /// Return true if successful.
    bool neighborTriangulation(VertexPointer i)
    {
      // get delaunay faces for more than 4 neighbors
      std::size_t nneigh = m_neighbors.size();
      // NB: I don't know why we need this, but we do -- it crashes otherwise
      m_tri.clear();
      if (nneigh >= 4)
      {
        m_tri = simple_delaunay_triangulation(simple_neighborhood_flattening(i->pos(), m_neighbors));
      }
      // For three neighbors, simply add the only triangle possible
      else if (nneigh == 3)
      {
        boost::array<std::size_t,3> tmp = { {0,1,2} };
        m_tri.push_back(tmp);
      }
      // For two points or less, do not add any face. Note that if this happens, it means that the mesh
      // is in really bad shape -- actually the following might crash.
      else 
      {
        m_error = invalid_neighborhood;
        std::cout << "w2!";
        return false;
      }
      // Check consistency between number of points and number of faces
      if (m_tri.size() != nneigh - 2)
      {
        m_error = triangulation_failure;
        std::cout << "wK!" << m_tri.size() << "-" << nneigh-2;
        return false;
      }
      return true;
    }
  
  
    /// Initialize the oriented graph of the neighbors of the neighbors of i, as they were
    /// before triangulation.
    void initializeOrientedEdges(VertexPointer i)
    {
      std::size_t nneigh = m_neighbors.size();
      m_convmap.resize(nneigh);
      for (std::size_t j = 0; j < nneigh; ++j)
      {
        // first pass: vertices
        {
          for (typename TVertexNode::VertexIndexCollection::const_iterator k = m_neighbors[j]->neighbors().begin(); k != m_neighbors[j]->neighbors().end(); ++k)
          {
            if (*k == i) continue;
            NeighborGraphNode node(*k);
            m_convmap[j][*k] = m_neighborGraph[j].insert(m_neighborGraph[j].end(), node);
          }
        }
        // second pass: edges
        {
          typename TVertexNode::VertexIndexCollection::const_iterator k = m_neighbors[j]->neighbors().begin();
          const_cyclic_iterator<typename TVertexNode::VertexIndexCollection> k2(m_neighbors[j]->neighbors(), m_neighbors[j]->neighbors().begin());
          ++k2;
          const_cyclic_iterator<typename TVertexNode::VertexIndexCollection> k3(m_neighbors[j]->neighbors(), m_neighbors[j]->neighbors().begin());
          ++k3;
          ++k3;
          for (; k != m_neighbors[j]->neighbors().end(); ++k, ++k2, ++k3)
          {
            if (*k2 == i) continue;
            if (*k != i)  m_convmap[j][*k2]->from.push_back(m_convmap[j][*k]);
            if (*k3 != i) m_convmap[j][*k2]->to.push_back(m_convmap[j][*k3]);
          }
        }
      }
    }
  
    /// Complete the neighborhood graph of i's neighbors by taking into account the triangulation
    /// of the whole left by removing i.
    void updateNeighborGraph()
    {
      for (std::list<boost::array<std::size_t,3> >::iterator newf = m_tri.begin(); newf != m_tri.end(); ++newf)
      {
        for (std::size_t n = 0; n < 3; ++n)
        {
          typename NeighborGraph::iterator p2, p3;
          
          p2 = std::find_if(m_neighborGraph[(*newf)[n]].begin(), m_neighborGraph[(*newf)[n]].end(), valueFinder(m_neighbors[(*newf)[(n+1)%3]]));
          p3 = std::find_if(m_neighborGraph[(*newf)[n]].begin(), m_neighborGraph[(*newf)[n]].end(), valueFinder(m_neighbors[(*newf)[(n+2)%3]]));
          if (p2 == m_neighborGraph[(*newf)[n]].end()) p2 = m_neighborGraph[(*newf)[n]].insert(p2, NeighborGraphNode(m_neighbors[(*newf)[(n+1)%3]]));
          if (p3 == m_neighborGraph[(*newf)[n]].end()) p3 = m_neighborGraph[(*newf)[n]].insert(p3, NeighborGraphNode(m_neighbors[(*newf)[(n+2)%3]]));
  
          p2->to.push_back(p3);
          p3->from.push_back(p2);
        } 
      }
    }
  
  
    /// Check that neighbors have at least three neighbors.
    bool checkNeighborsHaveThreeNeighbors()
    {
      for (std::size_t k = 0; k < m_neighborGraph.size(); ++k)
      {
        if (m_neighborGraph[k].size() < 3)
        {
          //std::cout << "WREJECT!(small-neighborhood)";
          return false;
        }
      }
      return true;
    }
  
    /// Check that neighbors have an adequate circular neighborhood.
    bool checkCorrectNeighborhood()
    {
      for (std::size_t k = 0; k < m_neighborGraph.size(); ++k)
      {
        //std::cout << k << "/" << m_neighborGraph.size() << std::endl;
        typename NeighborGraph::iterator p = m_neighborGraph[k].begin();
        typename NeighborGraph::iterator p0 = p;
        bool flag = false;
        for (;;)
        {
          //std::cout << ":" << &*p << " " << std::flush;
          if (flag)
          {
            //std::cout << "R" << std::flush;
            if (p == p0) break;
            //std::cout << "T" << std::flush;
          }
          else
          {
            //std::cout << "S" << std::flush;
            flag = true;
          }
          //std::cout << "P" << std::flush;
          if (p->to.size() != 1 || p->from.size() != 1)
          {
            //std::cout << "WREJECT!(incorrect-neighborhood)";
            return false;
          }
          //std::cout << "%" << std::flush;
          p = p->to.front();
          //std::cout << "*" << std::flush;
        }
      }
      return true;
    }
  
  
    void addFaces(VertexPointer i)
    {
      for (std::list<boost::array<std::size_t,3> >::iterator newf = m_tri.begin(); newf != m_tri.end(); ++newf)
      {
        // convert into graph faces
        TFaceNode f;
        for (std::size_t n = 0; n < 3; ++n)
        {
          f.face[n] = m_neighbors[(*newf)[n]];
        }
        // add face to graph
        typename std::list<TFaceNode>::iterator gf = m_graph_faces.insert(m_graph_faces.end(), f);
        numeric_array<float,3> normal = cross(
          i->faces().front()->face[0]->pos() - i->faces().front()->face[1]->pos(),
          i->faces().front()->face[0]->pos() - i->faces().front()->face[2]->pos()
        );
        // add face index to face points
        for (std::size_t n = 0; n < 3; ++n)
        {
          GroComp<VertexPointer> grocomp(f.face[0], f.face[1], f.face[2]);
          typename TVertexNode::FaceIndexCollection::iterator p = find_if(f.face[n]->faces().begin(), f.face[n]->faces().end(), grocomp);
          if (p != f.face[n]->faces().end())
          {
            std::cout << "FATAL: face already there: " << std::endl;
            std::cout << &*f.face[0] << " " << &*f.face[1] << " " << &*f.face[2] << std::endl;
            std::cout << &*(*p)->face[0] << " " << &*(*p)->face[1] << " " << &*(*p)->face[2] << std::endl;
            std::cout << "I'll continue as if nothing happened, but this is a bug and your algorithm is likely to go bananas :)" << std::endl;
          }
          
          // checking normal consistency
          // I removed this test because normality consistency is broken very often, and yet this is legal.
          // Actually, I know it is legal -- it is just surprising that this happens so often (like 1/10th of the 
          // cases!).
          /*
          {
            Vector<float,3> mynormal = cross(
              f.face[0]->pos() - f.face[1]->pos(),
              f.face[0]->pos() - f.face[2]->pos());
            if (dot(normal, mynormal) < 0) std::cout << "wY!";
          }
          */
          f.face[n]->faces().push_back(gf);
        }                
      }
    }
  
    void addNeighbors(VertexPointer i)
    {
      for (std::size_t j = 0; j < i->neighbors().size(); ++j)
      {
        typename TVertexNode::VertexIndexCollection ntmp = m_neighbors[j]->neighbors();
        
        // remove all neighbors
        m_neighbors[j]->neighbors().clear();
  
        typename NeighborGraph::iterator p = m_neighborGraph[j].begin();
        typename NeighborGraph::iterator p0 = p;
        bool flag = false;
        for (;;)
        {
          if (flag) { if (p == p0) break; }
          else flag = true;
          m_neighbors[j]->neighbors().push_back(p->value);
          if (p->to.size() != 1)
          {
            std::cout << "wW!";
            /*
            std::cout << "ori neighbors" << std::endl;
            for (typename std::vector<VertexPointer>::const_iterator i = neighbors.begin(); i != neighbors.end(); ++i)
            {
              std::cout << (*i)->pos() << std::endl;
            }
            std::cout << "end ori neighbors" << std::endl;
            
            std::cout << "my neighbor number is " << j << std::endl;
            
            std::cout << "neighbors" << std::endl;
            for (typename TVertexNode::VertexIndexCollection::const_iterator i = ntmp.begin(); i != ntmp.end(); ++i)
            {
              std::cout << (*i)->pos() << std::endl;
            }
            std::cout << "end neighbors" << std::endl;
            std::cout << "triangles" << std::endl;
            for (typename std::list<boost::array<std::size_t,3> >::const_iterator i = tri.begin(); i != tri.end(); ++i)
            {
              std::cout << (*i)[0] << " " << (*i)[1] << " " << (*i)[2] << std::endl;
            }
            //std::copy(tri.begin(), tri.end(), std::ostream_iterator<boost::array<std::size_t,3> >(std::cout, std::endl));
            std::cout << "end triangles" << std::endl;
            //std::copy(norms.begin(), norms.end(), std::ostream_iterator<double>(std::cout, " "));
            for (std::list<boost::array<std::size_t,3> >::iterator i = tri.begin(); i != tri.end(); ++i)
            {
              std::cout << (*i)[0] << " " << (*i)[1] << " " << (*i)[2] << std::endl;
            }
            
            exit(1);
            */
          }
          if (p->to.size() == 0) break;
          p = p->to.front();
        }
      }
    }  
  
  public: // functions
  
  
    /// Check that vertex i can be removed from the mesh
    bool isRemovable(VertexPointer i)
    {
      m_error = none;
      //std::cout << "A" << std::endl;
      // copy list of neighbors in a vector, for fast random indexing
      //std::cout << m_neighbors.size() << " " << i->neighbors().size() << " = " << std::flush;
      m_neighbors.resize(i->neighbors().size());
      //std::cout << m_neighbors.size() << std::endl;
      std::copy(i->neighbors().begin(), i->neighbors().end(), m_neighbors.begin());
      if (!this->neighborTriangulation(i)) return false;
      // initialize oriented edges
      //std::cout << m_neighborGraph.size() << " " << m_neighbors.size() << " = ";
      m_neighborGraph.clear();
      m_neighborGraph.resize(m_neighbors.size());
      //std::cout << m_neighborGraph.size() << std::endl;
      this->initializeOrientedEdges(i);
      // Updating neighbor graphs
      this->updateNeighborGraph();
      // Checking that all neighbors have at least 3 neighbors
      if (!this->checkNeighborsHaveThreeNeighbors()) return false;
      // Checking that neighbors have correct neighborhood
      // Incorrect neighborhood may appear when a neighbor was already sharing an edge with another neighbor of
      // i that is not on its side, i.e that is not n-1 or n+1.
      //std::cout << "H" << std::endl;
      if (!this->checkCorrectNeighborhood()) return false;
      //std::cout << "I" << std::endl;
      return true;
    }
  
    /// Remove a vertex that has already been nodded as removable by the 'isRemovable' method.
    /// Return a pointer to the vertex 'after' the vertex that has been removed, in the list sense.
    VertexPointer remove(VertexPointer i)
    {
      // adding faces
      this->addFaces(i);
      // adding neighbors
      this->addNeighbors(i);
      // remove point
      return remove_vertex(i, m_graph_vertices, m_graph_faces);
    }
  
    /*
    /// Checks if a vertex is removable, and if so, removes it. 
    /// Returns true if the operation is successful.
    bool operator()(VertexPointer & i)
    {
      //std::cout << "." << std::flush;
      if (!this->isRemovable(i)) return false;
      i = this->remove(i);
      return true;
    }
    */
    
  
    /*  
    bool operator()(VertexPointer & i)
    {
      std::cout << "." << std::flush;
    
      // save list of neighbors
      // TODO: is it still necessary to save the neighbors?
      std::vector<VertexPointer> neighbors(i->neighbors().size());
      std::copy(i->neighbors().begin(), i->neighbors().end(), neighbors.begin());
      
      // Add new faces
      {
        std::list<boost::array<std::size_t,3> > tri;
        // get delaunay faces for more than 4 neighbors
        if (neighbors.size() >= 4)
        {
          tri = simple_delaunay_triangulation(simple_neighborhood_flattening(i->pos(), neighbors));
        }
        // For three neighbors, simply add the only triangle possible
        else if (neighbors.size() == 3)
        {
          boost::array<std::size_t,3> tmp = {0,1,2};
          tri.push_back(tmp);
        }
        // For two points or less, do not add any face. Note that if this happens, it means that the mesh
        // is in really bad shape -- actually the following might crash.
        else std::cout << "w2!";
        
        // Check consistency between number of points and number of faces
        if (tri.size() != neighbors.size() - 2)
        {
          std::cout << "wK!" << tri.size() << "-" << neighbors.size()-2;
        }
    
        //initialize oriented edges
        std::vector<NeighborGraph> neighborGraph(i->neighbors().size());
        std::vector<std::map<VertexPointer, typename std::list<NeighborGraphNode>::iterator, ItComp<VertexPointer> > > convmap(i->neighbors().size());
        for (std::size_t j = 0; j < i->neighbors().size(); ++j)
        {
          // first pass: vertices
          {
            for (typename TVertexNode::VertexIndexCollection::const_iterator k = neighbors[j]->neighbors().begin(); k != neighbors[j]->neighbors().end(); ++k)
            {
              if (*k == i) continue;
              NeighborGraphNode node(*k);
              convmap[j][*k] = neighborGraph[j].insert(neighborGraph[j].end(), node);
            }
          }
          // second pass: edges
          {
            typename TVertexNode::VertexIndexCollection::const_iterator k = neighbors[j]->neighbors().begin();
            const_cyclic_iterator<typename TVertexNode::VertexIndexCollection> k2(neighbors[j]->neighbors(), neighbors[j]->neighbors().begin());
            ++k2;
            const_cyclic_iterator<typename TVertexNode::VertexIndexCollection> k3(neighbors[j]->neighbors(), neighbors[j]->neighbors().begin());
            ++k3;
            ++k3;
            for (; k != neighbors[j]->neighbors().end(); ++k, ++k2, ++k3)
            {
              if (*k2 == i) continue;
              if (*k != i)  convmap[j][*k2]->from.push_back(convmap[j][*k]);
              if (*k3 != i) convmap[j][*k2]->to.push_back(convmap[j][*k3]);
            }
          }
        }
    
        // Updating neighbor graphs
        {
          for (std::list<boost::array<std::size_t,3> >::iterator newf = tri.begin(); newf != tri.end(); ++newf)
          {
            for (std::size_t n = 0; n < 3; ++n)
            {
              typename NeighborGraph::iterator p2, p3;
              
              p2 = std::find_if(neighborGraph[(*newf)[n]].begin(), neighborGraph[(*newf)[n]].end(), valueFinder(neighbors[(*newf)[(n+1)%3]]));
              p3 = std::find_if(neighborGraph[(*newf)[n]].begin(), neighborGraph[(*newf)[n]].end(), valueFinder(neighbors[(*newf)[(n+2)%3]]));
              if (p2 == neighborGraph[(*newf)[n]].end()) p2 = neighborGraph[(*newf)[n]].insert(p2, NeighborGraphNode(neighbors[(*newf)[(n+1)%3]]));
              if (p3 == neighborGraph[(*newf)[n]].end()) p3 = neighborGraph[(*newf)[n]].insert(p3, NeighborGraphNode(neighbors[(*newf)[(n+2)%3]]));
    
              p2->to.push_back(p3);
              p3->from.push_back(p2);
            } 
          }
        }
        
        // Checking that neighbors have at least 3 neighbors
        for (std::size_t k = 0; k < neighborGraph.size(); ++k)
        {
          if (neighborGraph[k].size() < 3)
          {
            std::cout << "WREJECT!(small-neighborhood)";
            return false;
          }
        }
        
        // Checking that neighbors have correct neighborhood
        // Incorrect neighborhood may appear when a neighbor was already sharing an edge with another neighbor of
        // i that is not on its side, i.e that is not n-1 or n+1.
        for (std::size_t k = 0; k < neighborGraph.size(); ++k)
        {
          typename NeighborGraph::iterator p = neighborGraph[k].begin();
          typename NeighborGraph::iterator p0 = p;
          bool flag = false;
          for (;;)
          {
            if (flag) { if (p == p0) break; }
            else flag = true;
            if (p->to.size() != 1 || p->from.size() != 1)
            {
              std::cout << "WREJECT!(incorrect-neighborhood)";
              return false;
            }
            p = p->to.front();
          }
        }
    
    
        // adding faces
        {
          for (std::list<boost::array<std::size_t,3> >::iterator newf = tri.begin(); newf != tri.end(); ++newf)
          {
            // convert into graph faces
            TFaceNode f;
            for (std::size_t n = 0; n < 3; ++n)
            {
              f.face[n] = neighbors[(*newf)[n]];
            }
            // add face to graph
            typename std::list<TFaceNode>::iterator gf = m_graph_faces.insert(m_graph_faces.end(), f);
            Vector<float,3> normal = cross(
              i->faces().front()->face[0]->pos() - i->faces().front()->face[1]->pos(),
              i->faces().front()->face[0]->pos() - i->faces().front()->face[2]->pos()
            );
            // add face index to face points
            for (std::size_t n = 0; n < 3; ++n)
            {
              GroComp<VertexPointer> grocomp(f.face[0], f.face[1], f.face[2]);
              typename TVertexNode::FaceIndexCollection::iterator p = find_if(f.face[n]->faces().begin(), f.face[n]->faces().end(), grocomp);
              if (p != f.face[n]->faces().end())
              {
                std::cout << "FATAL: face already there: " << std::endl;
                std::cout << &*f.face[0] << " " << &*f.face[1] << " " << &*f.face[2] << std::endl;
                std::cout << &*(*p)->face[0] << " " << &*(*p)->face[1] << " " << &*(*p)->face[2] << std::endl;
                std::cout << "I'll continue as if nothing happened, but this is a bug and your algorithm is likely to go bananas :)" << std::endl;
              }          
              f.face[n]->faces().push_back(gf);
            }                
          }
        }
    
        // Adding neighbors
        for (std::size_t j = 0; j < i->neighbors().size(); ++j)
        {
          typename TVertexNode::VertexIndexCollection ntmp = neighbors[j]->neighbors();
          
          // remove all neighbors
          neighbors[j]->neighbors().clear();
    
          typename NeighborGraph::iterator p = neighborGraph[j].begin();
          typename NeighborGraph::iterator p0 = p;
          bool flag = false;
          for (;;)
          {
            if (flag) { if (p == p0) break; }
            else flag = true;
            neighbors[j]->neighbors().push_back(p->value);
            if (p->to.size() != 1)
            {
              std::cout << "wW!";
            }
            if (p->to.size() == 0) break;
            p = p->to.front();
          }
        }
      }
      // remove point
      i = remove_vertex(i, m_graph_vertices, m_graph_faces);
      return true;
    }
    */
    
    bool operator()(VertexPointer & i)
    {
      //std::cout << "." << std::flush;
    
      // save list of neighbors
      // TODO: is it still necessary to save the neighbors?
      
      m_neighbors.resize(i->neighbors().size());
      std::copy(i->neighbors().begin(), i->neighbors().end(), m_neighbors.begin());
      
      // Add new faces
      {
        // TODO: we need this -- dunno why
        m_tri.clear();
        // get delaunay faces for more than 4 neighbors
        if (m_neighbors.size() >= 4)
        {
          m_tri = simple_delaunay_triangulation(simple_neighborhood_flattening(i->pos(), m_neighbors));
        }
        // For three neighbors, simply add the only triangle possible
        else if (m_neighbors.size() == 3)
        {
          boost::array<std::size_t,3> tmp = { {0,1,2} };
          m_tri.push_back(tmp);
        }
        // For two points or less, do not add any face. Note that if this happens, it means that the mesh
        // is in really bad shape -- actually the following might crash.
        else std::cout << "w2!";
        
        // Check consistency between number of points and number of faces
        if (m_tri.size() != m_neighbors.size() - 2)
        {
          std::cout << "wK!" << m_tri.size() << "-" << m_neighbors.size() - 2;
        }
    
        //initialize oriented edges
        //std::vector<NeighborGraph> neighborGraph(i->neighbors().size());
        // TODO: we need this -- dunno why
        m_neighborGraph.clear();
        m_neighborGraph.resize(i->neighbors().size());
        //std::vector<std::map<VertexPointer, typename std::list<NeighborGraphNode>::iterator, ItComp<VertexPointer> > > convmap(i->neighbors().size());
        m_convmap.resize(i->neighbors().size());
        for (std::size_t j = 0; j < i->neighbors().size(); ++j)
        {
          // first pass: vertices
          {
            for (typename TVertexNode::VertexIndexCollection::const_iterator k = m_neighbors[j]->neighbors().begin(); k != m_neighbors[j]->neighbors().end(); ++k)
            {
              if (*k == i) continue;
              NeighborGraphNode node(*k);
              m_convmap[j][*k] = m_neighborGraph[j].insert(m_neighborGraph[j].end(), node);
            }
          }
          // second pass: edges
          {
            typename TVertexNode::VertexIndexCollection::const_iterator k = m_neighbors[j]->neighbors().begin();
            const_cyclic_iterator<typename TVertexNode::VertexIndexCollection> k2(m_neighbors[j]->neighbors(), m_neighbors[j]->neighbors().begin());
            ++k2;
            const_cyclic_iterator<typename TVertexNode::VertexIndexCollection> k3(m_neighbors[j]->neighbors(), m_neighbors[j]->neighbors().begin());
            ++k3;
            ++k3;
            for (; k != m_neighbors[j]->neighbors().end(); ++k, ++k2, ++k3)
            {
              if (*k2 == i) continue;
              if (*k != i)  m_convmap[j][*k2]->from.push_back(m_convmap[j][*k]);
              if (*k3 != i) m_convmap[j][*k2]->to.push_back(m_convmap[j][*k3]);
            }
          }
        }
    
        // Updating neighbor graphs
        {
          for (std::list<boost::array<std::size_t,3> >::iterator newf = m_tri.begin(); newf != m_tri.end(); ++newf)
          {
            for (std::size_t n = 0; n < 3; ++n)
            {
              typename NeighborGraph::iterator p2, p3;
              
              p2 = std::find_if(m_neighborGraph[(*newf)[n]].begin(), m_neighborGraph[(*newf)[n]].end(), valueFinder(m_neighbors[(*newf)[(n+1)%3]]));
              p3 = std::find_if(m_neighborGraph[(*newf)[n]].begin(), m_neighborGraph[(*newf)[n]].end(), valueFinder(m_neighbors[(*newf)[(n+2)%3]]));
              if (p2 == m_neighborGraph[(*newf)[n]].end()) p2 = m_neighborGraph[(*newf)[n]].insert(p2, NeighborGraphNode(m_neighbors[(*newf)[(n+1)%3]]));
              if (p3 == m_neighborGraph[(*newf)[n]].end()) p3 = m_neighborGraph[(*newf)[n]].insert(p3, NeighborGraphNode(m_neighbors[(*newf)[(n+2)%3]]));
    
              p2->to.push_back(p3);
              p3->from.push_back(p2);
            } 
          }
        }
        
        // Checking that neighbors have at least 3 neighbors
        for (std::size_t k = 0; k < m_neighborGraph.size(); ++k)
        {
          if (m_neighborGraph[k].size() < 3)
          {
            //std::cout << "WREJECT!(small-neighborhood)";
            return false;
          }
        }
        
        // Checking that neighbors have correct neighborhood
        // Incorrect neighborhood may appear when a neighbor was already sharing an edge with another neighbor of
        // i that is not on its side, i.e that is not n-1 or n+1.
        for (std::size_t k = 0; k < m_neighborGraph.size(); ++k)
        {
          typename NeighborGraph::iterator p = m_neighborGraph[k].begin();
          typename NeighborGraph::iterator p0 = p;
          bool flag = false;
          for (;;)
          {
            if (flag) { if (p == p0) break; }
            else flag = true;
            if (p->to.size() != 1 || p->from.size() != 1)
            {
              //std::cout << "WREJECT!(incorrect-neighborhood)";
              return false;
            }
            p = p->to.front();
          }
        }
    
    
        // adding faces
        {
          for (std::list<boost::array<std::size_t,3> >::iterator newf = m_tri.begin(); newf != m_tri.end(); ++newf)
          {
            // convert into graph faces
            TFaceNode f;
            for (std::size_t n = 0; n < 3; ++n)
            {
              f.face[n] = m_neighbors[(*newf)[n]];
            }
            // add face to graph
            typename std::list<TFaceNode>::iterator gf = m_graph_faces.insert(m_graph_faces.end(), f);
            /*numeric_array<float,3> normal = cross(
              i->faces().front()->face[0]->pos() - i->faces().front()->face[1]->pos(),
              i->faces().front()->face[0]->pos() - i->faces().front()->face[2]->pos()
            );*/
            // add face index to face points
            for (std::size_t n = 0; n < 3; ++n)
            {
              GroComp<VertexPointer> grocomp(f.face[0], f.face[1], f.face[2]);
              typename TVertexNode::FaceIndexCollection::iterator p = find_if(f.face[n]->faces().begin(), f.face[n]->faces().end(), grocomp);
              if (p != f.face[n]->faces().end())
              {
                std::cout << "FATAL: face already there: " << std::endl;
                std::cout << &*f.face[0] << " " << &*f.face[1] << " " << &*f.face[2] << std::endl;
                std::cout << &*(*p)->face[0] << " " << &*(*p)->face[1] << " " << &*(*p)->face[2] << std::endl;
                std::cout << "I'll continue as if nothing happened, but this is a bug and your algorithm is likely to go bananas :)" << std::endl;
              }            
              f.face[n]->faces().push_back(gf);
            }                
          }
        }
    
        // Adding neighbors
        for (std::size_t j = 0; j < i->neighbors().size(); ++j)
        {
          typename TVertexNode::VertexIndexCollection ntmp = m_neighbors[j]->neighbors();
          
          // remove all neighbors
          m_neighbors[j]->neighbors().clear();
    
          typename NeighborGraph::iterator p = m_neighborGraph[j].begin();
          typename NeighborGraph::iterator p0 = p;
          bool flag = false;
          for (;;)
          {
            if (flag) { if (p == p0) break; }
            else flag = true;
            m_neighbors[j]->neighbors().push_back(p->value);
            if (p->to.size() != 1)
            {
              std::cout << "wW!";
            }
            if (p->to.size() == 0) break;
            p = p->to.front();
          }
        }
      }
      // remove point
      i = remove_vertex(i, m_graph_vertices, m_graph_faces);
      return true;
    }
  
  
  private: // data, input
  
    std::list<TVertexNode> &  m_graph_vertices;
    std::list<TFaceNode> &    m_graph_faces;  
  
  private: // data, output
  
    Error m_error;
  
  private: // data, internal
  
    std::vector<VertexPointer> m_neighbors;
    std::vector<NeighborGraph> m_neighborGraph;
    std::vector<std::map<VertexPointer, typename std::list<NeighborGraphNode>::iterator, Lesser_PointeeAddress<VertexPointer> > > m_convmap;
    std::list<boost::array<std::size_t,3> > m_tri;
  };


  //---------------------------------------------------------------------------

    //---------------------------------//
   //  Remove_indexed_graph_vertices  //
  //---------------------------------//

  /// Remove vertices from a graph list mesh, given a set of index.
  /// This obviously assume that the graph nodes contains such an index.
  template < typename TIndexCollection >
  class Remove_indexed_graph_vertices
  {
  public: // typedefs
    typedef MeshVertexNodeX<std::size_t> VertexNode;
    typedef MeshFaceNodeX<VertexNode> FaceNode;
    typedef std::list<VertexNode> VertexNodeCollection;
    typedef std::list<FaceNode> FaceNodeCollection;
    
  public: // set & get
  
    /// Returns true if it has been able to remove all desired points.
    /// NB: it is very frequent that due to topology constraints on the mesh,
    /// some points could not have been removed.
    bool complete() const { return m_complete; }
    
    /// Number of iterations that has been necessary to remove all points or to
    /// reach a configuration where no more points can be removed.
    unsigned int niter() const { return m_iter; }

    /// Returns the number of points that have been effectively removed.
    unsigned int count() const { return m_count; }

    /// Returns a binary index of those points who have been removed.
    std::vector<unsigned char> const & removed() { return m_res; }

  public: // operators
 
    void operator()
    (
      std::list<VertexNode> & graph_vertices
    , std::list<FaceNode> & graph_faces
    , TIndexCollection const & removed            ///< [input] 1 if wished to be removed, 0 otherwise.
    )
    {
      std::size_t nVertices = graph_vertices.size();
      assert(removed.size() == nVertices);
      
      // set all vertices as unremoved (yet).
      m_res.resize(nVertices, 0);
      m_iter = 0;
      m_count = 0;
      bool done;
      Vertex_remover<VertexNode, FaceNode> vertexRemover(graph_vertices, graph_faces);
      // Sometimes, a vertex cannot be removed. But it may be that it can be removed after
      // other points are removed. That's why we have this loop where we iterate on points that
      // have not yet been able to be removed.
      do
      {
        m_complete = true;
        done = true;
        std::list<VertexNode>::iterator i = graph_vertices.begin();
        std::list<FaceNode>::iterator itmp;
        while (i != graph_vertices.end())
        {
          if (i->neighbors().size() != i->faces().size())
          {
            std::cout << "[SERIOUS TOPOLOGY ERROR] F!=N -- I most probably screwed up your mesh, please report error; " << std::flush;
          }
          // NB: we save i->attr to access it when i is deleted
          std::size_t iInd = i->attribute();
          // do we wish to remove current point?
          if (removed[iInd] == 1)
          {
            // Here, the 'done' boolean should be here.
            // Note that it means that to be robust, I should add a check that the number of points
            // has strictly decreased during two iterations.
            m_complete = false;
            if (vertexRemover(i))
            {
              // NB: i has moved, nothing can rely on i anymore because, especially since it might
              // now point to end().
              // TODO: the complete/false structure is wrong here, because a
              // final pass necessarily means that nothing is done, while we
              // could be smart and detect for the last pass if any point was
              // denied removal.
              done = false;
              // increase removed counter
              ++m_count;
              // mark vertex as removed in output array.
              m_res[iInd] = 1;
              // do *not* increment i, since it has a new value already
              continue;
            }
          }
          ++i;
        }
        ++m_iter;
      } while (!done);
    }

  private: // data, output
    bool m_complete;
    unsigned int m_iter;
    unsigned int m_count;
    std::vector<unsigned char> m_res;
  };

  //---------------------------------------------------------------------------

    //---------------------------//
   //  Remove_indexed_vertices  //
  //---------------------------//

  /// Remove vertices from a mesh, given a set of index.
  template < typename TVertexCollection, typename TFaceCollection, typename TCNeighborhoods, typename TIndexCollection >
  class Remove_indexed_vertices
  {
  public: // typedefs
    typedef typename Remove_indexed_graph_vertices<TIndexCollection>::VertexNodeCollection VertexNodeCollection;
    typedef typename Remove_indexed_graph_vertices<TIndexCollection>::FaceNodeCollection FaceNodeCollection;
  public: // set & get

    Remove_indexed_graph_vertices<TIndexCollection> & remover() { return m_remover; }

    /// Returns true if it has been able to remove all desired points.
    /// NB: it is very frequent that due to topology constraints on the mesh,
    /// some points could not have been removed.
    bool complete() { return m_remover.complete(); }

    /// Returns the number of points that have been effectively removed.
    unsigned int count() { return m_remover.count(); }

    /// Returns a binary index of those points who have been removed.
    std::vector<unsigned char> const & removed() { return m_remover.removed(); }

  public: // operators
    void operator()
    (
      TVertexCollection const & vertices
    , TFaceCollection const & faces
    , TCNeighborhoods const & neighc
    , TIndexCollection const & removed            ///< [input] 1 if wished to be removed, 0 otherwise.
    , TVertexCollection & verticesOut
    , TFaceCollection & facesOut
    )
    {
      // convert mesh into a graph
      VertexNodeCollection graph_vertices;
      FaceNodeCollection graph_faces;
      list2graph_mesh_conversion(vertices, faces, neighc, graph_vertices, graph_faces);     
      // add an attribute to each vertex corresponding to its initial rank
      {
        typename VertexNodeCollection::iterator it = graph_vertices.begin();
        std::size_t gvsize = graph_vertices.size();
        for (std::size_t i = 0; i < gvsize; ++i, ++it)
        {
          it->attribute() = i;
        }
      }
      // remove vertices in the graph
      m_remover(graph_vertices, graph_faces, removed);
      // convert graph back into mesh
      Graph2ListMeshConvertor2<TVertexCollection, TFaceCollection, VertexNodeCollection, FaceNodeCollection> g2l;
      g2l(graph_vertices, graph_faces, verticesOut, facesOut);
    }
  private: // data, internal
    Remove_indexed_graph_vertices<TIndexCollection> m_remover;
  };

  //---------------------------------------------------------------------------


  /// NB: vertices are not removed one by one but in groups via this kind of
  /// functions, because the static "vector" collection is incompatible with
  /// this kind of operation. Consequently, decimation needs to transform the
  /// mesh into a linked list type of graph, which is costly, and therefore
  /// which is not done each time.
  /// Of course, the best would be to stick with a linked list graph in the
  /// first place.
  // TODO:  this should obsiously go into a class.
  template < typename TVertexCollection, typename TFaceCollection, typename TCNeighborhoods, typename TIndexCollection >
  //std::pair<std::vector<unsigned char>, std::size_t>
  std::vector<unsigned char>
  remove_vertices
  (
    TVertexCollection       &  vertices, 
    TFaceCollection         &  faces,
    TCNeighborhoods   const &  neighc,
    TIndexCollection  const &  removed
  )
  {
    Remove_indexed_vertices< TVertexCollection, TFaceCollection, TCNeighborhoods, TIndexCollection >
      remover;
    remover(vertices, faces, neighc, removed, vertices, faces);
    //return std::make_pair(remover.removed(), remover.count());
    return remover.removed();
  }

  //---------------------------------------------------------------------------
  
  /// Attempt to label a vertex as removable, and if successfull, label its
  /// neighbors as unremovable.
  template < typename TVertexCollection, typename TNeighborhoodCollection >
  class Neigh_decimation_labeler
  {
  public: // typedefs
    typedef unsigned char Label;
    typedef std::vector<Label> LabelCollection;
  public: // constants
    /// Vertex labels.
    /// Unprocessed: no decision has been made yet.
    /// Remove: vertex labeled as to be removed.
    /// Keep: vertex labeled as unremovable.
    enum { UNPROCESSED = 0, KEEP, REMOVE }; 
  public: // constructors
    Neigh_decimation_labeler
    (
      TVertexCollection const & vertices
    , TNeighborhoodCollection const & neighs
    , LabelCollection & label
    )
      : m_vertices(vertices)
      , m_neighs(neighs)
      , m_label(label)
    {
      assert(vertices.size() == neighs.size());
      m_label.resize(vertices.size(), UNPROCESSED);
    }
  public: // set & get
    TVertexCollection const & vertices() { return m_vertices; }
    TNeighborhoodCollection const & neighs() { return m_neighs; }
    LabelCollection & label() { return m_label; }
  public: // operator
    /// Try to remove i-th point. Returns true if successful.
    bool operator()(std::size_t i)
    {
      // skip point if a decision has already been taken about it
      if (label()[i] != UNPROCESSED) return false;
      // label point as removed
      label()[i] = REMOVE;
      // mark all of its neighbors as kept
      for (std::size_t j = 0; j < neighs()[i].size(); ++j)
      {
        assert(neighs()[i][j] < label().size());
        assert(label()[neighs()[i][j]] != REMOVE);
        label()[neighs()[i][j]] = KEEP;
      }
      return true;
    }
  private: // data, input
    TVertexCollection const & m_vertices;
    TNeighborhoodCollection const & m_neighs;
  private: // data, internal
    LabelCollection & m_label;
  };
  
  //---------------------------------------------------------------------------

  /// Decimates a mesh until a specified number of vertices is reached.
  /// Do nothing if the desired size is actually larger than the current size.
  template < typename TVertexCollection, typename TFaceCollection >
  std::list<std::size_t> 
  quantizer
  (
    TVertexCollection & vertices  ///< [input/output] mesh vertices
  , TFaceCollection & faces       ///< [input/output] mesh faces
  , std::size_t newSize           ///< [input] desired number of vertices
  )
  {
    
    // initialize a list containing the index of the vertices.
    // at the end only the index of the remaining vertices will remain.
    std::list<std::size_t> res(vertices.size());
    std::generate(res.begin(), res.end(), Incrementor<std::size_t>(0));

    // do nothing if nothing to do
    if (newSize >= vertices.size()) return res;
  
    // main loop, which stops when enough points are removed
    do
    {
      std::size_t n = vertices.size();
      //std::vector<unsigned char> removed(n, 0);
      typedef std::vector<std::vector<std::size_t> > Neighborhoods;
      shared_ptr<Neighborhoods> neighc = circular_neighborhoods(vertices, faces);

      // Compute a criterion to order candidates to removal
      /*
      Mesh_curvature<MeshTraits<Mesh_N>::VertexCollection, std::vector<std::vector<std::size_t> >, float>
        mc(getVertices(mesh), *neighc);
      */
      MeshWaveletEnergy<TVertexCollection, Neighborhoods, float> mwe(vertices, *neighc);
      std::vector<std::pair<std::size_t, float> > sortedCurv(n);
      for (std::size_t i = 0; i < n; ++i)
      {
        /*
        mc.process(i);
        sortedCurv[i] = std::make_pair(i, min(1.0f, std::abs(mc.gaussianCurvature())));
        */
        sortedCurv[i] = std::make_pair(i, mwe(i));
      }
      // order points according to their curvatures, higher curvatures first.
      std::sort(sortedCurv.begin(), sortedCurv.end(), Greater_Pair2<std::size_t, float>());
      
      //for (std::size_t i = 0; i < n; ++i) cout << sortedCurv[i].second << " ";
      //std::cout << std::endl;

      // loop through all vertices, higher curvature first
      std::vector<unsigned char> removed;
      Neigh_decimation_labeler<TVertexCollection, Neighborhoods>
        declabeler(vertices, *neighc, removed);
      std::size_t count = 0;
      for (std::size_t i0 = 0; i0 < n; ++i0)
      {
        // get real index of current vertex
        std::size_t i = sortedCurv[i0].first;
        // Try to label i as removable, go to next vertex if unsuccessful
        if (!declabeler(i)) continue;
        // add a fraction of the displacement to neighbors
        for (std::size_t j = 0; j < (*neighc)[i].size(); ++j)
        { 
          vertices[(*neighc)[i][j]] += (vertices[i] - vertices[(*neighc)[i][j]]) * (1.0f / (*neighc)[i].size());
        }
        ++count;
        // stop here if enough points are labeled
        if (n - count < newSize) break;
      }
      // Put 'removed' in standard binary form
      for (std::size_t i = 0; i < n; ++i) removed[i] = (removed[i] == 2 ? 1 : 0);
      //std::pair<std::vector<unsigned char>, std::size_t> tmp
      std::vector<unsigned char> tmp
        = remove_vertices(vertices, faces, *neighc, removed);
      // remove indices of removed vertices from the list of remaining vertices
      {
        std::size_t i = 0;
        std::list<std::size_t>::iterator iRes = res.begin();
        for (; iRes != res.end(); ++i)
        {
          if (tmp[i])
          {
            iRes = res.erase(iRes);
          }
          else
          {
            ++iRes;
          }
        }
      }      
      //mesh.getNeighborIndices() = getNeighborIndices(mesh);
    //} while (0);
    } while (vertices.size() > newSize);
    
    return res;
  }

  //---------------------------------------------------------------------------

  /// Decimates a mesh until no more vertex with unremoved neighbors can be
  /// removed
  template < typename TVertexCollection, typename TFaceCollection >
  std::list<std::size_t> 
  quantizer_2
  (
    TVertexCollection & vertices  ///< [input/output] mesh vertices
  , TFaceCollection & faces       ///< [input/output] mesh faces
  )
  {
    std::size_t n = vertices.size();
    // initialize a list containing the index of the vertices.
    // at the end only the index of the remaining vertices will remain.
    std::list<std::size_t> res(n);
    std::generate(res.begin(), res.end(), Incrementor<std::size_t>(0));
  
    typedef std::vector<std::vector<std::size_t> > Neighborhoods;
    shared_ptr<Neighborhoods> neighc = circular_neighborhoods(vertices, faces);

    // Compute a criterion to order candidates to removal
    MeshWaveletEnergy<TVertexCollection, Neighborhoods, float> mwe(vertices, *neighc);
    std::vector<std::pair<std::size_t, float> > sortedCurv(n);
    for (std::size_t i = 0; i < n; ++i)
    {
      sortedCurv[i] = std::make_pair(i, mwe(i));
    }
    // order points according to their curvatures, higher curvatures first.
    std::sort(sortedCurv.begin(), sortedCurv.end(), Greater_Pair2<std::size_t, float>());
    
    // loop through all vertices, higher curvature first
    std::vector<unsigned char> removed;
    Neigh_decimation_labeler<TVertexCollection, Neighborhoods>
      declabeler(vertices, *neighc, removed);
    for (std::size_t i0 = 0; i0 < n; ++i0)
    {
      // get real index of current vertex
      std::size_t i = sortedCurv[i0].first;
      // Try to label i as removable, go to next vertex if unsuccessful
      if (!declabeler(i)) continue;
      // add a fraction of the displacement to neighbors
      for (std::size_t j = 0; j < (*neighc)[i].size(); ++j)
      { 
        vertices[(*neighc)[i][j]] += (vertices[i] - vertices[(*neighc)[i][j]]) * (1.0f / (*neighc)[i].size());
      }
    }
    // Put 'removed' in standard binary form
    for (std::size_t i = 0; i < n; ++i) removed[i] = (removed[i] == 2 ? 1 : 0);
    //std::pair<std::vector<unsigned char>, std::size_t> tmp
    std::vector<unsigned char> tmp
      = remove_vertices(vertices, faces, *neighc, removed);
    // remove indices of removed vertices from the list of remaining vertices
    {
      std::size_t i = 0;
      std::list<std::size_t>::iterator iRes = res.begin();
      for (; iRes != res.end(); ++i)
      {
        if (tmp[i])
        {
          iRes = res.erase(iRes);
        }
        else
        {
          ++iRes;
        }
      }
    }
    return res;
  }


  /*
  template < typename TIndex, typename TVertexXsr, typename TFaceXsr, typename TNeighborhoodXsr >
  std::list<std::size_t> 
  quantizer2
  (
    TIndex begin,
    TIndex end,
    TVertexXsr vertexXsr,
    TFaceXsr faceXsr,
    TNeighborhoodXsr neighXsr,
    std::size_t newSize
  )
  {
    typedef TIndex                                  index_type;
    typedef typename TNeighborhoodXsr::reference    NeighborhoodRef;
    typedef typename TNeighborhoodXsr::value_type   Neighborhood;
    
    std::size_t bigCount = 0;
    
    // Create a list holding the index number of the vertices.
    std::list<std::size_t> res(vertices.size());
    std::generate(res.begin(), res.end(), Incrementor<std::size_t>(0));
  
    do
    {
      std::size_t n = vertices.size();
      std::vector<unsigned char> removed(n, 0);
      // compute decimation criterion
      MeshWaveletEnergy2<TVertexCollection, Neighborhoods, float> mwe(vertexXsr, neighXsr);
      std::vector<std::pair<std::pair<TIndex, std::size_t>, float> > sortedCurv(n);

      std::size_t i = 0;  
      for (index_type ind = begin; ind != end; ++ind, ++i)
      {
        sortedCurv[i] = std::make_pair(std::make_pair(ind, i), mwe(ind));
      }
      // order points according to their curvatures, higher curvatures first.
      std::sort(sortedCurv.begin(), sortedCurv.end(), Greater_Pair2<std::size_t, float>());
      // loop through all vertices, higher curvature first
      std::size_t count = 0;
      for (std::size_t i0 = 0; i0 < n; ++i0)
      {
        // get real index of current vertex
        std::pair<TIndex, std::size_t> tmp = sortedCurv[i0].first;
        TIndex ind = tmp.first;
        std::size_t i = tmp.second;
        // skip point if not labeled as 'unprocessed'
        if (removed[i]) continue;
        // label point as removed
        removed[i] = 2;
        ++count;
        if (n - count < newSize)
        {
          break;
        }
        // mark all of its neighbors as kept
        {
          NeighborhoodRef nh = neighXsr(ind);
          for (typename Neighborhood::const_iterator iN = nh.begin(); iN != nh.end(); ++iN)
          {
            assert(removed(*iN) != 2);
            removed(*iN)  = 1;
          }
        }
        for (std::size_t j = 0; j < neighXsr() (*neighc)[i].size(); ++j)
        {
          if (removed[(*neighc)[i][j]] == 2)
          {
            std::cout << "(e!neighbor_marked_as_removed)";
          }
          removed[(*neighc)[i][j]] = 1;
        }
        }
        // add a fraction of the displacement to neighbors
        {
          for (std::size_t j = 0; j < (*neighc)[i].size(); ++j)
          {
            vertices[(*neighc)[i][j]] += (vertices[i] - vertices[(*neighc)[i][j]]) * (1.0f / (*neighc)[i].size());
          }
        }
      }
      // Put 'removed' in standard binary form
      for (std::size_t i = 0; i < n; ++i) removed[i] = (removed[i] == 2 ? 1 : 0) ;
      //std::cout << "Marked " << count << " points out of " << n << " for deletion" << std::endl;
      std::pair<std::vector<unsigned char>, std::size_t> tmp = remove_vertices(vertices, faces, *neighc, removed);
      //std::cout << count << "points could be effectively removed" << std::endl;
      bigCount += tmp.second;
      
      {
        std::size_t i = 0;
        std::list<std::size_t>::iterator iRes = res.begin();
        for (; iRes != res.end(); ++i)
        {
          if (tmp.first[i])
          {
            iRes = res.erase(iRes);
          }
          else
          {
            ++iRes;
          }
        }
      }
      
      //mesh.getNeighborIndices() = getNeighborIndices(mesh);
    //} while (0);
    //} while (bigCount < N);
    } while (vertices.size() > newSize);
    
    return res;
  }
  */
  
} // namespace til


//#include "mesh_decimation.tpp"

#endif

