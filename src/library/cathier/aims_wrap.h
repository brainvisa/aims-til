#ifndef AIMS_WRAP_H_
#define AIMS_WRAP_H_


/// \file Contains all code related to interfacing with the aims library.

// includes from AIMS
#include "aims/io/reader.h"
#include "aims/io/writer.h"
#include "aims/mesh/surface.h"
#include "aims/mesh/texture.h"
#include "aims/utility/anytype_reader.h"
#include "aims/vector/vector.h"

// includes from TIL
#include "til/cat2type.h"
#include "til/is_traits.h"
#include "til/numeric_array.h"
//#include "til/Point.h"

// includes from TIL
//#include "convert.h"
#include "globalTraits.h"
#include "MeshTraits.h"
#include "meshUtils.h"

#define TIL_COMMA , 

//-----------------------------------------------------------------------------

/// Lexicographic order.
/// This means, v1 < v2 iff there is an index 0<= i < D so that
/// v1[j] = v2[j], j<i, and v1[i] < v2[i]
/// NB: an order is necessary even just to use a std::set
template < typename T, int D >
inline
bool operator < (const AimsVector<T,D> &v1,
                 const AimsVector<T,D> &v2)
{
  for (int d = 0; d < D; ++d)
    {
      if (v1.item(d) < v2.item(d))
        return true;
      else if (v2.item(d) < v1.item(d))
        return false;
    }
  // v1 and v2 are equal: return false
  return false;
}

//-----------------------------------------------------------------------------

/// Hack around bad AimsApplication definition.
/// Very bad and unsafe! Don't use it for anything else than AimsApplication.
inline 
const char ** aims_const_hack(char** argv)
{
  return static_cast<const char**>(static_cast<void*>(argv));
  //const char * &temp = const_cast<char*&>(*argv);
  //return &temp;
}

//-----------------------------------------------------------------------------

namespace til
{
  //TODO: add an aims namespace for those functions that are not extending functionality of til but that are really
  // adding functions in the namespace for manual convertions.
  
  ////////////////////////////////// global ///////////////////////////////////////////
  
  /// Returns the number of elements in a AimsVector
  template < class T, int D >
  inline std::size_t size(const AimsVector<T,D> &)
  {
    return D;
  }
  
  template < typename T, int D >
  void
  add(AimsVector<T,D> & v1, const AimsVector<T,D> & v2)
  {
    v1 += v2;
  }
  template <typename T, int D>
  void mul(AimsVector<T,D> &v, T d)
  {
    v *= d;
  }
  
  template < typename T, int D >
  struct value_type_of<AimsVector<T,D> > { typedef T type; };
  
  ////////////////////////////////misc////////////////////////////////////
  
  
  
  
  
  ////////////////////////////////// istraits ///////////////////////////////////////////
  
  TIL_DECLARE_IS_SPEC_T(3DPoint, typename T, AimsVector<T TIL_COMMA 3>);
  
  ////////////////////////////////// convert ///////////////////////////////////////////
  
  
  namespace functor
  {
    template < >
    class CastTo< boost::array<size_t, 3>, AimsVector<uint, 3> >
     : public std::binary_function<boost::array<size_t, 3> &, const AimsVector<uint, 3> &, void>
    {
    public:
      void operator()(boost::array<size_t, 3> & y, const AimsVector<uint, 3> & x) const
      {
        y[0] = x[0]; y[1] = x[1]; y[2] = x[2];
      }
    };

    template < >
    class CastTo< AimsVector<uint, 3>, boost::array<size_t, 3> >
     : public std::binary_function<AimsVector<uint, 3> &, const boost::array<size_t, 3> &, void>
    {
    public:
      void 
      operator()(AimsVector<uint, 3> & y, const boost::array<size_t, 3> & x) const
      {
        y[0] = x[0]; y[1] = x[1]; y[2] = x[2];
      }
    };        

    template < >
    class CastTo< numeric_array<size_t, 3>, AimsVector<uint, 3> >
     : public std::binary_function<numeric_array<size_t, 3> &, const AimsVector<uint, 3> &, void>
    {
    public:
      void operator()(numeric_array<size_t, 3> & y, const AimsVector<uint, 3> & x) const
      {
        y[0] = x[0]; y[1] = x[1]; y[2] = x[2];
      }
    };

    template < >
    class CastTo< AimsVector<uint, 3>, numeric_array<size_t, 3> >
     : public std::binary_function<AimsVector<uint, 3> &, const numeric_array<size_t, 3> &, void>
    {
    public:
      void 
      operator()(AimsVector<uint, 3> & y, const numeric_array<size_t, 3> & x) const
      {
        y[0] = x[0]; y[1] = x[1]; y[2] = x[2];
      }
    };
    
    template < >
    class CastTo< numeric_array<float,3>, AimsVector<float,3> >
     : public std::binary_function<numeric_array<float,3> &, const AimsVector<float,3> &, void>
    {
    public:
      void 
      operator()(numeric_array<float,3> & y, const AimsVector<float,3> & x) const
      {
        y[0] = x[0]; y[1] = x[1]; y[2] = x[2];
      }
    };        

    template < >
    class CastTo< AimsVector<float,3>, numeric_array<float,3> >
     : public std::binary_function<AimsVector<float,3> &, const numeric_array<float,3> &, void>
    {
    public:
      void 
      operator()(AimsVector<float,3> & y, const numeric_array<float,3> & x) const
      {
        y[0] = x[0]; y[1] = x[1]; y[2] = x[2];
      }
    };
    /*
    template < >
    class CastTo< Point<float,3>, AimsVector<float,3> >
     : public std::binary_function<Point<float,3> &, const AimsVector<float,3> &, void>
    {
    public:
      void 
      operator()(Point<float,3> & y, const AimsVector<float,3> & x) const
      {
        y[0] = x[0]; y[1] = x[1]; y[2] = x[2];
      }
    };        

    template < >
    class CastTo< AimsVector<float,3>, Point<float,3> >
     : public std::binary_function<AimsVector<float,3> &, const Point<float,3> &, void>
    {
    public:
      void 
      operator()(AimsVector<float,3> & y, const Point<float,3> & x) const
      {
        y[0] = x[0]; y[1] = x[1]; y[2] = x[2];
      }
    };        
    */

    template <typename TContainer>
    class CastTo<Texture1d, TContainer>
     : public std::binary_function<Texture1d &, const TContainer &, void>
    {
    public:
      void
      operator()(Texture1d & t, const TContainer & c)
      {
        for (std::size_t i=0; i<size(c); ++i)
          t.item(i) = c[i];
      }
    };

  } // namespace functor

/*  
  template < typename T1, typename T2, int D >
  inline
  void convert(const AimsVector<T1,D> & x, boost::array<T2,D> & y)
  {
    detail::convert_fixedLoop<D>(x,y);
  }
  
  template < typename T1, typename T2, int D >
  inline
  void convert(const Vector<T1,D> & x, AimsVector<T2,D> & y)
  {
    detail::convert_fixedLoop<D>(x,y);
  }
  
  template < typename T1, typename T2, int D >
  inline
  void convert(const AimsVector<T1,D> & x, Point<T2,D> & y)
  {
    detail::convert_fixedLoop<D>(x,y);
  }
  */
  
  
  // This is assuming Texture1d is cleared
  template <typename TContainer>
  inline void convert(const TContainer &c, Texture1d & t)
  {
    for (std::size_t i=0; i<size(c); ++i)
      t.item(i) = c[i];
  }
  
  
  /*
  template < typename T1, typename T2>
  inline
  void convert(const AimsVector<T1,3> & x, Point<T2,3> & y)
  {
    detail::convert_fixedLoop<3>(x,y);
  }
  */
  
  
  /*
  template < typename T, int D >
  inline
  void convert(const AimsVector<T, D> & x, boost::array<T, D> & v)
  {
    for (int i=0; i<D; ++i) y[i]=x[i];
  }
  */
  /*
  /// conversion between vector and aimsvector
  template < typename T1, typename T2, int D >
  inline
  void convert(const AimsVector<T1,D> & x, Vector<T2,D> & y)
  {
    for (int i=0; i<D; ++i) convert(x[i], y[i]);
  }
  */
  
  ////////////////////////////////// meshtraits ///////////////////////////////////////////
  
  
  //template < int D, class T > class AimsTimeSurfaceFaceCollection;
  
  /// Mesh traits for AimsTimeSurface.
  template < int D, class T >
  //struct MeshTraits<AimsTimeSurface<D,T> >
  struct MeshTraits<AimsTimeSurface<D,T> > : public MeshTraits_default
  {
    // constants
    enum { has_faces_indices = 1 };
    
    // typedefs
    typedef Point3df                            Vertex;
    typedef std::vector<Vertex>                 VertexCollection;
    typedef AimsVector<uint, D>                 FaceIndex;
    typedef std::vector<FaceIndex>              FaceIndexCollection;
    //typedef AimsTimeSurfaceFaceCollection<D,T>  FaceCollection;
  };
  
  
  /*
  /// A class to wrap data in an AimsSurface as if it were a list
  /// of faces giving coordinates.
  /// Experimental.
  template < int D, class T >
  class AimsTimeSurfaceFaceCollection
  {
  public: // constructors & destructor
    AimsTimeSurfaceFaceCollection(AimsTimeSurface<D,T> &mesh) :
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
        Face(const AimsTimeSurface<D,T> & mesh, typename MeshTraits<AimsTimeSurface<D,T> >::FaceIndex &f)
          : m_mesh(mesh), m_f(f) {}
    
      public: // functions
        typename MeshTraits<AimsTimeSurface<D,T> >::Vertex &
          operator[](int i)
        {
          return getVertices(m_mesh)[m_f[i]];
        }
    
      private: // data
        const AimsTimeSurface<D,T> & m_mesh;
        typename MeshTraits<AimsTimeSurface<D,T> >::FaceIndex & m_f;
      };
    * /
    
    public: // constructors & destructor
      ConstIterator(AimsTimeSurface<D,T> & mesh, typename MeshTraits<AimsTimeSurface<D,T> >::FaceIndexCollection::iterator i)
        : m_vertices(getVertices(mesh)), m_i(i) {}
    
    public: // functions
      inline void operator++() { ++m_i; }
      inline const ConstIterator & operator*() const { return (*this); }
      inline bool operator!=(const ConstIterator &i) const { return m_i != i.m_i; }
      inline const typename MeshTraits<AimsTimeSurface<D,T> >::Vertex &
      operator[](int i) const
      {
        return m_vertices[(*m_i)[i]];
      }
          
    private: // data
      const typename MeshTraits<AimsTimeSurface<D,T> >::VertexCollection & m_vertices;
      //const AimsTimeSurface<D,T> & m_mesh;
      typename MeshTraits<AimsTimeSurface<D,T> >::FaceIndexCollection::iterator m_i;
      //Face m_face;
    };
    typedef ConstIterator const_iterator;
  
    ConstIterator begin()  { return ConstIterator(m_mesh, getFaceIndices(m_mesh).begin()); }
    ConstIterator end()    { return ConstIterator(m_mesh, getFaceIndices(m_mesh).end()); }
  
  private: // data
    AimsTimeSurface<D,T> & m_mesh;
  };
  */
  
  
  /// Get all edges in mesh.
  std::set<AimsVector<uint, 2> >
  getEdges(const AimsSurfaceTriangle *surf);
  
  //////////////////////////////////////////////////////////////////////////////
  
  
  //template < int D, class T > class AimsSurfaceFaceCollection;
  
  /// Mesh traits for AimsSurface.
  template < int D, class T >
  //struct MeshTraits<AimsTimeSurface<D,T> >
  struct MeshTraits<AimsSurface<D,T> > : public MeshTraits_default
  {
    // constants
    static const bool has_edges = false;
    static const bool has_faces = false;
    static const bool has_faces_indices = true;
    
    // typedefs
    typedef Point3df                        Vertex;
    typedef std::vector<Vertex>             VertexCollection;
    typedef AimsVector<uint, D>             FaceIndex;
    typedef std::vector<FaceIndex>          FaceIndexCollection;
    //typedef AimsSurfaceFaceCollection<D,T>  FaceCollection;
  };
  
  
  
  /*
  /// A class to wrap data in an AimsSurface as if it were a list
  /// of faces giving coordinates
  template < int D, class T >
  class AimsSurfaceFaceCollection
  {
  public: // constructors & destructor
    AimsSurfaceFaceCollection(AimsSurface<D,T> &mesh) :
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
        Face(const AimsTimeSurface<D,T> & mesh, typename MeshTraits<AimsTimeSurface<D,T> >::FaceIndex &f)
          : m_mesh(mesh), m_f(f) {}
    
      public: // functions
        typename MeshTraits<AimsTimeSurface<D,T> >::Vertex &
          operator[](int i)
        {
          return getVertices(m_mesh)[m_f[i]];
        }
    
      private: // data
        const AimsTimeSurface<D,T> & m_mesh;
        typename MeshTraits<AimsTimeSurface<D,T> >::FaceIndex & m_f;
      };
    * /
    
    public: // constructors & destructor
      ConstIterator(AimsSurface<D,T> & mesh, typename MeshTraits<AimsSurface<D,T> >::FaceIndexCollection::iterator i)
        : m_vertices(getVertices(mesh)), m_i(i) {}
    
    public: // functions
      inline void operator++() { ++m_i; }
      inline const ConstIterator & operator*() const { return (*this); }
      inline bool operator!=(const ConstIterator &i) const { return m_i != i.m_i; }
      inline const typename MeshTraits<AimsSurface<D,T> >::Vertex &
      operator[](int i) const
      {
        return m_vertices[(*m_i)[i]];
      }
          
    private: // data
      const typename MeshTraits<AimsSurface<D,T> >::VertexCollection & m_vertices;
      //const AimsSurface<D,T> & m_mesh;
      typename MeshTraits<AimsSurface<D,T> >::FaceIndexCollection::iterator m_i;
      //Face m_face;
    };
    typedef ConstIterator const_iterator;
  
    ConstIterator begin()  { return ConstIterator(m_mesh, getFaceIndices(m_mesh).begin()); }
    ConstIterator end()    { return ConstIterator(m_mesh, getFaceIndices(m_mesh).end()); }
  
  private: // data
    AimsSurface<D,T> & m_mesh;
  };
  */
  
  
  template < int D, class T >
  const typename MeshTraits<AimsSurface<D,T> >::VertexCollection &
  getVertices(const AimsSurface<D,T> &mesh)
  {
    return mesh.vertex();
  }
  
  template < int D, class T >
  typename MeshTraits<AimsSurface<D,T> >::VertexCollection &
  getVertices(AimsSurface<D,T> &mesh)
  {
    return mesh.vertex();
  }
  
  template < int D, class T >
  inline
  //const typename MeshTraits<AimsTimeSurface<D, T> >::Vertex &
  const Point3df &
  //getFaceVertex(const AimsTimeSurface<D,T> &mesh, const typename MeshTraits<AimsTimeSurface<D,T> >::FaceIndexCollection::const_iterator &iFC, int i)
  getFaceVertex
  (
   const AimsSurface<D,T> & mesh,
   const std::vector<AimsVector<unsigned int, 3> >::const_iterator & iFC,
   int i
  )
  {
    return mesh.vertex()[(*iFC)[i]];
  }
  
  
  template < int D, class T >
  inline
  const typename MeshTraits<AimsSurface<D,T> >::FaceIndexCollection &
  getFaceIndices(const AimsSurface<D,T> & mesh)
  {
    return mesh.polygon();
  }

  /*
  template < int D, class T >
  inline
  typename MeshTraits<AimsSurface<D,T> >::FaceIndexCollection &
  getFaceIndices_my(const AimsSurface<D,T> & mesh)
  {
    return mesh.polygon();
  }
  */
  
  template < int D, class T >
  inline
  typename MeshTraits<AimsSurface<D,T> >::FaceIndexCollection &
  getFaceIndices(AimsSurface<D,T> & mesh)
  {
    return mesh.polygon();
  }
  
  /*
  template < int D, class T >
  inline
  const typename MeshTraits<AimsSurface<D,T> >::FaceCollection
  getFaces(const AimsSurface<D,T> & mesh)
  {
    return AimsSurfaceFaceCollection<D,T>(mesh);
  }
  */
  /*
  template < int D, class T >
  inline
  //typename MeshTraits<AimsTimeSurface<D,T> >::FaceCollection
  typename MeshTraits<AimsSurface<D,T> >::FaceCollection
  getFaces(AimsSurface<D,T> & mesh)
  {
    return AimsSurfaceFaceCollection<D,T>(mesh);
  }
  */
  
  template < int D, class T >
  const typename MeshTraits<AimsTimeSurface<D,T> >::VertexCollection &
  getVertices(const AimsTimeSurface<D,T> &mesh)
  {
    return mesh.vertex();
  }
  
  template < int D, class T >
  typename MeshTraits<AimsTimeSurface<D,T> >::VertexCollection &
  getVertices(AimsTimeSurface<D,T> &mesh)
  {
    return mesh.vertex();
  }
  
  
  /// Returns the i-th vertex of a face of a mesh of type AimsTimeSurface.
  template < int D, class T >
  inline
  //const typename MeshTraits<AimsTimeSurface<D, T> >::Vertex &
  const Point3df &
  //getFaceVertex(const AimsTimeSurface<D,T> &mesh, const typename MeshTraits<AimsTimeSurface<D,T> >::FaceIndexCollection::const_iterator &iFC, int i)
  getFaceVertex
  (
   const AimsTimeSurface<D,T> & mesh,         ///< The mesh
   const std::vector<AimsVector<unsigned int, 3> >::const_iterator & iFC,   ///< An iterator pointing on the face
   int i                                      ///< The number of the face point
  )
  {
    return mesh.vertex()[(*iFC)[i]];
  }
  template < int D, class T >
  const typename MeshTraits<AimsTimeSurface<D,T> >::FaceIndexCollection &
  getFaceIndices(const AimsTimeSurface<D,T> & mesh)
  {
    return mesh.polygon();
  }
  
  /// Return the face indices of a mesh of type AimsTimeSurface.
  template < int D, class T >
  typename MeshTraits<AimsTimeSurface<D,T> >::FaceIndexCollection &
  getFaceIndices(AimsTimeSurface<D,T> & mesh)
  {
    return mesh.polygon();
  }
  
  /*
  template < int D, class T >
  const typename MeshTraits<AimsTimeSurface<D,T> >::FaceCollection
  getFaces(const AimsTimeSurface<D,T> & mesh)
  {
    return AimsTimeSurfaceFaceCollection<D,T>(mesh);
  }
  */
  
  /*
  template < int D, class T >
  //typename MeshTraits<AimsTimeSurface<D,T> >::FaceCollection
  typename MeshTraits<AimsTimeSurface<D,T> >::FaceCollection
  getFaces(AimsTimeSurface<D,T> & mesh)
  {
    return AimsTimeSurfaceFaceCollection<D,T>(mesh);
  }
  */

  /*
  class BundleLoader : public comist::BundleListener
  {
  public: // typedefs
  
    typedef numeric_array<float,3> Point;
    typedef std::vector<Point> Fiber;
    typedef std::vector<Fiber> Bundle;
      
  public: // constructors & destructor
  
    BundleLoader() : m_fibers(new Bundle) {}
    virtual ~BundleLoader() {}
  
  public: // set & get
  
    boost::shared_ptr<Bundle> getFibers() const { return m_fibers; }
  
  private: // functions
  
    void bundleStarted(const comist::BundleProducer &, const comist::BundleInfo &)
    {
    }
    void bundleTerminated(const comist::BundleProducer &, const comist::BundleInfo &)
    {
    }
    void fiberStarted(const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo &)
    {
      m_fibers->push_back(Fiber());
    }
    void fiberTerminated(const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo &)
    {
    }  
    void newFiberPoint(const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo &, const comist::FiberPoint & point)
    {
      m_fibers->back().push_back(til::convert<Point>(point));
    }
    void noMoreBundle(const comist::BundleProducer &)
    {
    }
    
  private: // data
  
    boost::shared_ptr<Bundle> m_fibers;
  };
  */


  //////////////////// MESH UTILS ////////////////////////////

  namespace functor
  {

    template < typename TParam, typename T >
    typename boost::enable_if<boost::is_pointer<typename TParam::FaceIndex::value_type> >::type
    convertAimsMesh(Mesh<TParam> & mesh, const T & aimsmesh)
    {
      detail::convert_mesh_2(aimsmesh, mesh);
    }

    template < typename TParam, typename T >
    typename boost::disable_if<boost::is_pointer<typename TParam::FaceIndex::value_type> >::type
    convertAimsMesh(Mesh<TParam> & mesh, const T & aimsmesh)
    {
////      detail::convert_mesh_1(aimsmesh, mesh);  //// in meshUtils.h
      detail::convert_aimsmeshTomesh1(aimsmesh, mesh);
    }

    template < typename TParam, typename T >
    typename boost::enable_if<boost::is_pointer<typename TParam::FaceIndex::value_type> >::type
    convertAimsMesh(T & aimsmesh, const Mesh<TParam> & mesh)
    {
      detail::convert_mesh_3(mesh, aimsmesh);
    }

    template < typename TParam, typename T >
    typename boost::disable_if<boost::is_pointer<typename TParam::FaceIndex::value_type> >::type
    convertAimsMesh(T & aimsmesh, const Mesh<TParam> & mesh)
    {
////      detail::convert_mesh_1(mesh, aimsmesh);  //// in meshUtils.h
      detail::convert_mesh1Toaimsmesh(mesh, aimsmesh);
    }

    template < typename TParam >
    class CastTo<Mesh<TParam>, AimsSurfaceTriangle >
     : public std::binary_function<Mesh<TParam> &, const AimsSurfaceTriangle &, void>
    {
    public:
      void
      operator()(Mesh<TParam> & mesh, const AimsSurfaceTriangle & aimsMesh) const
      {
        convertAimsMesh(mesh, aimsMesh);
      }
    };


    template < typename TParam >
    class CastTo<Mesh<TParam>, AimsSurface<3,Void> >
     : public std::binary_function<Mesh<TParam> &, const AimsSurface<3,Void> &, void>
    {
    public:
    
      void operator()(Mesh<TParam> & mesh, const AimsSurface<3,Void> & aimsMesh) const
      {
        convertAimsMesh(mesh, aimsMesh);
      }
    };

    template < typename T >
    class CastTo<AimsSurfaceTriangle, T >
     : public std::binary_function<AimsSurfaceTriangle &, const T &, void>
    {
    public:
      void
      operator()(AimsSurfaceTriangle & aimsMesh, const T & mesh) const
      {
        convertAimsMesh(aimsMesh, mesh);
      }
    };        


    /*
    template < typename TParam >
    class CastTo<AimsSurface<3,Void>, Mesh<TParam> >
     : public std::binary_function<AimsSurface<3,Void> &, const Mesh<TParam> &, void>
    {
    public:
      void
      operator()(AimsSurface<3,Void> & aimsMesh, const Mesh<TParam> & mesh) const
      {
        if (boost::is_pointer<typename TParam::FaceIndex::value_type>::value)
          detail::convert_mesh_3(mesh, aimsMesh);
        else
          detail::convert_mesh_1(mesh, aimsMesh);
      }
    };*/
    
    template < typename T >
    class CastTo<AimsSurface<3,Void>, T >
     : public std::binary_function<AimsSurface<3,Void> &, const T &, void>
    {
    public:
      void
      operator()(AimsSurface<3,Void> & aimsMesh, const T & mesh) const
      {
        convertAimsMesh(aimsMesh, mesh);
      }
    };
    
  } // namespace functor
  
  
  template < typename TVertex, typename TFace >
  void write_mesh(const std::vector<TVertex> & vertices, const std::vector<TFace> & faces, aims::Writer<AimsSurfaceTriangle> w)
  {
    AimsSurfaceTriangle s;
    Mesh1 tmp;
    getVertices(tmp) = vertices;
    getFaceIndices(tmp) = faces;
    convert(s, tmp);
    w.write(s);
  }

  template < typename TVertex, typename TFace >
  void write_mesh(const std::vector<TVertex> & vertices, const std::vector<TFace> & faces, const char * name)
  {
    write_mesh(vertices, faces, (aims::Writer<AimsSurfaceTriangle>(name)));
  }

  template < typename TMesh >
  void write_mesh(const TMesh & mesh, const char * name)
  {
    write_mesh(getVertices(mesh), getFaceIndices(mesh), name);
  }
  
  template < typename T >
  void writeTexture(const T & data, const char * name)
  {
    Texture1d t;
    t.reserve(data.size());
    for (std::size_t i = 0; i < data.size(); ++i) t.push_back(data[i]);
    aims::Writer<Texture1d> w(name);
    w.write(t);
  }
  
  template < typename T >
  void read_mesh(aims::Reader<T> & r, Mesh_N & mesh)
  {
    AimsSurfaceTriangle s;
    r.read(s);
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  
  template < typename T, typename TVertexCollection, typename TFaceCollection >
  void read_mesh(aims::Reader<T> const & r, TVertexCollection & vertices, TFaceCollection & faces)
  {
    T s;
    r.read(s);
    til::convert(vertices, til::getVertices(s));
    til::convert(faces, til::getFaceIndices(s));
  }

  template < typename T >
  void read_anyType(T & data, const std::string name)
  {
    aims::AnyTypeReader<AimsData<float> > r(name);
    r.read(data);
  }
  
  template < typename T >
  void aimswrite(const T & data, const std::string name)
  {
    aims::Writer<T> w(name);
    w.write(data);
  }
  
  template < typename TIterator1, typename TIterator2 >
  void convert_collection(TIterator1 begin1, TIterator1 end1, TIterator2 begin2)
  {
    for (; begin1 != end1; ++begin1, ++begin2)
    {
      convert2(*begin1).into(*begin2);
    }
  }
  
  template < typename VertexCollection, typename FaceCollection >
  void convert_mesh(AimsSurface<3, Void> const & mesh, VertexCollection & vertices, FaceCollection & faces)
  {
    vertices.resize(mesh.vertex().size());
    convert_collection(mesh.vertex().begin(), mesh.vertex().end(), vertices.begin());
    faces.resize(mesh.polygon().size());
    convert_collection(mesh.polygon().begin(), mesh.polygon().end(), faces.begin());
  } 

  template < typename VertexCollection, typename FaceCollection >
  void convert_mesh(VertexCollection const & vertices, FaceCollection const & faces, AimsSurface<3, Void> & mesh)
  {
    mesh.vertex().resize(vertices.size());
    convert_collection(vertices.begin(), vertices.end(), mesh.vertex().begin());
    mesh.polygon().resize(faces.size());
    convert_collection(faces.begin(), faces.end(), mesh.polygon().begin());
  }
  
} // namespace til

#undef TIL_COMMA

#endif /*AIMS_WRAP_H_*/
