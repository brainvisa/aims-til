#ifndef _MESH_H_
#define _MESH_H_

// includes from STL
#include <list>
#include <vector>

// includes from BOOST
#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/type_traits.hpp>
using boost::shared_ptr;

// includes from TIL
//#include "til/Point.h"
//#include "til/Vector3.h"
#include "til/numeric_array.h"

// includes from TIL
//#include "til_common.h"
#include "globalTraits.h"
#include "miscUtils.h"

// declarations
#include "til/til_declarations.h"


namespace til
{

  /// A class to represent a very basic mesh, consisting of a set of vertices
  /// and a set of edges, represented as vertex indices.
  template < typename TMeshParam >
  class Mesh : public TMeshParam
  {
  public: // typedefs
  
    // typedefs from TMeshParam
    typedef Mesh<TMeshParam>                             Self;
    typedef typename TMeshParam::Vertex                  Vertex;
    typedef typename TMeshParam::VertexCollection        VertexCollection;
    typedef typename TMeshParam::FaceIndex               FaceIndex;
    typedef typename TMeshParam::FaceIndexCollection     FaceIndexCollection;
  
  public: // constructors & destructor
  
    /// Default constructor.
    /// NB: we should initialize the containers so that we never have to
    /// check whether they have to be initialized or not. This is at the 
    /// price of a small performance loss (namely, the initialization to
    /// something that may not be used), but should be negligeable.
  
    Mesh() : m_vertices(new VertexCollection()), m_faceIndices(new FaceIndexCollection()) {}
  
    /*  
    Mesh(const Mesh & mesh) : 
        m_vertices(new VertexCollection(*(mesh.m_vertices)))
      , m_faceIndices(new FaceIndexCollection(*(mesh.m_faceIndices)))
    */
    
  public: // set & get
  
    const VertexCollection & getVertices() const { return *m_vertices; }
    VertexCollection & getVertices() { return *m_vertices; }
    
    const FaceIndexCollection & getFaceIndices() const { return *m_faceIndices; }
    FaceIndexCollection & getFaceIndices() { return *m_faceIndices; }
    
  public: // functions
  
    void deepCopy(const Mesh & mesh)
    {
      // NB: the til:: is necessary to disambiguate with the this->deepCopy member function
      //til::deepCopy(mesh.m_vertices, m_vertices);
      //til::deepCopy(mesh.m_faceIndices, m_faceIndices);
      //m_vertices->resize(size(*(mesh.m_vertices)));

      // TODO: this whole deepCopy crap sucks. First, if the new container is much smaller, are we sure we
      // gain memory? Second, it does not work for pointers (see code below). An elegant solution has yet
      // to be found for pointer indices. E.G. by letting this code be handled by an index collection class.

      *m_vertices = *mesh.m_vertices;
      *m_faceIndices = *mesh.m_faceIndices;
  
      /*  
      //if (boost::is_same<typename value_type<FaceIndex>::type, std::size_t>::value)
      if (boost::is_same<typename value_type_of<FaceIndex>::type, std::size_t>::value)
      {
        *m_faceIndices = *mesh.m_faceIndices;
      }
      else
      {
        m_faceIndices->resize(size(*(mesh.m_faceIndices)));
        // *m_faceIndices = *(mesh.m_faceIndices);
        for (std::size_t i = 0; i < size(*m_faceIndices); ++i)
        {
          for (std::size_t j = 0; j < size((*m_faceIndices)[i]); ++j)
          {
            (*m_faceIndices)[i][j] = (*mesh.m_faceIndices)[i][j]-&(*mesh.m_vertices)[0] + &(*m_vertices)[0];
          }
        }
      }
      */
    }
    
    shared_ptr<VertexCollection> & vertices() { return m_vertices; }
    const shared_ptr<VertexCollection> & vertices() const { return m_vertices; }
    shared_ptr<FaceIndexCollection> & faces() { return m_faceIndices; }
    const shared_ptr<FaceIndexCollection> & faces() const { return m_faceIndices; }
      
  private: //data
  
    // Vertices
    shared_ptr<VertexCollection>     m_vertices;
    
    // Faces as vertex indices
    shared_ptr<FaceIndexCollection>  m_faceIndices;
  };
  
  /*
  template < typename TVertex, typename TAttribute >
  class XVertex
  {
    
  private: //data
    TVertex     m_vertex;
    TAttribute  m_attributes;
    //TAttribute  * m_attributes;
  };
  
  template < typename TMesh, typename TAttribute >
  class AddVertexAttribute
  {
  };
  
  
  
  template < typename TMeshParam, typename TAttributeParam >
  class AttributeMesh : public Mesh
  {
  public: // typedefs
    typedef typename TAttributeParam::VertexAttributeCollection  VertexAttributeCollection;
    typedef typename TAttributeParam::FaceAttributeCollection    FaceAttributeCollection;
    typedef typename TAttributeParam::EdgeAttributeCollection    EdgeAttributeCollection;
    
  public: // constructors & destructors
  
  private: // data
  
    shared_ptr< VertexAttributeCollection >  m_vertexAttributes;
    shared_ptr< FaceAttributeCollection >    m_faceAttributes;
    shared_ptr< EdgeAttributeCollection >    m_edgeAttributes;
  };
  
  template < typename TVertexAttribute, typename TFaceAttribute, typename TEdgeAttribute, template <typename> TContainer = std::vector >
  struct MeshAttribute
  {
    typedef TVertexAttribute                VertexAttribute;
    typedef TFaceAttribute                  FaceAttribute;
    typedef TEdgeAttribute                  EdgeAttribute;
    typedef TContainer<VertexAttribute>     VertexAttributeCollection;
    typedef TContainer<FaceAttribute>       FaceAttributeCollection;
    typedef TContainer<EdgeAttribute>       EdgeAttributeCollection;
  };
  */ 
  
  
  
  
  /*
  class MyMesh
  {
  public: // typedefs
  
    // typedefs from TMeshParam
    typedef Point<float,3>                      Vertex;
    typedef std::vector<Vertex>                 VertexCollection;
    typedef boost::array<std::size_t, 3>        FaceIndex;
    typedef std::vector<FaceIndex>              FaceIndexCollection;  
  
  public: // constructors & destructor
  
    /// Default constructor.
    /// NB: we should initialize the containers so that we never have to
    /// check whether they have to be initialized or not. This is at the 
    /// price of a small performance loss (namely, the initialization to
    /// something that may not be used), but should be negligeable.
    //Mesh() : m_vertices(new VertexCollection()), m_faceIndices(new FaceIndexCollection()) {}
    MyMesh() {}
    
  public: // set & get
  
    / *
    const VertexCollection & getVertices() const { return *m_vertices; }
    VertexCollection & getVertices() { return *m_vertices; }
    
    const FaceIndexCollection & getFaceIndices() const { return *m_faceIndices; }
    FaceIndexCollection & getFaceIndices() { return *m_faceIndices; }
    * /
    
    const VertexCollection & getVertices() const { return m_vertices; }
    VertexCollection & getVertices() { return m_vertices; }
    
    const FaceIndexCollection & getFaceIndices() const { return m_faceIndices; }
    FaceIndexCollection & getFaceIndices() { return m_faceIndices; }
  
    
  private: //data
  
    // Vertices
    //boost::shared_ptr<VertexCollection>     m_vertices;
    VertexCollection     m_vertices;
    
    // Faces as vertex indices
    //boost::shared_ptr<FaceIndexCollection>  m_faceIndices;
    FaceIndexCollection  m_faceIndices;
  };
  */
  
  struct MeshParam1
  {
    //typedef Point<float,3>                      Vertex;
    typedef numeric_array<float,3>              Vertex;
    //typedef Point3df                            Vertex;
    //typedef Vector<float,3>                     Vertex;
    //typedef boost::array<float,3>                 Vertex;
    typedef std::vector<Vertex>                 VertexCollection;
    typedef numeric_array<std::size_t, 3>        FaceIndex;
    typedef std::vector<FaceIndex>              FaceIndexCollection;  
  };
  
  struct MeshParam2
  {
    typedef numeric_array<float,3>              Vertex;
    //typedef Point<float,3>                      Vertex;
    typedef std::vector<Vertex>                 VertexCollection;
    typedef boost::array<Vertex*, 3>            FaceIndex;
    typedef std::vector<FaceIndex>              FaceIndexCollection;  
  };
  
  typedef Mesh<MeshParam1> Mesh1;
  typedef Mesh<MeshParam2> Mesh2;
  
  
  /// The 'detail' namespace contains technical details of the library, which
  /// a library user should not need to use nor to know.
  namespace detail
  {
   
    /// A structure giving default template parameters for the Add_XXX_Attribute
    /// classes, for a given Mesh.
    template < class TParam >
    struct DefaultAttributes {};
   
    /// This specialization gives the default template arguments of the Add_XXX_Attribute
    /// classes for the library Mesh class.
    //TODO: right now, vector is used by default => should change to sth dependent
    // on TParam
    template <typename TParam>
    struct DefaultAttributes<Mesh<TParam> >
    {
      //typedef Vector<float,3>                                       Normal;
      typedef numeric_array<float,3>                                       Normal;
      typedef std::vector<Normal>                                   NormalCollection;
      typedef std::vector<typename TParam::FaceIndex::value_type>   NeighborIndex;
      typedef std::vector<NeighborIndex>                            NeighborIndexCollection;
    };
  
    /* 
    struct AddNormalAttributeParams1
    {
      typedef Vector<float,3>       Normal;
      typedef std::vector<Normal>   NormalCollection;
    };
    */
    
    /// This class enhance a mesh class with a normal vector attribute.
    /// The first template argument is the mesh class we want to upgrade.
    /// The second template argument is the parameters for this upgrade, giving
    /// for exemple the type of the normal vectors.
    template < class TMesh, typename TParam = DefaultAttributes<TMesh> >
    class AddNormalAttribute : public TMesh
    {
    public: // typedefs
    
      typedef typename TParam::Normal             Normal;
      typedef typename TParam::NormalCollaction   NormalCollection;
      
    public: // constructors & destructor
    
      AddNormalAttribute() : TMesh() {}
      // NB: We have to add explicit here, because otherwise it's a source of bugs which
      // are complicated to track. Without explicit, the mother can be passed silently
      // as a children during a function call. This kind of precaution should
      // always been taken, especially when a constructor exists with a
      // parent as the only argument.
      explicit AddNormalAttribute(const TMesh & mesh) : TMesh(mesh) {}
      AddNormalAttribute(const TMesh & mesh, 
        shared_ptr<NormalCollection> & ni) : TMesh(mesh), m_normals(ni) {}
  
    public: // set & get
    
      const NormalCollection & getNormals() const { return *m_normals; }
      NormalCollection & getNormals() { return *m_normals; }
      
    private: // data
    
      // Neighbors as vertex indices
      shared_ptr<NormalCollection> m_normals;
    };
    
    /*
    struct AddNeighborAttributeParam1
    {
      typedef std::vector<std::size_t>    NeighborIndex;
      typedef std::vector<NeighborIndex>  NeighborIndexCollection;
    };
    */
    
    /// This class enhance a mesh class with a neighbor index attribute.
    /// The first template argument is the mesh class we want to upgrade.
    /// The second template argument is the parameters for this upgrade, giving
    /// for exemple the type of the neigbor indices.
    template < class TMesh, typename TParam = DefaultAttributes<TMesh> >
    class AddNeighborIndexAttribute : public TMesh
    {
    public: // typedefs
    
      typedef typename TParam::NeighborIndex            NeighborIndex;
      typedef typename TParam::NeighborIndexCollection  NeighborIndexCollection;
      
    public: // constructors & destructor
    
      AddNeighborIndexAttribute() : TMesh() {}
      explicit AddNeighborIndexAttribute(const TMesh & mesh) : TMesh(mesh) {}
      AddNeighborIndexAttribute(const TMesh & mesh, 
        shared_ptr<NeighborIndexCollection> & ni) : TMesh(mesh), m_neighborIndices(ni) {}
  
    public: // set & get
    
      const NeighborIndexCollection & getNeighborIndices() const { return *m_neighborIndices; }
      NeighborIndexCollection & getNeighborIndices() { return *m_neighborIndices; }
      
    private: // data
    
      // Neighbors as vertex indices
      shared_ptr<NeighborIndexCollection> m_neighborIndices;
    };
  }
  
  typedef detail::AddNeighborIndexAttribute<Mesh1>  Mesh_N;
  typedef detail::AddNormalAttribute<Mesh1>         Mesh_No;
  typedef detail::AddNormalAttribute<Mesh_N>        Mesh_NNo;
  
  typedef detail::AddNeighborIndexAttribute<Mesh2>  Mesh2_N;
  typedef detail::AddNormalAttribute<Mesh2>         Mesh2_No;
  typedef detail::AddNormalAttribute<Mesh2_N>       Mesh2_NNo;



  /*
  template < typename TIterator, typename TAccessor >
  class subiterator
  {
  };

  template < typename TNodeIterator >
  class vertex_iterator
  {
  public: // constructors
    vertex_iterator() : m_i() {}
    vertex_iterator(TNodeIterator i) : m_i(i) {}
  public: // operatoooors
    void operator++() { ++m_i; }
    vertex operator*() { return m_i->pos(); }
  private: // data
    TNodeIterator m_i;
  };
  */


} // namespace til

#endif //_MESH_H_
