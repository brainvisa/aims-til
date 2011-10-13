#ifndef TRIANGLEMESHGEODESICMAP_H_
#define TRIANGLEMESHGEODESICMAP_H_

// includes from STL
#include <queue>
#include <vector>

// includes from TIL
#include "til/sparse_vector.h"

// includes from TIL
#include "cyclic_iterator.h"
#include "fraction.h"
#include "MeshTraits.h"
#include "meshUtils.h"
#include "miscUtils.h"
#include "poly_solver.h"

// declarations
#include "til/til_declarations.h"

namespace til
{
  // TODO: rename as plugin
  namespace ghost 
  {
  
    struct GMapStop_None
    {
      // Only one function instead of two (warn(I,M) and proceed(void)) so that we don't have
      // to store anything if not needed and then everything can be inlined.
      template < typename TGMap >
      bool proceed(std::size_t, TGMap &) { return true; }      
    };
  
    template < typename TPrec >
    struct GMapStop_AboveThreshold
    {
      //GMapStop_AboveThreshold() : m_threshold() {}
      GMapStop_AboveThreshold(TPrec threshold) : m_threshold(threshold) {}
      
      template < typename TGMap >
      typename boost::enable_if<boost::is_same<TPrec, typename TGMap::precision_type>, bool>::type
      proceed(std::size_t i, TGMap & gmap) const
      {
        return (*(gmap.distanceMap()))[i] < m_threshold;
      }
      
      TPrec m_threshold;
    };
    
    
    template < typename TIterator >
    struct GMapStop_MyPointsDone
    {
      GMapStop_MyPointsDone() : m_begin(), m_end() {}

      GMapStop_MyPointsDone(TIterator begin, TIterator end)
      {
        this->init(begin, end);
      }
      
      void init(TIterator begin, TIterator end)
      {
        m_begin = begin;
        m_end = end;
        m_count = std::distance(m_begin, m_end);
      }
            
      template < typename TGMap >
      bool proceed(std::size_t i, TGMap &)
      {
        if (std::find(m_begin, m_end, i) != m_end)
        {
          if (--m_count == 0) return false;
        }
        return true;
      }
      
      TIterator m_begin;
      TIterator m_end;
      int m_count;
    };
  }
  
  
  namespace policy 
  {
    /*
    template < template < typename, typename, typename, typename > class GM >
    struct GMapStop_DistAboveThreshold
    {
    };
    */
    
    /*
    template < typename T >
    struct VectorStorage
    {
      typename std::vector<T> Container;
      void init(Container & c, std::size_t size, T value)
      {
        c.resize(size);
        std::fill(c.begin(), c.end(), value);
      }
    };
    */
    
    //template < template <typename> class TContainer, typename TMesh, typename TPrec >
/*    template < template <typename> class TContainer, typename TPrec >
    struct GMap_DefaultStorage
    {
      typedef unsigned char      Label;
      typedef TContainer<Label>  LabelCollection;
      typedef TContainer<TPrec>  DistCollection;
    };
*/
    // created to replace GMap_DefaultStorage<til::sparse_vector, double> calls
    // not allowed in new gcc version
    struct GMap_DefaultStorage_sparse_vect_dbl
    {
      typedef unsigned char      Label;
      typedef double      TPrec;
      typedef til::sparse_vector<Label>  LabelCollection;
      typedef til::sparse_vector<TPrec>  DistCollection;
    };

  }

  
  
  struct DistanceMapLabels
  {
    enum {
      UNACTIVE = 0,
      DONE,
      ACTIVE
    };
  };
  
  //---------------------------------------------------------------------------
   
     //---------------------//
    //  Mesh_distance_map  //
   //---------------------//

  // TODO: Actually, this empty shell should better be named MeshGreedyGrowth, cause it's really what it is.
  //template < typename TMesh, typename TPrec, typename TStopGhost = ghost::GMapStop_None, typename TStoragePolicy = policy::GMap_DefaultStorage<std::vector, TMesh, TPrec> >
  // TODO: replace reference to containers with accessors.

    // til::policy::GMap_DefaultStorage<til::sparse_vector, TPrec > replaced by GMap_DefaultStorage_sparse_vect_dbl
  template < typename TVertices, typename TNeighborhoods, typename TPrec, typename TStopGhost = ghost::GMapStop_None, typename TStoragePolicy = policy::GMap_DefaultStorage_sparse_vect_dbl >
  class Mesh_distance_map
    : public DistanceMapLabels
  {
  private: // classes
   
   /*
    /// Comparison (operator>) of second member of pairs
    template < typename T1, typename T2 >
    struct PairComp
    {
      // is operator< on
      inline bool operator()(const std::pair<T1, T2> & p1, const std::pair<T1, T2> & p2) const
      //inline bool operator()(std::pair<T1, T2> p1, std::pair<T1, T2> p2) const
      {
        return p1.second > p2.second;
      }
    };
    */
  
  public: // typedefs
  
    typedef Mesh_distance_map<TVertices, TNeighborhoods, TPrec, TStopGhost, TStoragePolicy> Self;
  
    typedef TPrec precision_type;
  
    // import typdefs from policies
    typedef typename TStoragePolicy::DistCollection             DistCollection;
    typedef typename TStoragePolicy::Label                      Label;
    typedef typename TStoragePolicy::LabelCollection            LabelCollection;
  
    //typedef typename MeshTraits<TMesh>::FaceIndex::value_type   VertexIndex;
    //typedef typename MeshTraits<TMesh>::FaceIndex::value_type   VertexIndex;

    typedef typename TNeighborhoods::value_type                 Neighborhood;
    
    typedef std::size_t                                         VertexIndex;
    //typedef std::vector<VertexIndex>                            Neighborhood;
    //typedef std::vector<VertexIndex>                            Neighborhood;
    typedef std::priority_queue<
      std::pair<VertexIndex, TPrec>,
      std::vector<std::pair<VertexIndex, TPrec> >,
      Greater_Pair2<VertexIndex, TPrec> >                       Queue;
  
    
  public: // constructors & destructor
  
    Mesh_distance_map(const TVertices & vertices, const TNeighborhoods & neighc);
  
    Mesh_distance_map(const TVertices & vertices, const TNeighborhoods & neighc, const TStopGhost & sp);
      
    virtual ~Mesh_distance_map() {}
    
  public: // set & get
  
    TStopGhost & stopGhost() { return m_stopGhost; }
    const TStopGhost & stopGhost() const { return m_stopGhost; }
    
    /// Get computed distance map.
    /// Important warning: the raw computation table is returned. If you did not use the default stop
    /// policy, it is very likely that this raw distance map contains bad estimates at extra points, used
    /// as intermediate stage for internal computations. So, you should be careful that you collect
    /// also the point labels and consider the distance only of those points labeled as "DONE".
    shared_ptr<const DistCollection>  distanceMap() const { return m_pDist; }
    shared_ptr<const LabelCollection> labels() const { return m_pLabel; }
    shared_ptr<DistCollection>  distanceMap() { return m_pDist; }
    shared_ptr<LabelCollection> labels() { return m_pLabel; }
  
  public: // functions
  
    /// Initialize to compute the geodesic distance to a single point.
    void init(VertexIndex iVertex);
  
    /// Initialize to compute the geodesic distance to a list a points with predefined distance values.
    /// This is actually very helpful to approximate the distance to objects that are NOT mesh points.
    /// For exemple, the distance to a point that is inside a triangle, or even a line or a curve.
    // NB: we template on TVertexIndex to be able to tackle const Vertex* kind of stuff...
    template < typename TVertexIndex >
    void init(std::vector<TVertexIndex> & startPoints, std::vector<TPrec> & startDist);
  
    void process();
    
  private: // functions
  
    void _init();
  
  protected: // pure virtual functions
  
    virtual TPrec distanceEstimate(VertexIndex iVertex) = 0;
  
  protected: // data
    //const TMesh & m_mesh;
    const TVertices & m_vertices;
    //const TFaces & m_faces;
    const TNeighborhoods & m_neighc;
    // contains the status of vertices
    shared_ptr<LabelCollection> m_pLabel;
    // contains the geodesic distance, i.e. the result
    shared_ptr<DistCollection>  m_pDist;
    // Circular neighborhoods
    //shared_ptr<std::vector<Neighborhood> > m_pNeigh;
    TPrec F;  
    Queue m_queue;
    TStopGhost m_stopGhost;
    bool m_allDone;
  };

  //---------------------------------------------------------------------------
   
     //----------------------//
    //  Graph_distance_map  //
   //----------------------//

  //template < typename TMesh, typename TPrec, typename TStopGhost = ghost::GMapStop_None, typename TStoragePolicy = policy::GMap_DefaultStorage<std::vector, TMesh, TPrec> >
  // til::policy::GMap_DefaultStorage<til::sparse_vector, TPrec > replaced by GMap_DefaultStorage_sparse_vect_dbl
  template < typename TVertices, typename TNeighborhoods, typename TPrec, typename TStopGhost = ghost::GMapStop_None, typename TStoragePolicy = policy::GMap_DefaultStorage_sparse_vect_dbl >
  class Graph_distance_map
    //: public Mesh_distance_map<TMesh, TPrec, TStopGhost, TStoragePolicy>
    : public Mesh_distance_map<TVertices, TNeighborhoods, TPrec, TStopGhost, TStoragePolicy>
  {
  public: // typedefs

    //typedef Mesh_distance_map<TMesh, TPrec, TStopGhost, TStoragePolicy> Base;
    typedef Mesh_distance_map<TVertices, TNeighborhoods, TPrec, TStopGhost, TStoragePolicy> Base;
    typedef typename Base::Neighborhood Neighborhood;
    typedef typename Base::VertexIndex VertexIndex;
    
  public: // constructors & destructor
  
    Graph_distance_map(const TVertices & vertices, const TNeighborhoods & neighc)
      : Base(vertices, neighc) {}
    Graph_distance_map(const TVertices & vertices, const TNeighborhoods & neighc, const TStopGhost & sp)
      : Base(vertices, neighc, sp) {}
    
  private: // functions
    TPrec distanceEstimate(VertexIndex iVertex);
  };


  //---------------------------------------------------------------------------
   
     //------------------------------//
    //  Triangle_mesh_geodesic_map  //
   //------------------------------//

  // NB: the neighborhood have to be circular now. TODO: enforce that.
  //template < typename TMesh, typename TPrec, typename TStopGhost = ghost::GMapStop_None, typename TStoragePolicy = policy::GMap_DefaultStorage<std::vector, TMesh, TPrec> >
  // til::policy::GMap_DefaultStorage<til::sparse_vector, TPrec > replaced by GMap_DefaultStorage_sparse_vect_dbl
  template < typename TVertices, typename TCircularNeighborhood, typename TPrec, typename TStopGhost = ghost::GMapStop_None, typename TStoragePolicy = policy::GMap_DefaultStorage_sparse_vect_dbl >
  class Triangle_mesh_geodesic_map
    //: public Mesh_distance_map<TMesh, TPrec, TStopGhost, TStoragePolicy>
    : public Mesh_distance_map<TVertices, TCircularNeighborhood, TPrec, TStopGhost, TStoragePolicy>
  {
  public: // typedefs
  
    typedef Mesh_distance_map<TVertices, TCircularNeighborhood, TPrec, TStopGhost, TStoragePolicy> Base;
    typedef typename Base::Neighborhood Neighborhood;
    typedef typename Base::VertexIndex VertexIndex;
    
  public: // constructors & destructor
  
    Triangle_mesh_geodesic_map(const TVertices & vertices, const TCircularNeighborhood & neighc)
      : Base(vertices, neighc) {}
    Triangle_mesh_geodesic_map(const TVertices & vertices, const TCircularNeighborhood & neighc, const TStopGhost & sp)
      : Base(vertices, neighc, sp) {}
    
  private: // functions
    TPrec distanceEstimate(VertexIndex iVertex);
  };

  //---------------------------------------------------------------------------
      
     //--------------------------------------//
    //  Voronoi_triangle_mesh_geodesic_map  //
   //--------------------------------------//

  // TODO: this should be redone to avoid code duplication with Triangle_mesh_geodesic_map.
  // It shouldn't be too hard, provided the API of distance maps is modified. The trick though is that
  // I guess I still want to have only one priority queues, not as many as there are elements...
  //template < typename TMesh, typename TPrec, typename TStopGhost = ghost::GMapStop_None, typename TStoragePolicy = policy::GMap_DefaultStorage<std::vector, TMesh, TPrec> >
  // til::policy::GMap_DefaultStorage<til::sparse_vector, TPrec > replaced by GMap_DefaultStorage_sparse_vect_dbl
  template < typename TVertices, typename TCircularNeighborhoods, typename TPrec, typename TStopGhost = ghost::GMapStop_None, typename TStoragePolicy = policy::GMap_DefaultStorage_sparse_vect_dbl >
  class Voronoi_triangle_mesh_geodesic_map
    : public Mesh_distance_map<TVertices, TCircularNeighborhoods, TPrec, TStopGhost, TStoragePolicy>
  {
  public: // typedefs
  
    typedef Mesh_distance_map<TVertices, TCircularNeighborhoods, TPrec, TStopGhost, TStoragePolicy> Base;
    typedef typename Base::Neighborhood Neighborhood;
    typedef typename Base::VertexIndex VertexIndex;
    
  public: // constructors & destructor
  
    Voronoi_triangle_mesh_geodesic_map(const TVertices & vertices, const TCircularNeighborhoods & neighc)
      : Base(vertices, neighc)
      , m_clusterLabel(new std::vector<unsigned int>(vertices.size(),0))
    {}
    Voronoi_triangle_mesh_geodesic_map(const TVertices & vertices, const TCircularNeighborhoods & neighc, const TStopGhost & sp)
      : Base(vertices, sp)
      , m_clusterLabel(new std::vector<unsigned int>(vertices.size(),0))
    {}
    
  public: // functions

    shared_ptr<std::vector<unsigned int> > clusterLabels() { return m_clusterLabel; }
  
    template < typename TVertexIndex >
    void init(std::vector<TVertexIndex> & startPoints, std::vector<TPrec> & startDist)
    {
      unsigned int count = 1;
      for (typename std::vector<TVertexIndex>::iterator i = startPoints.begin(); i != startPoints.end(); ++i)
      {
        //(*m_clusterLabel)[getVertexNumber(Base::m_mesh, *i)] = count++;
        (*m_clusterLabel)[*i] = count++;
      }
      this->Base::init(startPoints, startDist);
    }
  
  private: // functions
    TPrec distanceEstimate(VertexIndex iVertex);

  private: // data
    shared_ptr<std::vector<unsigned int> > m_clusterLabel;
  };

  //---------------------------------------------------------------------------

     //--------------------//
    //  helper functions  //
   //--------------------//

  /// Returns neighbors that are below a certain geodesic distance, along with
  /// their distance.
  template < typename TVertices, typename TFaces, typename TPrec >
  void
  distance_to_neighbors
  (
    TVertices const & vertices
  , TFaces const & faces
  , TPrec distance
  , std::vector<til::sparse_vector<TPrec> > & res
  );
  
} // namespace til


#include "triangle_mesh_geodesic_map.tpp"

#endif /*TRIANGLEMESHGEODESICMAP_H_*/
