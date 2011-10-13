#ifndef TIL_TRIANGLE_MESH_GEODESIC_MAP_TPP_
#define TIL_TRIANGLE_MESH_GEODESIC_MAP_TPP_

namespace til
{
  //---------------------------------------------------------------------------

  template < typename TVertices, typename TNeighborhoods, typename TPrec, typename TStopGhost, typename TStoragePolicy >
  Mesh_distance_map< TVertices, TNeighborhoods, TPrec, TStopGhost, TStoragePolicy >::Mesh_distance_map(const TVertices & vertices, const TNeighborhoods & neighc)
    : m_vertices(vertices)
    //, m_faces(faces)
    , m_neighc(neighc)
    , m_pLabel()
    , m_pDist()
    //, m_pNeigh()
    , m_stopGhost()
  {
  }

  //---------------------------------------------------------------------------
  
  template < typename TVertices, typename TNeighborhoods, typename TPrec, typename TStopGhost, typename TStoragePolicy >
  Mesh_distance_map< TVertices, TNeighborhoods, TPrec, TStopGhost, TStoragePolicy >::Mesh_distance_map(const TVertices & vertices, const TNeighborhoods & neighc, const TStopGhost & sp)
    : m_vertices(vertices)
    //, m_faces(faces)
    , m_neighc(neighc)
    , m_pLabel()
    , m_pDist()
    //, m_pNeigh()
    , m_stopGhost(sp)
  {
  }

  //---------------------------------------------------------------------------
  
  template < typename TVertices, typename TNeighborhoods, typename TPrec, typename TStopGhost, typename TStoragePolicy >
  void Mesh_distance_map< TVertices, TNeighborhoods, TPrec, TStopGhost, TStoragePolicy >::init(VertexIndex iVertex)
  {
    std::vector<VertexIndex> startPoints(1, iVertex);
    std::vector<TPrec> startDist(1, TPrec(0));
    this->init(startPoints, startDist);
  }

  //---------------------------------------------------------------------------
  
  template < typename TVertices, typename TNeighborhoods, typename TPrec, typename TStopGhost, typename TStoragePolicy >
  void Mesh_distance_map< TVertices, TNeighborhoods, TPrec, TStopGhost, TStoragePolicy >::process()
  {
    // Main loop
    VertexIndex iVertex;
    std::size_t i, j;
    TPrec d;
    
    while (!m_queue.empty())
    {
      // Get the vertex with the lowest possible distance
      iVertex = m_queue.top().first;
      m_queue.pop();
      //i = getVertexNumber(m_mesh, iVertex);
      i = iVertex;
      // Skip point if it is done already
      // This happens because when an active neighbor is updated with a new (lower) distance estimate
      // it is pushed again in the queue instead of being updated. priority_queue indeed does not allow
      // that. The thing is that if we want both voxel access and distance access, one has to be linear
      // in time. Dunno which is better, in any case it's kind of a pain.
      if ((*m_pLabel)[i] == DONE)
      {
        continue;
      }

      if (!m_stopGhost.proceed(i, *this)) break;

      // Mark it as DONE
      (*m_pLabel)[i] = DONE;
  
      // Loop on its neighbors
      //for (typename Neighborhood::const_iterator iNeighbor = (*m_pNeigh)[i].begin();
      //     iNeighbor != (*m_pNeigh)[i].end(); ++iNeighbor)
      for (typename Neighborhood::const_iterator iNeighbor = m_neighc[i].begin();
           iNeighbor != m_neighc[i].end(); ++iNeighbor)
      {
        //j = getVertexNumber(m_mesh, *iNeighbor);
        j = *iNeighbor;
        // Skip DONE points
        if ((*m_pLabel)[j] == DONE)
        {
          continue;
        }          
        // Get new estimate
        d = distanceEstimate(*iNeighbor);
        (*m_pDist)[j] = d;
        // Add them in the priority queue
        m_queue.push(std::make_pair(*iNeighbor, d));
      }
    }
    m_allDone = true;
  }  
  
  //---------------------------------------------------------------------------
  
  template < typename TVertices, typename TNeighborhoods, typename TPrec, typename TStopGhost, typename TStoragePolicy >
  template < typename TVertexIndex >
  void Mesh_distance_map< TVertices, TNeighborhoods, TPrec, TStopGhost, TStoragePolicy >::init(std::vector<TVertexIndex> & startPoints, std::vector<TPrec> & startDist)
  {
    // Check that input sizes are the same
    assert(size(startPoints) == size(startDist));
    
    this->_init();
    
    typename std::vector<TVertexIndex>::const_iterator i = startPoints.begin();
    typename std::vector<TPrec>::const_iterator d = startDist.begin();
    for (; i < startPoints.end(); ++i, ++d)
    {
      //std::size_t index = getVertexNumber(m_mesh, *i);
      std::size_t index = *i;
      // Label input points as DONE
      (*m_pLabel)[index] = DONE;
      // Copy input distances in our results
      (*m_pDist)[index] = *d;
      // Process starting points neighbors
      //typename MeshTraits<TMesh>::NeighborIndex nic = getNeighborIndices(m_mesh)[getVertexNumber(m_mesh, *i)];
      const Neighborhood & nic = m_neighc[*i];
      //typename MeshTraits<TMesh>::NeighborIndex::const_iterator iNeighbor = nic.begin();
      typename Neighborhood::const_iterator iNeighbor = nic.begin();
      for (; iNeighbor != nic.end(); ++iNeighbor)
      {
        // Label neighbors of starting points as ACTIVE
        // label[getVertexNumber(mesh, *iNeighbor)] = ACTIVE;
        // Compute their geodesic distance
        TPrec d = distanceEstimate(*iNeighbor);
        
        //(*m_pDist)[getVertexNumber(m_mesh, *iNeighbor)] = d;
        (*m_pDist)[*iNeighbor] = d;
        // Add them in the priority queue
        m_queue.push(std::make_pair(*iNeighbor, d));
      }
    }
  }

  //---------------------------------------------------------------------------
  
  template < typename TVertices, typename TNeighborhoods, typename TPrec, typename TStopGhost, typename TStoragePolicy >
  void Mesh_distance_map<TVertices, TNeighborhoods, TPrec, TStopGhost, TStoragePolicy>::
  _init()
  {
    if (m_pLabel == 0)
      m_pLabel = shared_ptr<LabelCollection>(new LabelCollection(m_vertices.size()));
    if (m_pDist == 0)
      m_pDist = shared_ptr<DistCollection>(new DistCollection(m_vertices.size()));

    assert(size(*m_pLabel) == m_vertices.size());
    assert(size(*m_pDist) == m_vertices.size());

    // set all points as being unprocessed
    fill(*m_pLabel, Label(UNACTIVE));
    // set initial distance to infinity -- or close enough ;)
    fill(*m_pDist, std::numeric_limits<TPrec>::max());
    m_allDone = false;      
  }

  //---------------------------------------------------------------------------
  
  template < typename TVertices, typename TNeighborhoods, typename TPrec, typename TStopGhost, typename TStoragePolicy >
  TPrec Graph_distance_map< TVertices, TNeighborhoods, TPrec, TStopGhost, TStoragePolicy >::distanceEstimate(VertexIndex iVertex)
  {
    //std::size_t k = getVertexNumber(Base::m_mesh, iVertex);
    std::size_t k = iVertex;
    TPrec res = (*Base::m_pDist)[k];
    const Neighborhood & nh = Base::m_neighc[k];
    for (typename Neighborhood::const_iterator iN = nh.begin(); iN != nh.end(); ++iN)
    {
      //if ((*Base::m_pLabel)[getVertexNumber(Base::m_mesh, *iN)] == Base::DONE)
      if ((*Base::m_pLabel)[*iN] == Base::DONE)
      {
        //res = min(res, (*Base::m_pDist)[getVertexNumber(Base::m_mesh, *iN)] + dist(getVertex(iVertex, Base::m_mesh), getVertex(*iN, Base::m_mesh), prec<TPrec>()));
        res = min(res, (*Base::m_pDist)[*iN] + dist(Base::m_vertices[iVertex], Base::m_vertices[*iN], prec<TPrec>()));
      }
    }
    return res;
  }
  
  //---------------------------------------------------------------------------
  
  template < typename TVertices, typename TCircularNeighborhoods, typename TPrec, typename TStopGhost, typename TStoragePolicy >
  TPrec Triangle_mesh_geodesic_map<TVertices, TCircularNeighborhoods, TPrec, TStopGhost, TStoragePolicy>::distanceEstimate(VertexIndex iVertex)
  {
    const TPrec F = TPrec(1);
    //std::size_t k = getVertexNumber(Base::m_mesh, iVertex);
    std::size_t k = iVertex;
    //Neighborhood &nh = (*Base::m_pNeigh)[k];
    const Neighborhood & nh = Base::m_neighc[k];
    //typename til::const_cyclic_iterator<Neighborhood> iN2((*m_neigh)[getVertexNumber(m_mesh, iVertex)], iN+1);
    math::PolySolver_real<TPrec, math::policy::InfinitySolutions_None> solver;

    TPrec res = (*Base::m_pDist)[k];
        
    // For each neighbor of iVertex
    //for (typename Neighborhood::const_iterator iN = nh.begin(); iN != nh.end(); ++iN, ++iN2)
    for (typename Neighborhood::const_iterator iN = nh.begin(); iN != nh.end(); ++iN)
    {
      //if ((*Base::m_pLabel)[getVertexNumber(Base::m_mesh, *iN)] == Base::DONE)
      if ((*Base::m_pLabel)[*iN] == Base::DONE)
      {
        // NB: Notations of "Computing Geodesic Paths on Manifolds", Kimmel & Sethian, 97
        TPrec b2 = dist2(Base::m_vertices[iVertex], Base::m_vertices[*iN], prec<TPrec>());
        std::size_t iA = *iN;
        TPrec dA = (*Base::m_pDist)[iA];
        
        typename Neighborhood::const_iterator iN2 = cyclic_advance(iN, nh);

        if ((*Base::m_pLabel)[*iN2] == Base::DONE)
        {
          //std::size_t iB = getVertexNumber(Base::m_mesh, *iN2);
          std::size_t iB = *iN2;
          TPrec dB = (*Base::m_pDist)[iB];            
          //TPrec a2 = dist2(getVertex(iVertex, Base::m_mesh), getVertex(*iN2, Base::m_mesh), prec<TPrec>());
          //TPrec c2 = dist2(getVertex(*iN, Base::m_mesh), getVertex(*iN2, Base::m_mesh), prec<TPrec>());
          TPrec a2 = dist2(Base::m_vertices[iVertex], Base::m_vertices[*iN2], prec<TPrec>());
          TPrec c2 = dist2(Base::m_vertices[*iN], Base::m_vertices[*iN2], prec<TPrec>());
          
          if (dB < dA)
          {
            std::swap(a2,b2);
            std::swap(iA, iB);
            std::swap(dB, dA);
          }
          
          // We compute these after the previous test to avoid having to swap them as well
          TPrec u = dB - dA;
          TPrec a = std::sqrt(a2);
          TPrec b = std::sqrt(b2);
          
          // Question: what do we want to do when div by zero, given that this never happened?
          TPrec ctheta = fraction<policy::ZeroByZero_Zero, TPrec>(a2+b2-c2, 2*a*b);
          
          if      (ctheta >  1) ctheta =  1;
          else if (ctheta < -1) ctheta = -1;

          solver.solve(c2, 2*u*b*(a*ctheta-b), b2*(u*u-F*F*a2*(1-ctheta*ctheta)));
          
          // TODO: this has to change
          if (solver.nsols() == 0)
          {
            res = min(res, b*F + dA, a*F + dB);
          }
          else
          {
            for (int i = 0; i < solver.nsols(); ++i)
            {
              TPrec t = solver.sols()[i];
              if ((u < t) && 
                  (t*a*ctheta < b*(t-u)) && 
                  (ctheta <= 0  || b*(t-u)*ctheta < a*t))
              {
                res = min(res, t + dA);
              }
              else
              {
                res = min(res, b*F + dA, a*F + dB);
              }
            }
          }
          
          // TODO: Optimized version still in the boxes
          /*
          Vector<float,3> base = getVertexNumber(m_mesh, *iN2)    - getVertexNumber(m_mesh, *iN);
          Vector<float,3> side = getVertexNumber(m_mesh, iVertex) - getVertexNumber(m_mesh, *iN)
          math::Poly2solver_real p(
            norm2(base),
            -2*u*dot(base,side),
            u*u*norm2(side) - F^2*crossnorm2(base, side)
          );

          for (int i = 0; i < p.nsols(); ++i)
          {
            if (fabs(u) < t && 
          }
          */
        }
        else
        {
          res = min(res, std::sqrt(b2)*F + dA);
        }
      }
    }
    return res;
  }
  
  //---------------------------------------------------------------------------
  
  template < typename TVertices, typename TCircularNeighborhoods, typename TPrec, typename TStopGhost, typename TStoragePolicy >
  TPrec Voronoi_triangle_mesh_geodesic_map<TVertices, TCircularNeighborhoods, TPrec, TStopGhost, TStoragePolicy>::distanceEstimate(VertexIndex iVertex)
  {
    const TPrec F = TPrec(1);
    //std::size_t k = getVertexNumber(Base::m_mesh, iVertex);
    std::size_t k = iVertex;
    const Neighborhood &nh = Base::m_neighc[k];
    //typename til::const_cyclic_iterator<Neighborhood> iN2((*m_neigh)[getVertexNumber(m_mesh, iVertex)], iN+1);
    math::PolySolver_real<TPrec, math::policy::InfinitySolutions_None> solver;

    TPrec res = (*Base::m_pDist)[k];
        
    // For each neighbor of iVertex
    //for (typename Neighborhood::const_iterator iN = nh.begin(); iN != nh.end(); ++iN, ++iN2)
    for (typename Neighborhood::const_iterator iN = nh.begin(); iN != nh.end(); ++iN)
    {
      //std::size_t iA = getVertexNumber(Base::m_mesh, *iN);
      std::size_t iA = *iN;
      if ((*Base::m_pLabel)[iA] == Base::DONE)
      {
        // NB: Notations of "Computing Geodesic Paths on Manifolds", Kimmel & Sethian, 97
        //TPrec b2 = dist2(getVertex(iVertex, Base::m_mesh), getVertex(*iN, Base::m_mesh), prec<TPrec>());
        TPrec b2 = dist2(Base::m_vertices[iVertex], Base::m_vertices[*iN], prec<TPrec>());
        TPrec dA = (*Base::m_pDist)[iA];
        
        typename Neighborhood::const_iterator iN2 = cyclic_advance(iN, nh);

        //std::size_t iB = getVertexNumber(Base::m_mesh, *iN2);
        std::size_t iB = *iN2;
        if ((*Base::m_pLabel)[iB] == Base::DONE && (*m_clusterLabel)[iB] == (*m_clusterLabel)[iA])
        {
          TPrec dB = (*Base::m_pDist)[iB];            
          TPrec a2 = dist2(Base::m_vertices[iVertex], Base::m_vertices[*iN2], prec<TPrec>());
          TPrec c2 = dist2(Base::m_vertices[*iN], Base::m_vertices[*iN2], prec<TPrec>());
          
          if (dB < dA)
          {
            std::swap(a2,b2);
            std::swap(iA, iB);
            std::swap(dB, dA);
          }
          
          // We compute these after the previous test to avoid having to swap them as well
          TPrec u = dB - dA;
          TPrec a = std::sqrt(a2);
          TPrec b = std::sqrt(b2);
          
          // Question: what do we want to do when div by zero, given that this never happened?
          TPrec ctheta = fraction<policy::ZeroByZero_Zero, TPrec>(a2+b2-c2, 2*a*b);
          
          if      (ctheta >  1) ctheta =  1;
          else if (ctheta < -1) ctheta = -1;

          solver.solve(c2, 2*u*b*(a*ctheta-b), b2*(u*u-F*F*a2*(1-ctheta*ctheta)));
          
          // TODO: this has to change
          if (solver.nsols() == 0)
          {
            TPrec tmp = min(b*F + dA, a*F + dB);
            if (res > tmp)
            {
              res = tmp;
              (*m_clusterLabel)[k] = (*m_clusterLabel)[iA];
            }
          }
          else
          {
            for (int i = 0; i < solver.nsols(); ++i)
            {
              TPrec t = solver.sols()[i];
              if ((u < t) && 
                  (t*a*ctheta < b*(t-u)) && 
                  (ctheta <= 0  || b*(t-u)*ctheta < a*t))
              {
                TPrec tmp = t + dA;
                if (res > tmp)
                {
                  res = tmp;
                  (*m_clusterLabel)[k] = (*m_clusterLabel)[iA];
                }
              }
              else
              {
                TPrec tmp = min(b*F + dA, a*F + dB);
                if (res > tmp)
                {
                  res = tmp;
                  (*m_clusterLabel)[k] = (*m_clusterLabel)[iA];                    
                }
              }
            }
          }            
        }
        else
        {
          TPrec tmp = std::sqrt(b2)*F + dA;
          if (res > tmp)
          {
            res = tmp;
            (*m_clusterLabel)[k] = (*m_clusterLabel)[iA];              
          }
        }
      }
    }
    return res;
  }
  
  //---------------------------------------------------------------------------

  template < typename TVertices, typename TNeighbors, typename TPrec >
  void
  distance_to_neighbors
  (
    TVertices const & vertices
  , TNeighbors const & neighbors
  , TPrec distance
  , std::vector<til::sparse_vector<TPrec> > & res
  )
  {
    std::size_t n = vertices.size();
    typedef std::vector<std::vector<std::size_t> > CNeighborhoods;

    typedef til::ghost::GMapStop_AboveThreshold<TPrec> Stop_ghost;
    Stop_ghost stop_ghost(distance);

    // til::policy::GMap_DefaultStorage<til::sparse_vector, TPrec > replaced by GMap_DefaultStorage_sparse_vect_dbl
    til::Triangle_mesh_geodesic_map<TVertices, CNeighborhoods, TPrec, Stop_ghost, til::policy::GMap_DefaultStorage_sparse_vect_dbl >

    //til::Graph_distance_map<TMesh, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage<til::sparse_vector, TMesh, double > > 
    //til::Graph_distance_map<typename TMesh::VertexCollection, CNeighborhoods, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage<til::sparse_vector, double > >
    //til::Graph_distance_map<TMesh, double, til::ghost::GMapStop_AboveThreshold<double> > 
      geomap(vertices, neighbors, stop_ghost);
    std::vector<std::size_t> startPoints(1);
    std::vector<TPrec> dist(1, TPrec(0.0));
    res.resize(n);
    shared_ptr<til::sparse_vector<unsigned char> > labels;
    for (std::size_t i = 0; i < n; ++i)
    {
      startPoints[0] = i;
      geomap.init(startPoints, dist);
      geomap.process();
      res[i] = *(geomap.distanceMap());
      labels = geomap.labels();
      // this removed temporary, boundary points that are generaly there when
      // a stopping ghost is used.
      for (typename til::sparse_vector<TPrec>::Map::iterator j = res[i].getMap().begin();
           j != res[i].getMap().end();)
      {
        if ((*labels)[j->first] != (unsigned char)(1))
        {
          res[i].getMap().erase(j++);
        }
        else
        {
          ++j;
        }
      }
    }
  }
  
} // namespace til

#endif /*TRIANGLE_MESH_GEODESIC_MAP_TPP_*/
