

namespace til
{
  
  //---------------------------------------------------------------------------

  namespace detail
  {

    //-------------------------------------------------------------------------

    // Exported out of get_n_neighborhood_step to make sure that neigh_n is const so
    // that if we modify it, we still don't screw up this part: we should not modify
    // neigh_n while we're iterating on it, because this is a set!
    // TODO: this could actually be some kind of dilation on a graph.
    template < typename NeighborCollection >
    std::set<std::size_t>
    get_n_neighborhood_unitstep
    (
      std::set<std::size_t> const & neigh_n,
      NeighborCollection const & neighs
    )
    {
      typedef std::set<std::size_t> NeighN;
      NeighN new_neigh_n = neigh_n;
      // for all neighbors
      for (typename NeighN::const_iterator in = neigh_n.begin(); in != neigh_n.end(); ++in)
      //for (std::size_t j = 0; j < neigh_n.size(); ++j)
      {
        // for all neighbors of neighbors
        //for (std::size_t k = 0; k < neighs[neigh_n[j]].size(); ++k)
        for (std::size_t k = 0; k < neighs[*in].size(); ++k)
        {
          // add neighbor's neighbors to the neighborhood
          //new_neigh_n.insert(neighs[neigh_n[j]][k]);
          new_neigh_n.insert(neighs[*in][k]);
        }
      }
      return new_neigh_n;
    }
          
    //-------------------------------------------------------------------------
    
    template < typename NeighborCollection >
    void
    get_n_neighborhood_step
    (
      NeighborCollection const & neigh,
      std::vector<std::set<std::size_t> > & neigh_n
    )
    {
      assert(neigh.size() == neigh_n.size());
      // for all vertices
      for (std::size_t i = 0; i < neigh_n.size(); ++i)
      {
        neigh_n[i] = get_n_neighborhood_unitstep(neigh_n[i], neigh);
      }
    }

    //-------------------------------------------------------------------------

  } // namespace detail

  //---------------------------------------------------------------------------

  template < typename NeighborCollection, typename NeighborCollectionOut >
  void
  get_n_neighborhood
  (
    NeighborCollection const & neigh,     ///< [input] The 1-neighborhoods
    NeighborCollectionOut & neigh_n_out,  ///< [output] The computed n-neighborhood 
    unsigned int n                        ///< [input] the desired neighborhood width
  )
  {
    std::vector<std::set<std::size_t> > neigh_n;
    neigh_n.resize(neigh.size());
    // initialize n-neighborhood by including the point itself
    for (std::size_t i = 0; i < neigh_n.size(); ++i)
    {
      neigh_n[i].insert(i);
    }
    // iterate
    // TODO: clearly, this is not optimal for large neighborhoods.
    for (unsigned int i = 0; i < n; ++i)
    {
      detail::get_n_neighborhood_step(neigh, neigh_n);
    }
    // finally, remove point
    for (std::size_t i = 0; i < neigh_n.size(); ++i)
    {
      neigh_n[i].erase(i);
    }
    neigh_n_out.resize(neigh_n.size());
    for (std::size_t i = 0; i < neigh_n.size(); ++i)
    {
      neigh_n_out[i].resize(neigh_n[i].size());
      std::copy(neigh_n[i].begin(), neigh_n[i].end(), neigh_n_out[i].begin());
    }
  }

  //---------------------------------------------------------------------------

  template < typename TFaceCollection >
  shared_ptr<std::vector<std::vector<typename TFaceCollection::value_type::value_type> > >
  circular_neighborhoods(TFaceCollection const & faces, std::size_t nVertices)
  {
    typedef typename TFaceCollection::value_type::value_type VertexIndex;
    typedef std::list<VertexIndex> Neighborhood;
  
    // Allocate the result -- it should have as many elements as
    // the number of vertices of the mesh.
    std::vector<Neighborhood> res(nVertices);
    // To check if the circular neighbor of a point is finished
    // NB: we don't use a vector<bool> as it is apparently quite weird...
    std::vector<unsigned char> isFinished(nVertices, 0);
    // To check if everything is finished
    bool allFinished = false;
    
    // NB: the strategy here is quite dumb; basically, we loop through all triangles but will discard a neighbor
    // if it is not the one that fits at the end of our chain. Thus, we have to loop again and again through
    // this triangle list -- up to the maximum number of neighbors per point, that is, something like 14 times
    // for brain meshes. Plus, during the last pass, nothing should change, i.e. we'll do an extra pass for
    // nothing!
       
    // TODO: actually we could insert front and back to be hopefully twice as fast...  
    while (!allFinished)
    {
      allFinished = true;
  
      // for all faces
      for (typename TFaceCollection::const_iterator iFic = faces.begin(); iFic != faces.end(); ++iFic)
      {
        // for all couple of points on the face
        typename TFaceCollection::value_type::const_iterator iFi1 = iFic->begin();
        const_cyclic_iterator<typename TFaceCollection::value_type> iFi2(*iFic, iFi1+1);
  
        // NB: We decide that iFi2 is the main point, and iFi1 is the neighbor of iFi2 we want to add.
        // There is a reason why we do not also consider the reverse couple at the same time: this is to
        // ensure that all circular neighbors are turning in the same order, if the triangles do so.
        for (;iFi1 != iFic->end(); ++iFi1, ++iFi2)
        {
          //std::size_t i = getVertexNumber(mesh, *iFi2);
          std::size_t i = *iFi2;
  
          // Skip voxel if its neighborhood is already done
          if (isFinished[i])
          {
            continue;
          }
  
          // Check size of the neighborhood
          if (res[i].size() == 0)
          {
            // This is the first neighbor added: no need to check anything
            res[i].push_back(*iFi1);
            allFinished = false;
          }
          else
          {
            // We add the neighbor only if it is the "next" neighbor of the chain. That means
            // that the last point of the neighbor chain has to be in the current face as well.
            // Note also that we have to end as well: if this next neighbor also happens to be
            // the first of the chain, we mark this voxel neighborhood as done.
            
            // Loop through all face vertices
            VertexIndex lastPoint = res[i].back();
            typename TFaceCollection::value_type::const_iterator iFaceVertex = iFic->begin();
            for (; iFaceVertex != iFic->end(); ++iFaceVertex)
            {
              if (*iFaceVertex == lastPoint)
              {
                // Check that the point we are about to add is not actually the first point
                // of the circular neighborhood.
                if (*iFi1 == res[i].front())
                {
                  isFinished[i] = 1;
                }
                else
                {
                  res[i].push_back(*iFi1);
                  allFinished = false;
                }
                break;
              }
            }
          }        
        }
      }
    }  
    
    //convert the result into a vector
    shared_ptr<std::vector<std::vector<VertexIndex> > > res2(new std::vector<std::vector<VertexIndex> >);
    allocate_sameSize(res, *res2);
    //loop(*res2, res, Convert());

    {
      using namespace til::expr;
      for (std::size_t i = 0; i < size(res); ++i)
      {
        detail::loop_xx(castTo(*_1, *_2), (*res2)[i], res[i]);      
      }
      //detail::loop_xx(castTo(*_1, *_2), *res2, res);
    }
    return res2;
  }

  //---------------------------------------------------------------------------

  template < typename T >
  void
  neighbors2edges
  (
    std::vector<std::vector<T> >  const & neighbors,
    std::vector<std::pair<T, T> >       & edges
  )
  {
    typedef std::pair<std::size_t, std::size_t> Edge;
    std::set<Edge> myedges;
#ifndef NDEBUG
    std::size_t nNeighbors = 0;
#endif
    for (std::size_t i = 0; i < neighbors.size(); ++i)
    {
#ifndef NDEBUG
      nNeighbors += neighbors.size();
#endif
      for (std::size_t j = 0; j < neighbors[i].size(); ++j)
      {
        myedges.insert(make_sorted_pair(i, neighbors[i][j]));
      }
    }
    // works if neighbors is symmetric only...
    assert(nNeighbors == 2 * myedges.size());
    // converting back into vector
    edges.resize(myedges.size());
    std::copy(myedges.begin(), myedges.end(), edges.begin());
  }

  //---------------------------------------------------------------------------

  template < typename T >
  void
  Neighboring_faces<T>::operator()
  (
    std::vector<std::vector<T> >        const &   cneighs
  , std::vector<numeric_array<T, 3> >   const &   faces
  , std::vector<std::vector<T> >              &   neighborfaces
  )
  {
    typedef numeric_array<T, 3> Vec3D;
    typedef std::map<Vec3D, T, til::Lexicographical_compare<Vec3D> > Face_indices;
    Face_indices face_indices;
    for (std::size_t i = 0; i < faces.size(); ++i)
    {
      face_indices[sorted_vector(faces[i])] = i;
    }
    neighborfaces.resize(cneighs.size());
    for (std::size_t i = 0; i < cneighs.size(); ++i)
    {
      // NB: this test should never fail, actually a well behaved neighborhood
      // should have a size >= 3.
      assert(cneighs[i].size() >= 2);
      neighborfaces[i].resize(cneighs[i].size());
      std::size_t nNeighs = cneighs[i].size();
      for (std::size_t j = 0; j < nNeighs; ++j)
      {
        std::size_t j2 = (j+1) % nNeighs;
        typename Face_indices::iterator pos;
        pos = face_indices.find(sorted_vector<Vec3D>(i, cneighs[i][j], cneighs[i][j2]));
        if (pos != face_indices.end())
        {        
          neighborfaces[i][j] = pos->second;
        }
        else
        {
          throw InconsistentArguments();
        }
      }
    }
  }

  //---------------------------------------------------------------------------
  
} // namespace til
