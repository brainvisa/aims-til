#ifndef TIL_FORCES_H_
#define TIL_FORCES_H_

// includes from STL
#include <cmath>
#include <vector>

// includes from TIL
#include "til/til_declarations.h"
#include "til/is_traits.h"
#include "til/numeric_array.h"
#include "til/traits.h"

// includes from TIL
//#include "declarations.h"
#include "globalTraits.h"

namespace til
{
  namespace functor
  {  
    /// Label class for Mesh neighborhood functors.
    /// A mesh neighborhood functor should have an operator() that is called
    /// for each couple (vertex, neighbor). It should also have a
    /// value Functor::PerVertex, which, if non zero, indicates that it has 
    /// another function perVertex that is called for each vertex.
    struct MeshNeighborhoodFunctor_label { };
    
    template < typename TPoint3D, typename TVector3DRes >
    struct LaplacianForce : public MeshNeighborhoodFunctor_label
    {  
      enum { PerVertex = 1 };
  
      void
      operator()(const TPoint3D & v, const TPoint3D &, TVector3DRes & res) const
      {
        res += v;
        // add(res,v);
        //NB: don't use this, as this cannot be overloaded (say between types)
        //res += v;
      }
  
      void
      perVertex(const TPoint3D & v0, TVector3DRes & res, std::size_t n) const
      {
        res *= 1.0/n;
        res -= v0;
        //mul(res, 1.0/n);
        //sub(res, v0);
        //NB: don't use this, as this cannot be overloaded (say between types)
        //res *= 1.0/n;
        //res -= v0;
      }
      
      void foo() {};
    };
    
    struct SpringEnergy
    {
      template < typename TVector3D, typename TPrecision  >
      void operator()(const TVector3D & v, const TVector3D & v0, const TPrecision & l0, TPrecision & res) const
      {
        res += square(dist<TPrecision>(v, v0) - l0);
      }
    };
    
    /// Spring force functor.
    template < typename TPoint3D, typename TVector3DRes >
    struct SpringForce : public MeshNeighborhoodFunctor_label
    {
      // The precision of the computations (float, double...) is deduced fro,
      // the precision of the returning force vector type.
      typedef typename til::value_type_of<TVector3DRes>::type precision_t;
  
      enum { PerVertex = 0 };
      
      /// Compute -(v-v0)/||v-v0|| * (||v-v0|| - l0)^2, store result in res.
      /// This corresponds to -grad(E) where E = 1/2.(||v-v0|| - l0)^2.
  /*    template < typename TVector3D, typename TPrecision, typename TVector3DRes  >
      typename boost::enable_if_c<
        til::is_3DVector<TVector3D>::value
      >::type*/
      void operator()(const TPoint3D & v, const TPoint3D & v0, const precision_t & l0, TVector3DRes & res) const
      {
        precision_t k = static_cast<precision_t>(
        std::sqrt(norm2(v[0]-v0[0],
                        v[1]-v0[1],
                        v[2]-v0[2])));
  
        // v and v0 are too close to determine a force direction: return 0
        if (k < 128 * std::numeric_limits<precision_t>::epsilon())
        {
          return;
        }
        k = (k-l0) / k;
        //sub(v0, v, res);
        //mul(res, k);
        
        res[0] += (v0[0]-v[0]) * k;
        res[1] += (v0[1]-v[1]) * k;
        res[2] += (v0[2]-v[2]) * k;
      }
    };
  }
  
  
  /// A class to compute spring forces on a mesh.
  /// The first template argument corresponds to the mesh type, which should have
  /// vertex neighbor indices. The second template argument corresponds to the
  /// precision floating point operations should be performed on (typically float
  /// or double -- TODO: do we actually want to enforce TPrecision to be a numerical
  /// floating type?).
  template < class TMesh, typename TForceVector >
  class SpringForce
  {  
  public: // typedefs
    
    // The precision of the computations (float, double...) is deduced fro,
    // the precision of the returning force vector type.
    typedef typename til::value_type_of<TForceVector>::type precision_t;
    
  public: // functions
    
    /// Set the rest lengths of the springs as the length of the edges
    /// of the argument mesh.
    void initializeLengths(const TMesh & mesh)
    {
      getEdgeLengths(mesh, m_length0);
    }
  
    /// Computes the forces on the mesh, given the initial length already stored
    /// in the class, and store it in second argument.
    void getForces
    (
     const TMesh & mesh, 
     //std::vector<std::vector<typename ChangePrecision<typename MeshTraits<TMesh>::Vertex, TPrecision>::type > > & forces
     std::vector<TForceVector> & forces
    )
    {
      if (size(m_length0) == 0)
      {
        std::runtime_error("SpringForce::length0 uninitialized");
      }
      // for all vertices
      for_each_neighbors(mesh, m_length0, forces, functor::SpringForce<typename MeshTraits<TMesh>::Vertex, TForceVector>());
    }
  
  public: // set & get
  
    const std::vector<std::vector<precision_t> > & getLengths() const { return m_length0; }
  
  private: // check
  
    // Check that TMesh has neighbor indices
    typedef typename boost::enable_if_c<MeshTraits<TMesh>::has_neighbor_indices>::type CheckHasNeighborIndices;
  
  private: // data
  
    std::vector<std::vector<precision_t> > m_length0;
  };
  
  
  template < typename TForceVector, class TMesh >
  void
  laplacianForces(const TMesh & mesh, std::vector<TForceVector> & forces)
  {
    // The precision of the computations (float, double...) is deduced fro,
    // the precision of the returning force vector type.
    typedef typename til::precision<TForceVector>::type precision_t;
  
    // Allocate output with the same structure of neighbor indices of mesh
    // if the main container has the wrong size.
    if (size(forces) != size(getVertices(mesh)))
    {
      forces.resize(size(getVertices(mesh)));
    }
  
    //functor::LaplacianForce<typename MeshTraits<TMesh>::Vertex, TForceVector> lf;
    functor::LaplacianForce<numeric_array<precision_t,3>, numeric_array<precision_t,3> > lf;
    // for all vertices
    typename MeshTraits<TMesh>::VertexCollection                          v         = getVertices(mesh);
    typename MeshTraits<TMesh>::VertexCollection::const_iterator          iV        = v.begin();
    typename MeshTraits<TMesh>::NeighborIndexCollection::const_iterator   iNic      = getNeighborIndices(mesh).begin();
    typename std::vector<numeric_array<precision_t,3> >::iterator         iC2       = forces.begin();
    for (; iNic != getNeighborIndices(mesh).end(); ++iNic, ++iC2, ++iV)
    {
      // for all neighbors of current vertex
      for (typename MeshTraits<TMesh>::NeighborIndex::const_iterator iNi = iNic->begin(); iNi != iNic->end(); ++iNi)
      {
        lf(getVertexNeighbor(mesh, *iNi), *iV, *iC2);
      }
      lf.perVertex(*iV, *iC2, iNic->size());
    }
  }
} // namespace til

#endif //_FORCES_H_
