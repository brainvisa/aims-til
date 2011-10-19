
// includes from STL
#include <algorithm>

// includes from AIMS
#include "aims/mesh/surfacegen.h"

// includes from TIL
#include "til/std_wrap.h"

// includes from TIL
#include "aims_wrap.h"
#include "MeshTraits.h"
#include "meshUtils.h"
#include "miscUtils.h"

using namespace aims;




void makeSphere_oneSubdivisionStep(const AimsSurfaceTriangle * in, AimsSurfaceTriangle *out, const Point3df & center, float radius)
{
	// convenient typedefs
	typedef AimsVector<uint, 2> Edge;
	typedef AimsVector<uint, 3> Face;
	
	// Preallocate necessary space
	out->vertex().reserve(2*in->vertex().size());
	out->polygon().reserve(4*in->polygon().size());
	
	// Start copying all vertices
	out->vertex().resize(in->vertex().size());
	std::copy(in->vertex().begin(), in->vertex().end(), out->vertex().begin());
	
	// Collect edges
	std::set<Edge> edges = til::getEdges(in);
	
	// Split all edges into two equal parts
	std::set<Edge>::const_iterator iEdge;
	std::map<Edge, int> index;
	for (iEdge = edges.begin(); iEdge != edges.end(); ++iEdge)
	{
		// compute the center of the edge
		Point3df middle = in->vertex()[(*iEdge)[0]] + in->vertex()[(*iEdge)[1]];
		middle /= 2.0;
	
		// Project this point back on the sphere
		middle -= center;
		middle *= radius / middle.norm();
		middle += center;
		
		// Put this point in new mesh and keep track of its index
		index[*iEdge] = out->vertex().size();
		out->vertex().push_back(middle);
	}
	
	// Split triangles into four smaller triangles, making sure that
	// vertex order in the new triangles is compatible with the original
	// triangle
	std::vector<Face>::const_iterator iFace;
	for (iFace = in->polygon().begin(); iFace != in->polygon().end(); ++iFace)
	{
		out->polygon().push_back(Face((*iFace)[0],
		                              index[til::sortedVector<AimsVector<uint,2> >((*iFace)[0], (*iFace)[1])],
		                              index[til::sortedVector<AimsVector<uint,2> >((*iFace)[0], (*iFace)[2])]));
		
		out->polygon().push_back(Face((*iFace)[1],
		                              index[til::sortedVector<AimsVector<uint,2> >((*iFace)[1], (*iFace)[2])],
		                              index[til::sortedVector<AimsVector<uint,2> >((*iFace)[1], (*iFace)[0])]));
		
		out->polygon().push_back(Face((*iFace)[2],
		                              index[til::sortedVector<AimsVector<uint,2> >((*iFace)[2], (*iFace)[0])],
		                              index[til::sortedVector<AimsVector<uint,2> >((*iFace)[2], (*iFace)[1])]));
		
		out->polygon().push_back(Face(index[til::sortedVector<AimsVector<uint,2> >((*iFace)[0], (*iFace)[1])],
		                              index[til::sortedVector<AimsVector<uint,2> >((*iFace)[1], (*iFace)[2])],
		                              index[til::sortedVector<AimsVector<uint,2> >((*iFace)[2], (*iFace)[0])]));
	}	
}

/// Creates a triangulated sphere by iterative subdivision of 
/// an icosahedron.
/// The total number of triangles is 3*4^(n+1)
AimsSurfaceTriangle * 
makeSphere
(
 const Point3df &center, ///< the center of the sphere
 float radius,           ///< the radius of the spere
 int iter                ///< the number of subdivision iterations
 )
{	
	// Generate an icosahedron
	AimsSurfaceTriangle * res = 
		SurfaceGenerator::icosahedron(center, radius);  
  
	// Iterate subdivision
  AimsSurfaceTriangle *tmp;
	for (int i=0; i < iter; ++i)
	{
    tmp = new AimsSurfaceTriangle();
		makeSphere_oneSubdivisionStep(res, tmp, center, radius);
    delete res;
    res = tmp;
	}
	return res;
}

namespace til
{
  /// Get all edges in mesh.
  std::set<AimsVector<uint, 2> >
  getEdges(const AimsSurfaceTriangle *surf)
  {
    
    std::set<AimsVector<uint, 2> > edges;
    std::vector<AimsVector<uint,3> >::const_iterator iFace;
    for (iFace = surf->polygon().begin(); 
         iFace != surf->polygon().end(); 
         ++iFace)
    {
      
      // We sort vertex by index to input only edges (a,b) where
      // a < b. This is to ensure that the same edge is not
      // input twice as (a,b) and (b,a).
      edges.insert(sortedVector<AimsVector<uint, 2> >((*iFace)[0], (*iFace)[1]));
      edges.insert(sortedVector<AimsVector<uint, 2> >((*iFace)[1], (*iFace)[2]));
      edges.insert(sortedVector<AimsVector<uint, 2> >((*iFace)[2], (*iFace)[0]));      
    }
    return edges;
  }
}