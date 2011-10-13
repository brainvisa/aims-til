#ifndef TIL_NEIGHBORHOOD_TOOLS
#define TIL_NEIGHBORHOOD_TOOLS

// includes from STL
#include <vector>

// includes from TIL library
#include "til/Neighborhood.h"
#include "til/numeric_array.h"


// Namespace 
namespace til
{

  // Keep only 'past' neighbord in the neighborhood.
  // i.e. the connectivity used by a forward filter (say 
  // a connected component filter).
  TIL_API void forwardize(Neighborhood & nh);
  
  
  /// Pushes the neighbor coordinates in a container.
  /// NB: coordinates are pushed inside the container. The container is not cleared
  /// from whatever it might contain beforehand.
  TIL_API void getNeighbors(const Neighborhood &nh, std::vector<numeric_array<int,3> > & res);
  
  
  // Returns the connectivity of the neighborhood.
  TIL_API int connectivity(const Neighborhood & nh);
  
  
  // Creates a 2D neighborhood in 3D, orthogonal to the vector (nx, ny, nz)
  // The normal (nx, ny, nz) should be binary, i.e. each componant should
  // be either 0 or one
  // A certain connectivity can be required. If set to 0, then maximal connectivity
  // is assumed
  TIL_API void create2DNeighborhood(Neighborhood & nh, int nx, int ny, int nz, int connectivity = 0);


} // namespace til

#endif



