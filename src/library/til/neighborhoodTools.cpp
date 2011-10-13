
#include "til/neighborhoodTools.h"

// includes from STL
#include <cmath>


#define TIL_FOR_ALL_NEIGHBORS   \
for (i=-1; i<=1; ++i)           \
	for (j=-1; j<=1; ++j)         \
		for (k=-1; k<=1; ++k)       \



namespace til
{


	/// Tests whether neighbor (i,j,k) is a 'future' neighbor w.r.t raster scan order
	static inline bool isFutureNeighbor(int i, int j, int k)
	{
		return k>0 || (k==0 && j>0) || (k==0 && j==0 && i>0) ;
	}

/// Remove all 'future' neighbors of neighbor nh, in the raster scan sens
/// This makes it suitable for a 'forward' pass used in some algorithms such as 
/// connected componant labeling.
TIL_API void forwardize(Neighborhood &nh)
{
	int i, j, k;
	TIL_FOR_ALL_NEIGHBORS
	{
		if (isFutureNeighbor(i,j,k)) nh.remove(i,j,k);
	}
}

TIL_API void getNeighbors(const Neighborhood & nh, std::vector<numeric_array<int,3> > & res)
{
	int i, j, k;

	TIL_FOR_ALL_NEIGHBORS
	{
		if (nh.isNeighbor(i,j,k))
		{
			res.push_back(numeric_array<int,3>(i, j, k));
		}
	}
}

TIL_API int connectivity(const Neighborhood &nh)
{
	int res = 0;
	int i, j, k;

	// add one for any neighbor present
	TIL_FOR_ALL_NEIGHBORS
	{
		if (nh.isNeighbor(i,j,k))
		{
			++res;
		}
	}

	// subtract one if the center was there too
	if (nh.isNeighbor(0,0,0)) --res;

	return res;
}


static inline bool isCorner(int x, int y, int z)
{
	return std::abs(x)==1 &&
         std::abs(y)==1 &&
         std::abs(z)==1;
}


static inline bool isAValidDiscreteNormal(int x, int y, int z)
{
	return 
		std::abs(x) <= 1 &&
		std::abs(y) <= 1 &&
		std::abs(z) <= 1 &&
		std::abs(x) + std::abs(y) + std::abs(z) != 0;
}


TIL_API void create2DNeighborhood(Neighborhood &nh, int nx, int ny, int nz, int connectivity)
{

	if (!isAValidDiscreteNormal(nx, ny, nz))
	{
		throw std::invalid_argument("Invalid normal vector");
	}


	// initialize neighborhood
	nh.reset();


	int i, j, k;

	// Set to one voxels orthogonal to normal vector

	TIL_FOR_ALL_NEIGHBORS
	{
		// dot product between normal vector and current position
		if (i*nx+j*ny+k*nz == 0)
		{
			nh.set(i,j,k);
		}
	}

	// Remove the center;
	nh.remove(0,0,0);

	int count = til::connectivity(nh);

	// Check if we get the desired connectivity

	if ((connectivity == count) || (connectivity == 0))
	{
		// If so, we are done
		return;
	}

	// If not:
	// If the desired connectivity was greater, then that argument was invalid

	else if (connectivity > count)
	{
		throw std::invalid_argument("Connectivity too high");
	}

	// If it was smaller, we remove the corner points
	TIL_FOR_ALL_NEIGHBORS
	{
		if (isCorner(i, j, k))
		{
			nh.remove(i,j,k);
		}
	}


	// Check if we obtain the desired connectivity
	count = til::connectivity(nh);

	if (connectivity != count)
	{
		// If not, then that argument was invalid
		throw std::invalid_argument("Connectivity too high");
	}
}


} // namespace til

