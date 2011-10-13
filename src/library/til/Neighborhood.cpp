
#include "Neighborhood.h"

#define FOR_ALL_NEIGHBORS		\
for (i=-1; i<=1; ++i)			\
	for (j=-1; j<=1; ++j)		\
		for (k=-1; k<=1; ++k)	\




// Namespace 

namespace til {

Neighborhood create3DNeighborhood(double radius)
{
	Neighborhood res;
	int i, j, k;

	FOR_ALL_NEIGHBORS
	{
		if (norm(i,j,k) <= radius)
		{
			res.set(i,j,k);
		}
	}
	res.reset(0,0,0);

	return res;
}


// TODO: without the extern, VC complains!
// However, strange to have an extern here...!?

//extern TIL_API const Neighborhood N26(create3DNeighborhood(1.8));
extern TIL_API const Neighborhood N18(create3DNeighborhood(1.5));
extern TIL_API const Neighborhood N6(create3DNeighborhood(1.2));


} // namespace

