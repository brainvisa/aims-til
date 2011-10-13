#ifndef TIL_REGIONGROWING_H
#define TIL_REGIONGROWING_H


// Includes from TIL library
#include "til/til_common.h"
#include "til/functors.h"
#include "til/imageTools.h"
#include "til/Neighborhood.h"
#include "til/numeric_array.h"


//////////////////////////////////////////////////////////////////////////////
//
// Generic tools to do region growing.
//
// Region growing is split into two parts. The first part is the
// regionGrowing function, that does all the mechanics: finding neighbors
// (neighborhood is a parameter) and adding them to a segmentation image
// (another parameter) in a given color (parameter) starting from a
// list of given points (parameter).
//
// However, all the criteria for this mechanic, i.e. the logic of region
// growing, is put elsewhere in a Ghost class, which is the last parameter
// of the region growing function. Basically, the ghost class interacts with
// the region growing in three ways. (1), it tells the region growing whether
// to accept a point as being part of the region, via its test() member.
// (2), it tells the region growing whether to stop, via its stop() member.
// (3), it collects information from the region growing via its update
// member, which takes the position of any accepted point as an input.
//
// This separation allows minimal recoding. E.g. to do region growing for
// intensies below or above (or equal, or different, or...) from a certain
// threshold, one has only to code this tests in a region growing ghost,
// not to re-code an entire region-growing algorithm.
//
// Actually, for standard region growings, that provides only an acceptance
// criterion (no extra stopping criterion, no information collection), 
// a conveniance SimpleGhost is provided, that takes a PixelTest class as
// a parameter.
//
//
//////////////////////////////////////////////////////////////////////////////






namespace til
{

	typedef std::vector<numeric_array<int,3> > VoxelList;

template < class TExpr, class TImage >
class RegionGrowingExpr
{
public: // constructors & destructor

	RegionGrowingExpr(TImage & im, TExpr expr) :
		m_im(im), m_expr(expr) {};

public: // functions
	// nothing to update
	static void update(const numeric_array<int,3> & pos) {}
	// never stops
	static bool stop() { return false; }
	// use the pixel test to decide whether to add a value
	bool test(const numeric_array<int,3> & pos) const { return m_expr(&((*m_im)(pos))); }

private:
	
	TImage m_im;
	TExpr m_expr;
};

template < class PixelTest >
class SimpleGhost
{

public: // typedefs

	typedef typename PixelTest::TImage TImage;

public:
	
	SimpleGhost(TImage &im, const PixelTest & test)
	{
		m_im.shallowCopy(im);
		m_test = &test;
	}

	// nothing to update
	static void update(const numeric_array<int,3> &) {}

	// use the pixel test to decide whether to add a value
	bool test(const numeric_array<int,3> &pos) const { return m_test->compute(m_im, pos); }

	// never stops
	static bool stop() { return false; }

private:

	TImage m_im;
	const PixelTest *m_test;
};



template < class PixelTest, class Plugin >
class PluginGhost
{
public: // typedefs

	typedef typename PixelTest::TImage TImage;

public:
	
	PluginGhost(TImage & im, const PixelTest & test, Plugin & plugin)
	{
		m_im.shallowCopy(im);
		m_test = &test;
		m_plugin = &plugin;
	}

	// nothing to update
	void update(const numeric_array<int,3> &pos) { m_plugin->update(pos); }

	// use the pixel test to decide whether to add a value
	bool test(const numeric_array<int,3> &pos) const { return m_test->compute(m_im, pos); }

	// never stops
	static bool stop() { return false; }

private:

	TImage m_im;
	const PixelTest *m_test;
	Plugin *m_plugin;
};


/// If point at position 'pos' satisfy region-growing criteria, it is
/// added to the list of voxel 'vl' and also put in the image 'seg'
// TODO: Actually, should we not also mark points that are rejected,
// to avoid that they are selected again?
// TODO: change order: I think it is more natural if pos comes first
template < class TImage, class Ghost >
INLINE void
_addPoint
(
  VoxelList & vl,
  TImage & seg,
  const numeric_array<int,3> & pos,
  Ghost & ghost,
  typename TImage::value_type newColor
)
{
	if (contains(seg, pos) && seg(pos) != newColor && ghost.test(pos))
	{
		vl.push_back(pos);

		// NB: update could be directly linked in Ghost::test() and not
		// appear here...
		// NB: It is important that the update takes place before the seg value
		// is actually changed
		ghost.update(pos);

		seg(pos) = newColor;
	}
}

template < class TImage, class Ghost >
INLINE void _addPoint2
(
  VoxelList & vl,
  TImage & seg,
  const numeric_array<int,3> & pos,
  Ghost & ghost,
  typename TImage::value_type newColor
)
{
	if (seg(pos) != newColor && ghost.test(pos))
	{
		vl.push_back(pos);

		// NB: update could be directly linked in Ghost::test() and not
		// appear here...
		// NB: It is important that the update takes place before the seg value
		// is actually changed
		ghost.update(pos);

		seg(pos) = newColor;
	}
}

// TODO: there is a deep thinking to be done on this. E.g. I am not sure
// having the new color parameter here is right. It should probably be part
// of some policy.

template < class TImage, class Ghost >
std::auto_ptr<VoxelList>                                      ///< New boundary points
addNeighbors(TImage &seg,                                     ///< Segmentation image
             const std::vector<numeric_array<int,3> > &vl,    ///< Boundary points
             const std::vector<numeric_array<int,3> > &vnh,   ///< Neighborhood
             Ghost &ghost,                                    ///< Region growing criteria
             typename TImage::value_type newColor)            ///< New color
{
	std::auto_ptr<VoxelList> newVl(new VoxelList);
	// A factor of 3 is usually more than enough
	newVl->reserve(3 * vl.size());
	VoxelList::const_iterator iVl;

	std::vector<numeric_array<int,3> >::const_iterator iVnh;

	numeric_array<int,3> neighborCoord;

	for (iVl = vl.begin(); iVl != vl.end(); ++iVl)
	{
		for (iVnh = vnh.begin(); iVnh != vnh.end(); ++iVnh)
		{
			// Computes the coordinates of current neighbor
			neighborCoord = *iVl + *iVnh;
      //addTo(*iVl, *iVnh, neighborCoord);
			_addPoint(*(newVl.get()), seg, neighborCoord, ghost, newColor);
			
			if (ghost.stop()) goto endloop;
		}
	}

endloop:

	return newVl;
}

// Local macro to help writing function
// To be undefined after use
//	add_T<(i),(j),(k)>(*iPl, neighborCoord);								
#define ADD_NEIGHBORS2(i,j,k)                                                           \
if (nh.template isNeighbor<(i),(j),(k)>() && containsNeighbor<(i),(j),(k)>(iSeg))       \
{                                                                                       \
	addTo<(i),(j),(k)>(*iVl, neighborCoord);                                              \
	_addPoint2(*(newVl.get()), seg, neighborCoord, ghost, newColor);                      \
	if (ghost.stop()) break;                                                              \
}                                                                                       \


template < class TImage, class Ghost, class TNeighborhood >
std::auto_ptr<VoxelList>
addNeighbors2
(
  TImage & seg,
  const VoxelList & vl,
  const TNeighborhood & nh,
  Ghost & ghost,
  typename TImage::value_type newColor
)
{
	std::auto_ptr<VoxelList> newVl(new VoxelList);
	// A factor of 3 is usually more than enough
	newVl->reserve(3 * vl.size());
	VoxelList::const_iterator iVl;

	numeric_array<int,3> neighborCoord;

	typename Iterator<TImage>::Volumetric iSeg(seg);

	for (iVl = vl.begin(); iVl != vl.end(); ++iVl)
	{
		iSeg.setUnsafePos(*iVl);
		
		// First 6 neighbors
		ADD_NEIGHBORS2(-1, 0, 0)
		ADD_NEIGHBORS2( 0,-1, 0)
		ADD_NEIGHBORS2( 0, 0,-1)
		ADD_NEIGHBORS2(+1, 0, 0)
		ADD_NEIGHBORS2( 0,+1, 0)
		ADD_NEIGHBORS2( 0, 0,+1)

		// Next 12
		ADD_NEIGHBORS2(-1,-1, 0)
		ADD_NEIGHBORS2(-1, 0,-1)
		ADD_NEIGHBORS2( 0,-1,-1)
		ADD_NEIGHBORS2(+1,+1, 0)
		ADD_NEIGHBORS2(+1, 0,+1)
		ADD_NEIGHBORS2( 0,+1,+1)
		ADD_NEIGHBORS2(-1,+1, 0)
		ADD_NEIGHBORS2(-1, 0,+1)
		ADD_NEIGHBORS2( 0,-1,+1)
		ADD_NEIGHBORS2(+1,-1, 0)
		ADD_NEIGHBORS2(+1, 0,-1)
		ADD_NEIGHBORS2( 0,+1,-1)

		// Last 8
		ADD_NEIGHBORS2(-1,-1,-1)
		ADD_NEIGHBORS2(+1,+1,+1)
		ADD_NEIGHBORS2(+1,-1,-1)
		ADD_NEIGHBORS2(-1,+1,-1)
		ADD_NEIGHBORS2(-1,-1,+1)
		ADD_NEIGHBORS2(-1,+1,+1)
		ADD_NEIGHBORS2(+1,-1,+1)
		ADD_NEIGHBORS2(+1,+1,-1)
	}

	return newVl;
}

// Undefine local macro after use
#undef ADD_NEIGHBORS2


template < class TImage, class Ghost>
std::auto_ptr<VoxelList>
addSeeds
(
  TImage &seg,
  const VoxelList &vl,
  Ghost &ghost,
  typename TImage::value_type newColor
)
{
	std::auto_ptr<VoxelList> newVl(new VoxelList);
	VoxelList::const_iterator iVl;

	for (iVl = vl.begin(); iVl != vl.end(); ++iVl)
	{
		_addPoint(*(newVl.get()), seg, *iVl, ghost, newColor);
			
		if (ghost.stop()) break;
	}

	return newVl;
}


template < typename TImage, typename RegionGrowingGhost, typename TNeighborhood >
size_t regionGrowing2
(
  TImage & seg,
  const VoxelList & seeds,
  const TNeighborhood & nh,
  RegionGrowingGhost & ghost,
  typename TImage::value_type color
)
{
	// Initialize point list
	std::auto_ptr<VoxelList> vl;


	// Initialize the region growing
	// NB: Not all seeds are necesarily taken into account
	// They have to follow the same rules that apply to future
	// neighbors
	// Seeds that no not meet requirements are thrown out.
	// Seeds order might therefore be important, depending
	// on the ghost used.

	//seg.reset();
	vl = addSeeds(seg, seeds, ghost, color);

	// Counts the number of points in region
	size_t nPoints = vl->size();

	while (vl->size() > 0 && !ghost.stop())
	{
		vl = addNeighbors2(seg, *vl.get(), nh, ghost, color);
		nPoints += vl->size();
	}

	return nPoints;
}



// 
template < typename TImage, typename RegionGrowingGhost >
size_t regionGrowing
(
  TImage & seg,
  const VoxelList & seeds,
  const Neighborhood & nh,
  RegionGrowingGhost & ghost,
  typename TImage::value_type color
)
{
	// Initialize point list
	std::auto_ptr<VoxelList> vl;

	// Initialize the region growing
	// NB: Not all seeds are necesarily taken into account
	// They have to follow the same rules that apply to future
	// neighbors
	// Seeds that no not meet requirements are thrown out.
	// Seeds order might therefore be important, depending
	// on the ghost used.

	//seg.reset();
	vl = addSeeds(seg, seeds, ghost, color);

	// Initialize neighborhood
	std::vector<numeric_array<int,3> > vnh;
	getNeighbors(nh, vnh);


	// Counts the number of points in region
	size_t nPoints = vl->size();

	while (vl->size() > 0 && !ghost.stop())
	{
		vl = addNeighbors(seg, *(vl.get()), vnh, ghost, color);
		nPoints += vl->size();
	}

	return nPoints;
}


template < typename TImage, typename RegionGrowingGhost >
size_t regionGrowing
(
  TImage & seg,
  const numeric_array<int,3> & seed,
  const Neighborhood & nh,
  RegionGrowingGhost & ghost,
  typename TImage::value_type color
)
{
	VoxelList seeds;
	seeds.push_back(seed);
	return regionGrowing(seg, seeds, nh, ghost, color);
}


/*
// If seg = im, assumes that 'color' is not present in image

template < typename TImage, typename PixelTest, typename PlugIn = DoNothing >
int regionGrowing(Ptr<TImage> &seg,
				  ConstPtr<PointList<int> > &seeds,
				  const Neighborhood &nh,
				  const PixelTest &pixelTest,
				  typename TImage::value_type color,
				  PlugIn &plugin = DoNothing())
{
	// Initialize the segmentation
	seg.reset();
	dumpPointListInImage(seeds, seg, color);

	// Initialize neighborhood
	std::vector<numeric_array<int,3> > vnh;
	getNeighbors(nh, vnh);

	// Initialize point list
	Ptr<PointList<int> > pl = new PointList<int>(*seeds);

	// Counts the number of points in region
	int nPoints = pl->size();

	while (pl->size() > 0)
	{
		addNeighbors(seg, pl, vnh, pixelTest, color);
		plugin.update(pl);
		nPoints += pl->size();
	}

	return nPoints;
}
*/

} // namespace til

#endif

