#ifndef TIL_POINTLIST_H
#define TIL_POINTLIST_H


// Standard library includes

#include <vector>


// Local includes

#include "til_common.h"

#include "SmartObject.h"
#include "Vector3.h"


// Namespace 

namespace til {


// A class to contain a list of point
// It is really just an std::vector<Vector3<T> > with SmartObject feature
// so that it can be used with Ptr and ConstPtr
/*
template < typename T >
class PointList : public SmartObject, public std::vector<Vector3<T> >
{
public:
	explicit PointList() : SmartObject(), std::vector<Vector3<T> >() {};
	explicit PointList(int i) : SmartObject(), std::vector<Vector3<T> >(i) {};
	explicit PointList(const std::vector<Vector3<T> > &pl) : SmartObject(), std::vector<Vector3<T> >(pl) {};
};
*/
} // namespace

#endif

