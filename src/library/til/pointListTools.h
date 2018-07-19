#ifndef TIL_POINTLISTTOOLS_H
#define TIL_POINTLISTTOOLS_H

#include <cartobase/config/cartobase_config.h>

// Includes from STL
#include <memory>
#include <vector>

// Includes from TIL library
#include "til/til_common.h"


// Namespace 
namespace til {

template < typename T >
bool contains(const std::vector<numeric_array<T,3> > &pl, const numeric_array<T,3> &v)
{
	typename std::vector<numeric_array<T,3> >::const_iterator iPl;

	for (iPl = pl->begin(); iPl != pl->end(); ++iPl)
	{
		if ((*iPl) == v) return true;
	}
	return false;
}


/// Add all the vectors of pl to v.
/// NB: v is not initialized to zero (yet?).
// TODO: should be called sum
template < typename T1, typename T2 >
void add(const std::vector<numeric_array<T1,3> > &pl, numeric_array<T2,3> &v)
{
	typename std::vector<numeric_array<T1,3> >::const_iterator iPl;

	for (iPl = pl.begin(); iPl != pl.end(); ++iPl)
	{
		v += iPl->data();
	}
}


template < typename T1, typename T2 >
void add2ndOrder(const std::vector<numeric_array<T1,3> > &pl, SymMatrix3<T2> &mat)
{
	typename std::vector<numeric_array<T1,3> >::const_iterator iPl;

	// TODO: there is probably something better than using this temp variable
	SymMatrix3<T1> temp;

	for (iPl = pl->begin(); iPl != pl->end(); ++iPl)
	{
		tdot(*iPl, temp);
		add(temp, mat);
	}
}
/*
template < typename TImage >
class SetImageVoxel
{
public:
	SetImageVoxel(const Ptr<TImage> &im, typename TImage::value_type value) : m_im(im), m_value(value) {}
	void operator()(const numeric_array<int,3> &pos)
	{ 
		if (contains(m_im, *iPl)) (*m_im)(pos) = m_value;
	}
private:
	Ptr<TImage> m_im;
	typename TImage::value_type m_value;
};


/// This factory has mainly the role of using template argument deduction
template < class TImage >
SetImageVoxel<TImage> setImageVoxel(const Ptr<TImage> &im, typename TImage::value_type value)
{
	return SetImageVoxel<TImage>(im, value);
}
*/

/// std::for_each(pl.begin(), pl.end(), setImageVoxel(im, value))
/*
template < class TImage >
void dumpPointListInImage(const std::vector<numeric_array<int,3> > &pl, const Ptr<TImage > &im, typename TImage::value_type value)
{
	std::vector<numeric_array<int,3> >::const_iterator iPl;

	for (iPl = pl->begin(); iPl != pl->end(); ++iPl)
	{
		if (contains(im, *iPl))
		{
			im(*iPl) = value;
		}
	}
}
*/

// TODO: that's not worth a function, coz there is probably a grep out there...
template < class TImage1, class TImage2>
#if __cplusplus >= 201103L
std::unique_ptr<std::vector<numeric_array<int,3> > >
#else
std::auto_ptr<std::vector<numeric_array<int,3> > >
#endif
keepPointsAboveThreshold
(
  const std::vector<numeric_array<int,3> > &pl,
  const TImage1 &im,
  TImage2 &seg,
  typename TImage1::value_type threshold,
  typename TImage2::value_type background
)
{

#if __cplusplus >= 201103L
	std::unique_ptr<std::vector<numeric_array<int,3> > > newPl = new std::vector<numeric_array<int,3> >;
#else
	std::auto_ptr<std::vector<numeric_array<int,3> > > newPl = new std::vector<numeric_array<int,3> >;
#endif
	std::vector<numeric_array<int,3> >::const_iterator iPl;

	for (iPl = pl.begin(); iPl != pl.end(); ++iPl)
	{
		if (im.getValue(*iPl) > threshold)
		{
			newPl->push_back(*iPl);
		}
		else
		{
			seg(*iPl) = background;
		}
	}
	return newPl;
}

} // namespace

#endif
