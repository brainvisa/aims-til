#ifndef TIL_IMAGENCVOLUMETRICITERATOR_H
#define TIL_IMAGENCVOLUMETRICITERATOR_H

// includes from TIL library
#include "til/til_common.h"
#include "til/ConstImageNCVolumetricIterator.h"
#include "til/Range.h"

namespace til
{

// A class to iterate through all values of an image while knowing the
// image coordinates of the current pixel
// Allows to read and write to image
template < typename T >
class VolumetricIterator<ImageNC<T> > : public ConstVolumetricIterator<ImageNC<T> >
{

public: // typedefs

	typedef T & reference;
	//typedef T value_type;
	/// Type returned by operator* (needed for template expression)
	//typedef T& StarType;

public: // constructors & destructor

    /// default constructor
	VolumetricIterator<ImageNC<T> >() : ConstVolumetricIterator<ImageNC<T> >() {};

	/// constructor over the whole image range
	VolumetricIterator<ImageNC<T> >(ImageNC<T> &im) 
	: ConstVolumetricIterator<ImageNC<T> >(im) { }

	// constructor over a volume of interest
	VolumetricIterator<ImageNC<T> >(ImageNC<T> &im, Range<int,3> &voi)
    : ConstVolumetricIterator<ImageNC<T> >(im, voi) { }

public: // functions

	/// Get reference of the current pixel
	reference operator* () const { return *(this->m_index); }
};


} // namespace


#endif

