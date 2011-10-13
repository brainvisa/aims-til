#ifndef TIL_IMAGENCLINEARITERATOR_H
#define TIL_IMAGENCLINEARITERATOR_H

// includes from TIL library
#include "til/til_common.h"
#include "til/ConstImageNCLinearIterator.h"


// Namespace 
namespace til {


// A class to iterate through all values of an image without having to
// know the current position (hence the fastest possible)
// Allows to read and write to the image
template < typename T >
class LinearIterator<ImageNC<T> > : public ConstLinearIterator<ImageNC<T> >
{

public: // typedefs

	typedef T & reference;
	//typedef T value_type;
	/// Type returned by operator* (needed for template expression)
	//typedef T& StarType;

public: // constructors & destructor

	LinearIterator<ImageNC<T> >(ImageNC<T> &im)
		: ConstLinearIterator<ImageNC<T> >(im) {};

public: // functions

	/// Get reference of the current pixel
	reference operator* () { return *(this->m_index); }
};


} // namespace

#endif
