#ifndef TIL_IMAGECLINEARITERATOR_H
#define TIL_IMAGECLINEARITERATOR_H

// includes from TIL library
#include "til/ConstImageLinearIterator.h"


// Namespace 
namespace til {


/// A class to iterate through all values of an image.
/// Allows to read and write to the image.
template < typename T >
//class ImageCLinearIterator : public ConstLinearIterator<ImageC<T> >, public VolumetricImageIterator_label
class LinearIterator<ImageC<T> > : public ConstLinearIterator<ImageC<T> >, public VolumetricImageIterator_label
{

public: // typedef

	//typedef T	value_type;
	typedef T&	reference;
	/// Type returned by operator* (needed for template expression)
	//typedef T& StarType;


public: // constructors & destructor

	LinearIterator<ImageC<T> >(ImageC<T> &im)
		: ConstLinearIterator<ImageC<T> >(im) {};

public: // functions

	/// Get reference of the current pixel
	reference operator* () const { return *(this->getIndex()); }
};

} // namespace


#endif
