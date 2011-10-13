#ifndef TIL_IMAGECVOLUMETRICITERATOR_H
#define TIL_IMAGECVOLUMETRICITERATOR_H

// includes from TIL library
#include "til/ConstImageVolumetricIterator.h"
#include "til/Range.h"


// Namespace 
namespace til {


// A class to iterate through all values of an image while knowing the
// image coordinates of the current pixel
// Allows to read and write to image


template < typename T >
//class ImageCVolumetricIterator : public ConstImageCVolumetricIterator<T>
class VolumetricIterator<ImageC<T> > : public ConstVolumetricIterator<ImageC<T> >
{

public: // typedefs

	//typedef T value_type;
	typedef T& reference;
	/// Type returned by operator* (needed for template expression)
	// TODO: SHould'nt this be a trait?
	//typedef T& StarType;

public: // constructors & destructor

    VolumetricIterator<ImageC<T> >() : ConstVolumetricIterator<ImageC<T> >() {};
	VolumetricIterator<ImageC<T> >(ImageC<T> &im) 
	: ConstVolumetricIterator<ImageC<T> >(im) { }
	VolumetricIterator<ImageC<T> >(ImageC<T> &im, const Range<int,3> & box)
    : ConstVolumetricIterator<ImageC<T> >(im, box) { }

public: // functions

	// Get reference of the current pixel
	reference operator* () const { return *(this->getIndex()); }
};





} // namespace

#endif

