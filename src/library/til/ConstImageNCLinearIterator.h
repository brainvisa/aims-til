#ifndef TIL_CONSTIMAGENCLINEARITERATOR_H
#define TIL_CONSTIMAGENCLINEARITERATOR_H

// includes from TIL library
#include "til/til_common.h"
#include "til/imageIterator.h"
#include "til/labels.h"

// Namespace 
namespace til
{


// Predeclaration
template < typename T > class ImageNC;


// A class to iterate through all values of an image without having to
// know the current position (hence the fastest possible)
// Allows to read but not to write to image
template < typename T >
//class ConstImageNCLinearIterator : public ImageIterator_label
class ConstLinearIterator<ImageNC<T> > : public ImageIterator_label
{

public: // typedefs

	typedef T			value_type;
	typedef const T &	reference;
	typedef ImageNC<T>	TImage;
	/// Type returned by operator* (needed for template expression)
	//typedef T StarType;

public: // constructors & destructor

	//ConstLinearIterator<ImageNC<T> >(){ this->init(); }
	ConstLinearIterator<ImageNC<T> >(const ImageNC<T> &im) : m_im(const_cast<ImageNC<T>&>(im)) { this->init(); }
	virtual ~ConstLinearIterator<ImageNC<T> >() {};


public: // initialization

	void init();
	//void init(const ImageNC<T> &im);


public: // functions

	void operator++()
	{
		// Check whether current slice is done
		if (--m_nPixelsLeft == 0)
		{
			// We are done: go to next slice
			++m_sliceIndex;

			// Check whether all slices are already done
			if (m_sliceIndex == m_im.dim()[2]) 
			{
				// Yes: stop
				m_index = 0;
				return;
			}
			else
			{
				// No: initialize variables for the current slice
				this->initForSlice(m_sliceIndex);
			}
		}
		else
		{
			// We are still in the current slice
			// simply increment the pointer
			++m_index;
		}
	}

	bool isAtEnd() const { return (m_index == 0); }
	reference operator*() const { return *m_index; }

protected: // functions
	
	T * getIndex() { return m_index; }

private: // functions

	void initForSlice(int i)
	{
		m_nPixelsLeft = m_sliceSize;
		m_index = const_cast<T*>(m_im.getSlicePointer(i));		
	}


protected: // data

	// Pointer to current element
	T* m_index;


private: // data

	// The current image
	ImageNC<T> & m_im;

	// the number of pixels per slice
	int m_sliceSize;

	// the number of pixels that still haven't been run into
	// in the current slice
	int m_nPixelsLeft;

	// the current slice
	int m_sliceIndex;
};


































////////////////////////////////////////////////////
//
// CODE
//
////////////////////////////////////////////////////


template < typename T >
void ConstLinearIterator<ImageNC<T> >::init()
{
/*
	if (!isAllocated(im))
	{
		m_index = 0; 
		m_nPixelsLeft = 0;
		m_sliceSize = 0;
		m_sliceIndex = 0;
		return;
	}*/

	m_sliceSize = m_im.getSliceSize();

	// Start at slice 0
	m_sliceIndex = 0;
	this->initForSlice(m_sliceIndex);
}

} // namespace

#endif

