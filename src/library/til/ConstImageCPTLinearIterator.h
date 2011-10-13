#ifndef TIL_CONSTIMAGERLELINEARITERATOR_H
#define TIL_CONSTIMAGERLELINEARITERATOR_H

// includes from STL
#include <list>
#include <vector>

// includes from TIL library
#include "til/til_common.h"
#include "til/labels.h"


// Namespace 
namespace til {


// Predeclaration
template < typename T > class ImageRLE;


// A class to iterate through all values of an image without having to
// know the current position (hence the fastest possible)
// Allows to read but not to write to image


template < typename T >
class ConstLinearIterator<ImageRLE<T> > 
	: public ImageIterator_label
{
public: // typedefs

	typedef T value_type;
	typedef ImageRLE<T> TImage;
	// Actually, not a reference at all, just the return type of operator*
	// -- which is assumed to be the same thing in the STD... :/
	typedef T reference;

public: // constructors & destructor

//	ConstLinearIterator<ImageRLE<T> >() { this->init(); }
	ConstLinearIterator<ImageRLE<T> >(const ImageRLE<T> &im) : m_im(const_cast<ImageRLE<T>&>(im)) { this->init(); }
	virtual ~ConstLinearIterator<ImageRLE<T> >() {};

public: // initialization

	// Initialize the iterator
	void init();
	void init(const ImageRLE<T> &im);

public: // functions

	bool isAtEnd() const { return m_flagIsAtEnd; }

	reference operator*() const { return m_value.getValue(); }

	void operator++()
	{	
		--(m_value.getRepeat());

		// Check whether there are no more pixels remaining
		if (m_value.getRepeat() == 0)
		{
			// go to next pixel on the line
			++m_iList;

			// Check whether there are any pixel remaining
			if (m_iList == m_iLine->end())
			{
				// If not, go to next line
				++m_iLine;

				// Check whether there are indeed any line remaining
				if (m_iLine == m_im.m_data.end())
				{
					m_flagIsAtEnd = true;
					return;
				}
				m_iList = m_iLine->begin();
			}
			m_value = *m_iList;
		}
	}

	 const ImageRLE<T> & image() const { return m_im; }

private: // typedefs

	typedef typename ImageRLE<T>::Line Line;
	typedef typename ImageRLE<T>::Data Data;

private: // data

	// The current stack of point
	typename ImageRLE<T>::RepeatedValue m_value;

	// The current image
	ImageRLE<T> & m_im;

	// the current slice
	int m_lineIndex;

	// indexes of current processed point
	typename Line::const_iterator m_iList;
	typename Data::const_iterator m_iLine;

	//

	bool m_flagIsAtEnd;
};


































////////////////////////////////////////////////////
//
// CODE
//
////////////////////////////////////////////////////


template < typename T >
void ConstLinearIterator<ImageRLE<T> >::init()
{
	/*
	if (!isAllocated(im))
	{
		m_lineIndex = 0;
		m_flagIsAtEnd = true;
		return;
	}
	*/

	m_lineIndex = 0;
	m_iLine = m_im.m_data.begin();
	m_iList = (*m_iLine).begin();
	m_value = *m_iList;
	m_flagIsAtEnd = false;
}

} // namespace


#endif


