#ifndef TIL_IMAGERLELINEARITERATOR_H
#define TIL_IMAGERLELINEARITERATOR_H

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
class LinearIterator<ImageRLE<T> > : public ImageIterator_label
{
public: // predeclarations

	class VirtualVoxel;

public: // typedefs

	typedef T				value_type;
	typedef ImageRLE<T>		TImage;
	typedef VirtualVoxel &	reference;

public: // constructors & destructor

	//LinearIterator<ImageRLE<T> >() { this->init(); }
	LinearIterator<ImageRLE<T> >(ImageRLE<T> &im) : m_im(im) { this->init(); }
	virtual ~LinearIterator<ImageRLE<T> >() {};

public: // initialization

	// Initialize the iterator
	void init();
	//void init(ImageRLE<T> &im);
	
public: // classes

	friend class VirtualVoxel;

	class VirtualVoxel
	{
	public: // friends
		friend class LinearIterator<ImageRLE<T> >;
	public:
		VirtualVoxel() : m_it() {};
		//VirtualVoxel(LinearIterator<ImageRLE<T> > &it) : m_it(it) {}
		void operator=(const T &value)
		{
			if (value != m_it->getValue())
			{
				m_it->setUnsafeValue(value);
			}
		}

		template < typename X >
		void operator*=(const X & value)
		{
			m_it->setUnsafeValue(m_it->getValue() * value);
		}

		operator T () { return m_it->getValue(); }

	private: // set & get
		void set(LinearIterator<ImageRLE<T> > *it) { m_it = it; }
	private: // data
		LinearIterator<ImageRLE<T> > *m_it;
	};

public: // functions

	bool isAtEnd() { return m_isAtEnd; }

	const T & getValue() const { return m_iRepValue->getValue(); } //{ return m_value.getValue(); }

	void setUnsafeValue(T value)
	{
		m_im.setUnsafeValue(value, *m_iLine, m_iRepValue, m_count);
	}


	//m_VirtualVoxel operator*() { return m_VirtualVoxel(*this); }
	VirtualVoxel & operator*() { return m_virtualVoxel; }

	//T& operator* () { return m_value.value; }

	void operator++()
	{	
		// Check whether value has been changed
		/*
		if (m_savedValue != m_value.value)
		{
			// If so, insert new value
			m_im._setUnsafeValue(m_value.value, *m_iLine, m_iRepValue, m_value.repeat);
			m_savedValue = m_value.value = (*m_iRepValue).value;
		}
		*/
		--m_count;

		// Check whether there are no more pixels remaining
		if (m_count == 0)
		{
			// go to next pixel on the line
			++m_iRepValue;

			// Check whether there are indeed any pixel remaining
			if (m_iRepValue == m_iLine->end())
			{
				// If not, go to next line
				++m_iLine;

				// Check whether there are indeed any line remaining
				if (m_iLine == m_im.m_data.end())
				{
					m_isAtEnd = true;
					return;
				}
				m_iRepValue = m_iLine->begin();
			}
			m_count = m_iRepValue->getRepeat();
			//m_value = *m_iRepValue;
			//m_savedValue = m_value.value;
		}
	}


private: // typedefs

	typedef typename ImageRLE<T>::Line Line;
	typedef typename ImageRLE<T>::Data Data;

private: // functions


private: // data

	// The current stack of point
	//typename TImage::RepeatedValue m_value;

	// The current image
	TImage & m_im;

	// the current slice
	int m_lineIndex;
	
	// indexes of current processed point
	typename Line::iterator m_iRepValue;
	typename Data::iterator m_iLine;
	// This is a countdown of the number of voxels remaining in the current RepeatedValue.
	int m_count;	

	// Flag to indicate whether we already scanned the entire image
	bool m_isAtEnd;

	// The value we are currently pointing at
	//T m_savedValue;

	// We allocate a virtual voxel so that we don't have to create one every time operator*
	// is called
	VirtualVoxel m_virtualVoxel;
};


































////////////////////////////////////////////////////
//
// CODE
//
////////////////////////////////////////////////////


template < typename T >
void LinearIterator<ImageRLE<T> >::init()
{
	/*
	if (!isAllocated(im))
	{
		m_lineIndex = 0;
		m_isAtEnd = true;
		return;
	}
	*/

	m_lineIndex = 0;
	m_iLine = m_im.m_data.begin();
	m_iRepValue = m_iLine->begin();
	m_count = m_iRepValue->getRepeat();
	//m_value = *m_iRepValue;
	m_isAtEnd = false;
	//m_savedValue = m_value.getValue();
	m_virtualVoxel.set(this);
}

} // namespace


#endif

