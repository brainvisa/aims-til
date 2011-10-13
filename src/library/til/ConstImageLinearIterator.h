#ifndef TIL_CONSTIMAGECLINEARITERATOR_H
#define TIL_CONSTIMAGECLINEARITERATOR_H

// includes from TIL library
#include "til/labels.h"

// Namespace 
namespace til
{


// Predeclaration
template < typename T > class ImageC;


// A class to iterate through all values of an image without having to
// know the current position (hence the fastest possible)
// Allows to read but not to write to image

template < typename T >
class ConstLinearIterator<ImageC<T> > : public ImageIterator_label
{

public: // typedefs

	typedef T			value_type;
	typedef ImageC<T>	TImage;
	typedef const T &	reference;

public: // constructors & destructor

	ConstLinearIterator<ImageC<T> >() { m_index = 0; };
	ConstLinearIterator<ImageC<T> >(const ImageC<T> &im) { this->init(im); }
	virtual ~ConstLinearIterator<ImageC<T> >() {};

public: // initialization

	// Initialize the iterator 
	void init(const ImageC<T> &im);

public: // set & get

	void setEnd(const T * end)
	{
		if (end < m_index)
		{
			throw std::invalid_argument("end is after current position");
		}
		m_end = const_cast<T*>(end);
	}

public: // functions

	// NB: it appeared that _not_ being compliant with the standard, by returning
	// nothing (as opposed to returning a reference on the object) is 50%
	// faster (tested on release mode w/ profiling).
	// It was assumed here that this performance was preferable to enabling
	// such statement as iter1 = ++iter2;

	INLINE void operator++()
	{
		++m_index;
	}

	INLINE bool next() { if (m_index == m_end) return 0; ++m_index; return 1;}

	bool isAtEnd() const { return (m_index > m_end); }
	reference operator*() const { return *m_index; }

public:
//protected: // set & get

	// get pointer on the current voxel
	T* getIndex() const { return m_index; }


private: // data

	// Pointer to current element	
	T* m_index;
	// Pointer to last element
	// NB: this is not the (last+1) element (like ::std), because we want this
	// to work even for volumetric iterators. But, volumetric iterators can run
	// along y and z axis. So end cannot be simply the 'next after the end' because
	// the notion of 'next' depends on which direction we move on. So end
	// has to be the real end element.
	T* m_end;
};


































////////////////////////////////////////////////////
//
// CODE
//
////////////////////////////////////////////////////


template < typename T >
void ConstLinearIterator<ImageC<T> >::init(const ImageC<T> &im)
{
	if (!isAllocated(im))
	{
		m_index = 0;
		m_end = 0;
		return;
	}

	m_index = const_cast<T*>(im.getPointer());
	m_end = m_index + im.size() - 1;
}

} // namespace

#endif
