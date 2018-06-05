#ifndef TIL_CONSTIMAGENCVOLUMETRICITERATOR_H
#define TIL_CONSTIMAGENCVOLUMETRICITERATOR_H

// includes from TIL library
#include "til/til_common.h"
#include "til/labels.h"
#include "til/miscTools.h"
#include "til/numeric_array.h"
#include "til/Range.h"


namespace til
{
  
// Predeclaration
template < typename T > class ImageNC;


// A class to iterate through all values of an image while knowing the
// image coordinates of the current pixel
// Allows to read but not to write to image


template < typename T >
//class ConstImageNCVolumetricIterator : public VolumetricImageIterator_label
class ConstVolumetricIterator<ImageNC<T> > : public VolumetricImageIterator_label
{

public: // typedefs

	typedef T			value_type;
	typedef const T &	reference;
	typedef ImageNC<T>	TImage;
	/// Type returned by operator* (needed for template expression)
	//typedef T StarType;


public: // constuctors & destructors

	// Default constructor
  ConstVolumetricIterator<ImageNC<T> >() { m_index = 0; }

	// Iterates on the whole image
	ConstVolumetricIterator<ImageNC<T> >(const ImageNC<T> &im) : m_im(const_cast<ImageNC<T>&>(im)) { this->init(); }

	// Iterates only on the volume of interest specified by the box
	ConstVolumetricIterator<ImageNC<T> >(const ImageNC<T> &im, const Range<int,3> & roi) : m_im(im) { this->init(roi); }
	
	virtual ~ConstVolumetricIterator<ImageNC<T> >() {}


public: // initialization

	virtual void init();
	virtual void init(const Range<int,3> & box);


public: // set & get

	// Get current position
  /*
  int getX() const { return m_pos.getX(); }
	int getY() const { return m_pos.getY(); }
	int getZ() const { return m_pos.getZ(); }
	int getPos(int i) const { return m_pos.get(i); }
  */
	const numeric_array<int,3> & pos() const { return m_pos; }

	// Set current position
	//void setPos(int x, int y, int z);
	void set_pos(const numeric_array<int,3> & pos);// { this->setPos(EXPAND_VECTOR(pos)); }

	// Same as above, without range checking
	//void setUnsafePos(int x, int y, int z);
	void setUnsafePos(const numeric_array<int,3> &pos);// { this->setUnsafePos(EXPAND_VECTOR(pos)); }

	// Get the region of interest
	const Range<int,3> & getRoi() const { return m_roi;}

	// Get current image
	const ImageNC<T> & image() const { return m_im; }


public: // functions


	// Get value of a neighbor in an unsafe way
	// NB: the offset is passed, not the actual position
	INLINE T getUnsafeValue(int offsetx, int offsety, int offsetz) const
	{
		return m_im.getUnsafeValue(
			this->pos()[0] + offsetx,
			this->pos()[1] + offsety,
			this->pos()[2] + offsetz);
	}

	template < int offsetx, int offsety, int offsetz >
		INLINE T getUnsafeValue() const
	{
		// TODO: there might be something better to do here, e.g.
		// newplanebegin - oldplanebegin + m_index + offsetx + offsety * m_offset.getY()
		// because then the multiplication on getY() is not done for offsety = -1, 0, 1...
		if (offsetz) return m_im.getUnsafeValue( this->pos()[0] + offsetx, this->pos()[1] + offsety, this->pos()[2] + offsetz);
		else return *(m_index + offsetx + offsety * m_offset[1] );
	}

	template < class Extrapolator, int offsetx, int offsety, int offsetz >
		INLINE T getValue() const
	{
		if (containsNeighbor<offsetx, offsety, offsetz>(*this))
		{
			return this->getUnsafeValue<offsetx, offsety, offsetz>();
		}
		else
		{
			return Extrapolator::getExtrapolatedValue(m_im, this->pos() + 
        numeric_array<int,3>(offsetx, offsety, offsetz));
// TODO: Did I broke efficiency again???
      /*
				this->pos()[0] + offsetx,
				this->pos()[1] + offsety,
				this->pos()[2] + offsetz);
        */
		}
	}



	template < class Extrapolator >
	T getValue(int offsetx, int offsety, int offsetz) const
	{
		return Extrapolator::getValue(m_im, 
			this->pos()[0] + offsetx,
			this->pos()[1] + offsety,
			this->pos()[2] + offsetz);
	}

	template < class Extrapolator >
	T getValue(const numeric_array<int,3> & offset) const
	{
		return this->getValue<Extrapolator>(EXPAND_VECTOR(offset));
	}

	// Get value of a neighbor
	// NB: the offset is passed, not the actual position
	INLINE T operator()(int offsetx, int offsety, int offsetz) const
	{
		return m_im.getValue(
			this->pos()[0] + offsetx,
			this->pos()[1] + offsety,
			this->pos()[2] + offsetz);
	}

	INLINE T operator()(const numeric_array<int,3> &v) const
	{
		return this->operator()(EXPAND_VECTOR(v));
	}


	INLINE ConstVolumetricIterator<ImageNC<T> > & operator++();


	// Print some info about the object
	//friend void printInfo FRIEND_TEMPLATE_NO_ARG (ConstVolumetricIterator<ImageNC<T> >&);


	// Got to the next element along given axis
	INLINE void next(ImageAxis axis);
	bool isAtEnd() const { return (m_index == 0); }
	reference operator*() const { return *m_index; }

	const T * const getIndex() { return m_index; }

protected: // data

	// Current position in the volume
	numeric_array<int,3> m_pos;

	// Precomputed offsets to move in each direction
	numeric_array<int,3> m_offset;
	
	// Image
	ImageNC<T> & m_im;

	// Region of interest in the image
	Range<int,3> m_roi;

	// Pointer to current element
	T* m_index;
};

































////////////////////////////////////////////////////
//
// CODE
//
////////////////////////////////////////////////////


template < typename T >
void ConstVolumetricIterator<ImageNC<T> >::init(const Range<int,3> &box)
{
	/*
	if (!isAllocated(im))
	{
		m_index = 0;
		return;
	}
	*/
	
	if (!contains(getRange(m_im), box))
	{
		throw std::domain_error("ROI lies outside image range");
	}

	m_roi = box;
	m_pos = box.min_bounds();

	m_index = m_im.getUnsafePointerOf(m_pos);

	m_offset[0] = 1;
	m_offset[1] = m_im.dim()[0];
	m_offset[2] = m_im.dim()[0]*m_im.dim()[1];
}




template < typename T >
void ConstVolumetricIterator<ImageNC<T> >::init()
{
	/*
	if (!isAllocated(im))
	{
		m_index = 0;
		return;
	}
	*/
	this->init(getRange(m_im));
}


template < typename T >
INLINE ConstVolumetricIterator<ImageNC<T> > & 
ConstVolumetricIterator<ImageNC<T> >::operator++()
{
	// We hit a border

	if (++(m_pos[0]) > m_roi.max_bounds()[0])
	{
		m_pos[0] = m_roi.min_bounds()[0];

		if (++(m_pos[1]) > m_roi.max_bounds()[1])
		{
			m_pos[1] = m_roi.min_bounds()[1];
			
			if (++(m_pos[2]) > m_roi.max_bounds()[2])
			{
				m_index = 0;
				return *this;
			}
		}

	// TODO: to be faster, use offsets to jump in y and z directions
	m_index = m_im.getUnsafePointerOf(m_pos);
	}


	// We're still in the volume

	else
	{
		++m_index;
	}

	return *this;
}

/*
template < typename T >
void ConstVolumetricIterator<ImageNC<T> >::setPos(int x, int y, int z)
{
	if (!contains(m_roi, x, y, z))
	{
		throw std::out_of_range("Point does not lie within iterator range");
	}

	this->setUnsafePos(x, y, z);
}
*/

template < typename T >
void ConstVolumetricIterator<ImageNC<T> >::set_pos(const numeric_array<int,3> & pos)
{
  if (!contains(m_roi, pos))
    throw std::out_of_range("Point does not lie within iterator range");
  this->setUnsafePos(pos);
}

/*
template < typename T >
INLINE void ConstVolumetricIterator<ImageNC<T> >::setUnsafePos(int x, int y, int z)
{
	m_pos.setX(x);
	m_pos.setY(y);
	m_pos.setZ(z);

	m_index = m_im.getUnsafePointerOf(m_pos);
}
*/
template < typename T >
INLINE void ConstVolumetricIterator<ImageNC<T> >::setUnsafePos(const numeric_array<int,3> & pos)
{
  m_pos = pos;
	m_index = m_im.getUnsafePointerOf(m_pos);
}

template < typename T >
INLINE void
ConstVolumetricIterator<ImageNC<T> >::next(ImageAxis axis)
{

	if (++m_pos[axis] > m_roi.max_bounds()[0])
	{
		// We hit a border
		m_index = 0;
	}
	else
	{
		// We're still in the volume
		m_index += m_offset[axis];
	}
}


template < typename T >
void printInfo(ConstVolumetricIterator<ImageNC<T> > &it)
{
	std::cout << "Pointer: " << it.m_index << std::endl;
	std::cout << "Position: " << it.m_pos << std::endl;
	std::cout << "Offset: " << it.m_offset << std::endl;
	std::cout << "Range: " << it.m_roi.min_bounds() <<" "<< it.m_roi.max_bounds() << std::endl;
}

} // namespace

#endif
