#ifndef TIL_IMAGERLE_H
#define TIL_IMAGERLE_H

// includes from STL
#include <cassert>
//#include <float.h>
#include <limits>
#include <list>
#include <cstring>
#include <vector>
#include <iostream>


// includes from TIL library
#include "til/til_common.h"
#include "til/Box.h"
#include "til/ImageBase.h"
#include "til/image_common.h"
#include "til/numeric_array.h"
#include "til/templateTools.h"

#include "til/ConstImageCPTLinearIterator.h"
#include "til/ImageCPTLinearIterator.h"
//#include "til/ConstImageRLEVolumetricIterator.h"
//#include "til/ImageRLEVolumetricIterator.h"


// Namespace 
namespace til {


/// Image class using run-length encoded data
template < typename T >
//class ImageRLE : public SmartObject, public ImageBase
class ImageRLE : public ImageBase
{
private: // typedefs

	class RepeatedValue;
	typedef typename std::list<RepeatedValue>	Line;
	typedef typename std::vector<Line>			Data;


public: // typedefs

	typedef T value_type;
	typedef ImageRLE<T> Self;

public: // constructors & destructor
	
	ImageRLE();
	ImageRLE(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz);
	ImageRLE(const numeric_array<int,3> &dim, const numeric_array<t_voxsize,3> &vDim);
	ImageRLE(const ImageRLE<T> *);
	//template < typename TImage >
	//ImageRLE(const TImage *im, typename enable_if<is_Image<TImage> >::type);
	ImageRLE(const ImageParameter &param);
	
	/*
	/// Default constructor
	/// Buffer is not allocated
	/// Every other internal variables are set to zero
	static Self* New() { return new Self(); }
	/// Allocate an image of dimension (x, y, z) with voxel size (vx, vy, vz)
	static Self* New(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz) { return new Self(x,y,z,vx,vy,vz); }
	/// Same as above, using vectors
	static Self* New(const numeric_array<int,3> &dim, const numeric_array<t_voxsize,3> &vDim) { return new Self(dim, vDim); }
	/// Creates an image from an ImageParamter structure.
	static Self* New(const ImageParameter &param) {return new Self(param); }
	*/

public: // initialization

	// Initialization functions
	// Basically similar to constructors

	void init();
	void init(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz);
	void init(const numeric_array<int,3> &dim, const numeric_array<t_voxsize,3> &vDim);
	void init(const ImageRLE<T> *);
	void init(const ImageParameter &param);
	//template < typename TImage >
	//void init(const TImage *im, typename enable_if<is_Image<TImage> >::type);

public: // functions

	// Get value

	// Access (read-write) the value at point (i,j,k)
	// Range checking is done
	
	INLINE T operator()(int i, int j, int k);
	T operator()(const numeric_array<int,3> &v) { return this->operator()(EXPAND_VECTOR(v)); }

	/// Get the value at point (i,j,k) for read-only purpose
	// If point lies outside image, an exception is thrown
	INLINE const T getValue(int i, int j, int k) const;
	const T getValue(const numeric_array<int,3> &v) const { return this->getValue(EXPAND_VECTOR(v)); }
	
	INLINE void setValue(T value, int i, int j, int k);
	void setValue(T value, const numeric_array<int,3> &v) { this->setValue(value, EXPAND_VECTOR(v)); }
	INLINE void setUnsafeValue(T value, int i, int j, int k);
	
	/// Get the value at point (i,j,k) for fast read-only purpose
	// WARNING: No range checking!
	T getUnsafeValue(int i, int j, int k) const;
	T getUnsafeValue(const numeric_array<int,3> &v) const { return this->getUnsafeValue(EXPAND_VECTOR(v)); }
	
	/// Check whether data has been allocated or not
	bool isAllocated() const { return (!m_data.empty());}
	
	/// Set all values to zero
	void reset();
	
	void debugLine(const Line & line) const
	{
		for (typename Line::const_iterator iRepValue = line.begin(); iRepValue != line.end(); ++iRepValue)
		{
			std::cout << "Value : " <<  iRepValue->getValue() << "  ";
			std::cout << "Repeat : " << iRepValue->getRepeat() << ",  ";
		}
		std::cout << std::endl;
	}

	/// for debugging purposes
	void debug() const
	{
		for (typename Data::const_iterator iLine = m_data.begin(); iLine != m_data.end(); ++iLine)
		{
			this->debugLine(*iLine);
		}
	}

	void copy(const Self & im)
	{
		if (!((this->dim()[0]  == im.dim()[0]) &&
			(this->dim()[1]  == im.dim()[1]) &&
			(this->dim()[2]  == im.dim()[2]) &&
			(this->vdim()[0] == im.vdim()[0]) &&
			(this->vdim()[1] == im.vdim()[1]) &&
			(this->vdim()[2] == im.vdim()[2])))
		{
			throw std::invalid_argument("Incompatible images");
		}

		// TODO: check wheter a simple call to std::copy on the vector would not do, i.e.
		// get rid of the loop. Or carremently using operator=??
		typename Data::iterator myData = m_data.begin();
		typename Data::const_iterator imData = im.m_data.begin();
		for(; myData != m_data.end(); ++myData, ++imData)
		{
			//std::copy(imData->begin(), imData->end(), myData->begin());
			myData->assign(imData->begin(), imData->end());
		}
		//this->debugLine(this->getLine(this->dim()[0]-1, this->dim()[1]-1));
		//im.debugLine(this->getLine(this->dim()[0]-1, this->dim()[1]-1));
	}

public: // friends

	friend class ConstLinearIterator<ImageRLE<T> >;
	friend class LinearIterator<ImageRLE<T> >;
//	friend class ImageRLELinearIterator<T>::m_RLEPixel;

private: // constructors


private: // classes

	// The basic atom of run-length encoding: a structure that comprises
	// a value, and the number of time this value is repeated.
	// TODO: should repeat be coded on int, unsigned short,...?
	class RepeatedValue
	{
	public: // typedef
		typedef int t_repeat;
	public: // constructors & destructor
		RepeatedValue() : m_value(), m_repeat(0) {};
		RepeatedValue(const T &v, t_repeat r) : m_value(v), m_repeat(r) {};

		// I think a good reason not to return a t_repeat for const is that if we
		// do ++(x.getRepeat()) we can't be sure if the const has not actually been called!
		const t_repeat & getRepeat() const { return m_repeat; }
		t_repeat & getRepeat() { return m_repeat; }
		const T & getValue() const { return m_value; }
		T & getValue() { return m_value; }
	private: // data
		T		m_value;			///< A value
		t_repeat	m_repeat;			///< Number of contiguous occurences of the value
	};

private: // methods
		
	void allocateData()
	{
		// Not good: should not allocate that much!
		/*
		if ( this->dim()[0]<0 || this->dim()[1]<0 || this->dim()[2]<0)
	    {
			throw std::invalid_argument("Image size < 0");
	    }
		*/
		int nLine = this->dim()[1] * this->dim()[2];
		m_data.resize(nLine);
	}

	// Return a repeated value corresponding to an entire line filled with default value
	// (zeros being an abusive shorthand).
	RepeatedValue zeroLine()
	{
		RepeatedValue zeros;
		zeros.getValue() = T();
		zeros.getRepeat() = this->dim()[0];
		return zeros;
	}

	// Common operations done in init methods
	void _init(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz);

	void setUnsafeValue(
		T value,
		Line & line,
		typename Line::iterator & iLine,
		int & remaining // I used const just to check I didn't screw up
		);
/*
								   void _setUnsafeValue(
		T value, Line line, typename Line::iterator &iLine,
//		typename std::list<RepeatedValue> line,
//		typename std::list<RepeatedValue>::iterator & iLine,
		int &remaining);
*/
	// Copy constructor is intentionally left undefined
	// Non-explicit image copying is unwanted!
	// use copy instead
	ImageRLE(const ImageRLE<T> &);


	Line & getLine(int y, int z)
	{
		return m_data[y + z * this->dim()[1]];
	}
	const Line & getLine(int y, int z) const
	{
		return m_data[y + z * this->dim()[1]];
	}

	// for debuggin purposes: nb of elements on a line
	int lineLength(const Line & line)
	{
		int count = 0;
		typename Line::const_iterator iRepValue = line.begin();
		for (; iRepValue != line.end(); ++iRepValue)
		{
			count += iRepValue->getRepeat();
		}
		return count;
	}

private: // data
	
	// data
	Data m_data;
};


/// Iterator traits of ImageRLE
template < typename T >
struct Iterator<ImageRLE<T> >
{
	typedef ConstLinearIterator<ImageRLE<T> >	ConstLinear;
	typedef LinearIterator<ImageRLE<T> >		Linear;
};































////////////////////////////////////////////////////
//
// CODE
//
////////////////////////////////////////////////////


// Constructors

template < typename T >
ImageRLE<T>::ImageRLE() : ImageBase()
{
    this->init();
}

template < typename T >
void ImageRLE<T>::init()
{
}
/*
template < typename T >
template < typename TImage >
ImageRLE<T>::ImageRLE(const TImage *im, typename enable_if<is_Image<TImage> >::type) : ImageBase()
{
	this->init<TImage>(im);
}

template < typename T >
template < typename TImage >
void ImageRLE<T>::init(const TImage *im, typename enable_if<is_Image<TImage> >::type)
{
	return;
}
*/

template < typename T >
ImageRLE<T>::ImageRLE(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz) : ImageBase()
{
    this->init(x, y, z, vx, vy, vz);
}

template < typename T >
void ImageRLE<T>::init(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz)
{
	// set dim and vdim
	this->_init(x, y, z, vx, vy, vz);
	// allocate data
	this->allocateData();
	// set data to zero
    this->reset();
	//std::cout << "constr " << this->getLine(this->dim()[1]-1, this->dim()[2]-1).rbegin()->getRepeat() << std::endl;
}

template < typename T >
ImageRLE<T>::ImageRLE(const numeric_array<int,3> &dim, const numeric_array<t_voxsize,3> &vDim) : ImageBase()
{
	this->init(dim, vDim);
}

template < typename T >
void ImageRLE<T>::init(const numeric_array<int,3> &dim, const numeric_array<t_voxsize,3> &vDim)
{
	this->init(EXPAND_VECTOR(dim), EXPAND_VECTOR(vDim));
}

template < typename T >
void ImageRLE<T>::_init(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz)
{
  this->set_dim(til::numeric_array<int, 3>(x, y, z));
  this->set_vdim(til::numeric_array<float,3>(vx, vy, vz));
}

template < typename T >
ImageRLE<T>::ImageRLE(const ImageRLE<T> *pIm) : ImageBase()
{
    this->init(pIm);
}

template < typename T >
void ImageRLE<T>::init(const ImageParameter &param)
{
	this->init(param.m_dim, param.m_vDim);
	//this->setOrigin(param.m_ori);
}

template < typename T >
void ImageRLE<T>::init(const ImageRLE<T> *pIm)
{
    this->init(param(pIm));
	m_data = pIm.p_data;
}

template < typename T >
ImageRLE<T>::ImageRLE(const ImageParameter &param) : ImageBase(), m_data(0)
{
	this->init(param);
}


template < typename T >
void ImageRLE<T>::reset()
{
	// Make a structure corresponding to a zero repeated dimX times
	RepeatedValue zeros = this->zeroLine();
	
	// Wrap that in a list: you get an image line
	typename std::list<RepeatedValue> emptyLine(1, zeros);

	// Push that as many times as there are lines in the image
	int nLine = this->dim()[1] * this->dim()[2];
	// TODO: assign might be for vectors only, check if std::fill is what we want.
	// Also, what is the policy of assign wrt allocation?
	m_data.assign(nLine, emptyLine);
}

template < typename T >
INLINE T ImageRLE<T>::operator()(int i, int j, int k)
{
	return this->getValue(i,j,k);
}

//template < typename T >
//const T* ImageRLE<T>::getConstPointerOf(int i, int j, int k) const 
//{
//	return this->getPointerOf(i,j,k);
//}



template < typename T >
INLINE const T ImageRLE<T>::getValue(int i, int j, int k) const
{
    if (this->contains(numeric_array<int,3>(i,j,k)))
    {
        return this->getUnsafeValue(i,j,k);
    }
    else
    {
		throw std::out_of_range("Position out of image range");
    }
}


template < typename T >
T ImageRLE<T>::getUnsafeValue(int i, int j, int k) const
{
	//std::cout << "before: size " << this->getLine(this->dim()[1]-1, this->dim()[2]-1).size() << " " << this->getLine(this->dim()[1]-1, this->dim()[2]-1).rbegin()->getRepeat() << std::endl;
	// Get the (y,z) line
	const Line & line = this->getLine(j,k);
	// Scan the line, and stop when the count just went over our x coordinate
	typename std::list<RepeatedValue>::const_iterator iRepValue;
	int current_i = 0;
	for (iRepValue = line.begin(); iRepValue != line.end(); ++iRepValue)
	{
		current_i += iRepValue->getRepeat();
		if (current_i > i)
		{
			return iRepValue->getValue();
		}
	}
	// Well -- we scanned the line but the count never went that far!
	// This is an error probably due to a bug, let's print some debug info.
	std::cerr << "Problem in ImageRLE probably due to a bug\n";
	std::cerr << "Coordinate asked: " << i << " on line ( " << j << " , " << k << " ) \n";
	std::cerr << "Max coordinate reached: " << current_i << "\n";
	std::cerr << "Current line: \n";
	for (iRepValue = line.begin(); iRepValue != line.end(); ++iRepValue)
	{
		std::cerr << "Repeat: " << iRepValue->getRepeat() << ", value: " << iRepValue->getValue() << "\n";
	}
	std::cerr << std::endl;
	//std::cout << "after " << this->getLine(this->dim()[1]-1, this->dim()[2]-1).rbegin()->getRepeat() << std::endl;
	throw std::underflow_error("Incomplete image -- cannot return a value");
}



template < typename T >
INLINE void ImageRLE<T>::setValue(T value, int i, int j, int k)
{
    if (this->contains(numeric_array<int,3>(i,j,k)))
    {
        this->setUnsafeValue(value,i,j,k);
    }
    else
    {
		throw std::out_of_range("Index out of range");
    }
}

// NB: right now, the whole setValue setup is quite inefficient, as one would have to
// redo all the scanning to set another value, which might just be next to it. In order
// to have an efficient iterator, one should rather return the current pointer on the
// current Repeated value, as well as an update of the 'remaining' number. Right now
// this is not used anywhere (and foremost not in the iterator) so I commented out related shit).
// Also that would require the remaining passed by reference. Probably one should look at
// a cleaner way (?) to do this.

// 'remaining' gives the 'local coordinate' inside the repeated value pointed at by iLine.
// remaining = 1 if the point is exactly at the end of the RepeatedValue. It is equal
// to repeat if it is exactly at the beginning.
// NB: this assumes that we have already checked that the value we put is different
// from the one that is already there.
template < typename T >
INLINE void ImageRLE<T>::setUnsafeValue
(
 T value,							///< The new value we want to set
 Line & line,							///< The current line on which the pixel lies
 typename Line::iterator & iRepValue,	///< An iterator pointing on the current repeated value
 int & remaining						///< A count giving us the position of the pixel inside the repeated value  // I used const just to check I didn't screw up
 )
{
	// If the pixel already has this value, stop here
	if (value == iRepValue->getValue()) return;

	// Get pointers on previous and next repeated value, as well as flags indicating whether
	// we lie at the very beginning/end of the current line
	typename Line::iterator iRepValuePrevious = iRepValue;
	bool atBeginning = (iRepValue == line.begin());
	if (!atBeginning) --iRepValuePrevious;
	bool atEnd = (iRepValue == line.end());
	typename Line::iterator iRepValueNext = iRepValue;
	++iRepValueNext;
	

	// Is the current repeated sequence made of only one element?
	if (iRepValue->getRepeat() == 1)
	{
		/*
		typename Line::iterator iRepValuePrevious = iRepValue;
		if (iRepValue != line.begin()) --iRepValuePrevious;
		typename Line::iterator iRepValueNext = iRepValue;
		++iRepValueNext;
		*/

		// Is the next value the same?
		if (!atBeginning && iRepValueNext->getValue() == value)
		{
			remaining += iRepValueNext->getRepeat();
			// Is the previous value the same?
			if (!atEnd && iRepValuePrevious->getValue() == value)
			{
				// Merge all three values
				iRepValueNext->getRepeat() += iRepValuePrevious->getRepeat()+1;
				line.erase(iRepValuePrevious);
			}
			else
			{
				// Merge current and next
				++(iRepValueNext->getRepeat());
			}
			line.erase(iRepValue);
			iRepValue = iRepValueNext;
			return;
		}

		// Is the previous value the same?
		else if (!atBeginning && iRepValuePrevious->getValue() == value)
		{
			// Merge current and previous
			line.erase(iRepValue);
			++(iRepValuePrevious->getRepeat());
			
			iRepValue = iRepValuePrevious;
			//assert(remaining == 1);
			//remaining = 1;
			return; // Merge previous + 1;
		}
		// Then we remain a single point : simply update the value.
		else
		{
			// Simply overwrite the value
			iRepValue->getValue() = value;

			//assert(remaining == 1);
			//remaining = 1;
			return; // 1
		}
	}

	// Are we exactly at the end of current repeated sequence?
	else if (remaining == 1)
	{	
		//typename Line::iterator iRepValueNext = iRepValue;
		//++iRepValueNext;
		
		// Then, decrease current repeat number by one
		--(iRepValue->getRepeat());
		
		// and look at the next sequence
		// If there is no next sequence, or,
		// if the next sequence has a different value
		if (atEnd || iRepValueNext->getValue() != value)
		{
			// Insert a new repeat value
			//RepeatedValue rv;
			//rv.value = value;
			//rv.repeat = 1;

			iRepValue = 
				line.insert(iRepValueNext, RepeatedValue(value, 1));
			remaining = 1;
			return; // Split previous + 1
		}
		else
		{
			// If the next sequence has the same value, just increase
			// its repeat number by one
			remaining = ++(iRepValueNext->getRepeat());
			iRepValue = iRepValueNext;
			return; // Split previous and merge 1+next
		}
	}
	// Are we at the very first of current repeated value?
	else if (remaining == iRepValue->getRepeat())
	{
		//typename Line::iterator iRepValuePrevious = iRepValue;
		//if (iRepValue != line.begin()) --iRepValuePrevious;
		
		// If so, decrease current repeat number by one
		--(iRepValue->getRepeat());
		// Set remaining to one because anyway we'll be at the end of a repvalue
		remaining = 1;

		// Is the previous value the same?
		if (!atBeginning && iRepValuePrevious->getValue() == value)
		{
			// Then just increase its repeat number
			++(iRepValuePrevious->getRepeat());
			
			iRepValue = iRepValuePrevious;
			return; // Merge previous+1
		}
		else
		{
			// else insert a new repeat value
			//RepeatedValue rv;
			//rv.value = value;
			//rv.repeat = 1;

			iRepValue = 
				line.insert(iRepValue, RepeatedValue(value, 1));
			return; // Split 1+current
		}
	}
	// we are in the middle of a repeated value
	else
	{
		//assert(remaining > 1);
		// First half
		line.insert(iRepValue, RepeatedValue(iRepValue->getValue(), iRepValue->getRepeat()-remaining));
		// New element
		line.insert(iRepValue, RepeatedValue(value, 1));
		// last half
		iRepValue->getRepeat() = remaining-1;

		--iRepValue;
		remaining = 1;
		return; // split previous + 1 + next
	}
}




template < typename T >
INLINE void ImageRLE<T>::setUnsafeValue(T value, int i, int j, int k)
{
	// Get the (y,z) line
	Line & line = this->getLine(j,k);

	// Loop till we reach the x coordinate
	int current_i = 0;
	typename Line::iterator iRepValue;
	for (iRepValue = line.begin(); iRepValue != line.end(); ++iRepValue)
	{
		current_i += iRepValue->getRepeat();
		// Check whether we finally arrived
		if (current_i > i)
		{
			// We have to create this temporary variable 'coz the following function takes a reference
			int tmp = current_i-i;
			this->setUnsafeValue(value, line, iRepValue, tmp);
			// Debugging stuff. TODO: remove this
			if (this->dim()[0] != this->lineLength(line))
			{
				throw std::runtime_error("ImageRLE length not preserved");
			}
			return;
		}
	}

	// Well -- we couldn't reach the desired x coordinate! This is probably
	// due to a bug.
	throw std::underflow_error("Incomplete ImageRLE (probably a bug) -- cannot set value");
}

} // namespace


#endif

