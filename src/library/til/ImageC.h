#ifndef TIL_IMAGEC_H
#define TIL_IMAGEC_H

// Includes from STL
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

// Includes from BOOST
#include <boost/shared_array.hpp>

// Include from TIL library
#include "til/image_common.h"
#include "til/ImageBase.h"
#include "til/imageTools.h"
#include "til/numeric_array.h"

// Package includes
#include "til/imageIterator.h"
#include "til/ConstImageLinearIterator.h"
#include "til/ImageLinearIterator.h"
#include "til/ConstImageVolumetricIterator.h"
#include "til/ImageVolumetricIterator.h"


// Ignore specific warnings
#ifdef _MSC_VER
#pragma warning (push)
#pragma warning (disable:4244)
#endif

namespace til
{

/// Image class using contiguous memory.

template < typename T >
class ImageC : public ImageBase
{

public: // typedefs

	typedef ImageC<T> Self;
	typedef T value_type;


public: // constructors & destructor

	ImageC();
	ImageC(const numeric_array<int,3> & dim, const numeric_array<t_voxsize,3> & vdim);
	ImageC(const ImageParameter &param);
	ImageC(T *data, const numeric_array<int,3> & dim, const numeric_array<t_voxsize,3> & vdim);

	/// Shallow copy.
	ImageC(const Self & im)
	{
		this->shallowCopy(im);
	}

	/*
	/// Default constructor.
	/// Buffer is not allocated.
	/// Every other internal variables are set to zero.
	static Self* New() { return new Self(); }
	/// Allocate an image of dimension (x, y, z) with voxel size (vx, vy, vz).
	static Self* New(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz) { return new Self(x,y,z,vx,vy,vz); }
	/// Allocate an image of dimension 'dim' with voxel size 'vDim'.
	static Self* New(const numeric_array<int,3> &dim, const numeric_array<t_voxsize,3> &vDim) {return new Self(dim, vDim); }
	/// Creates an image from an ImageParameter structure.
	static Self* New(const ImageParameter &param) { return new Self(param); }
	/// Create an image of dimension (x, y, z) with voxel size (vx, vy, vz)
	/// with data 'data'.
	/// NB: the data is not copied, ImageC is just a wrapper around
	/// the contiguous buffer.
	static Self* New(value_type *data, int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz) { return new Self(data,x,y,z,vx,vy,vz); }
	*/

	/// Destructor.
	~ImageC();
	

public: // initialization
	
	// Initialization functions
	// Basically similar to constructors 
	//void init(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz);
	//void init(T *data, int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz);
  void init(T * data, const numeric_array<int,3> & dim, const numeric_array<t_voxsize,3> & vdim);
	void init(const numeric_array<int,3> & dim, const numeric_array<t_voxsize,3> & vdim);
	void init(const Self &im);
	void init(const ImageParameter &param);


public: // functions

	// Get value

	// Access (read-write) the value at point (i,j,k)
	// Range checking is done
	//T & operator()(int i, int j, int k) { return *(this->getPointerOf(i,j,k)); }
	T & operator()(const numeric_array<int,3> & p) { return *(this->getPointerOf(p)); }
  const T & operator()(const numeric_array<int,3> & p) const { return *(this->getPointerOf(p)); }

	// Set value
	//void setValue(const T & v, int i, int j, int k) { *(this->getPointerOf(i,j,k)) = v; }
	void setValue(const T & v, const numeric_array<int,3> &p) { (*this)(p)=v; }

	// Get the value at point (i,j,k) for read-only purpose
	// If point lies outside image, an exception is thrown
	//INLINE T getValue(int i, int j, int k) const;
	//T getValue(const numeric_array<int,3> &v) const { return this->getValue(EXPAND_VECTOR(v)); }
	
	// Get the value at point (i,j,k) for fast read-only purpose
	// WARNING: No range checking!
	//INLINE T getUnsafeValue(int i, int j, int k) const;
	T getUnsafeValue(const numeric_array<int,3> & p) const { return *(this->getUnsafePointerOf(p)); }


	// Get pointer to data
	// Use to quickly port existing code
	// For any other job, use operator(), getValue, or iterators

	// For very dirty jobs only
	T * getPointer() { return m_data.get(); }
	const T * getPointer() const { return m_data.get(); }
	// For dirty jobs that do not mess with the data
	//const value_type* getConstPointer() const { return m_data.get();}
	// For dirty jobs that want range checking
	//INLINE const value_type* getConstPointerOf(int i, int j, int k) const;
	
	/// Check whether data has been allocated or not.
	bool isAllocated() const { return (m_data!=0);}
	
	// TODO: myreset should become reset, reset should become clear
	void myreset()
	{
		this->init();
	}

	/// Set all values to default value (0 for numerical types)
	void reset();
	
	/// Copy an image.
	/// NB: it is put here, rather than in a separate function, for
	/// optimization reasons, that cannot be templated.
	/// Using memcopy is (supposedly) much faster than any iterator.
	/// But of course the use of memcopy depends on the implementation
	/// of the class.
	/// Also, there is intentionally no copy constructor or operator=
	/// in order that copy is always explicit.
	void copy(const Self & im);


	/// Shallow copy of an image.
	/// A shallow copy does not duplicate the buffer: it points to the same buffer
	/// as the original image
	void shallowCopy(const Self & im)
	{
		this->init(param(im));
		m_data = im.m_data;
	}

public: // operators


	/// Check whether images point on the same buffer
	bool operator!=(const Self & im) const
	{
		return m_data != im.m_data;
	}
	/// Check whether images point on the same buffer
	bool operator==(const Self & im) const
	{
		return m_data == im.m_data;
	}


public: // friends

	friend class ConstVolumetricIterator<Self>;
	// yields an internal compiler error under MSVC7.1
	/*
	template < class TImage>
	friend	typename enable_if_c<is_Image<TImage>::value, bool>::type
		operator==(const TImage &im1, const TImage &im2);
	template < class TImage>
	friend	typename enable_if_c<is_Image<TImage>::value, bool>::type
		operator!=(const TImage &im1, const TImage &im2);
	*/

private: // constructors
	
	/// Copy pointed image.
	ImageC(const Self *);	


private:	// methods

	void init();

	// Unsafe but fast access	
	//T * getUnsafePointerOf(int i, int j, int k) { return this->getPointer() + i + this->dim()[0] * ( j + this->dim()[1]*k); }
	T * getUnsafePointerOf(const numeric_array<int,3> & p)
  {
    return this->getPointer() + p[0] + this->dim()[0] * ( p[1] + this->dim()[1] * p[2]);
  }
	const T * getUnsafePointerOf(const numeric_array<int,3> & p) const
  {
    return this->getPointer() + p[0] + this->dim()[0] * ( p[1] + this->dim()[1] * p[2]);
  }
	T * getUnsafePointerOf(int n) { return this->getPointer() + n;}
	//const T * getUnsafePointerOf(int i, int j, int k) const { return this->getPointer() + i + this->dim()[0] * ( j + this->dim()[1]*k); }
	const T * getUnsafePointerOf(int n) const { return this->getPointer() + n;}
	
  /*
	T * getPointerOf(int i, int j, int k)
	{
		if (!this->contains(i,j,k))
		{
			return 0;
		}
		return this->getUnsafePointerOf(i,j,k);
	}
  */


  T * getPointerOf(const numeric_array<int,3> & p)
  {
    if (!this->contains(p)) return 0;
    return this->getUnsafePointerOf(p);
  }

  const T * getPointerOf(const numeric_array<int,3> & p) const
  {
    if (!this->contains(p)) return 0;
    return this->getUnsafePointerOf(p);
  }

	void setData(T * data)
	{
		//if (m_data == data) return;
		//if (m_data) delete [] m_data;
		//m_data = data;
		m_data.reset(data);
	}

	// Copy constructor is intentionally left undefined
	// Non-explicit image copying is unwanted!
	// use copy instead

private: // data
	
	// Pointer to data
	boost::shared_array<T> m_data;
};


/// Iterator traits of ImageC
template < typename T >
struct Iterator< ImageC<T> >
{
	typedef ConstLinearIterator< ImageC<T> >		ConstLinear;
	typedef ConstVolumetricIterator< ImageC<T> >	ConstVolumetric;
	typedef LinearIterator< ImageC<T> >				Linear;
	typedef VolumetricIterator< ImageC<T> >			Volumetric;
};



// Constructors

template < typename T >
ImageC<T>::ImageC() : ImageBase(), m_data(0)
{
    this->init();
}



template < typename T >
void ImageC<T>::init()
{
	// Null buffer, zero dimension, zero voxel size
  this->init(numeric_array<int,3>(0,0,0), numeric_array<t_voxsize,3>(0,0,0));
}


/*
template < typename T >
ImageC<T>::ImageC(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz) : ImageBase(), m_data(0)

{
    this->init(x, y, z, vx, vy, vz);
}
*/


template < typename T >
void ImageC<T>::init
(
 const numeric_array<int,3> & dim,
 const numeric_array<t_voxsize,3> & vdim
)
{
	// NB: no catch(std::bad_alloc& ex) is done here
    T * data = new T[dim[0]*dim[1]*dim[2]];
    this->init(data, dim, vdim);
    this->reset();
}


template < typename T >
ImageC<T>::ImageC(const numeric_array<int,3> &dim, const numeric_array<t_voxsize,3> &vDim) : ImageBase(), m_data(0)
{
	this->init(dim, vDim);
}


template < typename T >
ImageC<T>::ImageC(const ImageParameter &param) : ImageBase(), m_data(0)
{
	this->init(param);
}


template < typename T >
void ImageC<T>::init(const ImageParameter &param)
{
	this->init(param.m_dim, param.m_vDim);
	//this->setOrigin(param.m_ori);
}


template < typename T >
ImageC<T>::ImageC(T *data, const numeric_array<int,3> & dim, const numeric_array<t_voxsize,3> & vdim) : ImageBase(), m_data(0)
{
  this->init(data, dim, vdim);
}



template < typename T >
void ImageC<T>::init(T *data, const numeric_array<int,3> & dim, const numeric_array<t_voxsize,3> & vdim)
{
  this->set_dim(dim);
  this->set_vdim(vdim);
  this->setData(data);
}


/*
template < typename T >
ImageC<T>::ImageC(const ImageC<T> & im) : ImageBase(), m_data(0)
{
  this->init(im);
}
*/

// I comment this first to see if anybody is using this function.
// If not then I can safely change its content to shallowCopy
/*
template < typename T >
void ImageC<T>::init(const ImageC<T> & im)
{
  
  this->init(param(im));
  memcpy((void*)this->getPointer(), 
         (void*)im.getPointer(),
         im.size()*sizeof(T));
  
}
*/

template < typename T >
ImageC<T>::~ImageC()
{
  // NB: deleting a null pointer is also safe
  // delete[] m_data;
}



template < typename T >
void ImageC<T>::reset()
{
  // use a memset for numeric types only
  if (std::numeric_limits<T>::is_specialized)
  {
    memset((void*)this->getPointer(), 0, this->size()*sizeof(T));
  }
  // For other types, just use a regular loop with default constructor
  else
  {
    for (int i = 0; i < this->size(); ++i)
    {
      m_data[i] = T();
    }
  }
}


/*
template < typename T >
INLINE const T* ImageC<T>::getConstPointerOf(int i, int j, int k) const 
{
  return this->getPointerOf(i,j,k);
}
*/

/*
template < typename T >
INLINE T ImageC<T>::getValue(int i, int j, int k) const
{
  if (this->contains(i,j,k))
  {
    return this->getUnsafeValue(i,j,k);
  }
  else
  {
    throw std::out_of_range("Position out of image range");
  }
}
*/

/*
template < typename T >
INLINE T ImageC<T>::getUnsafeValue(int i, int j, int k) const
{
	return *(this->getUnsafePointerOf(i,j,k));
}
*/

template < typename T >
void ImageC<T>::copy(const ImageC<T> &im)
{
	
	if (!(this->isAllocated() && til::isAllocated(im)))
	{
		throw std::invalid_argument("Unallocated image");
	}

	if (!((this->dim()[0]  == im.dim()[0]) &&
		  (this->dim()[1]  == im.dim()[1]) &&
		  (this->dim()[2]  == im.dim()[2]) &&
		  (this->vdim()[0] == im.vdim()[0]) &&
		  (this->vdim()[1] == im.vdim()[1]) &&
		  (this->vdim()[2] == im.vdim()[2])))
	{
		throw std::invalid_argument("Incompatible images");
	}
	
	memcpy(
		(void*)this->getPointer(), 
		(void*)(const_cast<T*>(im.getPointer())),
		im.size()*sizeof(T));
}


} // namespace


#ifdef _MSC_VER
#pragma warning (pop)
#endif

#endif

