#ifndef TIL_IMAGENC_H
#define TIL_IMAGENC_H


// Includes from STL
#include <float.h>
#include <limits.h>
#include <string.h>
#include <vector>
#include <iostream>


// Includes from TIL library
#include "til/til_common.h"
#include "til/Box.h"
#include "til/ImageBase.h"
#include "til/image_common.h"

#include "til/ConstImageNCLinearIterator.h"
#include "til/ImageNCLinearIterator.h"
#include "til/ConstImageNCVolumetricIterator.h"
#include "til/ImageNCVolumetricIterator.h"


// Namespace 

namespace til {


/// Image class storing data slice-by-slice

template < typename T >
class ImageNC : public ImageBase
{
private: // typedefs

	typedef T** Buffer;

public: // typedefs

	typedef ImageNC<T> Self;
	typedef T value_type;


public: // constructors & destructor

	ImageNC();
	ImageNC(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz);
	ImageNC(const numeric_array<int,3> &dim, const numeric_array<t_voxsize,3> &vDim);
	ImageNC(T *data, int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz);	
	ImageNC(T **data, int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz);
	ImageNC(const std::vector<T*> & slices, int x, int y, t_voxsize vx, t_voxsize vy, t_voxsize vz);
	ImageNC(const ImageParameter & param);
	ImageNC(const ImageNC<T> *);
	
	/*
	/// Default constructor.
	/// Buffer is not allocated, every other variables are set to default value
	static Self* New() { return new Self(); }
	/// Allocate an image of dimension (x, y, z) with voxel size (vx, vy, vz)
	static Self* New(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz) { return new Self(x,y,z,vx,vy,vz); }
	/// Same as above, using vectors
	static Self* New(const Vector<int,3> &dim, const Vector<t_voxsize,3> &vDim) { return new Self(dim,vDim); }
	/// Create an image of dimension (x, y, z) with voxel size (vx, vy, vz)
	/// with data 'data'
	/// NB: the data is not copied, ImageNC is just a wrapper around
	/// the contiguous buffer (hence it is passed as a non-const pointer).
	static Self* New(T *data, int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz) { return new Self(data,x,y,z,vx,vy,vz); }
	static Self* New(T **data, int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz) {return new Self(data,x,y,z,vx,vy,vz); }
	/// Create an image from a list of slice pointers.
	/// The z dimension is deduced from the length of the list
	static Self* New(const std::vector<T*> & slices, int x, int y, t_voxsize vx, t_voxsize vy, t_voxsize vz) { return new Self(slices,x,y,vx,vy,vz); }
	/// Creates an image from an ImageParamter structure
	static Self* New(const ImageParameter & param) { return new Self(param); }
	*/

	/// Destructor
	~ImageNC();
	

public: // initialization

	// Initialization functions
	// Basically similar to constructors
	void init();
	void init(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz);
	void init(T *data, int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz);
	void init(T **data, int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz);
	void init(const numeric_array<int,3> &dim, const numeric_array<t_voxsize,3> &vDim);
	void init(const std::vector<T*> & slices, int x, int y, t_voxsize vx, t_voxsize vy, t_voxsize vz);
	void init(const ImageNC<T> *);
	void init(const ImageParameter &param);

public: // functions

	void setValue(const T & value, int i, int j, int k)
	{
		*(this->getUnsafePointerOf(i,j,k)) = value;
	}
  void setValue(const T & value, const numeric_array<int,3> & p)
  {
		*(this->getUnsafePointerOf(p)) = value;
  }

	// Get value

	// Access (read-write) the value at point (i,j,k)
	// Range checking is done
	INLINE T& operator()(int i, int j, int k);
	T& operator()(const numeric_array<int,3> &v) { return this->operator()(EXPAND_VECTOR(v)); }

	// Get the value at point (i,j,k) for read-only purpose
	// If point lies outside image, an exception is thrown
	INLINE const T getValue(int i, int j, int k) const;
	const T getValue(const numeric_array<int,3> &v) const { return this->getValue(EXPAND_VECTOR(v)); }
	
	// Get the value at point (i,j,k) for fast read-only purpose
	// WARNING: No range checking!
	INLINE const T getUnsafeValue(int i, int j, int k) const;
	const T getUnsafeValue(const numeric_array<int,3> &v) const { return this->getUnsafeValue(EXPAND_VECTOR(v)); }
	

public: // functions

	// Check whether data has been allocated or not
	bool isAllocated() const { return (m_data!=0);}
	
	// Set all values to zero
	void reset();

	void copy(const Self &im)
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

		// TODO: Replace this by iterators
		for (int i = 0; i < this->dim()[2]; ++i)
		{
			memcpy(
				(void*)this->getUnsafeSlicePointer(i),
				(void*)(const_cast<T*>(im.getUnsafeSlicePointer(i))),
				im.getSliceSize()*sizeof(T));
		}
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
	friend class ConstLinearIterator<Self>;


private:	// methods
		
	// Get pointer to data
	// Use to quickly port existing code
	// For any other job, use operator(), getValue, or iterators

	// For very dirty jobs only
	//T** getPointer() { return m_data;} 
	// For dirty jobs that do not mess with the data
	//const T* const * getConstPointer() const { return m_data;}
	// For dirty jobs that want range checking
	//const T* getConstPointerOf(int i, int j, int k) const;

	const T* getUnsafeSlicePointer(int i) const { return m_data[i]; }
	T*	      getUnsafeSlicePointer(int i)		 { return m_data[i]; }

	const T* getSlicePointer(int i) const
	{
		if (i < 0 || i >= this->dim()[2])
		{
			throw std::out_of_range("Slice number out of range");
		}
		return this->getUnsafeSlicePointer(i);
	}

	//TODO: find a way to optimize those duplicate const/non-const get functions
	T* getSlicePointer(int i)
	{
		if (i < 0 || i >= this->dim()[2])
		{
			throw std::out_of_range("Slice number out of range");
		}
		return this->getUnsafeSlicePointer(i);
	}


	// Unsafe but fast access
	T* getUnsafePointerOf(int i, int j, int k) const { return m_data[k] + i + this->dim()[0] * j; }
	T* getUnsafePointerOf(const numeric_array<int,3> &v) const { return this->getUnsafePointerOf(EXPAND_VECTOR(v)); }
			
	T* getPointerOf(int i, int j, int k) const
	{
		if (!this->contains(numeric_array<int,3>(i,j,k)))
		{
			return 0;
		}
		return this->getUnsafePointerOf(i,j,k);
	}
	T* getPointerOf(const numeric_array<int,3> &v) const { return this->getPointerOf(EXPAND_VECTOR(v)); }
		
	void setData(T** data)
	{
		if (m_data == data) return;
		if (m_data) this->deallocateData();
		m_data = data;
	}

	void setData(const std::vector<T*> & slices)
	{
		if (m_data) this->deallocateData();
		m_data = new T*[slices.size()];
		
		for (int i=0; i < slices.size(); ++i)
		{
			m_data[i] = slices[i];
		}
	}

	void allocateData()
	{
		if ( this->dim()[0]<0 || this->dim()[1]<0 || this->dim()[2]<0)
	    {
			throw std::invalid_argument("Image size < 0");
	    }
		m_data = new T*[this->dim()[2]];
		for (int i = 0; i < this->dim()[2]; ++i)
		{
			m_data[i] = new T[this->getSliceSize()];
		}
	}

	void deallocateData()
	{
		for (int i=0; i < this->dim()[2]; ++i)
		{
			delete [] m_data[i];
		}
		delete [] m_data;
		m_data = 0;
	}

	// Common operations done in init methods
	void _init(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz);
	
	// Get number of element in a slice
	int getSliceSize() const { return this->dim()[0] * this->dim()[1]; }

	// Copy constructor is intentionally left undefined
	// Non-explicit image copying is unwanted!
	// use copy instead
	ImageNC(const ImageNC<T> &);

private: // data
	
	// Double pointer to data
	Buffer m_data;
};


/// Iterator traits of ImageNC
template < typename T >
struct Iterator<ImageNC<T> >
{
	typedef ConstLinearIterator<ImageNC<T> >		ConstLinear;
	typedef ConstVolumetricIterator<ImageNC<T> >	ConstVolumetric;
	typedef LinearIterator<ImageNC<T> >				Linear;
	typedef VolumetricIterator<ImageNC<T> >			Volumetric;
};































////////////////////////////////////////////////////
//
// CODE
//
////////////////////////////////////////////////////


// Constructors

template < typename T >
ImageNC<T>::ImageNC() : ImageBase(), m_data(0)
{
    this->init();
}



template < typename T >
void ImageNC<T>::init()
{
    m_data = 0;
}



template < typename T >
ImageNC<T>::ImageNC(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz) : ImageBase(), m_data(0)
{
    this->init(x, y, z, vx, vy, vz);
}



template < typename T >
void ImageNC<T>::init(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz)
{
	this->_init(x, y, z, vx, vy, vz);
	this->allocateData();
    this->reset();
}



template < typename T >
ImageNC<T>::ImageNC(const numeric_array<int,3> &dim, const numeric_array<t_voxsize,3> &vDim) : ImageBase(), m_data(0)
{
	this->init(dim, vDim);
}



template < typename T >
void ImageNC<T>::init(const numeric_array<int,3> &dim, const numeric_array<t_voxsize,3> &vDim)
{
	this->init(EXPAND_VECTOR(dim), EXPAND_VECTOR(vDim));
}


template < typename T >
ImageNC<T>::ImageNC(const ImageParameter &param) : ImageBase(), m_data(0)
{
	this->init(param);
}


template < typename T >
void ImageNC<T>::init(const ImageParameter &param)
{
	this->init(param.m_dim, param.m_vDim);
	//this->setOrigin(param.m_ori);
}

template < typename T >
ImageNC<T>::ImageNC(T *data, int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz) : ImageBase(), m_data(0)
{
    this->init(data, x, y, z, vx, vy, vz);
}


template < typename T >
void ImageNC<T>::_init(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz)
{
  this->set_dim(til::numeric_array<int,3>(x, y, z));
  this->set_vdim(til::numeric_array<float,3>(vx, vy, vz));
}


template < typename T >
void ImageNC<T>::init(T *data, int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz)
{
	this->_init(x, y, z, vx, vy, vz);
	T** data2 = simple2DoublePointer(data, x, y, z);
	this->setData(data2);
}

template < typename T >
ImageNC<T>::ImageNC(T **data, int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz) : ImageBase(), m_data(0)
{
	this->init(data, x, y, z, vx, vy, vz);
}

template < typename T >
void ImageNC<T>::init(T **data, int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz)
{
	this->_init(x, y, z, vx, vy, vz);
    this->setData(data);
}

template < typename T >
ImageNC<T>::ImageNC(const std::vector<T*> & slices, int x, int y, t_voxsize vx, t_voxsize vy, t_voxsize vz) : ImageBase(), m_data(0)
{
	this->init(slices, x, y, vx, vy, vz);
}

template < typename T >
void ImageNC<T>::init(const std::vector<T*> & slices, int x, int y, t_voxsize vx, t_voxsize vy, t_voxsize vz)
{
	this->_init(x, y, slices.size(), vx, vy, vz);
	this->setData(slices);
}

template < typename T >
ImageNC<T>::ImageNC(const ImageNC<T> *pIm) : ImageBase(), m_data(0)
{
    this->init(pIm);
}

template < typename T >
void ImageNC<T>::init(const ImageNC<T> *pIm)
{
    this->init(param(pIm));
	
	for (int i = 0; i < this->getZ(); ++i)
	{
		memcpy(
			(void*)this->getUnsafeSlicePointer(i), 
			(void*)pIm.getUnsafeConstSlicePointer(i),
			pIm.getSliceSize()*sizeof(T));
	}
}


template < typename T >
ImageNC<T>::~ImageNC()
{
	this->deallocateData();
}


template < typename T >
void ImageNC<T>::reset()
{
	for (int i=0; i<this->dim()[2]; ++i)
	{
		memset((void*)this->getUnsafeSlicePointer(i), 0, this->getSliceSize()*sizeof(T));
	}
}



template < typename T >
INLINE T& ImageNC<T>::operator()(int i, int j, int k)
{
	return *(this->getPointerOf(i,j,k));
}



//template < typename T >
//const T* ImageNC<T>::getConstPointerOf(int i, int j, int k) const 
//{
//	return this->getPointerOf(i,j,k);
//}



template < typename T >
INLINE const T ImageNC<T>::getValue(int i, int j, int k) const
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
INLINE const T ImageNC<T>::getUnsafeValue(int i, int j, int k) const
{
	return *(this->getUnsafePointerOf(i,j,k));
}


} // namespace

#endif

