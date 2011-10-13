#ifndef TIL_IMAGEBP_H
#define TIL_IMAGEBP_H

#pragma warning(disable:4244)


// Standard library includes

#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>


// Local includes

#include "til_common.h"

#include "Box.h"
#include "image_common.h"
#include "Ptr.h"
#include "SmartObject.h"


#include "ConstImageBPLinearIterator.h"
#include "ImageBPLinearIterator.h"
#include "ConstImageBPVolumetricIterator.h"
#include "ImageBPVolumetricIterator.h"


// Namespace 

namespace til {


typedef uchar t_value;



template < int N >
class ImageBP : public SmartObject
{

public: // constructors & destructor
	

	// Default constructor
	// Buffer is not allocated
	// Every other internal variables are set to zero
	ImageBP();


	// Allocate an image of dimension (x, y, z) with voxel size (vx, vy, vz)
	ImageBP(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz);


	// Same as above, using vectors
	ImageBP(const numeric_array<int,3> &dim, const numeric_array<t_voxsize,3> &vDim);
	

	// Create an image of dimension (x, y, z) with voxel size (vx, vy, vz)
	// with data 'data'
	// NB: the data is not copied, Image is just a wrapper around
	// the contiguous buffer (hence it is passed as t_value*, not const t_value*).
	ImageBP(t_value *data, int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz);
	
	
	// Copy pointed image
	ImageBP(const ImageBP<N> *);


	// Copy constructor is intentionally left undefined
	// Non-explicit image copying is unwanted!
	// use ::copy instead
	//Image(const Image<N> &);
	
	
	// Destructor
	virtual ~ImageBP();
	


public: // initialization
	
	// Initialization functions
	// Basically similar to constructors for situations where
	// the image already exists

	void init(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz);
	void init(uchar *data, int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz);
	void init(const numeric_array<int,3> & dim, const numeric_array<t_voxsize,3> & vDim);
	void init(const ImageBP<N> *);


public: // set & get

	// Get image size
	
	int getDimX() const { return m_dim.getX();}
	int getDimY() const { return m_dim.getY();}
	int getDimZ() const { return m_dim.getZ();}
	int getDim(int i) const { return m_dim.get(i); }
	const numeric_array<int,3> & getDim() const { return m_dim; }
	
	
	// Get total number of elements
	
	int getSize() const { return this->getDimX()*this->getDimY()*this->getDimZ(); }
	
	
	// Get voxel size
	
	t_voxsize getVx() const { return m_vDim.getX(); }
	t_voxsize getVy() const { return m_vDim.getY(); }
	t_voxsize getVz() const { return m_vDim.getZ(); }
	t_voxsize getVDim(int i) const { return m_vDim.get(i); }
	const numeric_array<t_voxsize,3> & getVDim() const { return m_vDim; }
	

	
public: // functions

	// Get value

	// Access (read-write) the value at point (i,j,k)
	// Range checking is done
	
	t_value& operator()(int i, int j, int k) { return *(this->getPointerOf(i,j,k)); }
	t_value& operator()(const numeric_array<int,3> &v) { return this->operator()(EXPAND_VECTOR(v)); }

	
	// Get the value at point (i,j,k) for read-only purpose
	// If point lies outside image, an exception is thrown

	INLINE const t_value getValue(int i, int j, int k) const;
	const t_value getValue(const numeric_array<int,3> &v) const { return this->getValue(EXPAND_VECTOR(v)); }
	

	// Get the value at point (i,j,k) for fast read-only purpose
	// WARNING: No range checking!

	INLINE const t_value getUnsafeValue(int i, int j, int k) const;
	const t_value getUnsafeValue(const numeric_array<int,3> &v) const { return this->getUnsafeValue(EXPAND_VECTOR(v)); }


	// Get pointer to data
	// Use to quickly port existing code
	// For any other job, use operator(), getValue, or iterators

	// For very dirty jobs only
	t_value* getPointer() { return m_data;} 
	// For dirty jobs that do not mess with the data
	const t_value* getConstPointer() const { return m_data;}
	// For dirty jobs that want range checking
	inline const t_value* getConstPointerOf(int i, int j, int k) const;
	
	
	// Check whether data has been allocated or not
	
	bool isAllocated() const { return (m_data!=0);}
	
	
	// Set all values to zero

	void reset();
	
	// Copy an image
	// NB: it is put here, rather than in a separate function, for
	// optimization reasons, that cannot be templated.
	// using memcopy is (supposedly) much faster than any iterator.
	// but of course the use of memcopy depends on the implementation
	// of the class
	// Also, there is intentionally no copy constructor or operator=
	// in order that copy is always explicit.

	void copy(const ConstPtr<ImageBP<N> > &im);
	

public: // friends

	friend class ConstImageBPVolumetricIterator<N>;


public: // typedefs

	typedef t_value TPixel;
	typedef ConstImageBPLinearIterator<N> ConstLinearIterator;
	typedef ImageBPLinearIterator<N> LinearIterator;
	typedef ConstImageBPVolumetricIterator<N> ConstVolumetricIterator;
	typedef ImageBPVolumetricIterator<N> VolumetricIterator;


private:	// methods
		
	void init();

	// Unsafe but fast access
	
	t_value* getUnsafePointerOf(int i, int j, int k) const { return m_data + i + m_dim.getX() * ( j + m_dim.getY()*k); }
	t_value* getUnsafePointerOf(const numeric_array<int,3> &v) const { return this->getUnsafePointerOf(EXPAND_VECTOR(v)); }
	t_value* getUnsafePointerOf(int n) const { return m_data + n;}
	
	t_value* getPointerOf(int i, int j, int k) const
	{
		if (!this->contains(i,j,k))
		{
			return 0;
		}
		return this->getUnsafePointerOf(i,j,k);
	}
	t_value* getPointerOf(const numeric_array<int,3> &v) const { return this->getPointerOf(EXPAND_VECTOR(v)); }
	
	void setDim(int x, int y, int z)
	{
	    if ((x<0)||(y<0)||(z<0))
	    {
			throw std::invalid_argument("Image size < 0");
	    }
    	m_dim.set(x, y, z);
	}
	
	void setVDim(t_voxsize vx, t_voxsize vy, t_voxsize vz)
	{
        if ((vx<0)||(vy<0)||(vz<0))
	    {
			throw std::invalid_argument("Voxel size < 0");
	    }
		m_vDim.set(vx, vy, vz);
	}
	
	void setData(t_value* data)
	{
		if (m_data == data) return;
		if (m_data) delete [] m_data;
		m_data = data;
	}


	// Test whether these coordinates fall into image range
	
	bool contains(int i, int j, int k) const
	{
		return (i >= 0 && j >= 0 && k >= 0 && 
			i < m_dim.getX() && 
			j < m_dim.getY() &&
			k < m_dim.getZ());
	}


private: // data
	

	// Image dimensions
	numeric_array<int,3> m_dim;
	
	// Voxel size
	numeric_array<t_voxsize,3> m_vDim;
	
	// Origin
	numeric_array<t_voxsize,3> m_ori;

	// Pointer to data
	t_value* m_data;
};
































////////////////////////////////////////////////////
//
// CODE
//
////////////////////////////////////////////////////


// Constructors

template < int N >
ImageBP<N>::ImageBP() : SmartObject(), m_data(0)
{
    this->init();
}



template < int N >
void ImageBP<N>::init()
{
    m_data = 0;
}



template < int N >
ImageBP<N>::ImageBP(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz) : SmartObject(), m_data(0)
{
    this->init(x, y, z, vx, vy, vz);
}



template < int N >
void ImageBP<N>::init(int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz)
{
    t_value * data = new t_value[x*y*z];
    this->init(data, x, y, z, vx, vy, vz);
    this->reset();
}



template < int N >
ImageBP<N>::ImageBP(const numeric_array<int,3> &dim, const numeric_array<t_voxsize,3> &vDim) : SmartObject(), m_data(0)
{
	this->init(dim, vDim);
}



template < int N >
void ImageBP<N>::init(const numeric_array<int,3> &dim, const numeric_array<t_voxsize,3> &vDim)
{
	this->init(EXPAND_VECTOR(dim), EXPAND_VECTOR(vDim));
}



template < int N >
ImageBP<N>::ImageBP(t_value *data, int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz) : SmartObject(), m_data(0)
{
    this->init(data, x, y, z, vx, vy, vz);
}



template < int N >
void ImageBP<N>::init(t_value *data, int x, int y, int z, t_voxsize vx, t_voxsize vy, t_voxsize vz)
{
    this->setDim(x, y, z);
    this->setVDim(vx, vy, vz);
    this->setData(data);
}



template < int N >
ImageBP<N>::ImageBP(const ImageBP<N> *pIm) : SmartObject(), m_data(0)
{
    this->init(pIm);
}



template < int N >
void ImageBP<N>::init(const ImageBP<N> *pIm)
{
    this->init(param(pIm));
	memcpy(
		(void*)this->getPointer(), 
		(void*)pIm->getConstPointer(),
		pIm->getSize()*sizeof(t_value));
}



template < int N >
ImageBP<N>::~ImageBP()
{
	// NB: deleting a null pointer is also safe
	delete[] m_data;
}



template < int N >
void ImageBP<N>::reset()
{
	memset((void*)m_data, 0, this->getSize()*sizeof(t_value));
}



template < int N >
INLINE const t_value* ImageBP<N>::getConstPointerOf(int i, int j, int k) const 
{
	return this->getPointerOf(i,j,k);
}



template < int N >
INLINE const t_value ImageBP<N>::getValue(int i, int j, int k) const
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



template < int N >
INLINE const t_value ImageBP<N>::getUnsafeValue(int i, int j, int k) const
{
	return *(this->getUnsafePointerOf(i,j,k));
}

template < int N >
void ImageBP<N>::copy(const ConstPtr<ImageBP<N> > &im)
{
	
	if (!(this->isAllocated() && til::isAllocated(im)))
	{
		throw std::invalid_argument("Unallocated image");
	}

	if (!((this->getDimX()  == im->getDimX()) &&
		  (this->getDimY()  == im->getDimY()) &&
		  (this->getDimZ()  == im->getDimZ()) &&
		  (this->getVx() == im->getVx()) &&
		  (this->getVy() == im->getVy()) &&
		  (this->getVz() == im->getVz())))
	{
		throw std::invalid_argument("Incompatible images");
	}
	
	memcpy(
		(void*)this->getPointer(), 
		(void*)(const_cast<t_value*>(im->getConstPointer())),
		im->getSize()*sizeof(t_value));
}


) // namespace

#endif


