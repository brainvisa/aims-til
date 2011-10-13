#ifndef TIL_BASIC_PIXEL_TESTS_H
#define TIL_BASIC_PIXEL_TESTS_H


// Local includes
#include "til/til_common.h"
#include "til/numeric_array.h"

////////////////////////////////////////////////////////////////////////////////
//
// Classes that make simple tests on pixels, used in loops of some functions
// that enable to template over such classes (e.g. region growing).
//
// For consistency between iterator approaches and position approaches, the
// image has always to be passed as an argument, i.e. the pixel class 
// discovers every time the image it has to process. It comes to the price that
// no preprocessing can be done inside a test class. On the other hand, 
// iterators can be used without fear that the image being processed is 
// actually different from the image that the iterator is refering to.
//
////////////////////////////////////////////////////////////////////////////////


namespace til
{


// This class is not a test, it merely collects common code over a range of
// pixel test class.

/*
template < class TImage >
class PT_NoValue
{

public: // typedefs

	typedef typename TImage::ConstVolumetricIterator ConstVolumetricIterator;
	typedef typename TImage::value_type value_type;

public: // constructors & destructor

	PT_NoValue() {}

public: // set & get

	void setImage(ConstPtr<TImage> &im) { m_im = im; }
	const ConstPtr<TImage> & image() const { return m_im; }

public: // functions

	static void next() {}

private: // data

	ConstPtr<TImage> m_im;
};



template < class TImage >
class PT_IsOnImageBorder : public PT_NoValue<TImage>
{

public: // function

	bool compute(const ConstVolumetricIterator &iIm) const
	{
		return this->_compute(iIm.pos(), iIm.image()->getDim());
	}

	bool compute(const Vector<int,3> &pos) const
	{
		return this->_compute(this->image()->getValue(pos), m_im.size());
	}

private: // function

	bool _compute(const Vector<int,3> &pos, const Vector<int,3> &imDim) const
	{
		return (
			(pos.getX() == 0) ||
			(pos.getY() == 0) ||
			(pos.getZ() == 0) ||
			(pos.getX() == imDim.getX()-1) ||
			(pos.getY() == imDim.getY()-1) ||
			(pos.getZ() == imDim.getZ()-1));
	}
};
*/

// This class is not a test, it merely collects common code over a range of
// pixel test class.

template < class _TImage >
class PT_SingleValue
{

public: // typedefs

	typedef _TImage TImage;
	typedef typename TImage::value_type value_type;
	typedef typename Iterator<TImage>::ConstVolumetric ConstVolumetricIterator;

public: // constructors & destructor

	PT_SingleValue() { this->setValue(0); }
	PT_SingleValue(value_type value) { this->setValue(value); }

public: // set & get

	void setValue(value_type value) { m_value = value; }
	value_type getValue() const { return m_value; }

public: // functions

	static void next() {}

private: // data

	value_type m_value;
};



// Test is intensity is below some fixed threshold

template < class TImage >
class PT_IsBelow : public PT_SingleValue<TImage>
{
public: // typedefs

	typedef typename TImage::value_type value_type;
	typedef typename Iterator<TImage>::ConstVolumetric ConstVolumetricIterator;

public: // constructors & destructor

	PT_IsBelow() : PT_SingleValue<TImage>() {}
	PT_IsBelow(value_type value) : PT_SingleValue<TImage>(value) {}

public: // function

	bool compute(const ConstVolumetricIterator &iIm) const
	{
		return this->_compute(*iIm);
	}

	bool compute(const TImage &im, const numeric_array<int,3> &pos) const
	{
		return this->_compute(im.getValue(pos));
	}

private: // function

	bool _compute(value_type value) const
	{
		return value <= this->getValue();
	}
};


// Test if intensity is above some fixed threshold

template < class TImage >
class PT_IsAbove : public PT_SingleValue<TImage>
{
public: // typedefs

	typedef typename TImage::value_type value_type;
	typedef typename Iterator<TImage>::ConstVolumetric ConstVolumetricIterator;

public: // constructors & destructor

	PT_IsAbove() : PT_SingleValue<TImage>() {}
	PT_IsAbove(value_type value) : PT_SingleValue<TImage>(value) {}

public: // function

	bool compute(const ConstVolumetricIterator & iIm) const
	{
		return this->_compute(*iIm);
	}

	bool compute(const TImage & im, const numeric_array<int,3> & pos) const
	{
		return this->_compute(im(pos));
	}

private: // function

	bool _compute(value_type value) const
	{
		return value >= this->getValue();
	}
};


// Test if intensity is equal to some fixed value

template < class TImage >
class PT_IsEqual : public PT_SingleValue<TImage>
{
public: // typedefs

	typedef typename TImage::value_type value_type;
	typedef typename Iterator<TImage>::ConstVolumetric ConstVolumetricIterator;

public: // constructors & destructor
	
	PT_IsEqual() : PT_SingleValue<TImage>() {}
	PT_IsEqual(value_type value) : PT_SingleValue<TImage>(value) {}

public: // function

	bool compute(const ConstVolumetricIterator & iIm) const
	{
		return this->_compute(*iIm);
	}

	bool compute(const TImage &im, const numeric_array<int,3> & pos) const
	{
		return this->_compute(im(pos));
	}

private: // function

	bool _compute(value_type value) const
	{
		return value == this->getValue();
	}
};


// Test if intensity is different from a given fixed value
/*
template < class TImage >
class PT_IsNotEqual : public PT_SingleValue<TImage>
{

public: // constructors & destructor

	PT_IsNotEqual() : PT_SingleValue<TImage>() {}
	PT_IsNotEqual(value_type value) : PT_SingleValue<TImage>(value) {}

public: // function

	bool compute(const ConstVolumetricIterator &iIm) const
	{
		return this->_compute(*iIm);
	}

	bool compute(const ConstPtr<TImage> &im, const Vector<int,3> &pos) const
	{
		return this->_compute(im.getValue(pos));
	}

private: // function

	bool _compute(value_type value) const
	{
		return value != this->getValue();
	}
};
*/
} // namespace

#endif


