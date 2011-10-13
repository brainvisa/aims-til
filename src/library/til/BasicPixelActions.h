#ifndef TIL_BASIC_PIXEL_ACTIONS_H
#define TIL_BASIC_PIXEL_ACTIONS_H

// Local includes
#include "til/til_common.h"
#include "til/numeric_array.h"


// Namespace 

namespace til {

	// predeclaration
	template < class TImage > struct Iterator;

template < class _TImage >
class PA_SingleValue
{

public: // typedefs

    typedef _TImage TImage;
	typedef typename TImage::value_type value_type;
	typedef typename Iterator<TImage>::Volumetric VolumetricIterator;

public: // constuctors & destructor

	PA_SingleValue() { this->setValue(0); }
	PA_SingleValue(value_type value) { this->setValue(value); }

public: // set & get

	void setValue(value_type value) { m_value = value; }
	void setImage(const TImage &im) { m_im.shallowCopy(im); }

	value_type getValue() const { return m_value; }
	TImage & image() const { return m_im; }

public: // functions

	static void next() {}

private: // data

	TImage m_im;
	value_type m_value;
};




template < class TImage >
class PA_SetValue : public PA_SingleValue<TImage>
{
public:
	typedef typename TImage::value_type value_type;
	typedef typename Iterator<TImage>::Volumetric VolumetricIterator;

public: // constuctors & destructor

	PA_SetValue(value_type value) : PA_SingleValue<TImage>(value) {}


public: // function

	
	INLINE void apply(const VolumetricIterator &iIm) const
	{
		*iIm = this->getValue();
	}
	

	void apply(const numeric_array<int,3> &pos) const
	{
		(*(this->image()))(pos) = this->getValue();
	}
};


template < class TImage >
class PA_CopyImage
{

public: // typedefs

	typedef typename Iterator<TImage>::Volumetric VolumetricIterator;
	typedef typename TImage::value_type value_type;

public: // constuctors & destructor

	PA_CopyImage(TImage &im) : m_iIm(im) {}

public: // functions

	void next() { ++m_iIm; }

	void apply(const VolumetricIterator &iIm)
	{
		*iIm = *m_iIm;
	}

	void apply(const numeric_array<int,3> &pos)
	{
		(*(this->image()))(pos) = m_im.getValue(pos);
	}

private: // data

	TImage m_im;
	typename Iterator<TImage>::Linear m_iIm;
};


} // namespace


#endif

