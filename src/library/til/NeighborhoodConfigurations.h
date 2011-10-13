#ifndef TIL_NEIGHBORHOOD_CONFIGURATION_H
#define TIL_NEIGHBORHOOD_CONFIGURATION_H

// includes from TIL library
#include "til/til_common.h"
#include "til/labels.h"
#include "til/Neighborhood.h"
#include "til/templateTools.h"


// Namespace 
namespace til {


// This class only collects code common to most morphological assessment
// class. Basically, it handles neighborhoods and foreground/background
// intensity values

template <class _TImage, class _TNeighborhood>
class PT_NeighborhoodConfigurations
{

public: // typedefs

	typedef _TImage	TImage;
	typedef _TNeighborhood TNeighborhood;
	typedef typename TImage::value_type value_type;
	typedef typename Iterator<TImage>::ConstVolumetric ConstVolumetricIterator;


public: // constructors & destructor

	explicit PT_NeighborhoodConfigurations(const TNeighborhood &nh)
	{
		m_foreground = 1;
		m_background = 0;
		m_nh = nh;
		//this->setNeighborhood(nh);
	}


public: // set & get

	void setForeground(value_type foreground) { m_foreground = foreground; }
	void setBackground(value_type background) { m_background = background; }

	value_type getForeground() const { return m_foreground; }
	value_type getBackground() const { return m_background; }


public: // functions

	// This function is used in if_next functions
	// It says that nothing is to be updated when going to the next pixel
	static void next() {}


protected: // data

	// The neighborhood
	TNeighborhood m_nh;

private: // data
	
	// Foreground and background values
	value_type m_foreground;
	value_type m_background;
};


template <class TVolumetricIterator, class _TNeighborhood = Neighborhood>
class NeighborhoodConfigurations : public ImageFunctor_label
{
public: // typedefs

	typedef _TNeighborhood TNeighborhood;
	typedef typename TVolumetricIterator::value_type value_type;

public: // constructors & destructor


	NeighborhoodConfigurations()
		: m_foreground(1), m_background(0), m_nh()
	{}

	// TODO: when value_type is not numeric, it might be that this class won't
	// even compile because of the default values -- which is bad
	NeighborhoodConfigurations(const TNeighborhood &nh)
		: m_foreground(1), m_background(0), m_nh(nh) {};

	NeighborhoodConfigurations(value_type foreground, value_type background, const TNeighborhood &nh = TNeighborhood())
		: m_foreground(foreground), m_background(background), m_nh(nh) {};

public: // set & get

	//void setForeground(value_type foreground) { m_foreground = foreground; }
	//void setBackground(value_type background) { m_background = background; }

	value_type getForeground() const { return m_foreground; }
	value_type getBackground() const { return m_background; }


protected: // data

	// neighborhood
	TNeighborhood m_nh;

private: // data
	
	// foreground value for all morphological considerations
	const value_type m_foreground;
	// background value for all morphological considerations
	const value_type m_background;
};





// This class checks whether a point is interior, i.e. if all its neighbors
// are also foreground.

template < class _TImage, class _TNeighborhood = Neighborhood >
class PT_IsInterior : public PT_NeighborhoodConfigurations<_TImage,_TNeighborhood>
{
public: // typedefs

	typedef _TImage	TImage;
	typedef _TNeighborhood TNeighborhood;
	typedef typename TImage::value_type value_type;
	typedef typename Iterator<TImage>::ConstVolumetric ConstVolumetricIterator;

public: // constructors & destructor

	PT_IsInterior(const TNeighborhood &nh = _TNeighborhood()) 
		: PT_NeighborhoodConfigurations<TImage,TNeighborhood>(nh) {}

public:

	INLINE bool _compute(const ConstVolumetricIterator &iIm) const
	{
		if (
			// First 6 neighbors
			this->cond<-1, 0, 0>(iIm) &&
			this->cond< 0,-1, 0>(iIm) &&
			this->cond< 0, 0,-1>(iIm) &&
			this->cond<+1, 0, 0>(iIm) &&
			this->cond< 0,+1, 0>(iIm) &&
			this->cond< 0, 0,+1>(iIm) &&

			// Next 12
			this->cond<-1,-1, 0>(iIm) &&
			this->cond<-1, 0,-1>(iIm) &&
			this->cond< 0,-1,-1>(iIm) &&
			this->cond<+1,+1, 0>(iIm) &&
			this->cond<+1, 0,+1>(iIm) &&
			this->cond< 0,+1,+1>(iIm) &&
			this->cond<-1,+1, 0>(iIm) &&
			this->cond<-1, 0,+1>(iIm) &&
			this->cond< 0,-1,+1>(iIm) &&
			this->cond<+1,-1, 0>(iIm) &&
			this->cond<+1, 0,-1>(iIm) &&
			this->cond< 0,+1,-1>(iIm) &&

			// Last 8
			this->cond<-1,-1,-1>(iIm) &&
			this->cond<+1,+1,+1>(iIm) &&
			this->cond<+1,-1,-1>(iIm) &&
			this->cond<-1,+1,-1>(iIm) &&
			this->cond<-1,-1,+1>(iIm) &&
			this->cond<-1,+1,+1>(iIm) &&
			this->cond<+1,-1,+1>(iIm) &&
			this->cond<+1,+1,-1>(iIm)
			)
		{
			return true;
		}
		else
		{
			return false;
		}
	}


	INLINE bool compute(const ConstVolumetricIterator &iIm) const
	{
		if (*iIm != this->getForeground())
		{
			return false;
		}

		return _compute(iIm);
	}

private: // functions

	template <int i, int j, int k>
		bool cond(const ConstVolumetricIterator &iIm) const
	{
			return (
				!this->m_nh.template isNeighbor<i,j,k>() ||	// There is no neighbor
				!containsNeighbor<i,j,k>(iIm) ||	// Or it is not in image
				iIm.template getUnsafeValue<i,j,k>() == this->getForeground() // or it is foreground
				);
	}

};

// This class checks whether a point is interior, i.e. if all its neighbors
// are also foreground.

template < class _TImage, class _TNeighborhood = Neighborhood >
class PT_IsBoundary : public PT_NeighborhoodConfigurations<_TImage,_TNeighborhood>
{
public: // typedefs

	typedef _TImage	TImage;
	typedef _TNeighborhood TNeighborhood;
	typedef typename TImage::value_type value_type;
	typedef typename Iterator<TImage>::ConstVolumetric ConstVolumetricIterator;

public: // constructors & destructor

	PT_IsBoundary(const TNeighborhood &nh = _TNeighborhood()) 
		: PT_NeighborhoodConfigurations<TImage,TNeighborhood>(nh) {}

public:

	INLINE bool _compute(const ConstVolumetricIterator &iIm) const
	{
		if (
			// First 6 neighbors
			this->cond<-1, 0, 0>(iIm) ||
			this->cond< 0,-1, 0>(iIm) ||
			this->cond< 0, 0,-1>(iIm) ||
			this->cond<+1, 0, 0>(iIm) ||
			this->cond< 0,+1, 0>(iIm) ||
			this->cond< 0, 0,+1>(iIm) ||

			// Next 12
			this->cond<-1,-1, 0>(iIm) ||
			this->cond<-1, 0,-1>(iIm) ||
			this->cond< 0,-1,-1>(iIm) ||
			this->cond<+1,+1, 0>(iIm) ||
			this->cond<+1, 0,+1>(iIm) ||
			this->cond< 0,+1,+1>(iIm) ||
			this->cond<-1,+1, 0>(iIm) ||
			this->cond<-1, 0,+1>(iIm) ||
			this->cond< 0,-1,+1>(iIm) ||
			this->cond<+1,-1, 0>(iIm) ||
			this->cond<+1, 0,-1>(iIm) ||
			this->cond< 0,+1,-1>(iIm) ||

			// Last 8
			this->cond<-1,-1,-1>(iIm) ||
			this->cond<+1,+1,+1>(iIm) ||
			this->cond<+1,-1,-1>(iIm) ||
			this->cond<-1,+1,-1>(iIm) ||
			this->cond<-1,-1,+1>(iIm) ||
			this->cond<-1,+1,+1>(iIm) ||
			this->cond<+1,-1,+1>(iIm) ||
			this->cond<+1,+1,-1>(iIm)
			)
		{
			return true;
		}
		else
		{
			return false;
		}
	}


	INLINE bool compute(const ConstVolumetricIterator &iIm) const
	{
		if (*iIm != this->getForeground())
		{
			return false;
		}

		return _compute(iIm);
	}

private: // functions

  template <int i, int j, int k>
  bool cond(const ConstVolumetricIterator &iIm) const
  {
    return (
            this->m_nh.template isNeighbor<i,j,k>() &&                     // There is a neighbor...
            containsNeighbor<i,j,k>(iIm) &&                                // ...and it is inside the image...
            iIm.template getUnsafeValue<i,j,k>() != this->getForeground()  // ...and it is not foreground
           );
	}

};


template < class _TImage, class _TNeighborhood = Neighborhood >
class PT_IsIsolated : public PT_NeighborhoodConfigurations<_TImage, _TNeighborhood>
{
public: // typedefs

	typedef _TImage	TImage;
	typedef _TNeighborhood TNeighborhood;
	typedef typename TImage::value_type value_type;
	typedef typename Iterator<TImage>::ConstVolumetric ConstVolumetricIterator;

public: // constructors & destructor

	PT_IsIsolated(const TNeighborhood &nh = _TNeighborhood())
		: PT_NeighborhoodConfigurations<TImage,TNeighborhood>(nh) {}


public: // functions

	INLINE bool _compute(const ConstVolumetricIterator &iIm) const
	{
		if (

			// First 6 neighbors
			this->cond<-1, 0, 0>(iIm) &&
			this->cond< 0,-1, 0>(iIm) &&
			this->cond< 0, 0,-1>(iIm) &&
			this->cond<+1, 0, 0>(iIm) &&
			this->cond< 0,+1, 0>(iIm) &&
			this->cond< 0, 0,+1>(iIm) &&

			// Next 12
			this->cond<-1,-1, 0>(iIm) &&
			this->cond<-1, 0,-1>(iIm) &&
			this->cond< 0,-1,-1>(iIm) &&
			this->cond<+1,+1, 0>(iIm) &&
			this->cond<+1, 0,+1>(iIm) &&
			this->cond< 0,+1,+1>(iIm) &&
			this->cond<-1,+1, 0>(iIm) &&
			this->cond<-1, 0,+1>(iIm) &&
			this->cond< 0,-1,+1>(iIm) &&
			this->cond<+1,-1, 0>(iIm) &&
			this->cond<+1, 0,-1>(iIm) &&
			this->cond< 0,+1,-1>(iIm) &&

			// Last 8
			this->cond<-1,-1,-1>(iIm) &&
			this->cond<+1,+1,+1>(iIm) &&
			this->cond<+1,-1,-1>(iIm) &&
			this->cond<-1,+1,-1>(iIm) &&
			this->cond<-1,-1,+1>(iIm) &&
			this->cond<-1,+1,+1>(iIm) &&
			this->cond<+1,-1,+1>(iIm) &&
			this->cond<+1,+1,-1>(iIm)
			)
		{
			return true;
		}
		else
		{
			return false;
		}		
	}

	
	INLINE bool compute(const ConstVolumetricIterator &iIm) const
	{
		if (*iIm != this->getForeground())
		{
			return false;
		}

		return this->_compute(iIm);
	}

private: // functions


	template <int i, int j, int k>
		bool cond(const ConstVolumetricIterator &iIm) const
	{
		return (
			!this->m_nh.template isNeighbor<i,j,k>() ||                      // There is no neighbor
			!containsNeighbor<i,j,k>(iIm) ||                                 // or it is not in image
			iIm.template getUnsafeValue<i,j,k>() == this->getBackground()    // or it is background
			);
	}
};


template < class TVolumetricIterator, class _TNeighborhood = Neighborhood >
class IsIsolated :
	public NeighborhoodConfigurations<TVolumetricIterator, _TNeighborhood>, 
	public std::unary_function<TVolumetricIterator, bool>
{
public: // typedefs

	// The other type is empty, so we discard it as a "real" base and use
	// simple inheritance conventions.
	typedef NeighborhoodConfigurations<TVolumetricIterator, _TNeighborhood> Base;
	//typedef _TImage TImage;
	typedef _TNeighborhood TNeighborhood;
	typedef typename TVolumetricIterator::value_type value_type;
  // This redefinition is necessary, I think because of multiple inheritance.
  typedef bool result_type;
	//typedef typename TImage::value_type value_type;
	//typedef typename Iterator<TImage>::ConstVolumetric ConstVolumetricIterator;

public: // constructors & destructor

	IsIsolated() : Base() {}
	explicit IsIsolated(const TNeighborhood &nh) : Base(nh) {}
	IsIsolated(value_type foreground, value_type background, const TNeighborhood &nh = _TNeighborhood())
		: Base(foreground, background, nh) {}

public: // functions

	INLINE
//	typename enable_if<is_VolumetricImageIterator<VolumetricImageIterator>, bool>::type
	bool
	_compute(const TVolumetricIterator &iIm) const
	{
		if (
			// First 6 neighbors
			this->cond<-1, 0, 0>(iIm) &&
			this->cond< 0,-1, 0>(iIm) &&
			this->cond< 0, 0,-1>(iIm) &&
			this->cond<+1, 0, 0>(iIm) &&
			this->cond< 0,+1, 0>(iIm) &&
			this->cond< 0, 0,+1>(iIm) &&

			// Next 12
			this->cond<-1,-1, 0>(iIm) &&
			this->cond<-1, 0,-1>(iIm) &&
			this->cond< 0,-1,-1>(iIm) &&
			this->cond<+1,+1, 0>(iIm) &&
			this->cond<+1, 0,+1>(iIm) &&
			this->cond< 0,+1,+1>(iIm) &&
			this->cond<-1,+1, 0>(iIm) &&
			this->cond<-1, 0,+1>(iIm) &&
			this->cond< 0,-1,+1>(iIm) &&
			this->cond<+1,-1, 0>(iIm) &&
			this->cond<+1, 0,-1>(iIm) &&
			this->cond< 0,+1,-1>(iIm) &&

			// Last 8
			this->cond<-1,-1,-1>(iIm) &&
			this->cond<+1,+1,+1>(iIm) &&
			this->cond<+1,-1,-1>(iIm) &&
			this->cond<-1,+1,-1>(iIm) &&
			this->cond<-1,-1,+1>(iIm) &&
			this->cond<-1,+1,+1>(iIm) &&
			this->cond<+1,-1,+1>(iIm) &&
			this->cond<+1,+1,-1>(iIm)
			)
		{
			return true;
		}
		else
		{
			return false;
		}		
	}

	INLINE
	//typename enable_if<is_VolumetricImageIterator<VolumetricImageIterator>, bool>::type
	bool
	operator()(const TVolumetricIterator &iIm) const
	{
		if (*iIm != this->getForeground())
		{
			return false;
		}

		return this->_compute(iIm);
	}

private: // functions


	template <int i, int j, int k >
	bool cond(const TVolumetricIterator &iIm) const
	{
		return (
			!this->m_nh.template isNeighbor<i,j,k>() ||                    // There is no neighbor
			!containsNeighbor<i,j,k>(iIm) ||                               // or it is not in image
			iIm.template getUnsafeValue<i,j,k>() == this->getBackground()  // or it is background
			);
	}
};


template < class _TImage, class _TNeighborhood = Neighborhood >
class PT_HasNBackgroundNeighbors : public PT_NeighborhoodConfigurations<_TImage,_TNeighborhood>
{
public: // typedefs

	typedef _TImage	TImage;
	typedef _TNeighborhood TNeighborhood;
	typedef typename TImage::value_type value_type;
	typedef typename Iterator<TImage>::ConstVolumetric ConstVolumetricIterator;

public: // constuctors & destructor

	PT_HasNBackgroundNeighbors(int nNeighbors, const TNeighborhood &nh = _TNeighborhood())
		: PT_NeighborhoodConfigurations<TImage,TNeighborhood>(nh)
	{
		this->setNNeighbors(nNeighbors);
	}


public: // set & get

	void setNNeighbors(int nNeighbors) { m_nNeighbors = nNeighbors; }
	int getNNeighbors() { return m_nNeighbors; }


public: // functions

	INLINE bool _compute(const ConstVolumetricIterator &iIm) const
	{
		int count = 0;
		
		// First 6 neighbors
		if (this->cond<-1, 0, 0>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		if (this->cond< 0,-1, 0>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		if (this->cond< 0, 0,-1>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }

		if (this->cond<+1, 0, 0>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		if (this->cond< 0,+1, 0>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		if (this->cond< 0, 0,+1>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
				
		// Next 12
		if (this->cond<-1,-1, 0>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		if (this->cond<-1, 0,-1>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		if (this->cond< 0,-1,-1>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		
		if (this->cond<+1,+1, 0>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		if (this->cond<+1, 0,+1>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		if (this->cond< 0,+1,+1>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		
		if (this->cond<-1,+1, 0>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		if (this->cond<-1, 0,+1>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		if (this->cond< 0,-1,+1>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		
		if (this->cond<+1,-1, 0>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		if (this->cond<+1, 0,-1>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		if (this->cond< 0,+1,-1>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
				
		// Last 8
		if (this->cond<-1,-1,-1>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		if (this->cond<+1,+1,+1>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		
		if (this->cond<+1,-1,-1>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		if (this->cond<-1,+1,-1>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		if (this->cond<-1,-1,+1>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		
		if (this->cond<-1,+1,+1>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		if (this->cond<+1,-1,+1>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		if (this->cond<+1,+1,-1>(iIm)) { ++count; if (count >= m_nNeighbors) return true; }
		
		return false;
	}
	
	
	INLINE bool compute(const ConstVolumetricIterator &iIm) const
	{
		if (*iIm != this->getForeground())
		{
			return false;
		}

		return this->_compute(iIm);
	}

private: // functions

	template < int i, int j, int k >
	INLINE bool cond(const ConstVolumetricIterator &iIm) const
	{
		return (
			this->m_nh.template isNeighbor<i,j,k>() &&  // There is a neighbor
			containsNeighbor<i,j,k>(iIm) &&  // This neighbor is inside
			iIm.template getUnsafeValue<i,j,k>() == this->getBackground()  // and its value is background
			 );
	}

private:

	int m_nNeighbors;
};


} // namespace


#endif

