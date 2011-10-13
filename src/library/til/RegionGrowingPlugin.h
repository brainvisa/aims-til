#ifndef TIL_REGION_GROWING_PLUGIN_H
#define TIL_REGION_GROWING_PLUGIN_H

// includes from STL
#include <numeric>
#include <vector>

// includes from TIL library
#include "til/til_common.h"
#include "til/numeric_array.h"
#include "til/SymMatrix3.h"



// The region growing plugins are classes that can be used during
// the region growing process to compute some numbers, e.g. the
// size or the maximum intensity.

// Such numbers can also be computed after the region growing is
// done. But one might prefer to compute them on the fly.
// Especially when the region growing itself depends
// on these numbers (say, we stop region growing if a maximum
// size has been reached). This is after all the big feature
// of region growing over simple connected components algorithms.


namespace til
{

  template < typename TPoint >
  class RegionCentroid
  {
  public: // constructors
  
  	RegionCentroid() { this->init(); };
  
  public: // initializations
  
    void init()
    {
      m_nElem = 0;
      std::fill(m_sum.begin(), m_sum.end(), 0);
    }
  
  public: // set & get
  
    /// Get centroid.
    TPoint get() { return m_sum / m_nElem; }
  
    /*
    // TODO: this is really stupid, return v instead.
  	void getCentroid(numeric_array<double,3> & v)
  	{
      v = m_sum * (1.0/ m_nElem);
      / *
  		v[0] = m_sum[0] / m_nElem;
  		v[1] = m_sum[1] / m_nElem;
  		v[2] = m_sum[2] / m_nElem;
      * /
  	}
    */
  
  public: // functions
  
  
    template < typename TIterator >
    void update(TIterator begin, TIterator end)
    {
      for (; begin != end; ++begin)
      {
        m_sum += *begin;
      }
    }
    /*
  	void update(const std::vector<numeric_array<int,3> > & pl)
  	{
      //m_sum = std::accumulate(pl.begin(), pl.end(), m_sum, std::plus<numeric_array<double,3> >());
      //m_sum = std::accumulate(pl.begin(), pl.end(), m_sum, functor::Add<numeric_array<std::plus<numeric_array<double,3> >());
  		//add(pl, m_sum);
  		m_nElem += pl.size();
  	}
    */
  
  private: // data
  
  	//numeric_array<double,3> m_sum;
    TPoint m_sum;
  	std::size_t m_nElem;
  };



template < class TImage >
class RegionMaxInt
{
public: // typedefs

	typedef typename TImage::value_type value_type;

public: // constuctors & destructor

	RegionMaxInt() { m_max = std::numeric_limits<value_type>::min(); }

	RegionMaxInt(const TImage &im)
	{
		this->init(im);
	}

public: // set & get

	value_type getMaxInt() { return m_max; }

public: // functions

	void init(const TImage &im)
	{
		m_im.shallowCopy(im);
		m_max = std::numeric_limits<value_type>::min();
	}

	void update(const numeric_array<int,3> & v)
	{
		if (m_im(v) > m_max)
		{
			m_max = m_im(v);
		}
	}


private: // data

	TImage m_im;
	value_type m_max;
};


/*
class TIL_API RegionBoundingBox
{

public: // constructors & destructor

	RegionBoundingBox() { this->initialize(); }

public: // set & get

	int getMaxz() { return m_range.getZMax(); }

public: // functions
	
	void initialize()
	{
		m_minx = m_miny = m_minz = std::numeric_limits<int>::max();
		m_maxx = m_maxy = m_maxz = std::numeric_limits<int>::min();
	}

	void update(const ConstPtr<PointList<int> > &pl)
	{
		Range tmpRange;
		boundingBox(pl, tmpRange);
		include(tmpRange, m_range);
	}


private: // data

	Range m_range;
};
*/


// TODO: use Accumulator
class RegionMoments
{

public: // constructors & destructor

	RegionMoments() { this->initialize(); }

public: // set & get

	void getMoments(SymMatrix3<double> &mat)
	{

		// NB: we do not do m_xx /= m_nElem, etc. in case sb calls 
		// getMoments twice or more

    // TODO: replace this by an operator or something...

		mat(0,0) = (m_sum2(0,0) / m_nElem) - square(m_sum[0] / m_nElem);
		mat(1,1) = (m_sum2(1,1) / m_nElem) - square(m_sum[1] / m_nElem);
		mat(2,2) = (m_sum2(2,2) / m_nElem) - square(m_sum[2] / m_nElem);

		mat(0,1) = (m_sum2(0,1) / m_nElem) - (m_sum[0] / m_nElem) * (m_sum[1] / m_nElem);
		mat(0,2) = (m_sum2(0,2) / m_nElem) - (m_sum[0] / m_nElem) * (m_sum[2] / m_nElem);
		mat(1,2) = (m_sum2(1,2) / m_nElem) - (m_sum[2] / m_nElem) * (m_sum[1] / m_nElem);
	}


public: // functions

	void initialize()
	{
		m_nElem = 0;
	}

	void update(const numeric_array<int,3> &v)
	{
		static SymMatrix3<int> mat;
		m_sum += v;
		tdot(v, mat);
		add(mat, m_sum2);
		++m_nElem;
	}



private: // data

	// sum of coordinates
	numeric_array<double,3> m_sum;

	// sum of multiplication of coordinates
	SymMatrix3<double> m_sum2;

	// number of elements
	int m_nElem;
};



class DoNothing
{

public: // constuctors & destructor

	DoNothing() {};

public:

	void update(const std::vector<numeric_array<int,3> > &) {}
};


class RegionSize
{

public: // constructors & destructor

	RegionSize() { this->init(); };


public: // set & get

	size_t size()
	{
		return m_nElem;
	}

public: // functions

	void init() { m_nElem = 0; };

	void update(const std::vector<numeric_array<int,3> > & pl)
	{
		m_nElem += pl.size();
	}



private: // data

	size_t m_nElem;
};



} // namespace

#endif

