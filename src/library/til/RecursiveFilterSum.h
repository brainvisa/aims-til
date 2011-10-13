#ifndef TIL_RECURSIVE_FILTER_SUM_H
#define TIL_RECURSIVE_FILTER_SUM_H

// TODO: rename into bidirectionalFilter

// includes from STL
#include <cassert>
#include <vector>

// includes from BOOST
#include <boost/type_traits.hpp>

// includes from TIL library
#include "til/til_common.h"
#include "til/templateTools.h"


namespace til {

	/// Sum of two recursive filters, a forward and a backward one.
	// TODO: Could be templated w.r.t filter...
template < class TRecursiveFilter1, class TRecursiveFilter2 >
class RecursiveFilterSum
{
public: // typedef

	/// Type of the filter
	// NB: see also private rules
	typedef typename TRecursiveFilter1::T T;

public: // constuctors & destructor

	RecursiveFilterSum(
		const TRecursiveFilter1 & filter1,
		const TRecursiveFilter2 & filter2)
	{
		m_removeOverlap = true;
		m_filter1 = filter1;
		m_filter2 = filter2;
	}

public: // methods
	
	// Apply filter to 'in' (of given length), result outputed in 'out'
	// out has to be allocated with the same size as in
	void apply(const std::vector<T> &in, std::vector<T> &out) const
	{
		// Make sure that input and output vectors have same size
		assert(in.size() == out.size());

		if (m_removeOverlap)
		{
			// TODO: C'mon, I'm sure there is an STL algorithm to do that...
			typename std::vector<T>::const_iterator iIn = in.begin();
			typename std::vector<T>::iterator iOut = out.begin();
			for (; iIn != in.end(); ++iIn, ++iOut)
			{
				*iOut = - *iIn * this->bias();
			}
		}

		m_filter1.apply(in, out);
		m_filter2.apply(in, out);
	}
private: // rules

	// Check that types of both filters are equal
	typedef typename enable_if<boost::is_same<typename TRecursiveFilter1::T, typename TRecursiveFilter2::T> >::type Rule_both_filters_have_same_type;
	// TODO: also check that there are a forward and a backward filter

private: // methods

	// Pb: the value at zero is computed two times, once forward, once
	// backward. So it has to be removed.
	// Theoretically this value is m_mif[0] = m_mib[0]
	// But this equality may not hold (this is the case for
	// Gaussian 1st derivative for example).
	// So the mean is taken here.

	T bias() const { return (m_filter1.getMI()[0] + m_filter2.getMI()[0]) / 2.0; }

	void removeOverlap(bool flag) { m_removeOverlap = flag; }


private: // data

	bool m_removeOverlap;

	TRecursiveFilter1 m_filter1;
	TRecursiveFilter2 m_filter2;
};


}

#endif

