#ifndef TIL_RECURSIVE_FILTER1_H
#define TIL_RECURSIVE_FILTER1_H

// includes from STL
#include <cassert>
#include <limits> // numeric_limits
#include <vector>

// include from TIL library
#include "til/til_common.h"
#include "til/RecursiveFilter.h"


namespace til {


// 1st order recursive filter
	
// direction = +1 : forward filter
// direction = -1 : backward filter

template < typename T, int direction>
class RecursiveFilter<T, direction, 1>
{
public: // typedefs

	typedef RecursiveFilter<T, direction, 1> Self;

public: // types

	// Boundary conditions of the recursive filter
	// Zeros : values are supposed to be zero outside
	// Constant : boundary values are duplicated to +/- infinity
	// Mirror : image is assumed to have R/Z topology
	enum BoundaryConditions { Zero, Constant, Mirror };

public: // constuctors & destructor

	// Default constructor
	RecursiveFilter() { m_boundaryConditions = Self::Constant; }

public: // set & get

	// get filter coefficient for input
	T getMI() const { return m_mi0; }
	// get filter coefficient for output
	T getMO() const { return m_mo0; }
	
	// Set filter coefficients
	// The filter is output = mi0*input0 - mo0*output0
	void setFilter(T mi0, T mo0)
	{ 
		m_mi0 = mi0;
		m_mo0 = mo0;
	}

	// Set boundary conditions
	void setBoundaryConditions(BoundaryConditions bc) { m_boundaryConditions = bc; }


public: // methods

	// Apply filter to 'in' (of given length), result outputed in 'out'
	// out has to be allocated with the same size as in
	// NB: out and in cannot point to the same data
	void apply(const std::vector<T> &in, std::vector<T> &out) const;


private: // methods

	// Multiplicative coefficient used for 'Constant' boundary condition

	T borderFactor() const { return m_mi0 /(1 + m_mo0); }

	T mirrorBC(const std::vector<T> &in, size_t length) const;

	T apply(T in0, T out0) const
	{
		return m_mi0 * in0 - m_mo0 * out0;
	}

private: // data

	// Boundary conditions
	BoundaryConditions m_boundaryConditions;

	// Filter coefficients;
	T m_mi0;
	T m_mo0;
};


template < typename T, int direction >
void RecursiveFilter<T,direction,1>::apply(const std::vector<T> &in, std::vector<T> &out) const
{
	assert(in.size() == out.size());

	size_t length = in.size();

	typename std::vector<T>::const_iterator pIn;
	typename std::vector<T>::iterator pOut;

	if (direction == 1)
	{
		pIn = in.begin();
		pOut = out.begin();
	}
	else if (direction == -1)
	{
		pIn = in.begin() + length - 1;
		pOut = out.begin() + length - 1;
	}
	
	T x0;
	T y0 = 0; // useless initialization -- done just to avoid a warning.

	// Filter initialization on borders

	switch (m_boundaryConditions)
	{
	case Zero:
		y0 = 0;
		break;
	case Constant:
		y0 = *pIn * this->borderFactor();
	case Mirror:
		y0 = this->mirrorBC(in, length);
	}

	for(size_t count = 0; count < length; ++count)
	{
		x0 = *pIn;
		pIn  += direction;
		*pOut += y0 = this->apply(x0, y0);
		pOut += direction;
	}
}

template < typename T, int direction >
T RecursiveFilter<T, direction,1>::mirrorBC(const std::vector<T> &in, size_t length) const
{
	const T EPSILON = std::numeric_limits<T>::epsilon() * 128;
	typename std::vector<T>::const_iterator pIn;
	size_t count;

	if (direction == 1)
	{
		pIn = in.begin();
	}
	else if (direction == -1)
	{
		pIn == in.begin() + length - 1;
	}
	
	T res = T(0);
	T factor = m_mi0;


	// One way...

	for(count = 0; count < length; ++count)
	{
		res += (*pIn)*factor;
		factor *= m_mo0;
		if (factor <= EPSILON) return res;
		pIn  += direction;
	}

	// ...and back

	pIn -= direction;
	for(count = 1; count < length; ++count)
	{
		res += (*pIn)*factor;
		factor *= m_mo0;
		if (factor <= EPSILON) return res;
		pIn  -= direction;
	}

	return res / (1 - factor * m_mo0);
}

} // namespace

#endif

