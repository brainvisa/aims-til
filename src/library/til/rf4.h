#ifndef TIL_RECURSIVE_FILTER4_H
#define TIL_RECURSIVE_FILTER4_H

// includes from STL
#include <cassert>
#include <vector>

// includes from TIL library
#include "til/til_common.h"
#include "til/RecursiveFilter.h"

// Namespace 
namespace til {



// 4-th order recursive filtering
	
// direction = +1 : forward filter
// direction = -1 : backward filter

template < typename TT, int direction >
class RecursiveFilter<TT, direction, 4>
{
public: // typedefs

	typedef RecursiveFilter<TT, direction, 4> Self;
	typedef TT T;

public: // types

	// Boundary conditions of the recursive filter
	// Zeros : values are supposed to be zero outside
	// Constant : boundary values are duplicated to +/- infinity
	enum BoundaryConditions { Zero, Constant };

public: // constuctors & destructor

	// Default constructor
	RecursiveFilter() { m_boundaryConditions = Self::Constant; }

public: // set & get

	// get filter coefficient for input
	const T* getMI() const { return m_mi; }
	// get filter coefficient for output
	const T* getMO() const { return m_mo; }
	
	// Set filter coefficients
	// The filter is output = sum(mik*inputk) - sum(mok*outputk)

	void setFilter(T mi0, T mi1, T mi2, T mi3, T mo0, T mo1, T mo2, T mo3)
	{ 
		m_mi[0] = mi0; 
		m_mi[1] = mi1;
		m_mi[2] = mi2;
		m_mi[3] = mi3;
		m_mo[0] = mo0; 
		m_mo[1] = mo1;
		m_mo[2] = mo2;
		m_mo[3] = mo3;
	}

	// Set filter coefficients from impulse response formula
	// Filter formula:
	// ( a0 * cos(w0 * x/sigma) + a1 * cos(w0 * x/sigma) * exp( - b0 * x/sigma) ) +
	// ( c0 * cos(w1 * x/sigma) + c1 * cos(w1 * x/sigma) * exp( - b1 * x/sigma) )

	void setFilterFromFormula(T a0, T a1, T c0, T c1, T w0, T w1, T b0, T b1, T sigma);

	// Set boundary conditions

	void setBoundaryConditions(BoundaryConditions bc) { m_boundaryConditions = bc; }


	// TODO: To be exported outside class
	void normalize(T norm)
	{
		norm = 2.0*norm - m_mi[0];

		m_mi[0] /= norm;
		m_mi[1] /= norm;
		m_mi[2] /= norm;
		m_mi[3] /= norm;
	}


public: // methods

	// Apply filter to 'in' (of given length), result outputed in 'out'
	// out has to be allocated with the same size as in
	// NB: out and in cannot point to the same data
	void apply(const std::vector<T> &in, std::vector<T> &out) const;


private: // methods

	// Multiplicative coefficient used for 'Constant' boundary condition

	T borderFactor() const { return 
		(    m_mi[0] + m_mi[1] + m_mi[2] + m_mi[3]) /
		(1 + m_mo[0] + m_mo[1] + m_mo[2] + m_mo[3]); }

	T apply(T in0, T in1, T in2, T in3, T out0, T out1, T out2, T out3) const
	{
		return (
			m_mi[0] * in0 + 
			m_mi[1] * in1 + 
			m_mi[2] * in2 + 
			m_mi[3] * in3 -

			m_mo[0] * out0 - 
			m_mo[1] * out1 - 
			m_mo[2] * out2 - 
			m_mo[3] * out3);
	}

private: // data

	// Boundary conditions
	BoundaryConditions m_boundaryConditions;

	// 4 multiplicative coefficients of the input
	T m_mi[4];

	// 4 multiplicative coefficients of the output
	T m_mo[4];
};





template < typename TT, int direction >
void RecursiveFilter<TT,direction,4>::apply(const std::vector<T> &in, std::vector<T> &out) const
{
	typename std::vector<T>::const_iterator pIn;
	typename std::vector<T>::iterator pOut;

	size_t length = in.size();
	assert(in.size() == out.size());

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

	size_t count = length;
	
	T x0, x1 = 0, x2 = 0, x3 = 0; // useless initializations -- done just to avoid a warning.
	T y0 = 0, y1 = 0, y2 = 0, y3 = 0; // useless initializations -- done just to avoid a warning.

	// Filter initialization on borders

	switch (m_boundaryConditions)
	{
	case Zero:
		x1 = x2 = x3 = 0;
		y0 = y1 = y2 = y3 = 0;
		break;
	case Constant:
		x1 = x2 = x3 = *pIn;
		y0 = y1 = y2 = y3 = *pIn * this->borderFactor();
	}

	for(;;)
	{	
		x0 = *pIn;
		*pOut += y3 = this->apply(x0, x1, x2, x3, y0, y1, y2, y3);
		if (--count == 0) break;
		pOut += direction;
		pIn  += direction;
		
		x3 = *pIn;
		*pOut += y2 = this->apply(x3, x0, x1, x2, y3, y0, y1, y2);
		if (--count == 0) break;
		pOut += direction;
		pIn  += direction;
		
		x2 = *pIn;
		*pOut += y1 = this->apply(x2, x3, x0, x1, y2, y3, y0, y1);
		if (--count == 0) break;
		pOut += direction;
		pIn  += direction;
		
		x1 = *pIn;
		*pOut += y0 = this->apply(x1, x2, x3, x0, y1, y2, y3, y0);		
		if (--count == 0) break;
		pOut += direction;
		pIn  += direction;
	}
}

template < typename TT, int direction >
void RecursiveFilter<TT,direction,4>::setFilterFromFormula(T a0, T a1, T c0, T c1, T w0, T w1, T b0, T b1, T sigma)
{

	// Standard shortcuts variables

	T cw0 = cos(w0/sigma);
	T cw1 = cos(w1/sigma);
	T sw0 = sin(w0/sigma);
	T sw1 = sin(w1/sigma);

	T eb0 = exp(-b0/sigma);
	T eb1 = exp(-b1/sigma);
	T e2b0 = exp(-2.0*b0/sigma);
	T e2b1 = exp(-2.0*b1/sigma);
	T eb0b1 = exp(-(b0+b1)/sigma);
	T e2b0b1 = exp(-(2.0*b0+b1)/sigma);
	T eb02b1 = exp(-(b0+2.0*b1)/sigma);


	// Deduce iterative coefficients from previous filter

	// coefficient for input

	m_mi[0] = a0 + c0;
	m_mi[1] =
		eb1 * (sw1 * c1 - 
		cw1 * (c0 + 2*a0)) +
		eb0 * (sw0 * a1 -
		cw0 * (a0 + 2*c0));
	m_mi[2] =
		2.0 * eb0b1 * ((a0+c0)*cw0*cw1 - a1*sw0*cw1 - c1*sw1*cw0)
		+ c0 * e2b0 + a0 * e2b1;
	m_mi[3] = e2b0b1 * (c1*sw1 - c0*cw1) + eb02b1 * (a1*sw0 - a0*cw0);


	// coefficient for output

	m_mo[0] = -2.0*(eb1*cw1 + eb0*cw0);
	m_mo[1] = 4*cw1*cw0*eb0b1 + e2b0 + e2b1;
	m_mo[2] = -2.0*(cw0*eb02b1 + cw1*e2b0b1);
	m_mo[3] = exp(-2.0*(b0+b1)/sigma);
}





// 4th order recursive filter

	/*
template < typename T >
class RF4
{
public: // constuctors & destructor

	RF4() {}

public: // set & get

	// get forward filter coefficient for input
	const T* getMIF() const { return m_mif; }
	// get forward filter coefficient for output
	const T* getMOF() const { return m_mof; }
	// get backward filter coefficient for input
	const T* getMIB() const { return m_mib; }
	// get backward filter coefficient for input
	const T* getMOB() const { return m_mob; }
	
	// Set filter coefficients

	void setForwardFilter(T mi0, T mi1, T mi2, T mi3,
						  T mo0, T mo1, T mo2, T mo3)
	{ this->setFilter(mi0, mi1, mi2, mi3, mo0, mo1, mo2, mo3, m_mif, m_mof); }

	void setBackwardFilter(T mi0, T mi1, T mi2, T mi3,
						   T mo0, T mo1, T mo2, T mo3)
	{ this->setFilter(mi0, mi1, mi2, mi3, mo0, mo1, mo2, mo3, m_mib, m_mob); }


	// Set filter coefficients from impulse response formula

	// Filter formula:
	// ( a0 * cos(w0 * x/sigma) + a1 * cos(w0 * x/sigma) * exp( - b0 * x/sigma) ) +
	// ( c0 * cos(w1 * x/sigma) + c1 * cos(w1 * x/sigma) * exp( - b1 * x/sigma) )

	void setForwardFormula(T a0, T a1, T c0, T c1, T w0, T w1, T b0, T b1, T sigma)
	{
		this->formula2filter(a0, a1, c0, c1, w0, w1, b0, b1, sigma, m_mif, m_mof);
	}

	void setBackwardFormula(T a0, T a1, T c0, T c1, T w0, T w1, T b0, T b1, T sigma)
	{
		this->formula2filter(a0, a1, c0, c1, w0, w1, b0, b1, sigma, m_mib, m_mob);
	}


public: // methods

	// Apply filter to 'in' (of given length), result outputed in 'out'
	// out has to be allocated with the same size as in
	void apply(const T *in, T *out, int length) const;


	// Normalize filter
	// TODO: This is probably to be exported outside class
	void normalize(T norm)
	{
		this->normalize(norm, m_mif);
		this->normalize(norm, m_mib);
	}

private: // methods

	// Multiplicative coefficient used for RF initialization

	T borderFactorForward() const { return borderFactor(m_mif, m_mof); }
	T borderFactorBackward() const { return borderFactor(m_mib, m_mob); }
	T borderFactor(const T mi[4], const T mo[4]) const { return (mi[0]+mi[1]+mi[2]+mi[3])/(1+mo[0]+mo[1]+mo[2]+mo[3]); }

	T applyForward(T in0, T in1, T in2, T in3,
				   T out0, T out1, T out2, T out3) const
	{
		return (
			m_mif[0] * in0 + 
			m_mif[1] * in1 + 
			m_mif[2] * in2 + 
			m_mif[3] * in3 -

			m_mof[0] * out0 - 
			m_mof[1] * out1 - 
			m_mof[2] * out2 - 
			m_mof[3] * out3);
	}

	T applyBackward(T in0, T in1, T in2, T in3,
				    T out0, T out1, T out2, T out3) const
	{
		return (
			m_mib[0] * in0 + 
			m_mib[1] * in1 + 
			m_mib[2] * in2 + 
			m_mib[3] * in3 -

			m_mob[0] * out0 - 
			m_mob[1] * out1 - 
			m_mob[2] * out2 - 
			m_mob[3] * out3);
	}


	void normalize(T norm, T mi[4])
	{
		norm = 2.0*norm - mi[0];
		
		mi[0] /= norm;
		mi[1] /= norm;
		mi[2] /= norm;
		mi[3] /= norm;
	}

	void setFilter(T mi0, T mi1, T mi2, T mi3, T mo0, T mo1, T mo2, T mo3, T mi[4], T mo[4])
	{
		mi[0] = mi0; 
		mi[1] = mi1;
		mi[2] = mi2;
		mi[3] = mi3;
		mo[0] = mo0; 
		mo[1] = mo1;
		mo[2] = mo2;
		mo[3] = mo3;
	}

	void formula2filter(T a0, T a1, T c0, T c1, T w0, T w1, T b0, T b1, T sigma, T mi[4], T mo[4])
	{
		//std::cout << "FORMULA2FILTER" << std::endl;
		//std::cout << "ABC " << a0 __ a1 __ c0 __ c1 __ w0 __ w1 __ b0 __ b1 << std::endl;

		// Standard shortcuts variables
		
		T cw0 = cos(w0/sigma);
		T cw1 = cos(w1/sigma);
		T sw0 = sin(w0/sigma);
		T sw1 = sin(w1/sigma);
		
		T eb0 = exp(-b0/sigma);
		T eb1 = exp(-b1/sigma);
		T e2b0 = exp(-2.0*b0/sigma);
		T e2b1 = exp(-2.0*b1/sigma);
		T eb0b1 = exp(-(b0+b1)/sigma);
		T e2b0b1 = exp(-(2.0*b0+b1)/sigma);
		T eb02b1 = exp(-(b0+2.0*b1)/sigma);
		

		//std::cout << "CW " << cw0 __ cw1 __ sw0 __ sw1 __ eb0 __ eb1 __ e2b0 __ e2b1 __ eb0b1 __ e2b0b1 __ eb02b1 << std::endl;

		// Deduce iterative coefficients from previous filter
		
		// coefficient for input

		mi[0] = a0 + c0;
		mi[1] =
			eb1 * (sw1 * c1 - 
			cw1 * (c0 + 2*a0)) +
			eb0 * (sw0 * a1 -
			cw0 * (a0 + 2*c0));
		mi[2] =
			2.0 * eb0b1 * ((a0+c0)*cw0*cw1 - a1*sw0*cw1 - c1*sw1*cw0)
			+ c0 * e2b0 + a0 * e2b1;
		mi[3] = e2b0b1 * (c1*sw1 - c0*cw1) + eb02b1 * (a1*sw0 - a0*cw0);
		

		// coefficient for output
		
		mo[0] = -2.0*(eb1*cw1 + eb0*cw0);
		mo[1] = 4*cw1*cw0*eb0b1 + e2b0 + e2b1;
		mo[2] = -2.0*(cw0*eb02b1 + cw1*e2b0b1);
		mo[3] = exp(-2.0*(b0+b1)/sigma);

			
		//std::cout << mi[0] __ mi[1] __ mi[2] __ mi[3] << std::endl;
		//std::cout << mo[0] __ mo[1] __ mo[2] __ mo[3] << std::endl;

	}


	// Pb: the value at zero is computed two times, once forward, once
	// backward. So it has to be removed.
	// Theoretically this value is m_mif[0] = m_mib[0]
	// But this equality may not hold (this is the case for
	// Gaussian 1st derivative for example).
	// So the mean is taken here.

	T bias() const { return (m_mif[0] + m_mib[0]) / 2.0; }



private: // data

	// FORWARDS coefficients

	// 4 multiplicative coefficients of the input
	T m_mif[4];

	// 4 multiplicative coefficients of the output
	T m_mof[4];


	// BACKWARDS coefficients

	// 4 multiplicative coefficients of the input
	T m_mib[4];

	// 4 multiplicative coefficients of the output
	T m_mob[4];


};





template < typename T >
void RF4<T>::apply(const T *in, T *out, int nSeq) const
{

	const T *pIn = in;
	T *pOut = out;

	// Initialization of output
	
	for (int i = 0; i < nSeq; ++i)
	{
		out[i] = - in[i] * this->bias();
	}

	int count = nSeq;
	
	T x0, x1, x2, x3;
	T y0, y1, y2, y3;

	// Filter initialization on borders

	x1 = x2 = x3 = *pIn;
	y0 = y1 = y2 = y3 = *pIn * this->borderFactorForward();

	for(;;)
	{	
		x0 = *pIn;
		(*pOut) += y3 = this->applyForward(x0, x1, x2, x3, y0, y1, y2, y3);
		if (--count == 0) break;
		++pOut;
		++pIn;
		
		x3 = *pIn;
		(*pOut) += y2 = this->applyForward(x3, x0, x1, x2, y3, y0, y1, y2);
		if (--count == 0) break;
		++pOut;
		++pIn;
		
		x2 = *pIn;
		(*pOut) += y1 = this->applyForward(x2, x3, x0, x1, y2, y3, y0, y1);
		if (--count == 0) break;
		++pOut;
		++pIn;
		
		x1 = *pIn;
		(*pOut) += y0 = this->applyForward(x1, x2, x3, x0, y1, y2, y3, y0);		
		if (--count == 0) break;
		++pOut;
		++pIn;
	}
	

	pIn = in + nSeq - 1;
	pOut = out + nSeq - 1;

	x1 = x2 = x3 = *pIn;
	y0 = y1 = y2 = y3 = *pIn * this->borderFactorBackward();

	count = nSeq;

	for(;;)
	{	
		x0 = *pIn;
		(*pOut) += y3 = this->applyBackward(x0, x1, x2, x3, y0, y1, y2, y3);
		if (--count == 0) break;
		--pOut;
		--pIn;
		
		x3 = *pIn;		
		(*pOut) += y2 = this->applyBackward(x3, x0, x1, x2, y3, y0, y1, y2);
		if (--count == 0) break;
		--pOut;
		--pIn;
		
		x2 = *pIn;
		(*pOut) += y1 = this->applyBackward(x2, x3, x0, x1, y2, y3, y0, y1);		
		if (--count == 0) break;
		--pOut;
		--pIn;
		
		x1 = *pIn;
		(*pOut) += y0 = this->applyBackward(x1, x2, x3, x0, y1, y2, y3, y0);
		if (--count == 0) break;
		--pOut;
		--pIn;
	}
}
*/




} // namespace


#endif

