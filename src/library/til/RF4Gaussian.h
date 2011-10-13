#ifndef TIL_RF4_GAUSSIAN_H
#define TIL_RF4_GAUSSIAN_H



// Local includes

#include "til/til_common.h"

#include "til/rf4.h"
#include "til/recursiveFilters.h"
#include "til/RecursiveFilterSum.h"

// Namespace 

namespace til {


template < typename T >
void setRecursiveFilterCoefficientsToGaussian
(T &a0, T &a1, 
 T &c0, T &c1, 
 T &w0, T &w1,
 T &b0, T &b1) {
	
	// Deriche's Gaussian approximation:
	//  ( 1.68 * cos(0.6318 * x / sigma) 
	// + 3.735 * sin(0.6318 * x / sigma) ) * exp (-1.783 * x / sigma) -
	// (0.6803 * cos(1.997 * x / sigma) 
	//+ 0.2598 * sin(1.997 * x / sigma) ) * exp (-1.723 * x / sigma)
	
	a0 = (T) 1.68;
	a1 = (T) 3.735;
	c0 = (T) -0.6803;
	c1 = (T) -0.2598;
	w0 = (T) 0.6318;
	w1 = (T) 1.997;
	b0 = (T) 1.783;
	b1 = (T) 1.723;
}

template < typename T >
void setRecursiveFilterCoefficientsToDerivativeGaussian
(T &a0, T &a1, 
 T &c0, T &c1, 
 T &w0, T &w1,
 T &b0, T &b1) {
	
	// Deriche's First Gaussian Derivative approximation:
	// (-0.6472 * cos(0.6719 * x / sigma) 
	// - 4.531  * sin(0.6719 * x / sigma) ) * exp (-1.527 * x / sigma) +
	// (0.6494  * cos(2.072  * x / sigma)
	// + 0.9557 * sin(2.072  * x / sigma) ) * exp (-1.516 * x / sigma)
	
	a0 = (T) -0.6472;
	a1 = (T) -4.531;
	c0 = (T) 0.6494;
	c1 = (T) 0.9557;
	w0 = (T) 0.6719;
	w1 = (T) 2.072;
	b0 = (T) 1.527;
	b1 = (T) 1.516;
}


template < typename T >
void setRecursiveFilterCoefficientsToSecondDerivativeGaussian
(T &a0, T &a1, 
 T &c0, T &c1, 
 T &w0, T &w1,
 T &b0, T &b1) {

	a0 = (T) -1.331;
	a1 = (T) 3.661;
	c0 = (T) 0.3225;
	c1 = (T) -1.738;
	w0 = (T) 0.748;
	w1 = (T) 2.166;
	b0 = (T) 1.24;
	b1 = (T) 1.314;
}



template < typename T >
T recursiveGaussianNorm
(T a0, T a1,
 T c0, T c1,
 T cw0, T sw0,
 T cw1, T sw1,
 T eb0, T eb1)
{
	return 
		eb0 * (a0 * (cw0 - eb0) - a1*sw0)
		/ (T(2) * cw0 * eb0 - eb0*eb0 - T(1)) +
		eb1 * (c0 * (cw1 - eb1) - c1*sw1)
		/ (T(2) * cw1 * eb1 - eb1*eb1 - T(1));
}

template < typename T >
T recursiveGaussianDerivativeNorm
(T a0, T a1,
 T c0, T c1,
 T cw0, T sw0,
 T cw1, T sw1,
 T eb0, T eb1)
{
	return 
		eb0 * 
		( cw0*a0 + eb0*eb0*a0*cw0 + eb0*eb0*a1*sw0 - 2.0*a0*eb0 -a1*sw0 ) / 
		(4*square(cw0*eb0)-4*cw0*eb0-4*cw0*eb0*eb0*eb0+2.0*eb0*eb0+square(eb0*eb0)+1.0) + 
		eb1 * 
		( cw1*c0 + eb1*eb1*c0*cw1 + eb1*eb1*c1*sw1 - 2.0*c0*eb1 -c1*sw1 ) / 
		(4*square(cw1*eb1)-4*cw1*eb1-4*cw1*eb1*eb1*eb1+2.0*eb1*eb1+square(eb1*eb1)+1.0);
}

template < typename T >
T recursiveGaussianSecondDerivativeNorm
(T a0, T a1,
 T c0, T c1,
 T cw0, T sw0,
 T cw1, T sw1,
 T eb0, T eb1)
{
	return
		(-2*cube(cw0)*cube(eb0)*a1
		-2*cube(cw0)*eb0*a1
		-square(cw0)*square(square(eb0))*a1
		+2*square(cw0)*cube(eb0)*a0*sw0
		+6*square(eb0)*square(cw0)*a1
		-2*square(cw0)*eb0*a0*sw0
		-square(cw0)*a1
		+cw0*square(square(eb0))*a0*sw0
		+2*cube(eb0)*a1*cw0
		+2*cw0*a1*eb0
		-cw0*sw0*a0
		+a1*square(square(eb0))
		-4*cube(eb0)*a0*sw0
		-6*a1*square(eb0)
		+4*a0*eb0*sw0+a1)*eb0/(
		(4*square(cw0)*square(eb0)
		-4*cw0*cube(eb0)
		-4*eb0*cw0
		+square(square(eb0))
		+2*square(eb0)+1)*sw0*
		(-2*eb0*cw0+square(eb0)+1))
		+
		(
		-2*cube(cw1)*cube(eb1)*c1
		-2*cube(cw1)*eb1*c1
		-square(square(eb1))*square(cw1)*c1
		+square(square(eb1))*c1
		+square(square(eb1))*sw1*c0*cw1
		-4*cube(eb1)*sw1*c0
		+2*square(cw1)*cube(eb1)*sw1*c0
		+2*cube(eb1)*cw1*c1
		-6*c1*square(eb1)
		+6*square(cw1)*square(eb1)*c1
		-2*square(cw1)*eb1*c0*sw1
		+4*c0*eb1*sw1
		+2*cw1*c1*eb1
		+c1
		-cw1*c0*sw1
		-square(cw1)*c1)*eb1/
		((-4*cube(eb1)*cw1
		+square(square(eb1))
		+2*square(eb1)
		+4*square(cw1)*square(eb1)
		-4*eb1*cw1+1)*sw1*
		(-2*eb1*cw1+square(eb1)+1));
}



template < typename T >
RecursiveFilterSum < RecursiveFilter<T, +1, 4>,
					 RecursiveFilter<T, -1, 4> > RF4Gaussian(T sigma)
{
	typedef RecursiveFilterSum< RecursiveFilter<T, +1, 4>, RecursiveFilter<T, -1, 4> > RecGaussian;

	const T EPSILON = 128 * std::numeric_limits<T>::epsilon();

	RecursiveFilter<T,+1,4> resForward;
	RecursiveFilter<T,-1,4> resBackward;


	// If sigma is too small (division following), set recursive
	// filter to identity
	if (sigma <= EPSILON)
	{
		resForward.setFilter(1,0,0,0,0,0,0,0);
		resBackward.setFilter(1,0,0,0,0,0,0,0);
		return RecGaussian(resForward, resBackward);
	}

	T a0, a1, c0, c1, b0, b1, w0, w1;
	setRecursiveFilterCoefficientsToGaussian(a0,a1,c0,c1,w0,w1,b0,b1);

	resForward.setFilterFromFormula(a0, a1, c0, c1, w0, w1, b0, b1, sigma);
	resBackward.setFilterFromFormula(a0, a1, c0, c1, w0, w1, b0, b1, sigma);

	T cw0 = cos(w0/sigma);
	T cw1 = cos(w1/sigma);
	T sw0 = sin(w0/sigma);
	T sw1 = sin(w1/sigma);
	T norm = recursiveGaussianNorm(a0,a1,c0,c1,cw0,sw0,cw1,sw1,(T)exp(b0/sigma),(T)exp(b1/sigma));

	resForward.normalize(norm);
	resBackward.normalize(norm);

	return RecGaussian(resForward, resBackward);
}


template < typename T >
RecursiveFilterSum< RecursiveFilter<T, +1, 4>,
					RecursiveFilter<T, -1, 4> > RF4GaussianDerivative(T sigma, T stepSize = 1)
{
	typedef RecursiveFilterSum< RecursiveFilter<T, +1, 4>, RecursiveFilter<T, -1, 4> > RecGaussian;

	const T EPSILON = 128 * std::numeric_limits<T>::epsilon();

	if (sigma <= EPSILON)
	{
		// WHY: division by sigma later on
		throw std::domain_error("Sigma is too small");
	}

	T a0, a1, c0, c1, b0, b1, w0, w1;

	setRecursiveFilterCoefficientsToDerivativeGaussian(a0,a1,c0,c1,w0,w1,b0,b1);

	RecursiveFilter<T, +1, 4> resForward;
	RecursiveFilter<T, -1, 4> resBackward;
	resForward.setFilterFromFormula(-a0, -a1, -c0, -c1, w0, w1, b0, b1, sigma);
	resBackward.setFilterFromFormula(a0, a1, c0, c1, w0, w1, b0, b1, sigma);

	T cw0 = cos(w0/sigma);
	T cw1 = cos(w1/sigma);
	T sw0 = sin(w0/sigma);
	T sw1 = sin(w1/sigma);

	T norm = recursiveGaussianDerivativeNorm(a0,a1,c0,c1,cw0,sw0,cw1,sw1,(T)exp(b0/sigma),(T)exp(b1/sigma));

	norm *= stepSize;

	resForward.normalize(norm);
	resBackward.normalize(norm);

	return RecGaussian(resForward, resBackward);
}


template < typename T >
RecursiveFilterSum< RecursiveFilter<T, +1, 4>,
					RecursiveFilter<T, -1, 4> > RF4GaussianSecondDerivative(T sigma)
{
	typedef RecursiveFilterSum< RecursiveFilter<T, +1, 4>, RecursiveFilter<T, -1, 4> > RecGaussian;


	if (sigma < 0.7)
	{
		throw std::invalid_argument("Second order derivative computation using recursive filtering are unreliable for sigma < 0.7");
	}


	T a0, a1, c0, c1, b0, b1, w0, w1;
	setRecursiveFilterCoefficientsToSecondDerivativeGaussian(a0,a1,c0,c1,w0,w1,b0,b1);

	RecursiveFilter<T, +1, 4> resForward;
	RecursiveFilter<T, -1, 4> resBackward;
	resForward.setFilterFromFormula(a0, a1, c0, c1, w0, w1, b0, b1, sigma);
	resBackward.setFilterFromFormula(a0, a1, c0, c1, w0, w1, b0, b1, sigma);

	T cw0 = cos(w0/sigma);
	T cw1 = cos(w1/sigma);
	T sw0 = sin(w0/sigma);
	T sw1 = sin(w1/sigma);

	T norm = recursiveGaussianSecondDerivativeNorm(a0,a1,c0,c1,cw0,sw0,cw1,sw1,(T)exp(b0/sigma),(T)exp(b1/sigma));

	resForward.normalize(norm);
	resBackward.normalize(norm);

	return RecGaussian(resForward, resBackward);
}

} // namespace

#endif

