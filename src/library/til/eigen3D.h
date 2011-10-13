#ifndef TIL_EIGEN3D_H
#define TIL_EIGEN3D_H

// includes from STL
//#include <math.h>
#include <cmath>
#include <limits>
//#include <complex> // Warnings at compil time: can't use in DLL?

// includes from TIL library
#include "til/til_common.h"


// Namespace 
namespace til
{
  template <class T>
  void myComplexPow(T numReal, T numIm, T ex, T &resReal, T &resIm)
  {
	  // Value under which no computation is done
    const T EPSILON = 64 * (std::numeric_limits<T>::epsilon());

	  T numNorm = (T) std::sqrt(double(numReal)*numReal + double(numIm)*numIm);
    
    if (numNorm <= EPSILON)
    {
        resReal = resIm = T(0);
        return;
    }
      
    //T numNorm = (T) std::sqrt(numNorm2);

    T resNorm = (T) std::pow(numNorm, ex);

	  double angle;

	  if (numReal < 0)
	  {
		  if (numIm < 0)
		  {
        angle = (-M_PI - 2.0*std::atan(numIm / (numNorm - numReal)))*ex;
		  }
		  else
		  {
        angle = (M_PI - 2.0*std::atan(numIm / (numNorm - numReal)))*ex;
		  }
	  }
	  else 
	  {
      angle = 2.0*std::atan(numIm / (numNorm + numReal))*ex;
	  }

    resReal = T(std::cos(angle)*resNorm);
    resIm   = T(std::sin(angle)*resNorm);
  }

  template <class T>
  void myComplexSqrt(T num, T &resReal, T &resIm)
  {
	  if (num >= T(0))
	  {
		  resReal = (T) std::sqrt(num);
		  resIm   = (T) 0;
	  }
	  else
	  {
		  resReal = (T) 0;
		  resIm   = (T) std::sqrt(-num);
	  }
  }


  template <class T>
  void eigen3D(T A11, T A12, T A13, T A22, T A23, T A33,
			  T &ev1, T &ev2, T &ev3)
  {

	  T A2 = -(A11+A22+A33);
	  T A1 = -(A12*A12+A13*A13+A23*A23-A11*A22-A11*A33-A22*A33);
	  T A0 = -(A11*A22*A33+2*A12*A13*A23-A12*A12*A33-A11*A23*A23-A13*A13*A22);

	  T R = (9*A2*A1-27*A0-2*cube(A2))/T(54);
	  T D = (-square(A1*A2) + 27*square(A0)-18*A0*A1*A2 + 4*A0*cube(A2) + 4*cube(A1)) / 108;


	  T sqrtDReal;
	  T sqrtDIm;

	  myComplexSqrt(D, sqrtDReal, sqrtDIm);

    if (std::abs(sqrtDIm) < 128 * (std::numeric_limits<T>::epsilon()) && R < sqrtDReal )
	  {
		  sqrtDIm = 128 * (std::numeric_limits<T>::epsilon());
	  }

	  T SReal, SIm;
	  myComplexPow(R+sqrtDReal, sqrtDIm, T(1.0/3.0), SReal, SIm);
  	

	  T TReal, TIm;
	  myComplexPow(R-sqrtDReal, -sqrtDIm, T(1.0/3.0), TReal, TIm);

	  ev1 = T(-(1.0/3.0)*A2 + (SReal+TReal));
	  ev2 = T(-(1.0/3.0)*A2 - 0.5*(SReal+TReal) - 0.5*std::sqrt(3.0)*(SIm-TIm));
	  ev3 = T(-(1.0/3.0)*A2 - 0.5*(SReal+TReal) + 0.5*std::sqrt(3.0)*(SIm-TIm));
  }

} // namespace til

#endif

