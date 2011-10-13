#ifndef TIL_STRUCTURETENSOR_H
#define TIL_STRUCTURETENSOR_H


// Standard library includes

#include <vector>


// Local includes

#include "til/til_common.h"


// Namespace 

namespace til {


typedef double t_structtensprec;

template <class TImage1, class TImage2>
class StructureTensor
{
	private:
		
		TImage2 m_kernel;

		// Image derivatives

		TImage1 m_imdx;
		TImage1 m_imdy;
		TImage1 m_imdz;
		

	public:

		// Constructors and destructor

		StructureTensor() {};
		virtual ~StructureTensor() {};

		void setImages(
			TImage1 &imdx, 
			TImage1 &imdy,
			TImage1 &imdz, 
			TImage2 &kernel)
		{
			similarityCheck(imdx, kernel);
			similarityCheck(imdy, kernel);
			similarityCheck(imdz, kernel);

			m_imdx = imdx;
			m_imdy = imdy;
			m_imdz = imdz;

			m_kernel = kernel;
		}

		std::vector<t_structtensprec> compute();
};


template <class TImage1, class TImage2>
std::vector<t_structtensprec> StructureTensor<TImage1,TImage2>::compute()
{

	t_structtensprec dxx, dxy, dxz, dyy, dyz, dzz;

	dxx = dxy = dxz = dyy = dyz = dzz = 0.0;


	typename Iterator<TImage2>::ConstLinear ikernel(m_kernel);
	typename Iterator<TImage1>::ConstLinear idx(m_imdx);
	typename Iterator<TImage1>::ConstLinear idy(m_imdy);
	typename Iterator<TImage1>::ConstLinear idz(m_imdz);

	for (; !ikernel.isAtEnd(); ++ikernel, ++idx, ++idy, ++idz)
	{
		dxx += *ikernel**idx**idx;
		dxy += *ikernel**idx**idy;
		dxz += *ikernel**idx**idz;
		dyy += *ikernel**idy**idy;
		dyz += *ikernel**idy**idz;
		dzz += *ikernel**idz**idz;		
	}

	std::vector<t_structtensprec> res;

	res.push_back(dxx);
	res.push_back(dxy);
	res.push_back(dxz);
	res.push_back(dyy);
	res.push_back(dyz);
	res.push_back(dzz);

	return res;
}


} // namespace

#endif

