#ifndef TIL_GAUSSIANNORMALIZEDCORRELATION_H
#define TIL_GAUSSIANNORMALIZEDCORRELATION_H


// includes from TIL library
#include "til/til_common.h"
#include "til/imageTools.h"
#include "til/NormalizedCorrelation.h"

// Namespace 

namespace til {


template <class TImage1, class TImage2 >
class GaussianNormalizedCorrelation : public NormalizedCorrelation<TImage1, TImage2>
{
public:
	
	// Constructors and destructor
	
	GaussianNormalizedCorrelation() : NormalizedCorrelation<TImage1,TImage2>()
	{ 
		m_sigma[0] = m_sigma[1] = m_sigma[2] = (typename TImage2::value_type)(-1);
	}

	virtual ~GaussianNormalizedCorrelation() {};
	
	
	void setKernelParameters(typename TImage2::value_type sigmaX,
		typename TImage2::value_type sigmaY,
		typename TImage2::value_type sigmaZ,
		t_voxsize vx, t_voxsize vy, t_voxsize vz,
		int sx=-1, int sy=-1, int sz=-1)
	{
		if (sigmaX < 0 || sigmaY < 0 || sigmaZ < 0)
		{
			//throw Exception("Sigma cannot be < 0");
		}
		
		// If same parameters, don't do anything
		// TODO: change to a float equality test?
		if (m_sigma[0] == sigmaX &&
			m_sigma[1] == sigmaY &&
			m_sigma[2] == sigmaZ &&
			
			vx == this->getSecondImage().vdim()[0] &&
			vy == this->getSecondImage().vdim()[1] &&
			vz == this->getSecondImage().vdim()[2] )
		{
			return;
		}
		
		TImage2 tmp;
		generateGaussianKernel(tmp,
			sigmaX, sigmaY, sigmaZ,
			vx, vy, vz,
			sx, sy, sz);
		
		this->setSecondImage(tmp);
	}
	
	
private:
	
	typename TImage2::value_type m_sigma[3];
	
};

} // namespace

#endif

