#ifndef TIL_NORMALIZEDCORRELATION_H
#define TIL_NORMALIZEDCORRELATION_H


// Standard library includes

#include <math.h>
#include <limits>


// Local includes

#include "til/til_common.h"

#include "til/imageBasicStats.h"



// Namespace 

namespace til {



typedef double t_corprec;

// returns the normalized correlation between image 1 and image 2
// im1 and im2 should have the same size


template <class TImage1, class TImage2>
class NormalizedCorrelation
{
	private:

		TImage1 m_im1;
		TImage2 m_im2;

		t_corprec m_meanIm1;
		t_corprec m_meanIm2;

		t_corprec m_stdIm1;
		t_corprec m_stdIm2;


	public:

		// constructors and destructors

		NormalizedCorrelation()
		{
			m_meanIm1 = m_meanIm2 = t_corprec(0);
			m_stdIm1 = m_stdIm2 = t_corprec(0);
		};

		virtual ~NormalizedCorrelation() {};


		// Set up parameters

		void setFirstImage(TImage1 &im1)
		{
			allocationCheck(im1);

			if (m_im1 != im1)
			{
				m_im1 = im1;
				m_meanIm1 = (t_corprec) mean(m_im1);
				m_stdIm1 = (t_corprec) ::sqrt(var(m_im1));
			}
		}

		void setSecondImage(TImage2 &im2)
		{
			allocationCheck(im2);

			if (m_im2 != im2)
			{
				m_im2 = im2;
				m_meanIm2 = (t_corprec) mean(m_im2);
				m_stdIm2 = (t_corprec) ::sqrt(var(m_im2));
			}
		}

		TImage1 & getFirstImage()
		{
			return m_im1;
		}

		TImage2 & getSecondImage()
		{
			return m_im2;
		}

		void setImages(TImage1 &im1, TImage2 &im2)
		{
			allocationCheck(im1);
			allocationCheck(im2);
			similarityCheck(im1,im2);

			this->setFirstImage(im1);
			this->setSecondImage(im2);
		}

		t_corprec compute();
};































////////////////////////////////////////////////////
//
// CODE
//
////////////////////////////////////////////////////


template <class TImage1, class TImage2>
t_corprec NormalizedCorrelation<TImage1, TImage2>::compute()
{

	const t_corprec EPSILON = 4*std::numeric_limits<t_corprec>::epsilon() ;


	// If one of the standard deviation is too small
	// we will not be able to normalize and have to make
	// some decisions

	if (m_stdIm1 <= EPSILON || m_stdIm2 <= EPSILON)
	{
		if (m_stdIm1 <= EPSILON && m_stdIm2 <= EPSILON)
		{
			// If both images are very uniform,
			// we throw an exception
			// (chosing a numerical value to return is quite
			// hard in that case)

			throw std::runtime_error("Division by zero");
		}

		// otherwise, we return 0
		// This is not mathematically correct but works
		// for most applications.
		// TODO: maybe actually throwing a special exception
		// is better, and let the user decide if it wants 
		// a zero or not.

		return t_corprec(0.0);
	}


	allocationCheck(m_im1);
	allocationCheck(m_im2);
	similarityCheck(m_im1, m_im2);

	t_corprec res = (t_corprec) 0;

	typename Iterator<TImage1>::ConstLinear i1(m_im1);
	typename Iterator<TImage2>::ConstLinear i2(m_im2);


	for (; !i1.isAtEnd(); ++i1, ++i2)
	{
		res += t_corprec(*i1 - m_meanIm1) * (*i2 - m_meanIm2);
	}



	return res / (m_im1.size() * m_stdIm1 * m_stdIm2);
	
}


} // namespace

#endif

