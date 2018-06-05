#ifndef TIL_RF4_TOOLS_H
#define TIL_RF4_TOOLS_H


// includes from STL
#include <iostream>
#include <vector>

// includes from TIL library
#include "til/til_common.h"
#include "til/imageTools.h"
#include "til/rf4.h"



// Namespace 

namespace til {


// NB: out can be the same image as in

template < typename TImage, typename TLineFilter >
void filterAlongAxis
(const TImage &im, TImage &out, ImageAxis axis, const TLineFilter &filter)
{
	typedef typename TImage::value_type value_type;

	int barLength = im.dim()[axis];

	std::vector<value_type> barIn(barLength);
	std::vector<value_type> barOut(barLength);

  numeric_array <int,3> p;
	for (p[0]=0; p[0]<(axis==X_AXIS?1:im.dim()[0]); ++p[0])
	{
		for (p[1]=0; p[1]<(axis==Y_AXIS?1:im.dim()[1]); ++p[1])
		{
			for (p[2]=0; p[2]<(axis==Z_AXIS?1:im.dim()[2]); ++p[2])
			{
				extractBar(im, p, axis, barIn);
				filter.apply(barIn, barOut);
				insertBar(out, p, axis, barOut);
			}
		}
	}
}



/*
template < typename T >
void print(RF4<T> rf4)
{
	std::cout << "Forward Filter:" << std::endl;
	std::cout << rf4.getMIF()[0] <<" "<< rf4.getMIF()[1] <<" "<< rf4.getMIF()[2] <<" "<< rf4.getMIF()[3] << std::endl;
	std::cout << rf4.getMOF()[0] <<" "<< rf4.getMOF()[1] <<" "<< rf4.getMOF()[2] <<" "<< rf4.getMOF()[3] << std::endl;

	std::cout << "Backward Filter:" << std::endl;
	std::cout << rf4.getMIB()[0] <<" "<< rf4.getMIB()[1] <<" "<< rf4.getMIB()[2] <<" "<< rf4.getMIB()[3] << std::endl;
	std::cout << rf4.getMOB()[0] <<" "<< rf4.getMOB()[1] <<" "<< rf4.getMOB()[2] <<" "<< rf4.getMOB()[3] << std::endl;

}
*/

} // namespace

#endif

