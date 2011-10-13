#ifndef TIL_READTXTFILE_H
#define TIL_READTXTFILE_H

#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

// includes from TIL
#include "til/numeric_array.h"


template < class T >
void readTXT3DPoints(const char *filename, std::vector<numeric_array<T,3> > &v)
{

	ifstream f (filename); // ios::nocreate not defined???
	
	if (!f)
	{
		throw std::runtime_error("Unable to open file");
	}
	
	
	// Read points from file
	
	int count = 0;
	numeric_array<T,3> tv;
	double x, y, z;

	for (;;)
	{
		const int SIZE_BUFF = 100;
		char buff[SIZE_BUFF];

		f >> x >> y >> z;
		if (!f.good()) break;

		// Convert to vector's type

		tv[0] = castValue<double, T>(x);
		tv[1] = castValue<double, T>(y);
		tv[2] = castValue<double, T>(z);

		v.push_back(tv);
		++count;

		// Skip the remaining of the line
		
		f.getline(buff, SIZE_BUFF);
		if (!f.good()) break;
	}

}

#endif