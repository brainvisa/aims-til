#include <iostream>

#if 0
// includes from DICOMM library
#include "dicomm.h"

// includes from TIL library
#include "til/readDICOM.h"

namespace til {

	TIL_API void readDICOM(const char *fileName, ImageC<unsigned short> &im)
	{
		unsigned short *buff = 0;
		int xs, ys, zs;
		float vx, vy, vz;
		float *trans = 0;
		DicommRead(const_cast<char*>(fileName), &buff, &xs, &ys, &zs, &vx, &vy, &vz, 0);
		std::cout << "(ok" << buff << ")..." << std::flush;
		im.init(buff, xs, ys, zs, vx, vy, vz);
	}

} // namespace

#endif
