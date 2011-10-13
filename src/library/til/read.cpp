
// includes from TIL library
#include "til/read.h"
#include "til/CTImage.h"
#include "til/RAWImage.h"
#include "til/imageTools.h"



// Namespace 

namespace til {

#define readCT_DEFINITION(type,function)                                                  \
TIL_API void readCT(const char *fileName, ImageC< type > &im)                             \
{                                                                                         \
  int xs, ys, zs;                                                                         \
  float vx, vy, vz;                                                                       \
  type * data = function (const_cast<char*>(fileName), &xs, &ys, &zs, &vx, &vy, &vz);     \
  im.init(data, numeric_array<int,3>(xs, ys, zs), numeric_array<float,3>(vx, vy, vz));    \
}


TIL_API void readRAW(const char *fileName, ImageC<uchar> &im, int xs, int ys, int zs, float vx, float vy, float vz)
{
	uchar * data = readRawImage8 (const_cast<char*>(fileName), xs, ys, zs);
	im.init(data, numeric_array<int,3>(xs, ys, zs), numeric_array<float,3>(vx, vy, vz));
}


TIL_API void readRAW(const char *fileName, ImageC<ushort> &im, int xs, int ys, float vx, float vy, float vz)
{
	int zs;
	ushort * data = readRawImage16 (const_cast<char*>(fileName), xs, ys, 2, &zs);
	im.init(data, numeric_array<int,3>(xs, ys, zs), numeric_array<float,3>(vx, vy, vz));
}



readCT_DEFINITION(double, readCTImageDouble);
readCT_DEFINITION(float, readCTImageFloat);
readCT_DEFINITION(ushort, readCTImage16);
readCT_DEFINITION(uchar, readCTImage8);


} // namespace



