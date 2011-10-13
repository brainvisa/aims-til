#ifndef TIL_READ_H
#define TIL_READ_H

// includes from TIL
#include "til/til_common.h"
#include "til/ImageC.h"

// Namespace 

namespace til {

TIL_API void readCT(const char *fileName, ImageC<double> &);
TIL_API void readCT(const char *fileName, ImageC<float > &);
TIL_API void readCT(const char *fileName, ImageC<ushort> &);
TIL_API void readCT(const char *fileName, ImageC<uchar > &);

TIL_API void readRAW(const char *fileName, ImageC<ushort> &im, int xs, int ys, float vx, float vy, float vz);
TIL_API void readRAW(const char *fileName, ImageC<uchar>  &im, int xs, int ys, int zs, float vx = 1, float vy = 1, float vz = 1);

} // namespace

#endif

