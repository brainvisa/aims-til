#ifndef TIL_WRITE_H
#define TIL_WRITE_H

// includes from TIL
#include "til/ImageC.h"


// Namespace 

namespace til {


TIL_API void writeCT(const char *fileName, const ImageC<double> &im);
TIL_API void writeCT(const char *fileName, const ImageC<float> &im);
TIL_API void writeCT(const char *fileName, const ImageC<ushort> &im);
TIL_API void writeCT(const char *fileName, const ImageC<uchar> &im);

} // namespace

#endif

