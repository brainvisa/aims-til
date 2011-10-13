#ifndef TIL_READDICOM_H
#define TIL_READDICOM_H

// includes from TIL library
#include "til/til_common.h"
#include "til/ImageC.h"

namespace til {

	TIL_API void readDICOM(const char *fileName, ImageC<unsigned short> &);

} // namespace

#endif

