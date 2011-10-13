
// includes from TIL library
#include "til/write.h"
#include "til/CTImage.h"



// Namespace 

namespace til {

TIL_API void writeCT(const char *fileName, const ImageC<double> &im)
{
	writeCTImageDouble(const_cast<char*>(fileName), const_cast<double*>(im.getPointer()),
		im.dim()[0], im.dim()[1], im.dim()[2],
		(float) im.vdim()[0], (float) im.vdim()[1], (float) im.vdim()[2]);

}

TIL_API void writeCT(const char *fileName, const ImageC<float> &im)
{
	writeCTImageFloat(const_cast<char*>(fileName), const_cast<float*>(im.getPointer()),
		im.dim()[0], im.dim()[1], im.dim()[2],
		(float) im.vdim()[0], (float) im.vdim()[1], (float) im.vdim()[2]);

}

TIL_API void writeCT(const char *fileName, const ImageC<ushort> &im)
{
	writeCTImage16(const_cast<char*>(fileName), const_cast<ushort*>(im.getPointer()),
		im.dim()[0], im.dim()[1], im.dim()[2],
		(float)im.vdim()[0], (float)im.vdim()[1], (float)im.vdim()[2]);

}


TIL_API void writeCT(const char *fileName, const ImageC<uchar> &im)
{
	writeCTImage8(const_cast<char*>(fileName), const_cast<uchar*>(im.getPointer()),
		im.dim()[0], im.dim()[1], im.dim()[2],
		(float)im.vdim()[0], (float)im.vdim()[1], (float)im.vdim()[2]);
}


} // namespace

