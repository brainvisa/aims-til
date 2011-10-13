#ifndef TIL_CTIMAGE_H
#define TIL_CTIMAGE_H

/*
#ifdef __cplusplus
extern "C"
{
#endif
*/

#include <stdio.h>

#include "til/til_common.h"


#define CT_TYPE_UCHAR  1
#define CT_TYPE_INT8   1
#define CT_TYPE_UINT16 2
#define CT_TYPE_USHORT 2
#define CT_TYPE_UINT32 3
#define CT_TYPE_FLOAT  4
#define CT_TYPE_DOUBLE 5

/* Common functions */
TIL_API FILE *readCTImageHeader( char *filename, int *xs, int *ys, int *zs, 
						int *ctType, float *voxelSizeX, float *voxelSizeY, float *voxelSizeZ );
TIL_API void closeCTImage( FILE *fp, float voxelSizeX, float voxelSizeY, float voxelSizeZ );

/* Int 8 */
TIL_API unsigned char *readCTImage8( char *filename, int *xs, int *ys, int *zs,
                        float *voxelSizeX, float *voxelSizeY, float *voxelSizeZ );
TIL_API unsigned char *readCTImage8ByNumber( FILE *fp, int slice, int xs, int ys );
TIL_API unsigned char *readCTImage8Sequentially( FILE *fp, int xs, int ys );

TIL_API int writeCTImage8( char *filename, 
                  unsigned char *sPtr, int xs, int ys, int zs,
                  float voxelSizeX, float voxelSizeY, float voxelSizeZ );
TIL_API FILE *writeCTImage8Header( char *filename, 
						unsigned char *sPtr, int xs, int ys, int zs );
TIL_API int writeCTImage8Sequentially( FILE *fp, 
							  unsigned char *cPtr, int xs, int ys );

/* Int16 */
TIL_API unsigned short *readCTImage16( char *filename, int *xs, int *ys, int *zs,
                        float *voxelSizeX, float *voxelSizeY, float *voxelSizeZ );
TIL_API unsigned short *readCTImage16ByNumber( FILE *fp, int slice, int xs, int ys );
TIL_API unsigned short *readCTImage16Sequentially( FILE *fp, int xs, int ys );

TIL_API int writeCTImage16( char *filename, 
                   unsigned short *sPtr, int xs, int ys, int zs,
                   float voxelSizeX, float voxelSizeY, float voxelSizeZ );
TIL_API FILE *writeCTImage16Header( char *filename, unsigned short *sPtr, int xs, int ys, int zs );
TIL_API int writeCTImage16Sequentially( FILE *fp, 
							   unsigned short *sPtr, int xs, int ys );

/* Int32 */
TIL_API unsigned int *readCTImage32( char *filename, int *xs, int *ys, int *zs,
                        float *voxelSizeX, float *voxelSizeY, float *voxelSizeZ );
TIL_API unsigned int *readCTImage32ByNumber( FILE *fp, int slice, int xs, int ys );
TIL_API unsigned int *readCTImage32Sequentially( FILE *fp, int xs, int ys );

TIL_API int writeCTImage32( char *filename, 
                   unsigned int *iPtr, int xs, int ys, int zs,
                   float voxelSizeX, float voxelSizeY, float voxelSizeZ );
TIL_API FILE *writeCTImage32Header( char *filename, 
						unsigned int *iPtr, int xs, int ys, int zs );
TIL_API int writeCTImage32Sequentially( FILE *fp, 
							  unsigned int *iPtr, int xs, int ys );

/* Float */
TIL_API float *readCTImageFloat( char *filename, int *xs, int *ys, int *zs,
                        float *voxelSizeX, float *voxelSizeY, float *voxelSizeZ );
TIL_API float *readCTImageFloatByNumber( FILE *fp, int slice, int xs, int ys );
TIL_API float *readCTImageFloatSequentially( FILE *fp, int xs, int ys );

TIL_API int writeCTImageFloat( char *filename, 
                      float *fPtr, int xs, int ys, int zs,
                      float voxelSizeX, float voxelSizeY, float voxelSizeZ );
TIL_API FILE *writeCTImageFloatHeader( char *filename, 
						float *fPtr, int xs, int ys, int zs );
TIL_API int writeCTImageFloatSequentially( FILE *fp, 
							  float *fPtr, int xs, int ys );

/* Double */
TIL_API double *readCTImageDouble( char *filename, int *xs, int *ys, int *zs,
                        float *voxelSizeX, float *voxelSizeY, float *voxelSizeZ );
TIL_API double *readCTImageDoubleByNumber( FILE *fp, int slice, int xs, int ys );
TIL_API double *readCTImageDoubleSequentially( FILE *fp, int xs, int ys );

TIL_API int writeCTImageDouble( char *filename, 
                       double *dPtr, int xs, int ys, int zs,
                       float voxelSizeX, float voxelSizeY, float voxelSizeZ );
TIL_API FILE *writeCTImageDoubleHeader( char *filename, 
						double *dPtr, int xs, int ys, int zs );
TIL_API int writeCTImageDoubleSequentially( FILE *fp, 
							  double *dPtr, int xs, int ys );							  
							  
							  /*
#ifdef __cplusplus
}
#endif*/

#endif

