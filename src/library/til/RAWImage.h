#ifndef TIL_RAWIMAGE_H
#define TIL_RAWIMAGE_H

#include "til/til_common.h"

#include <stdio.h>

FILE *openRawImage( char *filename, int xs, int ys, int bytesPerPixel, int *zs );
unsigned short *readRawImage16( char *filename, int xs, int ys, 
							   int bytesPerPixel, int *zs );
unsigned short *readRawImage16ByNumber( FILE *fp, int slice, int xs, int ys,
									   int bytesPerPixel );
unsigned short *readRawImage16Sequentially( FILE *fp, int xs, int ys,
										   int bytesPerPixel );
short *readRawImageS16Sequentially( FILE *fp, int xs, int ys,
										   int bytesPerPixel );
unsigned char *readRawImage8( char *filename, int xs, int ys, int zs );
unsigned char*readRawImage8ByNumber( FILE *fp, int slice, int xs, int ys );
unsigned char *readRawImage8Sequentially( FILE *fp, int xs, int ys );
int writeRawImage16( char *filename, 
					unsigned short *sPtr, int xs, int ys, int zs );
int writeRawImage16Sequentially( FILE *fp, 
								unsigned short *sPtr, int xs, int ys );
int writeRawImageS16SequentiallyBO( FILE *fp, 
								short *sPtr, int xs, int ys, int bpp );
int writeRawImage16BO( char *filename, 
					unsigned short *sPtr, int xs, int ys, int zs, int bpp );
int writeRawImage16SequentiallyBO( FILE *fp, 
								unsigned short *sPtr, int xs, int ys, int bpp );
int writeRawImage8( char *filename, 
				   unsigned char *cPtr, int xs, int ys, int zs );
int writeRawImage8Sequentially( FILE *fp, 
							   unsigned char *cPtr, int xs, int ys );

#endif

