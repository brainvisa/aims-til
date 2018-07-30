
// inludes from STL
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // _swab (gcc)
#include <unistd.h>

// includes from TIL library
#include "til/RAWImage.h"

FILE *openRawImage( char *filename, int xs, int ys, int bytesPerPixel, int *zs )
{
    FILE *fp;
	long fileLen;

    fp = fopen( filename, "rb" );
    if ( !fp )
    {
        fprintf( stderr, "File %s not found\n", filename );
        exit( 0 );
    }

    /* Determine the length of the image */
	fseek( fp, 0L, SEEK_END );
    fileLen = ftell( fp );

    fseek( fp, 0L, SEEK_SET );

	*zs = fileLen / abs(bytesPerPixel) / xs / ys;

    return( fp );
}

unsigned short *readRawImage16( char *filename, int xs, int ys, 
							   int bytesPerPixel, int *zs )
{
    FILE *fp;
    unsigned short *imgPtr;
	unsigned char *iPtr, *s;
    int numgot;

	fp = openRawImage( filename, xs, ys, bytesPerPixel, zs );

    imgPtr = (unsigned short *)malloc( xs *ys **zs * sizeof(short) + 1 );
    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
	if ( bytesPerPixel == -2 )
	{
		/* Wrong byte order: Swap low and high bytes */
		numgot = fread( (char *)imgPtr + 1, sizeof(short), xs *ys **zs, fp );

		iPtr = (unsigned char *)imgPtr;
        s = (unsigned char *)imgPtr + 2 * xs * ys **zs;
        while ( iPtr < s )
		{
            *iPtr = *( iPtr + 2 );
            iPtr += 2;
		}

	}
	else
	{
		numgot = fread( (char *)imgPtr, sizeof(short), xs * ys **zs, fp );
	}

    fclose( fp );

    if ( numgot != ( xs *ys ) )
    {
        fprintf( stderr, "Warning: unexpected number of elements read\n" );
    }

    return( imgPtr );
}

unsigned short *readRawImage16ByNumber( FILE *fp, int slice, int xs, int ys,
									   int bytesPerPixel )
{
    unsigned short *imgPtr;
    unsigned char *iPtr, *s;
    int numgot;

    if ( slice < 0 )
    {
        fprintf( stderr, "Incorrect slice number\n" );
        exit( 0 );
    }

    fseek( fp, 0L, SEEK_SET );
    fseek( fp, (long)slice * xs * ys * sizeof(short), SEEK_SET );

    imgPtr = (unsigned short *)malloc( xs * ys * sizeof(short) + 1 );

    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
	if ( bytesPerPixel == -2 )
	{
		/* Wrong byte order: Swap low and high bytes */
		numgot = fread( (char *)imgPtr + 1, sizeof(short), xs * ys, fp );

		iPtr = (unsigned char *)imgPtr;
        s = (unsigned char *)imgPtr + 2 * xs * ys;
        while ( iPtr < s )
		{
            *iPtr = *( iPtr + 2 );
            iPtr += 2;
		}

	}
	else
	{
		numgot = fread( (char *)imgPtr, sizeof(short), xs *ys, fp );
	}

    if ( numgot != ( xs *ys ) )
    {
        fprintf( stderr, "Warning: unexpected number of elements read\n" );
    }

    return( imgPtr );
}

unsigned short *readRawImage16Sequentially( FILE *fp, int xs, int ys,
										   int bytesPerPixel )
{
    unsigned short *imgPtr;
    unsigned char *iPtr, *s;
    int numgot;

    imgPtr = (unsigned short *)malloc( xs * ys * sizeof(short) + 1 );

    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
	if ( bytesPerPixel == -2 )
	{
		/* Wrong byte order: Swap low and high bytes */
		numgot = fread( (char *)imgPtr + 1, sizeof(short), xs * ys, fp );

		iPtr = (unsigned char *)imgPtr;
        s = (unsigned char *)imgPtr + 2 * xs * ys;
        while ( iPtr < s )
		{
            *iPtr = *( iPtr + 2 );
            iPtr += 2;
		}
	}
	else
	{
		numgot = fread( (char *)imgPtr, sizeof(short), xs * ys, fp );
	}

	if (numgot != xs * ys)
		fprintf( stderr, "Warning: Red %d bytes out of %d\n", numgot, xs * ys );


    return( imgPtr );
}

short *readRawImageS16Sequentially( FILE *fp, int xs, int ys,
										   int bytesPerPixel )
{
    short *imgPtr;
    unsigned char *iPtr, *s;
    int numgot;

    imgPtr = (short *)malloc( xs * ys * sizeof(short) + 1 );

    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
	if ( bytesPerPixel == -2 )
	{
		/* Wrong byte order: Swap low and high bytes */
		numgot = fread( (char *)imgPtr + 1, sizeof(short), xs * ys, fp );

		iPtr = (unsigned char *)imgPtr;
        s = (unsigned char *)imgPtr + 2 * xs * ys;
        while ( iPtr < s )
		{
            *iPtr = *( iPtr + 2 );
            iPtr += 2;
		}
	}
	else
	{
		numgot = fread( (char *)imgPtr, sizeof(short), xs * ys, fp );
	}

	if (numgot != xs * ys)
		fprintf( stderr, "Warning: Red %d bytes out of %d\n", numgot, xs * ys );		

    return( imgPtr );
}

unsigned char *readRawImage8( char *filename, int xs, int ys, int zs )
{
    FILE *fp;
    unsigned char *imgPtr;
    int numgot;

    fp = fopen( filename, "rb" );
    if ( !fp )
    {
        fprintf( stderr, "File %s not found\n", filename );
        exit( 0 );
    }

    imgPtr = (unsigned char *)malloc( xs *ys *zs * sizeof(unsigned char) );
    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
    numgot = fread( (char *)imgPtr, sizeof(char), xs *ys *zs, fp );
    fclose( fp );

    if ( numgot != ( xs *ys ) )
    {
        fprintf( stderr, "Warning: unexpected number of elements read\n" );
    }

    return( imgPtr );
}

unsigned char *readRawImage8ByNumber( FILE *fp, int slice, int xs, int ys )
{
    unsigned char *imgPtr;
    int numgot;

    if ( slice < 0 )
    {
        fprintf( stderr, "Incorrect slice number\n" );
        exit( 0 );
    }

    fseek( fp, 0L, SEEK_SET );
    fseek( fp, (long)slice * xs * ys * sizeof(char), SEEK_SET );

    imgPtr = (unsigned char *)malloc( xs * ys * sizeof(unsigned char));

    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
    numgot = fread( (char *)imgPtr, sizeof(char), xs * ys, fp );

    if ( numgot != ( xs *ys ) )
    {
        fprintf( stderr, "Warning: unexpected number of elements read\n" );
    }

    return( imgPtr );
}

unsigned char *readRawImage8Sequentially( FILE *fp, int xs, int ys )
{
    unsigned char *imgPtr;
    int numgot;

    imgPtr = (unsigned char *)malloc( xs * ys * sizeof(unsigned char));

    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
    /* Wrong byte order: Swap low and high bytes */
    numgot = fread( (char *)imgPtr, sizeof(char), xs * ys, fp );

    if ( numgot != ( xs *ys ) )
    {
        fprintf( stderr, "Warning: unexpected number of elements read\n" );
    }

    return( imgPtr );
}

int writeRawImage16( char *filename, unsigned short *sPtr, int xs, int ys, int zs )
{
    FILE *fp;

    fp = fopen( filename, "wb" );
    if ( !fp )
    {
        fprintf( stderr, "Cannot open %s\n", filename );
        exit( 0 );
    }

    fwrite( sPtr, xs * ys * zs, sizeof(short), fp );

    fclose( fp );

    return( 1 );
}

int writeRawImage16Sequentially( FILE *fp, unsigned short *sPtr, int xs, int ys )
{
    fwrite( sPtr, xs * ys, sizeof(short), fp );

    return( 1 );
}

int writeRawImage16BO( char *filename, unsigned short *sPtr, 
					  int xs, int ys, int zs, int bpp )
{
    FILE *fp;
	unsigned short *tmpPtr;

    fp = fopen( filename, "wb" );
    if ( !fp )
    {
        fprintf( stderr, "Cannot open %s\n", filename );
        exit( 0 );
    }

	if ( bpp < 0 )
	{
       tmpPtr = (unsigned short *)malloc( xs * ys * zs * sizeof(unsigned short) );

       if ( !tmpPtr )
	   {
          fprintf( stderr, "Memory allocation failed\n" );
           exit( 0 );
	   }

       //_swab( (char *)sPtr, (char *)tmpPtr, xs * ys * zs * sizeof(short) );
       swab( (char *)sPtr, (char *)tmpPtr, xs * ys * zs * sizeof(short) );
       fwrite( tmpPtr, xs * ys * zs, sizeof(short), fp );

	   free( (void *)tmpPtr );
	}
	else
	{
		fwrite( sPtr, xs * ys * zs, sizeof(short), fp );
	}

    fclose( fp );

    return( 1 );
}

int writeRawImage16SequentiallyBO( FILE *fp, unsigned short *sPtr, 
								  int xs, int ys, int bpp )
{
	unsigned short *tmpPtr;

	if (bpp > 0 )
	{    
		fwrite( sPtr, xs * ys, sizeof(short), fp );
	}
	else
	{
       tmpPtr = (unsigned short *)malloc( xs * ys * sizeof(unsigned short) );

       if ( !tmpPtr )
	   {
          fprintf( stderr, "Memory allocation failed\n" );
          exit( 0 );
	   }

       swab( (char *)sPtr, (char *)tmpPtr, xs * ys * sizeof(short) );
       fwrite( tmpPtr, xs * ys, sizeof(short), fp );

	   free( (void *)tmpPtr );
	}

    return( 1 );
}

int writeRawImageS16SequentiallyBO( FILE *fp, short *sPtr, 
								  int xs, int ys, int bpp )
{
	short *tmpPtr;

	if (bpp > 0 )
	{    
		fwrite( sPtr, xs * ys, sizeof(short), fp );
	}
	else
	{
       tmpPtr = (short *)malloc( xs * ys * sizeof(short) );

       if ( !tmpPtr )
	   {
          fprintf( stderr, "Memory allocation failed\n" );
          exit( 0 );
	   }

       swab( (char *)sPtr, (char *)tmpPtr, xs * ys * sizeof(short) );
       fwrite( tmpPtr, xs * ys, sizeof(short), fp );

	   free( (void *)tmpPtr );
	}


    return( 1 );
}


int writeRawImage8( char *filename, unsigned char *cPtr, int xs, int ys, int zs )
{
    FILE *fp;

    fp = fopen( filename, "wb" );
    if ( !fp )
    {
        fprintf( stderr, "Cannot open %s\n", filename );
        exit( 0 );
    }

    fwrite( cPtr, xs * ys * zs, sizeof(char), fp );

    fclose( fp );

    return( 1 );
}
\
int writeRawImage8Sequentially( FILE *fp, unsigned char *cPtr, int xs, int ys )
{
    fwrite( cPtr, xs * ys, sizeof(char), fp );

    return( 1 );
}
