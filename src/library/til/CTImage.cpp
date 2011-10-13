
// includes from STL
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// includes from TIL
#include "til/CTImage.h"



TIL_API FILE *readCTImageHeader( char *filename, int *xs, int *ys, int *zs, int *ctType,
                        float *voxelSizeX, float *voxelSizeY, float *voxelSizeZ )
{
    FILE *fp;
    unsigned char dummy[2];
    int numgot, siz = 0;
    unsigned char ctt;
    float fdummy;

    fp = fopen( filename, "rb" );
    if ( !fp )
    {
        fprintf( stderr, "File %s not found\n", filename );
        exit( 0 );
    }

    /* Read first 6 bytes */
    fseek( fp, 0L, SEEK_SET );
    fread( &dummy, 2, 1, fp );
    *zs = dummy[0] * 256 + dummy[1];
    fread( &dummy, 2, 1, fp );
    *xs = dummy[0] * 256 + dummy[1];
    fread( &dummy, 2, 1, fp );
    *ys = dummy[0] * 256 + dummy[1];

    /* Determine the type of the image 8 or 16 bit */
    fread( &dummy, 2, 1, fp );
    ctt = dummy[1];

    switch (ctt) {
    case 241: *ctType = CT_TYPE_DOUBLE;
              siz = 8;    
        break;
    case 242: *ctType = CT_TYPE_FLOAT;
              siz = 4;
        break;
    case 243: *ctType = CT_TYPE_UINT32;
              siz = 4;
        break;
    case 254: *ctType = CT_TYPE_UINT16;
              siz = 2;
        break;
    case 255: *ctType = CT_TYPE_UCHAR;
              siz = 1;
    }

    /* Move to the end of the data block */
    fseek( fp, (long)*xs * *ys * *zs * siz, SEEK_CUR );

    /* Determine the voxelSize of the image data */
    numgot = fread( &fdummy, sizeof(float), 1, fp );
    if (numgot==1)
        *voxelSizeX = fdummy;
    else
        *voxelSizeX = 1.0f;
    numgot = fread( &fdummy, sizeof(float), 1, fp );
    if (numgot==1)
        *voxelSizeY = fdummy;
    else
        *voxelSizeY = 1.0f;
    numgot = fread( &fdummy, sizeof(float), 1, fp );
    if (numgot==1)
        *voxelSizeZ = fdummy;
    else
        *voxelSizeZ = 1.0f;


    return( fp );
}

TIL_API unsigned short *readCTImage16( char *filename, int *xs, int *ys, int *zs,
                        float *voxelSizeX, float *voxelSizeY, float *voxelSizeZ )
{
    FILE *fp;
    unsigned char dummy[2];
    float fdummy;
    unsigned short *imgPtr;
    unsigned char *iPtr, *s;
    int numgot;

    fp = fopen( filename, "rb" );
    if ( !fp )
    {
        fprintf( stderr, "File %s not found\n", filename );
        exit( 0 );
    }

    /* Read first 6 bytes */
    fseek( fp, 0L, SEEK_SET );
    fread( &dummy, 2, 1, fp );
    *zs = dummy[0] * 256 + dummy[1];
    fread( &dummy, 2, 1, fp );
    *xs = dummy[0] * 256 + dummy[1];
    fread( &dummy, 2, 1, fp );
    *ys = dummy[0] * 256 + dummy[1];

    /* Test the type of the image */
    fread( &dummy, 2, 1, fp );
    if ( dummy[1] != 254 )
    {
        fprintf( stderr, "Unknown/wrong data type\n" );
        exit( -1 );
    }

    imgPtr = (unsigned short *)malloc( *xs **ys **zs * sizeof(short) + 1 );
    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
    /* Wrong byte order: Swap low and high bytes */
    numgot = fread( (char *)imgPtr + 1, sizeof(short), *xs **ys **zs, fp );

    iPtr = (unsigned char *)imgPtr;
    s = (unsigned char *)imgPtr + 2 **xs **ys **zs;
    while ( iPtr < s )
    {
        *iPtr = *( iPtr + 2 );
        iPtr += 2;
    }

    /* Determine the voxelSize of the image data */
    numgot = fread( &fdummy, sizeof(float), 1, fp );
    if (numgot==1)
        *voxelSizeX = fdummy;
    else
        *voxelSizeX = 1.0f;
    numgot = fread( &fdummy, sizeof(float), 1, fp );
    if (numgot==1)
        *voxelSizeY = fdummy;
    else
        *voxelSizeY = 1.0f;
    numgot = fread( &fdummy, sizeof(float), 1, fp );
    if (numgot==1)
        *voxelSizeZ = fdummy;
    else
        *voxelSizeZ = 1.0f;

    fclose( fp );

    return( imgPtr );
}

TIL_API unsigned short *readCTImage16ByNumber( FILE *fp, int slice, int xs, int ys )
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
    fseek( fp, (long)slice * xs * ys * sizeof(short) + 8L, SEEK_SET );

    imgPtr = (unsigned short *)malloc( xs * ys * sizeof(short) + 1 );

    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
    /* Wrong byte order: Swap low and high bytes */
    numgot = fread( (char *)imgPtr + 1, sizeof(short), xs * ys, fp );

    iPtr = (unsigned char *)imgPtr;
    s = (unsigned char *)imgPtr + 2 * xs * ys;
    while ( iPtr < s )
    {
        *iPtr = *( iPtr + 2 );
        iPtr += 2;
    }

    return( imgPtr );
}

TIL_API unsigned short *readCTImage16Sequentially( FILE *fp, int xs, int ys )
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
    /* Wrong byte order: Swap low and high bytes */
    numgot = fread( (char *)imgPtr + 1, sizeof(short), xs * ys, fp );

    iPtr = (unsigned char *)imgPtr;
    s = (unsigned char *)imgPtr + 2 * xs * ys;
    while ( iPtr < s )
    {
        *iPtr = *( iPtr + 2 );
        iPtr += 2;
    }

    return( imgPtr );
}

TIL_API unsigned char *readCTImage8( char *filename, int *xs, int *ys, int *zs,
                        float *voxelSizeX, float *voxelSizeY, float *voxelSizeZ )
{
    FILE *fp;
    unsigned char dummy[2];
    float fdummy;
    unsigned char *imgPtr;
    int numgot;

    fp = fopen( filename, "rb" );
    if ( !fp )
    {
        fprintf( stderr, "File %s not found\n", filename );
        exit( 0 );
    }

    /* Read first 6 bytes */
    fseek( fp, 0L, SEEK_SET );
    fread( &dummy, 2, 1, fp );
    *zs = dummy[0] * 256 + dummy[1];
    fread( &dummy, 2, 1, fp );
    *xs = dummy[0] * 256 + dummy[1];
    fread( &dummy, 2, 1, fp );
    *ys = dummy[0] * 256 + dummy[1];

    /* Determine the type of the image 8 or 16 bit */
    fread( &dummy, 2, 1, fp );
    if ( dummy[1] != 255 )
    {
        fprintf( stderr, "Unknown/wrong data type\n" );
        exit( -1 );
    }

    imgPtr = (unsigned char *)malloc( *xs **ys **zs * sizeof(char) );
    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
    /* Wrong byte order: Swap low and high bytes */
    numgot = fread( (char *)imgPtr, sizeof(char), *xs **ys **zs, fp );
    /* Determine the voxelSize of the image data */

    numgot = fread( &fdummy, sizeof(float), 1, fp );
    if (numgot==1)
        *voxelSizeX = fdummy;
    else
        *voxelSizeX = 1.0f;
    numgot = fread( &fdummy, sizeof(float), 1, fp );
    if (numgot==1)
        *voxelSizeY = fdummy;
    else
        *voxelSizeY = 1.0f;
    numgot = fread( &fdummy, sizeof(float), 1, fp );
    if (numgot==1)
        *voxelSizeZ = fdummy;
    else
        *voxelSizeZ = 1.0f;

    fclose( fp );

    return( imgPtr );
}

TIL_API unsigned char *readCTImage8ByNumber( FILE *fp, int slice, int xs, int ys )
{
    unsigned char *imgPtr;
    int numgot;

    if ( slice < 0 )
    {
        fprintf( stderr, "Incorrect slice number\n" );
        exit( 0 );
    }

    fseek( fp, 0L, SEEK_SET );
    fseek( fp, (long)slice * xs * ys * sizeof(char) + 8L, SEEK_SET );

    imgPtr = (unsigned char *)malloc( xs * ys * sizeof(char) );

    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
    /* Wrong byte order: Swap low and high bytes */
    numgot = fread( (char *)imgPtr, sizeof(char), xs * ys, fp );

    return( imgPtr );
}

TIL_API unsigned char *readCTImage8Sequentially( FILE *fp, int xs, int ys )
{
    unsigned char *imgPtr;
    int numgot;

    imgPtr = (unsigned char *)malloc( xs * ys * sizeof(char) );

    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
    /* Wrong byte order: Swap low and high bytes */
    numgot = fread( (char *)imgPtr, sizeof(char), xs * ys, fp );

    return( imgPtr );
}

TIL_API unsigned int *readCTImage32( char *filename, int *xs, int *ys, int *zs,
                        float *voxelSizeX, float *voxelSizeY, float *voxelSizeZ )
{
    FILE *fp;
    unsigned char dummy[2];
    float fdummy;
    unsigned int *imgPtr;
    int numgot;

    fp = fopen( filename, "rb" );
    if ( !fp )
    {
        fprintf( stderr, "File %s not found\n", filename );
        exit( 0 );
    }

    /* Read first 6 bytes */
    fseek( fp, 0L, SEEK_SET );
    fread( &dummy, 2, 1, fp );
    *zs = dummy[0] * 256 + dummy[1];
    fread( &dummy, 2, 1, fp );
    *xs = dummy[0] * 256 + dummy[1];
    fread( &dummy, 2, 1, fp );
    *ys = dummy[0] * 256 + dummy[1];

    /* Test the type of the image */
    fread( &dummy, 2, 1, fp );
    if ( dummy[1] != 243 )
    {
        fprintf( stderr, "Unknown/wrong data type\n" );
        exit( -1 );
    }

    imgPtr = (unsigned int *)malloc( *xs **ys **zs * sizeof(int) + 1 );
    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
    numgot = fread( (void *)imgPtr, sizeof(int), *xs **ys **zs, fp );

    /* Determine the voxelSize of the image data */
    numgot = fread( &fdummy, sizeof(float), 1, fp );
    if (numgot==1)
        *voxelSizeX = fdummy;
    else
        *voxelSizeX = 1.0f;
    numgot = fread( &fdummy, sizeof(float), 1, fp );
    if (numgot==1)
        *voxelSizeY = fdummy;
    else
        *voxelSizeY = 1.0f;
    numgot = fread( &fdummy, sizeof(float), 1, fp );
    if (numgot==1)
        *voxelSizeZ = fdummy;
    else
        *voxelSizeZ = 1.0f;


    fclose( fp );


    return( imgPtr );
}

TIL_API unsigned int *readCTImage32ByNumber( FILE *fp, int slice, int xs, int ys )
{
    unsigned int *imgPtr;
    int numgot;

    if ( slice < 0 )
    {
        fprintf( stderr, "Incorrect slice number\n" );
        exit( 0 );
    }

    fseek( fp, 0L, SEEK_SET );
    fseek( fp, (long)slice * xs * ys * sizeof(int) + 8L, SEEK_SET );

    imgPtr = (unsigned int *)malloc( xs * ys * sizeof(int) );

    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
    numgot = fread( (void *)imgPtr, sizeof(int), xs * ys, fp );

    return( imgPtr );
}

TIL_API unsigned int *readCTImage32Sequentially( FILE *fp, int xs, int ys )
{
    unsigned int *imgPtr;
    int numgot;

    imgPtr = (unsigned int *)malloc( xs * ys * sizeof(int) );

    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
    numgot = fread( (void *)imgPtr, sizeof(int), xs * ys, fp );

    return( imgPtr );
}


TIL_API float *readCTImageFloat( char *filename, int *xs, int *ys, int *zs,
                        float *voxelSizeX, float *voxelSizeY, float *voxelSizeZ )
{
    FILE *fp;
    unsigned char dummy[2];
    float fdummy;
    float *imgPtr;
    int numgot;

    fp = fopen( filename, "rb" );
    if ( !fp )
    {
        fprintf( stderr, "File %s not found\n", filename );
        exit( 0 );
    }

    /* Read first 6 bytes */
    fseek( fp, 0L, SEEK_SET );
    fread( &dummy, 2, 1, fp );
    *zs = dummy[0] * 256 + dummy[1];
    fread( &dummy, 2, 1, fp );
    *xs = dummy[0] * 256 + dummy[1];
    fread( &dummy, 2, 1, fp );
    *ys = dummy[0] * 256 + dummy[1];

    /* Test the type of the image */
    fread( &dummy, 2, 1, fp );
    if ( dummy[1] != 242 )
    {
        fprintf( stderr, "Unknown/wrong data type\n" );
        exit( -1 );
    }

    //imgPtr = (float *)malloc( *xs **ys **zs * sizeof(float) + 1 );
	imgPtr = new float [*xs * *ys * *zs + 1];
    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
    numgot = fread( (void *)imgPtr, sizeof(float), *xs **ys **zs, fp );

    /* Determine the voxelSize of the image data */
    numgot = fread( &fdummy, sizeof(float), 1, fp );
    if (numgot==1)
        *voxelSizeX = fdummy;
    else
        *voxelSizeX = 1.0f;
    numgot = fread( &fdummy, sizeof(float), 1, fp );
    if (numgot==1)
        *voxelSizeY = fdummy;
    else
        *voxelSizeY = 1.0f;
    numgot = fread( &fdummy, sizeof(float), 1, fp );
    if (numgot==1)
        *voxelSizeZ = fdummy;
    else
        *voxelSizeZ = 1.0f;


    fclose( fp );


    return( imgPtr );
}

TIL_API float *readCTImageFloatByNumber( FILE *fp, int slice, int xs, int ys )
{
    float *imgPtr;
    int numgot;

    if ( slice < 0 )
    {
        fprintf( stderr, "Incorrect slice number\n" );
        exit( 0 );
    }

    fseek( fp, 0L, SEEK_SET );
    fseek( fp, (long)slice * xs * ys * sizeof(float) + 8L, SEEK_SET );

    imgPtr = (float *)malloc( xs * ys * sizeof(float) );

    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
    numgot = fread( (void *)imgPtr, sizeof(float), xs * ys, fp );

    return( imgPtr );
}

TIL_API float *readCTImageFloatSequentially( FILE *fp, int xs, int ys )
{
    float *imgPtr;
    int numgot;

    imgPtr = (float *)malloc( xs * ys * sizeof(float) );

    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
    numgot = fread( (void *)imgPtr, sizeof(float), xs * ys, fp );

    return( imgPtr );
}

TIL_API double *readCTImageDouble( char *filename, int *xs, int *ys, int *zs,
                        float *voxelSizeX, float *voxelSizeY, float *voxelSizeZ )
{
    FILE *fp;
    unsigned char dummy[2];
    float fdummy;
    double *imgPtr;
    int numgot;

    fp = fopen( filename, "rb" );
    if ( !fp )
    {
        fprintf( stderr, "File %s not found\n", filename );
        exit( 0 );
    }

    /* Read first 6 bytes */
    fseek( fp, 0L, SEEK_SET );
    fread( &dummy, 2, 1, fp );
    *zs = dummy[0] * 256 + dummy[1];
    fread( &dummy, 2, 1, fp );
    *xs = dummy[0] * 256 + dummy[1];
    fread( &dummy, 2, 1, fp );
    *ys = dummy[0] * 256 + dummy[1];

    /* Test the type of the image */
    fread( &dummy, 2, 1, fp );
    if ( dummy[1] != 241 )
    {
        fprintf( stderr, "Unknown/wrong data type\n" );
        exit( -1 );
    }

    //imgPtr = (double *)malloc( *xs **ys **zs * sizeof(double) + 1 );

	imgPtr = new double[*xs **ys **zs+1];

    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
    numgot = fread( (void *)imgPtr, sizeof(double), *xs **ys **zs, fp );

    /* Determine the voxelSize of the image data */
    numgot = fread( &fdummy, sizeof(float), 1, fp );
    if (numgot==1)
        *voxelSizeX = fdummy;
    else
        *voxelSizeX = 1.0f;
    numgot = fread( &fdummy, sizeof(float), 1, fp );
    if (numgot==1)
        *voxelSizeY = fdummy;
    else
        *voxelSizeY = 1.0f;
    numgot = fread( &fdummy, sizeof(float), 1, fp );
    if (numgot==1)
        *voxelSizeZ = fdummy;
    else
        *voxelSizeZ = 1.0f;


    fclose( fp );


    return( imgPtr );
}

TIL_API double *readCTImageDoubleByNumber( FILE *fp, int slice, int xs, int ys )
{
    double *imgPtr;
    int numgot;

    if ( slice < 0 )
    {
        fprintf( stderr, "Incorrect slice number\n" );
        exit( 0 );
    }

    fseek( fp, 0L, SEEK_SET );
    fseek( fp, (long)slice * xs * ys * sizeof(double) + 8L, SEEK_SET );

    imgPtr = (double *)malloc( xs * ys * sizeof(double) );

    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
    numgot = fread( (void *)imgPtr, sizeof(double), xs * ys, fp );

    return( imgPtr );
}

TIL_API double *readCTImageDoubleSequentially( FILE *fp, int xs, int ys )
{
    double *imgPtr;
    int numgot;

    imgPtr = (double *)malloc( xs * ys * sizeof(double) );

    if ( !imgPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Read image data */
    numgot = fread( (void *)imgPtr, sizeof(double), xs * ys, fp );

    return( imgPtr );
}


TIL_API FILE *writeCTImageDoubleHeader( char *filename, double * dPtr, int xs, int ys, int zs )
{
    FILE *fp;
    unsigned char dummy[2];

    fp = fopen( filename, "wb" );
    if ( !fp )
    {
        fprintf( stderr, "Cannot open %s\n", filename );
        exit( 0 );
    }

    /* Write first 8 bytes */
    dummy[0] = (unsigned char)(zs / 256);
    dummy[1] = (unsigned char)(zs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(xs / 256);
    dummy[1] = (unsigned char)(xs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(ys / 256);
    dummy[1] = (unsigned char)(ys % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = 255;
    dummy[1] = 241;
    fwrite( dummy, 2, 1, fp );

    return( fp );
}

TIL_API int writeCTImageDouble( char *filename, double *dPtr, int xs, int ys, int zs,
                   float voxelSizeX, float voxelSizeY, float voxelSizeZ )
{
    FILE *fp;
    unsigned char dummy[2];
    int rv;

    fp = fopen( filename, "wb" );
    if ( !fp )
    {
        fprintf( stderr, "Cannot open %s\n", filename );
        exit( 0 );
    }

    /* Write first 8 bytes */
    dummy[0] = (unsigned char)(zs / 256);
    dummy[1] = (unsigned char)(zs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(xs / 256);
    dummy[1] = (unsigned char)(xs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(ys / 256);
    dummy[1] = (unsigned char)(ys % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = 255;
    dummy[1] = 241;
    fwrite( dummy, 2, 1, fp );

    rv = fwrite( dPtr, xs * ys * zs, sizeof(double), fp );

    closeCTImage( fp, voxelSizeX, voxelSizeY, voxelSizeZ );

    return( rv==sizeof(double) );
}

TIL_API int writeCTImageDoubleSequentially( FILE *fp, double *dPtr, int xs, int ys )
{
    int rv;


    rv = fwrite( dPtr, xs * ys, sizeof(double), fp );

    return( rv==sizeof(double) );
}


TIL_API FILE *writeCTImageFloatHeader( char *filename, float *fPtr, int xs, int ys, int zs )
{
    FILE *fp;
    unsigned char dummy[2];

    fp = fopen( filename, "wb" );
    if ( !fp )
    {
        fprintf( stderr, "Cannot open %s\n", filename );
        exit( 0 );
    }

    /* Write first 8 bytes */
    dummy[0] = (unsigned char)(zs / 256);
    dummy[1] = (unsigned char)(zs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(xs / 256);
    dummy[1] = (unsigned char)(xs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(ys / 256);
    dummy[1] = (unsigned char)(ys % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = 255;
    dummy[1] = 242;
    fwrite( dummy, 2, 1, fp );

    return( fp );
}


TIL_API int writeCTImageFloatSequentially( FILE *fp, float *fPtr, int xs, int ys )
{
    int rv;


    rv = fwrite( fPtr, xs * ys, sizeof(float), fp );

    return( rv==sizeof(float) );
}


TIL_API int writeCTImageFloat( char *filename, float *fPtr, int xs, int ys, int zs,
                   float voxelSizeX, float voxelSizeY, float voxelSizeZ )
{
    FILE *fp;
    unsigned char dummy[2];
    int rv;

    fp = fopen( filename, "wb" );
    if ( !fp )
    {
        fprintf( stderr, "Cannot open %s\n", filename );
        exit( 0 );
    }

    /* Write first 8 bytes */
    dummy[0] = (unsigned char)(zs / 256);
    dummy[1] = (unsigned char)(zs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(xs / 256);
    dummy[1] = (unsigned char)(xs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(ys / 256);
    dummy[1] = (unsigned char)(ys % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = 255;
    dummy[1] = 242;
    fwrite( dummy, 2, 1, fp );

    rv = fwrite( fPtr, xs * ys * zs, sizeof(float), fp );

    closeCTImage( fp, voxelSizeX, voxelSizeY, voxelSizeZ );

    return( rv==sizeof(float) );
}


TIL_API FILE *writeCTImage32Header( char *filename, unsigned int *iPtr, int xs, int ys, int zs )
{
    FILE *fp;
    unsigned char dummy[2];

    fp = fopen( filename, "wb" );
    if ( !fp )
    {
        fprintf( stderr, "Cannot open %s\n", filename );
        exit( 0 );
    }

    /* Write first 8 bytes */
    dummy[0] = (unsigned char)(zs / 256);
    dummy[1] = (unsigned char)(zs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(xs / 256);
    dummy[1] = (unsigned char)(xs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(ys / 256);
    dummy[1] = (unsigned char)(ys % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = 255;
    dummy[1] = 243;
    fwrite( dummy, 2, 1, fp );

    return( fp );
}

TIL_API int writeCTImage32Sequentially( FILE *fp, unsigned int *iPtr, int xs, int ys )
{
    int rv;


    rv = fwrite( iPtr, xs * ys, sizeof(int), fp );

    return( rv==sizeof(int) );
}


TIL_API int writeCTImage32( char *filename, unsigned int *iPtr, int xs, int ys, int zs,
                   float voxelSizeX, float voxelSizeY, float voxelSizeZ )
{
    FILE *fp;
    unsigned char dummy[2];
    int rv;

    fp = fopen( filename, "wb" );
    if ( !fp )
    {
        fprintf( stderr, "Cannot open %s\n", filename );
        exit( 0 );
    }

    /* Write first 8 bytes */
    dummy[0] = (unsigned char)(zs / 256);
    dummy[1] = (unsigned char)(zs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(xs / 256);
    dummy[1] = (unsigned char)(xs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(ys / 256);
    dummy[1] = (unsigned char)(ys % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = 255;
    dummy[1] = 243;
    fwrite( dummy, 2, 1, fp );

    rv = fwrite( iPtr, xs * ys * zs, sizeof(int), fp );

    closeCTImage( fp, voxelSizeX, voxelSizeY, voxelSizeZ );

    return( rv==sizeof(int) );
}


TIL_API int writeCTImage16( char *filename, unsigned short *sPtr, int xs, int ys, int zs,
                   float voxelSizeX, float voxelSizeY, float voxelSizeZ )
{
    FILE *fp;
    unsigned char dummy[2];
    unsigned short *tmpPtr;
    int rv;

    fp = fopen( filename, "wb" );
    if ( !fp )
    {
        fprintf( stderr, "Cannot open %s\n", filename );
        exit( 0 );
    }

    tmpPtr = (unsigned short *)malloc( xs * ys * zs * sizeof(short) );
    if ( !tmpPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    /* Write first 8 bytes */
    dummy[0] = (unsigned char)(zs / 256);
    dummy[1] = (unsigned char)(zs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(xs / 256);
    dummy[1] = (unsigned char)(xs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(ys / 256);
    dummy[1] = (unsigned char)(ys % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = 255;
    dummy[1] = 254;
    fwrite( dummy, 2, 1, fp );

    //_swab( (char *)sPtr, (char *)tmpPtr, xs * ys * zs * sizeof(short) );
    swab( (char *)sPtr, (char *)tmpPtr, xs * ys * zs * sizeof(short) );
    rv = fwrite( tmpPtr, xs * ys * zs, sizeof(short), fp );

    closeCTImage( fp, voxelSizeX, voxelSizeY, voxelSizeZ );

    free( (void *)tmpPtr );

    return( rv==sizeof(short) );
}

TIL_API int writeCTImage16Sequentially( FILE *fp, unsigned short *sPtr, int xs, int ys )
{
    unsigned short *tmpPtr;
    int rv;

    tmpPtr = (unsigned short *)malloc( xs * ys * sizeof(short) );
    if ( !tmpPtr )
    {
        fprintf( stderr, "Memory allocation failed\n" );
        exit( 0 );
    }

    //_swab( (char *)sPtr, (char *)tmpPtr, xs * ys * sizeof(short) );
    swab( (char *)sPtr, (char *)tmpPtr, xs * ys * sizeof(short) );
    rv = fwrite( tmpPtr, xs * ys, sizeof(short), fp );

    free( (void *)tmpPtr );

    return( rv==sizeof(short) );
}

TIL_API FILE *writeCTImage16Header( char *filename, unsigned short *sPtr, int xs, int ys, int zs )
{
    FILE *fp;
    unsigned char dummy[2];

    fp = fopen( filename, "wb" );
    if ( !fp )
    {
        fprintf( stderr, "Cannot open %s\n", filename );
        exit( 0 );
    }

    /* Write first 8 bytes */
    dummy[0] = (unsigned char)(zs / 256);
    dummy[1] = (unsigned char)(zs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(xs / 256);
    dummy[1] = (unsigned char)(xs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(ys / 256);
    dummy[1] = (unsigned char)(ys % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = 255;
    dummy[1] = 254;
    fwrite( dummy, 2, 1, fp );

    return( fp );
}

TIL_API FILE *writeCTImage8Header( char *filename, unsigned char *sPtr, int xs, int ys, int zs )
{
    FILE *fp;
    unsigned char dummy[2];

    fp = fopen( filename, "wb" );
    if ( !fp )
    {
        fprintf( stderr, "Cannot open %s\n", filename );
        exit( 0 );
    }

    /* Write first 8 bytes */
    dummy[0] = (unsigned char)(zs / 256);
    dummy[1] = (unsigned char)(zs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(xs / 256);
    dummy[1] = (unsigned char)(xs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(ys / 256);
    dummy[1] = (unsigned char)(ys % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = 255;
    dummy[1] = 255;
    fwrite( dummy, 2, 1, fp );

    return( fp );
}

TIL_API int writeCTImage8( char *filename, unsigned char *sPtr, int xs, int ys, int zs,
                   float voxelSizeX, float voxelSizeY, float voxelSizeZ )
{
    FILE *fp;
    unsigned char dummy[2];
    int rv;

    fp = fopen( filename, "wb" );
    if ( !fp )
    {
        fprintf( stderr, "Cannot open %s\n", filename );
        exit( 0 );
    }

    /* Write first 8 bytes */
    dummy[0] = (unsigned char)(zs / 256);
    dummy[1] = (unsigned char)(zs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(xs / 256);
    dummy[1] = (unsigned char)(xs % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = (unsigned char)(ys / 256);
    dummy[1] = (unsigned char)(ys % 256);
    fwrite( dummy, 2, 1, fp );
    dummy[0] = 255;
    dummy[1] = 255;
    fwrite( dummy, 2, 1, fp );

    rv = fwrite( sPtr, xs * ys * zs, sizeof(char), fp );
    closeCTImage( fp, voxelSizeX, voxelSizeY, voxelSizeZ );

    return( rv==sizeof(char) );
}

TIL_API int writeCTImage8Sequentially( FILE *fp, unsigned char *cPtr, int xs, int ys )
{
    int rv;


    rv = fwrite( cPtr, xs * ys, sizeof(char), fp );

    return( rv==sizeof(char) );
}

TIL_API void closeCTImage( FILE *fp, float voxelSizeX, float voxelSizeY, float voxelSizeZ )
{
    fwrite( &voxelSizeX, sizeof(float), 1, fp );
    fwrite( &voxelSizeY, sizeof(float), 1, fp );
    fwrite( &voxelSizeZ, sizeof(float), 1, fp );

    fclose( fp );

    return;
}

