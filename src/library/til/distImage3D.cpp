#include "distImage3D.h"




// Namespace 

namespace til {


TIL_API void distance3D( Ptr<ImageC<unsigned char> > &img )
{
    int changed;
    int i, j, k, l, offset[26], count;
    int kEnd;
    int tmpy, tmpz, imgSize;
    int xSize, ySize, zSize;
    unsigned char *bPtr;


    //
    // Assumption: The foreground is set to 255!!!!!
    //

    xSize = img->getDimX();
    ySize = img->getDimY();
    zSize = img->getDimZ();

    imgSize = xSize * ySize;
    bPtr = img->getPointer();

    // Position of the neighbor with respect to the reference point in the middle
    offset[0] = -xSize;
    offset[1] = -1;
    offset[2] = +1;
    offset[3] = xSize;
    offset[4] = -imgSize;
    offset[5] = imgSize;

    count = 0;
    do
    {
        changed = 0;
        count++;
        if (count > 255)
        {
            // Fatal error
            // throw Exception("Fatal Error");
            memset( img->getPointer(), 0, xSize * ySize * zSize );
            return;
        }

        for ( l = 1; l < zSize - 1; ++l )
        {
            tmpz = l * imgSize;
            for ( j = 1; j < (ySize - 1); j++ )
            {
                tmpy = tmpz + j * xSize;

                kEnd = 6;

                for ( i = 1; i < (xSize - 1); i++ )
                {
                    if ( *(bPtr+tmpy+i) == distImage3D_FOREGROUND_VALUE )
                    {
                        for ( k = 0; k < kEnd; k++ )
                        {
                            if ( *(bPtr+tmpy+i+offset[k]) < count )
                            {
                                *(bPtr+tmpy+i) = count;
                                changed = 1;
                            }
                        }   // end for(6-topology)
                    }       // end if(point SET)
                }           // for i
            }               // for j
        }                   // for l
    } while ( changed == 1 );

    return;
}

} // namespace

