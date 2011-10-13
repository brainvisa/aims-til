
#include "connectedComponents3D.h"


int correctLabels(unsigned int n, int *equList);
int findLastLabel(int *list, int i);
int fillEquivList(int *lst, int label);
int incEquivList(unsigned int n, int *lst);


TIL_API void connectedComponents3D(Ptr<Image<unsigned short> > &img, unsigned short thres )
{
    int i;
    unsigned int j;
    int row, col, slice, imageSize, tmp, tmpz;
    int xSize, ySize, zSize;
    int tmpy, tmpyy, x, y, z;
    int ll, lu, lf;
    unsigned short *labelImg;
    unsigned int labelCounter;
    unsigned short c, f, l, u;
    int *lut;
    unsigned short *image;

    // Assumption: There is at least one free column/row at each side of
    // a slice. The first and the last slices are set to background!!!
    //
	// Search for all voxel which are higher or equal than the given threshold
    // Background is assumed to have the value 0

    /* A maximal number of 65536 volumes is possible */
    lut = (int *) malloc( 65536 * sizeof(int));
    if ( lut == NULL )
    {
        // throw Exception("Out of memory");
    }

    // Do not allocate new memory. Overwrite the original image.
    //
    xSize = img->getX();
    ySize = img->getY();
    zSize = img->getZ();

    labelImg = img->getPointer();

    labelCounter = 0;
    lut[0] = 0;

    image = img->getPointer();
    imageSize = xSize * ySize;

    // Set the first slice to 0
    memset((char *) image, 0, imageSize * sizeof(short));
    for (slice = 1; slice < zSize; ++slice)
    {
        // Set the left margin to 0
        tmpz = slice * imageSize;
        for (i = 1; i < ySize; image[tmpz + xSize * i++] = 0)
            ;
    }

    /* Calculate Label image */
    for (slice = 1; slice < zSize - 1; ++slice)
    {
		//printf("Completed %i %%\r", (int)((100.0 * slice) / (zSize-1) + 0.5));
        tmpz = slice * imageSize;
        for (row = 1; row < ySize - 1; ++row)
        {
            tmp = row * xSize + tmpz;
            for (col = 1; col < xSize; ++col)
            {
                i = tmp + col;

                // the current voxel
                c = *(image + i);

                if (c >= thres)
                {
                    // Check if the pixels above, left and in front
                    // do have a label already.
                    // The pixel above
                    u = *(labelImg+i-xSize);

                    // on the left side
                    l = *(labelImg+i-1);

                    // and one slice before
                    f = *(labelImg+i-imageSize);

                    if ( (!l) && (!u) && (!f) )
                    {
                        // None of the neighbors has a volume:
                        // Assign a new label to the volume
                        ++labelCounter;

                        // We can only handle 65535 labels.
                        // If we are running out of free labels, try to
                        // get additional labels by making all labels unique.
                        if (labelCounter == 65536)
                        {
                            labelCounter = correctLabels( 65535, lut );

                            if ( labelCounter >= 65534 )
                            {
                                // throw Exception("Too many labels");
                            }

                            for ( z = 0; z < zSize; ++z )
                            {
                                tmpy = z * imageSize;
                                for ( y = 0; y < ySize; ++y )
                                {
                                    tmpyy = y * xSize + tmpy;
                                    for ( x = 0; x < xSize; ++x )
                                    {
                                        if ( *(labelImg+tmpyy+x) )
                                            *(labelImg+tmpyy+x) = lut[*(labelImg+tmpyy+x)];
                                    }
                                }   
                            }   

                            ++labelCounter;

                            for ( j = 0; j < labelCounter; ++j )
                                lut[j] = (int)j;

                            for ( j = labelCounter; j < 65536; ++j )
                                lut[j] = 0;
                        }

                        labelImg[i] = (unsigned short) labelCounter;
                        lut[labelCounter] = (int) labelCounter;

                        continue;
                    }

                    // Only one of the three neighbors is set
                    if ((u) && (!l) && (!f))
                    {
                        labelImg[i] = labelImg[i - xSize];
                    }
                    else 
                        if ((l) && (!u) && (!f))
                        {
                            labelImg[i] = labelImg[i - 1];
                        } 
                        else 
                            if ((f) && (!u) && (!l))
                            {
                                labelImg[i] = labelImg[i - imageSize];
                            } 
                            else 
                                // Two of the neighbors are set
                                if ((l) && (u) && (!f))
                                {
                                    ll = labelImg[i - 1];
                                    lu = labelImg[i - xSize];

                                    labelImg[i] = (unsigned short) lu;

                                    if (lu != ll)
                                    {
                                        ll = findLastLabel(lut, ll);
                                        lu = findLastLabel(lut, lu);

                                        lut[ll] = lu;
                                    }
                                } 
                                else 
                                    if ((l) && (f) && (!u))
                                    {
                                        ll = labelImg[i - 1];
                                        lf = labelImg[i - imageSize];

                                        labelImg[i] = (unsigned short) lf;

                                        if (lf != ll)
                                        {
                                            ll = findLastLabel(lut, ll);
                                            lf = findLastLabel(lut, lf);

                                            lut[ll] = lf;
                                        }
                                    } 
                                    else 
                                        if ((u) && (f) && (!l))
                                        {
                                            lu = labelImg[i - xSize];
                                            lf = labelImg[i - imageSize];

                                            labelImg[i] = (unsigned short) lf;

                                            if (lf != lu)
                                            {
                                                lu = findLastLabel(lut, lu);
                                                lf = findLastLabel(lut, lf);

                                                lut[lu] = lf;
                                            }
                                        } 
                                        else 
                                            if ((u) && (l) && (f))
                                            {
                                                ll = labelImg[i - 1];
                                                lu = labelImg[i - xSize];
                                                lf = labelImg[i - imageSize];

                                                if ((lu == ll) && (ll == lf))
                                                {
                                                    // The 3 neighbors have the same label
                                                    labelImg[i] = (unsigned short) lu;
                                                } 
                                                else if ((lu == ll) && (lu != lf) && (ll != lf))
                                                {
                                                    labelImg[i] = lu;

                                                    lu = findLastLabel(lut, lu);
                                                    lf = findLastLabel(lut, lf);

                                                    lut[lf] = lu;
                                                } 
                                                else if ((lu == lf) && (lu != ll) && (lf != ll))
                                                {
                                                    labelImg[i] = lu;

                                                    lu = findLastLabel(lut, lu);
                                                    ll = findLastLabel(lut, ll);

                                                    lut[ll] = lu;
                                                } 
                                                else if ((lf == ll) && (lf != lu) && (ll != lu))
                                                {
                                                    labelImg[i] = lf;

                                                    lf = findLastLabel(lut, lf);
                                                    lu = findLastLabel(lut, lu);

                                                    lut[lu] = lf;
                                                } 
                                                else if ((lf != ll) && (lf != lu) && (ll != lu))
                                                {
                                                    labelImg[i] = lf;

                                                    lf = findLastLabel(lut, lf);
                                                    lu = findLastLabel(lut, lu);
                                                    ll = findLastLabel(lut, ll);

                                                    lut[lu] = lf;
                                                    lut[ll] = lf;
                                                }
                                            } 
                } // if (c>=thres)
            } // for col
        } // for row
    } // for slice

    correctLabels(labelCounter, lut );

    // Update the label image
    for ( i = 0; i < imageSize * zSize; ++i )
        *(labelImg+i) = lut[*(labelImg+i)];

    free((void *) lut);

    return;
}



TIL_API void connectedComponents2D(Ptr<Image<unsigned short> > &img, unsigned short thres )
{
    int i;
    unsigned int j;
    int row, col, slice, imageSize, tmp, tmpz;
    int xSize, ySize, zSize;
    int tmpy, tmpyy, x, y, z;
    //int ll, lu, lf;
    int ll, lu;
    unsigned short *labelImg;
    unsigned int labelCounter;
    //unsigned short c, f, l, u;
    unsigned short c, l, u;
    int *lut;
    unsigned short *image;

    // Assumption: There is at least one free column/row at each side of
    // a slice. The first and the last slices are set to background!!!
    //
	// Search for all voxel which are higher or equal than the given threshold
    // Background is assumed to have the value 0

    /* A maximal number of 65536 volumes is possible */
    lut = (int *) malloc( 65536 * sizeof(int));
    if ( lut == NULL )
    {
        // throw Exception("Out of memory");
    }

    // Do not allocate new memory. Overwrite the original image.
    //
    xSize = img->getX();
    ySize = img->getY();
    zSize = img->getZ();

    labelImg = img->getPointer();

    labelCounter = 0;
    lut[0] = 0;

    image = img->getPointer();
    imageSize = xSize * ySize;

    // Set the first slice to 0
    memset((char *) image, 0, imageSize * sizeof(short));
    for (slice = 1; slice < zSize; ++slice)
    {
        // Set the left margin to 0
        tmpz = slice * imageSize;
        for (i = 1; i < ySize; image[tmpz + xSize * i++] = 0)
            ;
    }

    /* Calculate Label image */
    for (slice = 1; slice < zSize - 1; ++slice)
    {
		printf("Completed %i %%\r", (int)((100.0 * slice) / (zSize-1) + 0.5));
        tmpz = slice * imageSize;
        for (row = 1; row < ySize - 1; ++row)
        {
            tmp = row * xSize + tmpz;
            for (col = 1; col < xSize; ++col)
            {
                i = tmp + col;

                // the current voxel
                c = *(image + i);

                if (c >= thres)
                {
                    // Check if the pixels above, left and in front
                    // do have a label already.
                    // The pixel above
                    u = *(labelImg+i-xSize);

                    // on the left side
                    l = *(labelImg+i-1);

                    // and one slice before
                    // f = *(labelImg+i-imageSize);

                    if ( (!l) && (!u) )//&& (!f) )
                    {
                        // None of the neighbors has a volume:
                        // Assign a new label to the volume
                        ++labelCounter;

                        // We can only handle 65535 labels.
                        // If we are running out of free labels, try to
                        // get additional labels by making all labels unique.
                        if (labelCounter == 65536)
                        {
                            labelCounter = correctLabels( 65535, lut );

                            if ( labelCounter >= 65534 )
                            {
                                // throw Exception("Too many labels");
                            }

                            for ( z = 0; z < zSize; ++z )
                            {
                                tmpy = z * imageSize;
                                for ( y = 0; y < ySize; ++y )
                                {
                                    tmpyy = y * xSize + tmpy;
                                    for ( x = 0; x < xSize; ++x )
                                    {
                                        if ( *(labelImg+tmpyy+x) )
                                            *(labelImg+tmpyy+x) = lut[*(labelImg+tmpyy+x)];
                                    }
                                }   
                            }   

                            ++labelCounter;

                            for ( j = 0; j < labelCounter; ++j )
                                lut[j] = (int)j;

                            for ( j = labelCounter; j < 65536; ++j )
                                lut[j] = 0;
                        }

                        labelImg[i] = (unsigned short) labelCounter;
                        lut[labelCounter] = (int) labelCounter;

                        continue;
                    }

                    // Only one of the three neighbors is set
                    if ((u) && (!l))// && (!f))
                    {
                        labelImg[i] = labelImg[i - xSize];
                    }
                    else 
                        if ((l) && (!u))// && (!f))
                        {
                            labelImg[i] = labelImg[i - 1];
                        } 
                        else /*
                            if ((f) && (!u) && (!l))
                            {
                               labelImg[i] = labelImg[i - imageSize];
                            } 
                            else */
                                // Two of the neighbors are set
                                if ((l) && (u))// && (!f))
                                {
                                    ll = labelImg[i - 1];
                                    lu = labelImg[i - xSize];

                                    labelImg[i] = (unsigned short) lu;

                                    if (lu != ll)
                                    {
                                        ll = findLastLabel(lut, ll);
                                        lu = findLastLabel(lut, lu);

                                        lut[ll] = lu;
                                    }
                                } /*
                                else
                                    if ((l) && (f) && (!u))
                                    {
                                        ll = labelImg[i - 1];
                                        lf = labelImg[i - imageSize];

                                        labelImg[i] = (unsigned short) lf;

                                        if (lf != ll)
                                        {
                                            ll = findLastLabel(lut, ll);
                                            lf = findLastLabel(lut, lf);

                                            lut[ll] = lf;
                                        }
                                    } 
                                    else *//*
                                        if ((u) && (f) && (!l))
                                        {
                                            lu = labelImg[i - xSize];
                                            lf = labelImg[i - imageSize];

                                            labelImg[i] = (unsigned short) lf;

                                            if (lf != lu)
                                            {
                                                lu = findLastLabel(lut, lu);
                                                lf = findLastLabel(lut, lf);

                                                lut[lu] = lf;
                                            }
                                        } *//*
                                        else 
                                            if ((u) && (l) && (f))
                                            {
                                                ll = labelImg[i - 1];
                                                lu = labelImg[i - xSize];
                                                lf = labelImg[i - imageSize];

                                                if ((lu == ll) && (ll == lf))
                                                {
                                                    // The 3 neighbors have the same label
                                                    labelImg[i] = (unsigned short) lu;
                                                } 
                                                else if ((lu == ll) && (lu != lf) && (ll != lf))
                                                {
                                                    labelImg[i] = lu;

                                                    lu = findLastLabel(lut, lu);
                                                    lf = findLastLabel(lut, lf);

                                                    lut[lf] = lu;
                                                } 
                                                else if ((lu == lf) && (lu != ll) && (lf != ll))
                                                {
                                                    labelImg[i] = lu;

                                                    lu = findLastLabel(lut, lu);
                                                    ll = findLastLabel(lut, ll);

                                                    lut[ll] = lu;
                                                } 
                                                else if ((lf == ll) && (lf != lu) && (ll != lu))
                                                {
                                                    labelImg[i] = lf;

                                                    lf = findLastLabel(lut, lf);
                                                    lu = findLastLabel(lut, lu);

                                                    lut[lu] = lf;
                                                } 
                                                else if ((lf != ll) && (lf != lu) && (ll != lu))
                                                {
                                                    labelImg[i] = lf;

                                                    lf = findLastLabel(lut, lf);
                                                    lu = findLastLabel(lut, lu);
                                                    ll = findLastLabel(lut, ll);

                                                    lut[lu] = lf;
                                                    lut[ll] = lf;
                                                }
                                            }*/
                } // if (c>=thres)
            } // for col
        } // for row
    } // for slice

    correctLabels(labelCounter, lut );

    // Update the label image
    for ( i = 0; i < imageSize * zSize; ++i )
        *(labelImg+i) = lut[*(labelImg+i)];

    free((void *) lut);

    return;
}

// correctLabel makes the label of a volume unique by merging neighbored volumes 
// with different label
int correctLabels(unsigned int n, int *equList)
{
    int i;
    int numVOI;

    // n: number of regions/volumes 
    // equList: look up table with all valid labels
	// !!!! Does not update the label image !!!! */
    for (i = 1; i <= (int) n; ++i)
        equList[i] = fillEquivList(equList, i);

    numVOI = incEquivList(n, equList);

    return (numVOI);
}


int findLastLabel(int *list, int i)
{
    if (i == list[i])
        return (i);
    else
        return (findLastLabel(list, list[i]));
}

// Termination criterion: Index == region/volume number */
int fillEquivList(int *lst, int label)
{
    int equ;

    equ = lst[label];
    if (equ == lst[equ])
        return (equ);
    else
        return (lst[equ] = fillEquivList(lst, equ));
}


// Renumber all labels so that they are continously and do start with 1
int incEquivList(unsigned int n, int *lst)
{
    int i, j, alt, neu;

    neu = -1;
    for (i = 1; i <= (int) n; ++i)
    {
        alt = lst[i];
        if (alt >= 0)
        {
            for (j = i; j <= (int) n; ++j)
                if (lst[j] == alt)
                {
                    lst[j] = neu;
                }

            --neu;
        }
    }

    for (i = 1; i <= (int) n; ++i)
    {
        lst[i] = -lst[i];
    }

    return (-neu);
}


