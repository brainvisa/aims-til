
#pragma warning (disable: 4786)

#include <iostream>



// NB: It is important to include all header files in there, for at least
// two reasons:
// (1) Non-template classes that are defined in the header only (such as
//     region growing plugins) would not be compiled otherwise
// (2) It makes it possible to instantiate and debug template classes

/*
#include "til/til_common.h"
#include "til/BasicPixelActions.h"
#include "til/BasicPixelTests.h"
#include "til/Box.h"
#include "til/connectedComponents3D.h"
#include "til/convert.h"
#include "til/CTimage.h"
#include "til/eigen3D.h"
#include "til/GaussianNormalizedCorrelation.h"
#include "til/if_then_else.h"
#include "til/ImageC.h"
#include "til/imageArith.h"
#include "til/imageBasicStats.h"
#include "til/imageNC.h"
#include "til/imagetools.h"
#include "til/morpho.h"
#include "til/neighborhoodconfigurations.h"
#include "til/NormalizedCorrelation.h"
#include "til/pointListTools.h"
#include "til/Ptr.h"
#include "til/Range.h"
#include "til/read.h"
#include "til/recursiveFilters.h"
#include "til/regionGrowing.h"
#include "til/regionGrowingPlugin.h"
#include "til/RF4.h"
//#include "til/SmartObject.h"
#include "til/StructureTensor.h"
#include "til/SubImageExtractor.h"
#include "til/TExprOperators.h"
#include "til/write.h"


#include "til/ConstImageLinearIterator.h"
#include "til/ConstImageVolumetricIterator.h"
#include "til/ImageLinearIterator.h"
#include "til/ImageVolumetricIterator.h"

#include "til/ConstImageNCLinearIterator.h"
#include "til/ConstImageNCVolumetricIterator.h"
#include "til/ImageNCLinearIterator.h"
*/


#if 0

template Box<double>;
template FeatureBox<unsigned char>;
template FeatureDiamond<unsigned char>;
template FeatureEnclosingBox<unsigned char, unsigned short>;
template FeatureGaussianCorrelation<Image<int>, Image<double> >;
template FeatureLocalIsotropy<Image<int>, Image<double> >;
template FeatureWallAttached<Image<double> >;
template GaussianNormalizedCorrelation<Image<int>, Image<double> >;
template Image<double>;
template ImageNC<double>;
template NormalizedCorrelation<Image<int>, Image<double> >;
template RF4<double>;
template StructureTensor<Image<int>, Image<double> >;
template SubImageExtractor<Image<int>,Image<float> >;
//template SymMatrix3<double>;

template ConstImageLinearIterator<int>;
template ImageLinearIterator<int>;
template ConstImageVolumetricIterator<int>;
template ImageVolumetricIterator<int>;
template < class T > inline T norm2(T a, T b) { return sqrt(a*a+b*b);}
template < class T > inline T norm2(T a, T b, T c) { return sqrt(a*a+b*b+c*c);}




// Gives the surface of the trace of a cube in the plane determined
// by normal n


float cubeSize(float vx, float vy, float vz, int nx, int ny, int nz)
{
	const char from[] = __FILE__ "cubeSize";

	float cubeSize;

	if (abs(nx) > 1 || abs(ny) > 1 || abs(nz) > 1)
	{
		throw std::invalid_argument("Invalid normal");
	}

	if (abs(nx)+abs(ny)+abs(nz) == 0)
	{
		throw std::invalid_argument("Invalid normal");
	}
	else if (abs(nx) + abs(ny) + abs(nz) == 1)
	{
		cubeSize = (nx?1:vx)*(ny?1:vy)*(nz?1:vz);
	}
	else if (abs(nx) + abs(ny) + abs(nz) == 2)
	{
		if (!nx)
		{
			cubeSize = vx * norm2(vy, vz);
		}
		else if (!ny)
		{
			cubeSize = vy * norm2(vx, vz);
		}
		else if (!nz)
		{
			cubeSize = vz * norm2(vx, vy);
		}
	}
	else
	{

		// Let us take a value lesser than heron

		cubeSize = max(vx * norm2(vy, vz), vy * norm2(vx, vz), vz * norm2(vx, vy));

		/*
		float l1, l2, l3;

		l1 = norm2(vx, vy);
		l2 = norm2(vx, vz);
		l3 = norm2(vy, vz);

		cubeSize = 2 * heron(l1, l2, l3);
		*/
	}

	return cubeSize;
}



template < class T >
T modifyAreaAccordingToVoxelSize
(T area,
 float vx, float vy, float vz,
 int nx, int ny, int nz)
{
	const char from[] = __FILE__ "modifyAreasAccordingToVoxelSize";

	float baseCubeSize = vx*vy;

	float cubeSize;

	if (abs(nx) > 1 || abs(ny) > 1 || abs(nz) > 1)
	{
		throw std::invalid_argument("Invalid normal");
	}

	cubeSize = cubeSize(vx, vy, vz, nx, ny, nz);

	return area * baseCubeSize / cubeSize;
}



template < class T >
T modifyAreaAccordingToVoxelSize
(T area,
 Vector<float, 3> voxSize,
 int nx, int ny, int nz)
{
	return modifyAreaAccordingToVoxelSize(area, EXPAND_VECTOR(voxSize), nx, ny, nz);
}



template < class TImage >
class MyGhost
{

public: // typenames

	typedef typename TImage::value_type T;


public:

	MyGhost(T value1, T value2, T value3, int maxSize)
		: t1(value1), t2(value2), t3(value3)
	{
		m_maxSize = maxSize;
		m_size = 0;
	}

public:

	void update(const Vector<int,3> &pos) { ++m_size; }

	bool test(const Vector<int,3> &pos) const
	{
		return t1.compute(pos) || t2.compute(pos) || t3.compute(pos);
	}
	
	bool stop() const { return (m_size >= m_maxSize); }


private:

	PT_IsEqual<TImage> t1;
	PT_IsEqual<TImage> t2;
	PT_IsEqual<TImage> t3;

	int m_size;
	int m_maxSize;
};




int keepConnectedComponentWithinSizeRangeFromBorderNoBorder
(
 Image<uchar> &seg,
 int nx, int ny, int nz,
 uchar foreground, uchar background, uchar border,
 int sizemax,
 int sizemin,
 int connectivity)
{
	uchar colorFound = findValueOtherThan(foreground, background, border);
	uchar colorDiscarded = findValueOtherThan(foreground, background, border, colorFound);
	uchar color = findValueOtherThan(foreground, background, border, colorFound, colorDiscarded);

	int count = 0, count2 = 0;

	Neighborhood nh;
	create2DNeighborhood(nh, nx, ny, nz, connectivity);

	int ccSize;
		
	Image<uchar>::VolumetricIterator iSeg(seg);

	for (; !iSeg.isAtEnd(); ++iSeg)
	{
		if (*iSeg == border)
		{
			{
				MyGhost<Image<uchar> > ghost(border, colorFound, colorDiscarded, sizemax);
				Ptr<PointList<int> > seed;
				seed->push_back(iSeg.pos());
				ccSize = regionGrowing(seg, seed, nh, ghost, color);
			}
			
			if ((ccSize <= sizemax) && (ccSize >= sizemin))
			{
				++count;
				{
					PT_IsEqual<Image<uchar> > test(color);
					SimpleGhost<PT_IsEqual<Image<uchar> > > ghost(test);
					Ptr<PointList<int> > seed = new PointList<int>(1);
					(*seed)[0] = iSeg.pos();
					regionGrowing(seg, seed, nh, ghost, colorFound);
				}
			}
			else
			{
				++count2;
				{
					PT_IsEqual<Image<uchar> > test(color);
					SimpleGhost<PT_IsEqual<Image<uchar> > > ghost(test);
					Ptr<PointList<int> > seed = new PointList<int>(1);
					(*seed)[0] = iSeg.pos();
					regionGrowing(seg, seed, nh, ghost, colorFound);
				}

			}
		}
	}
		

	// We return a result following the same foreground/background
	// convention as the input

	{
		PT_IsEqual<Image<uchar> > pt(colorFound);
		PA_SetValue<Image<uchar> > acIf(foreground);
		PA_SetValue<Image<uchar> > acElse(background);
		if_then_else(seg, pt, acIf, acElse);
	}

	//IFTHENELSE(seg[n] == colorFound, seg[n] = foreground, seg[n] = background);

	return count;
}



int removeNonCircularComponentsNoBorder
(
 Ptr<Image<uchar> > &seg,								
 uchar foreground, uchar background,								
 double distMax, const Neighborhood &nh)
{
    uchar color;
	for (color = 0; color==foreground || color == background; ++color);

	double dist;

	int count = 0;

	Image<uchar>::VolumetricIterator iSeg(seg);

	for (; !iSeg.isAtEnd(); ++iSeg)
	{
		if (*iSeg == foreground)
		{
			PT_IsEqual<Image<uchar> > test(foreground);
			RegionMoments plugin;

			PluginGhost<PT_IsEqual<Image<uchar> >, RegionMoments> ghost(test, plugin);

			Ptr<PointList<int> > seed = new PointList<int>(1);
			(*seed)[0] = iSeg.pos();

			regionGrowing(seg, seed, nh, ghost, color);

			SymMatrix3<double> moments;
			plugin.getMoments(moments);

			double e0, e1, e2;
			eigen3D(moments, e0, e1, e2);

			sort(e0, e1, e2);

			dist = 1-max(0.0,e1)/e0;
			
			if (dist > distMax)
			{
				++count;
				PT_IsEqual<Image<uchar> > test(foreground);
				SimpleGhost<PT_IsEqual<Image<uchar> > > ghost(test);
				Ptr<PointList<int> > seed = new PointList<int>(1);
				(*seed)[0] = iSeg.pos();
				regionGrowing(seg, seed, nh, ghost, color);
			}
		}
	}

	{
		PT_IsEqual<Image<uchar> > test(color);
		PA_SetValue<Image<uchar> > action(foreground);
		if_then(seg, test, action);
	}

	return count;
}




int removeNonCircularComponentsNoBorder(Ptr<Image<uchar> > &seg,
								uchar foreground, uchar background,
								double distMax,
								int nx, int ny, int nz, int connectivity)
{
	Neighborhood nh;
	create2DNeighborhood(nh, nx, ny, nz, connectivity);

	return removeNonCircularComponentsNoBorder(seg, foreground, background, distMax, nh);
}


						 
void detectDiscsInVolumeCutsFromBorder
(Ptr<Image<uchar> > &seg,
 Ptr<Image<uchar> > &output,
 uchar foreground, uchar background,
 int sizeMax, int sizeMin,
 float minCircular, bool useAnisotropy, bool newBorder)
{

	int noBorderTest = 1;

	Image<uchar> temp(param(seg));

	uchar border;

	
	// Give the object border pixel a special intensity

	border = findValueOtherThan(foreground, background);		
	if (!newBorder)
	{
		TIL_PT_HasNNeighbors<Image<uchar> > pt(N6, 3);
		PA_SetValue<Image<uchar> > ac(border);
		if_then(seg, pt, ac);
	}


	// Set the borders as background so region growing will not grow into it

	{
		PT_IsOnImageBorder<Image<uchar> > pt;
		PA_SetValue<Image<uchar> > ac(background);
		if_then(seg, pt, ac);
	}

	int norm;
	int nx, ny, nz;
	int sizeMaxModified;
	int sizeMinModified;


	// For all planes

	for (nx=0; nx<=1; ++nx)
		for (ny=(nx?-1:0); ny<=1; ++ny)
			for (nz=(nx||ny?-1:0); nz<=1; ++nz)
			{
				if ((norm = abs(nx)+abs(ny)+abs(nz)) != 0)
				{
					std::cout << '(' __ nx __ ',' __ ny __ ',' __ nz __ ')' << std::endl;

					copy(seg, temp);
					
					if (useAnisotropy)
					{
						sizeMinModified = (int) floor(modifyAreaAccordingToVoxelSize<float>(sizeMin, seg.vdim(), nx, ny, nz));
						sizeMaxModified = (int) ceil(modifyAreaAccordingToVoxelSize<float>(sizeMax, seg.vdim(), nx, ny, nz));
					}
					else
					{
						sizeMinModified = sizeMin;
						sizeMaxModified = sizeMax;
					}

					keepConnectedComponentWithinSizeRangeFromBorderNoBorder(temp,
						nx, ny, nz, foreground, background, border, sizeMaxModified, sizeMinModified, (norm == 3? 6: 8));
					
					if (minCircular > 0)
					{
						removeNonCircularComponentsNoBorder(temp,foreground,background,minCircular,nx,ny,nz,(norm == 3? 6: 8));
					}

					// Add result to output image

					{
						PT_IsNotEqual<Image<uchar> > test(background);
						PA_CopyImage<Image<uchar> > action(temp);
					}
				}
			}

	// We clear the borders to return the segmentation unchanged

	if (!newBorder)
	{
		PT_IsEqual<Image<uchar> > test(border);
		PA_SetValue<Image<uchar> > action(foreground);
		if_then(seg, test, action);
	}
}



void segmentColonNewCC
(					   
 Ptr<Image<ushort> > &im,
 ushort foreground,
 ushort background)
{

	// A color different from the background color
	ushort newColor = background+1;
	// A color different from all other colors
	ushort colorKeep = max(maxImPtr(im) + 1, newColor + 1);


	// Fill background regions that are connected to image X-Y borders


	PT_IsEqual<Image<ushort> > pt;
	PA_SetValue<Image<ushort> > ac(background);

	Image<ushort>::VolumetricIterator iIm(im);

	for (; !iIm.isAtEnd(); ++iIm)
	{
		if (
			(iIm.pos()[0] == 0) || (iIm.pos()[0] == im.dim()[0] - 1) ||
			(iIm.pos()[1] == 0) || (iIm.pos()[1] == im.dim()[1] - 1) )
		{
			if (*iIm != background)
			{
				pt.setValue(*iIm);
				if_then(im, im, pt, ac);
			}
		}
	}


	// Fill background regions that are connected to top border, but whose
	// z-dimension is not too large.
	// NB: the fact that this comes after the previous filling is important
	// because actually we want to do this test for regions that are inside
	// the patient. Regions outside the patient have already been filled by the
	// previous loop.


	//writeCTImage16("c:\\temp.ct", seg, xs, ys, zs, 1, 1, 1);
	int count = 0;
	ushort brot;


	const int Z_DEEP = int(im.getZ() * 0.38);

	Image<ushort>::VolumetricIterator iImSlice(im, Range<int,3>(0,0,Z_DEEP,im.getX()-1, im.getY()-1, Z_DEEP));


	for (; !iIm.isAtEnd(); ++iIm)
	{
		if (iIm.getZ() == 1)
		{
			brot = *iIm;
			if ((brot != background) && (brot != colorKeep))
			{
				bool flagDeep = false;
				{

					iImSlice.init(im, Range<int,3>(0, 0, Z_DEEP, iIm.getX()-1, iIm.getY()-1, Z_DEEP));
					
					for ( ; !iImSlice.isAtEnd(); ++iImSlice)
					{
						if (*iIm == brot) 
						{
							flagDeep = true;
							goto endloop;
						}
						
					}
				}
endloop:
				if (flagDeep)
				{
					pt.setValue(brot);
					ac.setValue(colorKeep);
					if_then(im, im, pt, ac);
				}
				else
				{
					pt.setValue(brot);
					ac.setValue(background);
					if_then(im, im, pt, ac);
				}
			}
		}
	}


	// Replace all non background pixels with foreground pixels
	{
		PT_IsNotEqual<Image<ushort> > pt(background);
		PA_SetValue<Image<ushort> > ac(foreground);
		if_then(im, pt, ac);
	}
}


void testColon(int argc, char* argv[])
{
	typedef Image<ushort> TInputImage;
	typedef Image<uchar> TSegImage;
	typedef Image<uchar> TOutImage;
	
	char * inputName = argv[1];
	char * outputName = argv[2];
	int minSize = atoi(argv[3]);
	float maxCirc = (float) atof(argv[4]);
	int maxSize = atoi(argv[5]);
	char *textname = argv[6];
	char *whichAlgo = argv[7];
	char *whichBorder = argv[8];
	ushort segThreshold = (ushort) atoi(argv[9]);

	bool useAniso = (strcmp(whichAlgo, "aniso") == 0);
	bool newBorder = (strcmp(whichBorder, "newBorder") == 0);


	Ptr<TInputImage> im;
	readCT(argv[1], im);

	// Threshold
	{
		PT_IsAbove<TInputImage> pt(segThreshold);
		PA_SetValue<TInputImage> acIf(1);
		PA_SetValue<TInputImage> acElse(0);
		if_then_else(im, im, pt, acIf, acElse);
	}


	// Remove isolated pixels
	{
		PT_IsIsolated<TInputImage> pt(N6);
		PA_SetValue<TInputImage> ac(0);
		if_then(im, im, pt, ac);
	}

	// Connected components

	connectedComponents(im, N26);

	// Remove connected components
	segmentColonNewCC(im, 1, 0);


	// Put result into an uchar image
	Ptr<TSegImage> seg = TSegImage::New(param(im));
	{
		PT_IsEqual<TInputImage> pt(0);
		PA_SetValue<TSegImage> acIf(1);
		PA_SetValue<TSegImage> acElse(0);
		if_then_else(im, seg, pt, acIf, acElse);
	}


	Ptr<TOutImage> out = TOutImage::New(param(im));
	detectDiscsInVolumeCutsFromBorder(seg, out, 1, 0, maxSize, minSize, maxCirc, useAniso, newBorder);
}


int main(int argc, char *argv[])
{
    std::cout << "START MAIN DLL" << std::endl;

	char *inputname = argv[1];
	char *outputname = argv[2];
	ushort segThreshold = (ushort) atoi(argv[3]);
	cout << "Threshold : " << segThreshold << endl;
	
	Ptr<Image<ushort> > seg;
	std::cout << "Reading image..." << std::flush;
	readCT(inputname, seg);
	std::cout << "OK" << std::endl;

	std::cout << "Thresholding..." << std::flush;
	Image<ushort>::LinearIterator iSeg(seg);
	for (; !iSeg.isAtEnd(); ++iSeg)
	{
		if (*iSeg >= segThreshold)
		{
			*iSeg = 0;
		}
		else
		{
			*iSeg = 1;
		}
	}
	std::cout << "OK" << std::endl;

	Neighborhood nh(N6);

	std::cout << "Connected components..." << std::flush;
	connectedComponents(seg, nh);
	std::cout << "OK" << std::endl;

	std::cout << "Writing image..." << std::flush;
	writeCT("concomp.ct", seg);
	std::cout << "OK" << std::endl;

    return 0;
}

#else
#if 0



template < typename T >
void removeInteriorPointsSlow(T * seg, int xs, int ys, int zs, T background)
{
	int i, j, k;
	int xy = xs*ys;
	
	for (k=0; k<zs; ++k)
	for (j=0; j<ys; ++j)
	for (i=0; i<xs; ++i)
	{
		if ((*seg != background) &&
			(!((i-1>=0)&&(i-1<=xs-1)&&(j>=0)&&(j<=ys-1)&&(k>=0)&&(k<=zs-1)) || *(seg-1)  != background) &&
			(!((i+1>=0)&&(i+1<=xs-1)&&(j>=0)&&(j<=ys-1)&&(k>=0)&&(k<=zs-1)) || *(seg+1)  != background) &&
			(!((i>=0)&&(i<=xs-1)&&(j-1>=0)&&(j-1<=ys-1)&&(k>=0)&&(k<=zs-1)) || *(seg-xs) != background) &&
			(!((i>=0)&&(i<=xs-1)&&(j+1>=0)&&(j+1<=ys-1)&&(k>=0)&&(k<=zs-1)) || *(seg+xs) != background) &&
			(!((i>=0)&&(i<=xs-1)&&(j>=0)&&(j<=ys-1)&&(k-1>=0)&&(k-1<=zs-1)) || *(seg-xy) != background) &&
			(!((i>=0)&&(i<=xs-1)&&(j>=0)&&(j<=ys-1)&&(k+1>=0)&&(k+1<=zs-1)) || *(seg+xy) != background))
		{
			*seg = background;
		}			
		++seg;
	}
}


template < typename T >
void removeInteriorPoints(T * seg, int xs, int ys, int zs, T background)
{
	int i, j, k;
	int xy = xs*ys;
	
	for (k=0; k<zs; ++k)
	for (j=0; j<ys; ++j)
	for (i=0; i<xs; ++i)
	{
		if ((*seg != background) &&
			(i == 0		|| *(seg-1)  != background) &&
			(i == xs-1	|| *(seg+1)  != background) &&
			(j == 0		|| *(seg-xs) != background) &&
			(j == ys-1	|| *(seg+xs) != background) &&
			(k == 0		|| *(seg-xy) != background) &&
			(k == zs-1	|| *(seg+xy) != background))
		{
			*seg = background;
		}			
		++seg;
	}
}

template < typename T >
void removeInteriorPoints_inInterior(T * seg2, int xs, int ys, int zs, T background)
{
	ImageC<T> im;
	int i, j, k;
	int xy = xs*ys;
	T* seg;
	for (k=1; k<zs-1; ++k)
	for (j=1; j<ys-1; ++j)
	for (i=1; i<xs-1; ++i)
	{
		seg = seg2 + i + xs * (j + ys*k);

		if ((*seg != background) &&
			(*(seg-1)  != background) &&
			(*(seg+1)  != background) &&
			(*(seg-xs) != background) &&
			(*(seg+xs) != background) &&
			(*(seg-xy) != background) &&
			(*(seg+xy) != background))
		{
			*seg = background;
		}
		++seg;
	}
}

template < typename T >
void next(typename ImageC<T>::VolumetricIterator &iIm)
{
	iIm.m_pos[0] = iIm.m_roi.min_bounds()[0];
	
	if (++iIm.m_pos[1] > iIm.m_roi.max_bounds()[1])
	{
		iIm.m_pos[1] = iIm.m_roi.min_bounds()[1];
		
		if (++iIm.m_pos[2] > iIm.m_roi.max_bounds()[2])
		{
			iIm.m_index = 0;
			return;
		}
	}
	
	// TODO: to be faster, use offsets to jump in y and z directions
	iIm.m_index = iIm.m_im.getUnsafePointerOf(iIm.m_pos);
}

template < typename T >
void removeInteriorPointsIterator1(T * seg, int xs, int ys, int zs, T background)
{
	int xy = xs*ys;
	
	writeCTImage16("c:\\temp\\temp1.ct", seg, xs, ys, zs, 1, 1, 1);
	
	Ptr<ImageC<T> > im = ImageC<T>::New(seg, xs, ys, zs, 1, 1, 1);
	typename ImageC<T>::VolumetricIterator iIm(im);

	int &i = (iIm.m_pos.m_values[0]);
	int &j = (iIm.m_pos.m_values[1]);
	int &k = (iIm.m_pos.m_values[2]);

	//int xs = iIm.m_roi.m_posMax.m_values[0];

	
	//for (; iIm.m_index != 0; ++iIm)
	for(;;)
	{

		/*
		if ((*seg != background) &&
			(i == 0		|| *(seg-1)  != background) &&
			(i == xs-1	|| *(seg+1)  != background) &&
			(j == 0		|| *(seg-xs) != background) &&
			(j == ys-1	|| *(seg+xs) != background) &&
			(k == 0		|| *(seg-xy) != background) &&
			(k == zs-1	|| *(seg+xy) != background))
		{
			//*seg = background;
		}
		*/
		++seg;
		

		if (++i >= xs)
		{
			next(iIm);
			if (iIm.m_index == 0) break;
		}
		
		
		// We're still in the volume
		
		else
		{
			++iIm.m_index;
		}
	}
	writeCTImage16("c:\\temp\\temp.ct", seg, xs, ys, zs, 1, 1, 1);
}

#define PRINT_TIME(comment, line)	\
{									\
	clock_t start, finish;			\
	start = clock();				\
	line;							\
	finish = clock();				\
	std::cout << comment << " : " << (double)(finish - start) << std::endl;	\
}									\


void testPerf()
{
	typedef unsigned short value_type;
	typedef ImageC<value_type> Image;
	int size = 300;

	{
		Ptr<Image> im = Image::New(size, size, size, 1, 1, 1);
		add(im, 1);
		PT_IsInterior<Image> test(N6);
		PA_SetValue<Image> action(0);
		PRINT_TIME("C++", if_then(im, test, action));
		//writeCT("c:\\temp\\interior1.ct", im);
	}

	{
		Ptr<Image> im = Image::New(size, size, size, 1, 1, 1);
		add(im, 1);
		PT_IsInteriorTN_nocheck<Image,TN6> test(N6);
		PA_SetValue<Image> action(0);
		PRINT_TIME("Template Neighbor, no check", if_then_inInterior(im, test, action));
		//writeCT("c:\\temp\\interior8.ct", im);
	}

	{
		Ptr<Image> im = Image::New(size, size, size, 1, 1, 1);
		add(im, 1);
		PT_IsInterior_nocheck<Image> test(N6);
		PA_SetValue<Image> action(0);
		PRINT_TIME("C++, No check", if_then_inInterior(im, test, action));
		//writeCT("c:\\temp\\interior7.ct", im);
	}

	{
		Ptr<Image> im = Image::New(size, size, size, 1, 1, 1);
		add(im, 1);
		PT_IsInteriorTN<Image, TN6> test(N6);
		PA_SetValue<Image> action(0);
		PRINT_TIME("Template Neigbor", if_then(im, test, action));
		//writeCT("c:\\temp\\interior6.ct", im);
	}

	{
		Ptr<Image> im = Image::New(size, size, size, 1, 1, 1);
		add(im, 1);
		PT_IsInteriorN6<Image> test(N6);
		PA_SetValue<Image> action(0);
		PRINT_TIME("C++ N6", if_then(im, test, action));
		//writeCT("c:\\temp\\interior4.ct", im);
	}

	{
		Ptr<Image> im = Image::New(size, size, size, 1, 1, 1);
		add(im, 1);
		PRINT_TIME("C, no check", removeInteriorPoints_inInterior<value_type>(im.getPointer(), im.dim()[0], im.dim()[1], im.dim()[2], 0));
		//writeCT("c:\\temp\\interior9.ct", im);
	}

	{
		Ptr<Image> im = Image::New(size, size, size, 1, 1, 1);
		add(im, 1);
		PRINT_TIME("C", removeInteriorPoints<value_type>(im.getPointer(), im.dim()[0], im.dim()[1], im.dim()[2], 0));
		//writeCT("c:\\temp\\interior2.ct", im);
	}
/*
	{
		Ptr<Image> im = Image::New(size, size, size, 1, 1, 1);
		add(im, 1);
		PRINT_TIME("Interior 5", removeInteriorPointsIterator1<value_type>(im.getPointer(), im.dim()[0], im.dim()[1], im.dim()[2], 0));
		writeCT("c:\\temp\\interior5.ct", im);
	}
*/
	{
		Ptr<Image> im = Image::New(size, size, size, 1, 1, 1);
		add(im, 1);
		PRINT_TIME("C, extensive check", removeInteriorPointsSlow<value_type>(im.getPointer(), im.dim()[0], im.dim()[1], im.dim()[2], 0));
		//writeCT("c:\\temp\\interior3.ct", im);
	}
}

int main(int argc, char *argv[])
{
    std::cout << "START MAIN DLL 2" << std::endl;
		testPerf();
		return 0;
	return 0;
}

#endif

#include <windows.h>
BOOL APIENTRY DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
					 )
{
    switch (ul_reason_for_call)
	{
		case DLL_PROCESS_ATTACH:
			std::cout << "TIL library loaded" << std::endl;
			break;
		case DLL_THREAD_ATTACH:
			std::cout << "TIL library loaded" << std::endl;
			break;
		case DLL_THREAD_DETACH:
			std::cout << "TIL library dismissed" << std::endl;
			break;
		case DLL_PROCESS_DETACH:
			std::cout << "TIL library dismissed" << std::endl;
			break;
    }
    return TRUE;
}


#endif


