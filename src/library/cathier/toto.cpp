
// includes from STL
//#include <math.h>
#include <cmath>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <numeric>
#ifndef _MSC_VER
#  include <sys/times.h>
#  include <sys/time.h>
#endif
#include <vector>

// includes from BOOST
//#include <boost/utility/result_of.hpp>
#include <boost/random/uniform_real.hpp>

// includes from MTL
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
// do not include MTL for GCC 3.4.4 as it does not compile
#ifndef _MSC_VER
#if GCC_VERSION < 030404
#include <mtl/matrix.h>
#include <mtl/mtl.h>
#include <mtl/utils.h>
#endif
#endif

// includes from CARTO
//#include "carto/Converter.h"

// includes from AIMS
//#ifdef USE_AIMS
#include "connectomist/fibertracking/bundles.h"
#include "aims/getopt/getopt2.h"
#include "aims/io/reader.h"
#include "aims/io/writer.h"
#include "aims/mesh/inflate.h"
#include "aims/mesh/surface.h"
#include "aims/mesh/surfacegen.h"
#include "aims/mesh/texture.h"
#include "aims/resampling/linearInterpolator.h"
#include "aims/roi/maskIterator.h"
#include "aims/utility/converter_volume.h"
#include "aims/vector/vector.h"
//#endif

// includes from TIL
#include "til/til_common.h"
#include "til/TExpr.h"


// includes from TIL
#include "aims_wrap.h"
//#include "BinaryTree.h"
#include "binary_tree.h"
#if GCC_VERSION < 030404
#include "connectivityMatrix.h"
#endif
#include "containers.h"
#include "Forces.h"
#include "func_iterator.h"
#include "globalTraits.h"
#include "histogram.h"
#include "kdtree.h"
#include "MeshTraits.h"
#include "meshUtils.h"
#include "minTools.h"
#include "miscUtils.h"
//#include "TriangleMeshGeodesicMap.h"
#include "triangle_mesh_geodesic_map.h"
#include "SparseVector.h"
#include "Vector.h"


# define PRINT_TIME3(comment, line)                                       \
{                                                                         \
  struct timeval start, finish;                                           \
  struct timezone tz;                                                     \
  gettimeofday(&start, &tz);                                              \
  line;                                                                   \
  gettimeofday(&finish, &tz);                                             \
  std::cout << comment << " : " << (double)((finish.tv_sec-start.tv_sec) * 1000000L + (finish.tv_usec-start.tv_usec)) << std::endl; \
}

# define DEFINE_TIMER_FUNCTION(comment, line)                             \
{                                                                         \
  struct timeval start, finish;                                           \
  struct timezone tz;                                                     \
  gettimeofday(&start, &tz);                                              \
  line;                                                                   \
  gettimeofday(&finish, &tz);                                             \
  std::cout << comment << " : " << (double)((finish.tv_sec-start.tv_sec) * 1000000L + (finish.tv_usec-start.tv_usec)) << std::endl; \
}


# define PRINT_TIME2(comment, line)                                       \
{                                                                         \
  struct tms start, finish;                                               \
  times(&start);                                                          \
  line;                                                                   \
  times(&finish);                                                         \
  std::cout << comment << " : " << (double)(finish.tms_utime - start.tms_utime) << std::endl; \
}

# define PRINT_TIME(comment, line)                                        \
{                                                                         \
  std::clock_t start, finish;                                             \
  start = std::clock();                                                   \
  std::cout << "Start: " << double(start) << std::endl;                   \
  line;                                                                   \
  finish = std::clock();                                                  \
  std::cout << comment << " : " << (double)(finish - start) << std::endl; \
}

inline unsigned long long int rdtsc()
{
  unsigned long long int x;
  __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
  return x;
}

using namespace aims;
using namespace carto;
using namespace std;

template < class TMatrix >
void workaround_init_matrix(TMatrix &mat)
{
  if (mat.nrows() != mat.ncols())
    throw std::invalid_argument("Argument is not a symmetric matrix");
  for (unsigned int i = 0; i < mat.ncols(); ++i) mat(i,i) = 0;
}



#if GCC_VERSION < 030404
void testMTL()
{
	std::cout << "testMTL begin" << std::endl;
	typedef mtl::matrix<double, mtl::symmetric<mtl::lower>, mtl::compressed<>, mtl::row_major>::type CMat;
	CMat B(10000);
	B(0,0) = 1;
	B(0,2) = 2;
	B(1,1) = 3;
	B(2,1) = 4;
	B(2,2) = 4;
  B(3,4) = 9;
  B(4,3) = 1;
  B(5,4) = 9;
  B(4,5) = 1;
	std::cout << B(1,2) << " " << B(9000,9000) << " " << B(1,1) << std::endl;
  std::cout << B(3,4) << std::endl;
  std::cout << B(4,3) << std::endl;
  std::cout << B(5,4) << std::endl;
  std::cout << B(4,5) << std::endl;
}

/// Sum matrix columns.
/// Columns of mat are summed and added to vec.
/// NB: vec should already be allocated.
/// NB: vec is not reset to zero here. Do that prior to calling this function
/// if that is what you wish.
template < class TMatrix >
void sumColumns(const TMatrix &mat, mtl::dense1D<double> &vec)
{
  for (int i = 0; i < mat.nrows(); ++i)
  {
    vec[i] += mtl::sum(mat[i]);
  }
}

template < class TMatrix >
void sumColumns2(const TMatrix &mat, mtl::dense1D<double> &vec)
{
  for (int i = 0; i < mat.nrows(); ++i)
  {
    for (int j = 0; j < mat.ncols(); ++j)
    {
      vec[i] += mat(i,j);
    }
  }
}

void testMul()
{
  const int N = 1000;
  //mtl::matrix<double, mtl::symmetric<mtl::upper>, mtl::compressed<>, mtl::row_major>::type mat(N,N);
  //mtl::matrix<double, mtl::symmetric<mtl::lower>, mtl::array< mtl::compressed<> >, mtl::column_major>::type mat(N,N);
  mtl::matrix<double, mtl::symmetric<mtl::lower>, mtl::array< mtl::compressed<> >, mtl::column_major>::type mat(N);
  workaround_init_matrix(mat);
  //mtl::matrix<double, mtl::rectangle<>, mtl::compressed<>, mtl::row_major>::type  mat(N,N);
  mtl::dense1D<double> vec(N);
  mtl::dense1D<double> res(N);
  
  mtl::set_value(mat, 0.0);
  mtl::set_value(vec, 0.0);
  mtl::set_value(res, 0.0);
  
  std::cout << "Attention mesdames et messieurs" << std::endl;
  
  mtl::mult(mat, vec, res);
}


void testConnectivityMatrix(int argc, char * argv[])
{
  aims::Reader<AimsSurfaceTriangle> reader;
//  aims::Write<aims::AimsSurfaceTriangle> writer;
  AimsApplication app( argc, aims_const_hack(argv), "testConnectivityMatrix" );
  app.addOption(reader, "-i", "input mesh" );
//  app.addOption(writer, "-o", "output mesh" );
  app.initialize();
  
  AimsSurfaceTriangle mesh;
  std::cout << "Reading mesh..."<< std::flush;
  reader.read(mesh);
  std::cout << "Done" << std::endl;
  std::cout << "Read mesh with" << size(getVertices(mesh)) << " vertices" << std::endl;
  
  std::cout << "Building connectivity matrix..." << std::flush;
  BuildConnectivityMatrix::SuggestedType mat(size(getVertices(mesh)), size(getVertices(mesh)));
  buildConnectivityMatrix(mesh, mat);
  std::cout << "OK" << std::endl;
  
  int n = (size(getVertices(mesh))-1) - 5;
  std::cout << "Printing neighbors of " << n << std::endl;
  for (int i = 0; i < mat.ncols(); ++i)
  {
    if (mat(n,i)) std::cout << i << "(" << mat(n,i) << ")"__ std::flush;
  }
  std::cout << std::endl;
  std::cout << "Printing neighbors (2) of " << n << std::endl;
  for (int i = 0; i < mat.ncols(); ++i)
  {
    if (mat(i,n)) std::cout << i << "(" << mat(i,n) << ")"__ std::flush;
  }
  std::cout << std::endl;

  mtl::dense1D<double> vec(mat.nrows());
  mtl::set_value(vec, 0.0);
  sumColumns(mat, vec);
  mtl::print_vector(vec);

  {
    mtl::dense1D<double> tmp(mat.ncols());
    mtl::set_value(tmp, 1.0);
    mtl::set_value(vec, 0.0);
    mtl::mult(mat,tmp,vec);
  }
  std::cout << "OK" << std::endl;
  mtl::print_vector(vec);
}

#endif


/*
//void foo(Writer<AimsSurfaceTriangle> &r2)
void foo()
{
  AimsSurface<3,Void>  t;
  vector<Point3df>			& vert = t.vertex();
  vector< AimsVector<uint,3> >		& poly = t.polygon();
  

  poly.push_back( AimsVector<uint,3> (0,1,2) );
  //poly[0].push_back( 1 );
  //poly[0].push_back( 2 );

  vert.push_back( Point3df(1,1,1) ) ;
  vert.push_back( Point3df(-1,1,1) );
  vert.push_back( Point3df(1,-1,1) ); 

  t.updateNormals();

  / *
  AimsSurfaceTriangle t;
  //t.push_back(pair<int,AimsSurface<3,Void> >(0,AimsSurface<3,Void>()));
  t[0] = AimsSurface<3,Void>();
  {
    std::vector<Point3df> v = t.vertex();
    v.push_back(Point3df( 1, 1, 1));
    v.push_back(Point3df(-1, 1, 1));
    v.push_back(Point3df( 1,-1, 1));
  }
  // Set faces
  {
    std::vector<AimsVector<uint,3> > p = t.polygon();
    p.push_back(AimsVector<uint,3>(0,1,2));
  }
  t.updateNormals();
  * /

  AimsSurfaceTriangle surface;
  surface[0] = t;

  Writer<AimsSurfaceTriangle> wt("bouzna");
  wt << surface;
  //r2.write(t);

 }
*/


template < typename TMesh >
void _testCenttemp(const TMesh & mesh)
{
  {
    PRINT_TIME3("standard4", ({
    typename MeshTraits<TMesh>::Vertex tmp;
    typename MeshTraits<TMesh>::FaceIndexCollection::const_iterator iFaceIndex;
    for (iFaceIndex = getFaceIndices(mesh).begin(); iFaceIndex != getFaceIndices(mesh).end(); ++iFaceIndex)
    {
      tmp = getFaceVertex(mesh, iFaceIndex, 0);
    }
    }));
  }
  {
    PRINT_TIME3("standard4", ({
    typename MeshTraits<TMesh>::Vertex tmp;
    typename MeshTraits<TMesh>::FaceIndexCollection::const_iterator iFaceIndex;
    for (iFaceIndex = getFaceIndices(mesh).begin(); iFaceIndex != getFaceIndices(mesh).end(); ++iFaceIndex)
    {
      tmp = getFaceVertex(mesh, iFaceIndex, 0);
    }
    }));
  }
  {
    PRINT_TIME3("superfast?", ({
    typename MeshTraits<TMesh>::Vertex tmp;
    typename MeshTraits<TMesh>::FaceIndexCollection::const_iterator iFaceIndex = getFaceIndices(mesh).begin();
    typename MeshTraits<TMesh>::FaceIndexCollection::const_iterator iFaceIndexEnd = getFaceIndices(mesh).end();
    for (; iFaceIndex != iFaceIndexEnd; ++iFaceIndex)
    {
      tmp = getFaceVertex(mesh, iFaceIndex, 0);
    }
    }));
  }
}


void _testCentroid(AimsSurfaceTriangle * s1)
{
  //AimsSurfaceTriangle &s3 = *s1;
  std::cout << "oki" << std::endl;
  
  AimsSurface<3, Void> * s = &(*s1)[0];
  AimsSurface<3, Void> & s2 = *s;

  std::cout << "Initializing writer" << std::endl;
  Writer<AimsSurfaceTriangle> w("bobo");
  //std::cout << "Writing mesh" << std::endl;
  //w.write(*s);
  
/////////////////////////////////////////////////////////////////////////
// AIMS
/////////////////////////////////////////////////////////////////////////

  std::cout << "AIMS" << std::endl;
  {
    PRINT_TIME3("standard", ({
    Point3df tmp;
    std::vector<AimsVector<unsigned int, 3> >::iterator iFaceIndex;
    for (iFaceIndex = s->polygon().begin(); iFaceIndex != s->polygon().end(); ++iFaceIndex)
    {
      tmp = s->vertex()[(*iFaceIndex)[0]];
      //std::cout << tmp << std::endl;
    }
    }));
  }  
  

  {
    PRINT_TIME3("standard2", ({
    Point3df tmp;
    MeshTraits<AimsSurfaceTriangle>::FaceIndexCollection::const_iterator iFaceIndex;
    //for (iFaceIndex = getFaceIndices(*s).begin(); iFaceIndex != getFaceIndices(*s).end(); ++iFaceIndex)
    for (iFaceIndex = getFaceIndices(*s).begin(); iFaceIndex != s->polygon().end(); ++iFaceIndex)
    {
      tmp = getVertices(s2)[(*iFaceIndex)[0]];
      //tmp = getVertices(*s)[(*iFaceIndex)[0]];
      //tmp = s->vertex()[(*iFaceIndex)[0]];
      //std::cout << tmp << std::endl;
    }
    }));
  }  

  {
    PRINT_TIME3("standard3=1", ({
    Point3df tmp;
    std::vector<AimsVector<unsigned int, 3> >::iterator iFaceIndex;
    for (iFaceIndex = s->polygon().begin(); iFaceIndex != s->polygon().end(); ++iFaceIndex)
    {
      tmp = s->vertex()[(*iFaceIndex)[0]];
      //std::cout << tmp << std::endl;
    }
    }));
  }  

  {
    PRINT_TIME3("superfast?", ({
    Point3df tmp;
    std::vector<AimsVector<unsigned int, 3> >::iterator iFaceIndex;
    std::vector<AimsVector<unsigned int, 3> > & polygon = s->polygon();
    std::vector<Point3df> & vertex = s->vertex();
    for (iFaceIndex = polygon.begin(); iFaceIndex != polygon.end(); ++iFaceIndex)
    {
      tmp = vertex[(*iFaceIndex)[0]];
      //std::cout << tmp << std::endl;
    }
    }));
  }  

  {
    PRINT_TIME3("standard4", ({
    Point3df tmp;
    std::vector<AimsVector<unsigned int, 3> >::iterator iFaceIndex;
    for (iFaceIndex = s->polygon().begin(); iFaceIndex != s->polygon().end(); ++iFaceIndex)
    {
      tmp = getFaceVertex(s2, iFaceIndex, 0);
      //std::cout << tmp << std::endl;
    }
    }));
  }

  {
    MeshTraits<AimsSurface<3,Void> >::Vertex tmp;
    PRINT_TIME3("traits", ({    
    MeshTraits<AimsSurface<3,Void> >::FaceCollection::const_iterator iFaces = getFaces(*s).begin();
    for ( ; iFaces != getFaces(*s).end(); ++iFaces)
    {
      tmp = (*iFaces)[0];      
      //std::cout << tmp << std::endl;
    }
    }));
  }

/////////////////////////////////////////////////////////////////////////
// MESH1
/////////////////////////////////////////////////////////////////////////


  std::cout << "Mesh1" << std::endl;

  {
    Mesh1 mesh;
    convert(*s, mesh);
  
    {
      PRINT_TIME3("standard", ({
      MeshTraits<Mesh1>::Vertex tmp;
      std::vector<boost::array<std::size_t, 3> >::iterator iFaceIndex;
      for (iFaceIndex = mesh.getFaceIndices().begin(); iFaceIndex != mesh.getFaceIndices().end(); ++iFaceIndex)
      {
        tmp = mesh.getVertices()[(*iFaceIndex)[0]];
        //std::cout << tmp << std::endl;
      }
      }));
    }  
  
    {
      PRINT_TIME3("standard2", ({
      MeshTraits<Mesh1>::Vertex tmp;
      MeshTraits<Mesh1>::FaceIndexCollection::const_iterator iFaceIndex;
      //for (iFaceIndex = getFaceIndices(*s).begin(); iFaceIndex != getFaceIndices(*s).end(); ++iFaceIndex)
      for (iFaceIndex = getFaceIndices(mesh).begin(); iFaceIndex != mesh.getFaceIndices().end(); ++iFaceIndex)
      {
        tmp = getVertices(mesh)[(*iFaceIndex)[0]];
        //tmp = getVertices(*s)[(*iFaceIndex)[0]];
        //tmp = s->vertex()[(*iFaceIndex)[0]];
        //std::cout << tmp << std::endl;
      }
      }));
    }  
  
    {
      PRINT_TIME3("standard3=1", ({
      MeshTraits<Mesh1>::Vertex tmp;
      std::vector<boost::array<std::size_t, 3> >::iterator iFaceIndex;
      for (iFaceIndex = mesh.getFaceIndices().begin(); iFaceIndex != mesh.getFaceIndices().end(); ++iFaceIndex)
      {
        tmp = mesh.getVertices()[(*iFaceIndex)[0]];
        //std::cout << tmp << std::endl;
      }
      }));
    }  
  
    {
      PRINT_TIME3("superfast?", ({
      MeshTraits<Mesh1>::Vertex tmp;
      MeshTraits<Mesh1>::FaceIndexCollection::iterator iFaceIndex;
      MeshTraits<Mesh1>::FaceIndexCollection & face = getFaceIndices(mesh);
      MeshTraits<Mesh1>::VertexCollection & vertex = getVertices(mesh);
      for (iFaceIndex = face.begin(); iFaceIndex != face.end(); ++iFaceIndex)
      {
        tmp = vertex[(*iFaceIndex)[0]];
        //std::cout << tmp << std::endl;
      }
      }));
    }  
    {
      PRINT_TIME3("standard4", ({
      MeshTraits<Mesh1>::Vertex tmp;
      std::vector<boost::array<std::size_t, 3> >::iterator iFaceIndex;
      for (iFaceIndex = mesh.getFaceIndices().begin(); iFaceIndex != mesh.getFaceIndices().end(); ++iFaceIndex)
      {
        tmp = getFaceVertex(mesh, iFaceIndex, 0);
        //std::cout << tmp << std::endl;
      }
      }));
    }
  
    {
      MeshTraits<Mesh1>::Vertex tmp;
      PRINT_TIME3("traits", {    
      MeshTraits<Mesh1>::FaceCollection::const_iterator iFaces = getFaces(mesh).begin();
      for ( ; iFaces != getFaces(mesh).end(); ++iFaces)
      {
        tmp = (*iFaces)[0];
        //std::cout << tmp << std::endl;
      }
      });
    }
    
    {
      const int N = 10000;
      Point<double, N> p;
      Vector<double, N> v;
      std::fill(p.begin(), p.end(), 1.23);
      std::fill(v.begin(), v.end(), 3.45);
      PRINT_TIME3("addto", add(p, v));
      PRINT_TIME3("direct", ({
      Point<double,N>::iterator iP = p.begin();
      Vector<double,N>::const_iterator iV = v.begin();
      for (; iP != p.end(); ++iP, ++iV)
      {
        *iP += *iV;
      }
      }));
    }
  }


/////////////////////////////////////////////////////////////////////////
// MESH2
/////////////////////////////////////////////////////////////////////////

  std::cout << "Mesh2" << std::endl;

  {
    Mesh2 mesh;
    convert(*s, mesh);
  
    {
      PRINT_TIME3("standard", ({
      MeshTraits<Mesh2>::Vertex tmp;
      MeshTraits<Mesh2>::FaceIndexCollection::iterator iFaceIndex;
      for (iFaceIndex = mesh.getFaceIndices().begin(); iFaceIndex != mesh.getFaceIndices().end(); ++iFaceIndex)
      {
        tmp = *((*iFaceIndex)[0]);
        //std::cout << tmp << std::endl;
      }
      }));
    }  
  
    {
      PRINT_TIME3("standard2", ({
      MeshTraits<Mesh2>::Vertex tmp;
      MeshTraits<Mesh2>::FaceIndexCollection::const_iterator iFaceIndex;
      //for (iFaceIndex = getFaceIndices(*s).begin(); iFaceIndex != getFaceIndices(*s).end(); ++iFaceIndex)
      for (iFaceIndex = getFaceIndices(mesh).begin(); iFaceIndex != mesh.getFaceIndices().end(); ++iFaceIndex)
      {
        tmp = *((*iFaceIndex)[0]);
        //tmp = getVertices(*s)[(*iFaceIndex)[0]];
        //tmp = s->vertex()[(*iFaceIndex)[0]];
        //std::cout << tmp << std::endl;
      }
      }));
    }  
  
    {
      PRINT_TIME3("standard3=1", ({
      MeshTraits<Mesh2>::Vertex tmp;
      MeshTraits<Mesh2>::FaceIndexCollection::iterator iFaceIndex;
      for (iFaceIndex = mesh.getFaceIndices().begin(); iFaceIndex != mesh.getFaceIndices().end(); ++iFaceIndex)
      {
        tmp = *((*iFaceIndex)[0]);
        //std::cout << tmp << std::endl;
      }
      }));
    }  
  
    {
      PRINT_TIME3("superfast?", ({
      MeshTraits<Mesh2>::Vertex tmp;
      MeshTraits<Mesh2>::FaceIndexCollection::iterator iFaceIndex;
      MeshTraits<Mesh2>::FaceIndexCollection & face = getFaceIndices(mesh);
      //MeshTraits<Mesh2>::VertexCollection & vertex = getVertices(mesh);
      for (iFaceIndex = face.begin(); iFaceIndex != face.end(); ++iFaceIndex)
      {
        tmp = *((*iFaceIndex)[0]);
        //std::cout << tmp << std::endl;
      }
      }));
    }  

    {
      PRINT_TIME3("standard4", ({
      MeshTraits<Mesh2>::Vertex tmp;
      MeshTraits<Mesh2>::FaceIndexCollection::iterator iFaceIndex;
      for (iFaceIndex = mesh.getFaceIndices().begin(); iFaceIndex != mesh.getFaceIndices().end(); ++iFaceIndex)
      {
        tmp = getFaceVertex(mesh, iFaceIndex, 0);
        //std::cout << tmp << std::endl;
      }
      }));
    }
  
    {
      MeshTraits<Mesh2>::Vertex tmp;
      PRINT_TIME3("traits", {    
      MeshTraits<Mesh2>::FaceCollection::const_iterator iFaces = getFaces(mesh).begin();
      for ( ; iFaces != getFaces(mesh).end(); ++iFaces)
      {
        tmp = (*iFaces)[0];
        //std::cout << tmp << std::endl;
      }
      });
    }
    
    {
      const int N = 10000;
      Point<double, N> p;
      Vector<double, N> v;
      std::fill(p.begin(), p.end(), 1.23);
      std::fill(v.begin(), v.end(), 3.45);
      PRINT_TIME3("addto", add(p, v));
      PRINT_TIME3("direct", ({
      Point<double,N>::iterator iP = p.begin();
      Vector<double,N>::const_iterator iV = v.begin();
      for (; iP != p.end(); ++iP, ++iV)
      {
        *iP += *iV;
      }
      }));
  }
  }
}

void testCentroid()
{

  std::cout << "Starting here" << std::endl;
  AimsSurfaceTriangle *s1 = makeSphere(Point3df(1.0, 2.0, 3.0), 0.7, 7);
  _testCentroid(s1);
}

template < typename T1, typename T2, typename T3 >
void _gra(const T1 &s, const T2 &mesh, const T3 &mesh2)
{
 _testCenttemp(s);
 _testCenttemp(mesh);
 _testCenttemp(mesh2);
}

void testCentroid3()
{
  AimsSurfaceTriangle *s1 = makeSphere(Point3df(1.0, 2.0, 3.0), 0.7, 7);
  AimsSurface<3, Void> & s = (*s1)[0];
  Mesh1 mesh;
  convert(s, mesh);
  Mesh2 mesh2;
  convert(s, mesh2);
  _gra(s,mesh,mesh2);
}


void testCentroid2()
{
  Mesh1 mesh;
  {
    std::cout << "Starting here" << std::endl;
    AimsSurfaceTriangle *s = makeSphere(Point3df(1.0, 2.0, 3.0), 0.7, 0);
    std::cout << "zero" << std::endl;
    convert(*s, mesh);
  }

  std::cout << "one" << std::endl;

  MeshTraits<Mesh1>::Vertex centr;
  centroid(getVertices(mesh), centr);

  std::cout << "two" << std::endl;

  std::cout << "Centroid : " << centr << " (should be 1.0, 2.0, 3.0 )" << std::endl;
  
  Mesh_N mesh2 = addNeighborsToMesh(mesh);

  std::cout << "three" << std::endl;
  
  std::cout << getVertices(mesh2)[0] << std::endl;
  std::cout << getVertices(mesh)[0] << std::endl;
  std::cout << &getVertices(mesh2) << std::endl;
  std::cout << &getVertices(mesh) << std::endl;

  MeshTraits<Mesh_N>::NeighborIndexCollection::const_iterator iNic = getNeighborIndices(mesh2).begin();
  MeshTraits<Mesh_N>::NeighborIndex::const_iterator iNi = iNic->begin();
  std::cout << "Voici" << std::endl;
  for (; iNi != iNic->end(); ++iNi)
  {
    std::cout << *iNi << std::endl;
  }
  SpringForce<Mesh_N, Vector<float,3> > sf;
  sf.initializeLengths(mesh2);
  std::vector<Vector<float,3> > f;
  sf.getForces(mesh2, f);
  std::cout << f[0] << std::endl;
}

template < typename T >
bool foo(T t)
{
  return boost::is_same<typename til::value_type<T>::type, std::size_t>::value;
}

void meshStats(Mesh2_N & mesh)
{
  int nVertices = size(getVertices(mesh));
  std::cout << "Number of vertices : " << nVertices << std::endl;
  int nFaces = size(getFaceIndices(mesh));
  std::cout << "Number of faces: " << nFaces << std::endl;

  int nEdges;
  
  {
    std::set<array2D<Mesh2_N::FaceIndex::value_type> >
      e = getEdges(mesh);
    nEdges = size(e);
    std::cout << "Number of edges: " << nEdges << std::endl;
  }
  
  std::cout << "Euler number: " << nVertices + nFaces - nEdges << std::endl;
  
  {
    std::cout << "Histogram: number of neighbors" << std::endl;
    //std::map<std::size_t, int> h = histogram(getNeighborIndices(mesh), functor::Size<MeshTraits<Mesh2_N>::NeighborIndexCollection>());
    til::Histogram<std::size_t> h;
    foo(til::func_it<functor::Size>(getNeighborIndices(mesh).begin()));
    h.accumulate(til::func_it<functor::Size>(getNeighborIndices(mesh).begin()), 
                 til::func_it<functor::Size>(getNeighborIndices(mesh).end()));
    //std::map<std::size_t, unsigned int> h = h.get();
    //std::for_each(h.begin(), h.end(), std::cout << _1 << " " << std::endl);
    /*
    std::map<std::size_t, unsigned int>::const_iterator iH = h.begin();
    for (; iH != h.end(); ++iH)
    {
      std::cout << iH->first << " : " << iH->second << std::endl;
    }
    */
    print(h);
  }
}


void myinflate(int argc, char* argv[])
{
  AimsApplication app( argc, aims_const_hack(argv), "myinflate" );
  aims::Reader<AimsSurfaceTriangle> reader;
  aims::Writer<AimsSurfaceTriangle> writer;
  aims::Writer<Texture1d> wt;
  typedef Mesh2_N MeshType;
  MeshType mesh;
  float kSpring = 0.2;
  float kCentroid = 0.05;
  //float kLapl = 0.01;
  int niter = 10000;
  int dump = 50;

  AimsSurfaceTriangle wmesh;
  // Read and convert input mesh
  {
    //Mesh1 mesh0;
    Mesh2 mesh0;
    app.addOption(reader, "-i", "input mesh" );
    app.addOption(writer, "-o", "output mesh" );
    app.addOption(wt, "-t", "");
    app.addOption(dump, "-dump", "");
    app.addOption(kSpring, "-kspring", "");
    app.addOption(kCentroid, "-kcentroid", "");
    app.addOption(niter, "-niter", "");
    app.initialize();
    
    //AimsSurfaceTriangle *s = makeSphere(Point3df(0,0,0), 50, 4);
    //AimsSurfaceTriangle aimsmesh = *s;
    //delete s;
    AimsSurfaceTriangle aimsmesh;
    std::cout << "Reading mesh..."<< std::flush;
    reader.read(aimsmesh);
    std::cout << "Done" << std::endl;
    std::cout << "Read mesh with" << size(getVertices(aimsmesh)) << " vertices" << std::endl;
    convert(aimsmesh, mesh0);
    mesh = addNeighborsToMesh(mesh0);
    std::cout << getVertices(aimsmesh)[0] << std::endl;
  }
  
  meshStats(mesh);
  
  Texture1d tex;
  
  // write texture
  /*
  {
    tex.reserve(size(getVertices(mesh)));
    for( unsigned i=0; i<size(getVertices(mesh)); ++i )
      tex.push_back( i==0 );
    wt.write( tex );
  }
  */
  
  // initialize spring force functor with initial lengths
  SpringForce<MeshType, Vector<float,3> > sf;
  //SpringForce<Mesh_N, float> sf;
  sf.initializeLengths(mesh);
  
  // Main loop
  MeshTraits<MeshType>::Vertex center;
  //MeshTraits<Mesh_N>::Vertex center;
  std::vector<Vector<float,3> > centroidForces(size(getVertices(mesh)));
  std::vector<Vector<float,3> > totalForces(size(getVertices(mesh)));
  std::vector<Vector<float,3> > springForces(size(getVertices(mesh)));
  std::vector<Vector<float,3> > lForces(size(getVertices(mesh)));
  int c = 0;
  int itex = 0;
  for (int i = 0; i < niter; ++i)
  {
    if (i % (niter/dump) == 0)
    {
      // write mesh
      {
        std::cout << "copying" << std::cout;
        AimsSurface<3,Void> tmp;
        convert(mesh, tmp);
        wmesh[c] = tmp;
        ++c;
      }
      // write texture
      {
        Texture<float> tmp(size(getVertices(mesh)));
        for_each_neighbors(mesh, sf.getLengths(), tmp.data(), functor::SpringEnergy());
        tex[itex] = tmp;
        ++itex;
      }
    }
    std::cout << "Iteration " << i << std::endl;
    
    // compute centroid forces
    centroid(getVertices(mesh), center);
    {    
      MeshTraits<MeshType>::Vertex tmp;
      MeshTraits<MeshType>::VertexCollection::const_iterator iVc = getVertices(mesh).begin();
      std::vector<Vector<float,3> >::iterator iC = centroidForces.begin();
      for (; iVc != getVertices(mesh).end(); ++iVc, ++iC)
      {
        sub(*iVc, center, tmp);
        mul(tmp, 1.0f/norm<float>(tmp), *iC);
      }
    }

    // compute spring forces
    zeros(springForces);
    sf.getForces(mesh, springForces);
    //std::cout << "sf " << springForces[0] << std::endl;

    // compute laplacian forces
    //laplacianForces(mesh, lForces);

    /*
    std::cout << size(getNeighborIndices(mesh)) << "*" << std::flush;
    std::cout << size(getVertices(mesh)) << "*" << std::flush;
    for (int j = 0; j < size(getNeighborIndices(mesh)); ++j)
    {
      std::cout << size(getNeighborIndices(mesh)[j]) << " " << std::flush;
    }
    */
    
    /*
    
    // hand-compute spring force for first voxel
    std::vector<Vector<float,3> > f(size(getVertices(mesh)));
    MeshTraits<Mesh_N>::VertexCollection v = getVertices(mesh);
    for (int j = 0; j < size(getVertices(mesh)); ++j)
    {
      for (std::size_t k = 0; k < size(getNeighborIndices(mesh)[j]); ++k)
      {
        //std::cout << v[getNeighborIndices(mesh)[0][i]] << std::endl;
        //std::cout << "dv " << v[0] - v[getNeighborIndices(mesh)[0][i]] << std::endl;
        //std::cout << "dn " << (norm<float>(v[0] - v[getNeighborIndices(mesh)[0][i]]) - sf.getLengths()[0][i]) << std::endl;
        //std::cout << "n " << norm<float>(v[0] - v[getNeighborIndices(mesh)[0][i]]) << std::endl;
        f[j] += (v[getNeighborIndices(mesh)[j][k]] - v[j]) * 
          ((norm<float>(v[j] - v[getNeighborIndices(mesh)[j][k]]) - sf.getLengths()[j][k])
           / norm<float>(v[j] - v[getNeighborIndices(mesh)[j][k]]));

        //std::cout << "f " << f << std::endl;
      }
    }
    */
    
    //std::cout << "my " << f << std::endl;
    // apply displacement to mesh
    {
      MeshTraits<MeshType>::VertexCollection::iterator iVc = getVertices(mesh).begin();
      std::vector<Vector<float,3> >::const_iterator iC = centroidForces.begin();
      std::vector<Vector<float,3> >::const_iterator iS = springForces.begin();
      //std::vector<Vector<float,3> >::const_iterator iL = lForces.begin();
      //for (; iVc != getVertices(mesh).end(); ++iVc, ++iC, ++iS, ++iL)
      float f = exp(-i/1000.0);
      for (; iVc != getVertices(mesh).end(); ++iVc, ++iC, ++iS)
      //std::vector<Vector<float,3> >::const_iterator iF = f.begin();
      //for (; iVc != getVertices(mesh).end(); ++iVc, ++iC, ++iF)
      {
        //(*iVc)[0] += kCentroid*(*iC)[0] + kSpring*(*iF)[0];
        //(*iVc)[1] += kCentroid*(*iC)[1] + kSpring*(*iF)[1];
        //(*iVc)[2] += kCentroid*(*iC)[2] + kSpring*(*iF)[2];
        //(*iVc)[0] += kCentroid*(*iC)[0] + kSpring*(*iS)[0] + kLapl*(*iL)[0];
        //(*iVc)[1] += kCentroid*(*iC)[1] + kSpring*(*iS)[1] + kLapl*(*iL)[1];
        //(*iVc)[2] += kCentroid*(*iC)[2] + kSpring*(*iS)[2] + kLapl*(*iL)[2];
/*        (*iVc)[0] += f * (kCentroid*(*iC)[0] + kSpring*(*iS)[0]);
        (*iVc)[1] += f * (kCentroid*(*iC)[1] + kSpring*(*iS)[1]);
        (*iVc)[2] += f * (kCentroid*(*iC)[2] + kSpring*(*iS)[2]);*/
        (*iVc)[0] += f * kCentroid*(*iC)[0] + kSpring*(*iS)[0];
        (*iVc)[1] += f * kCentroid*(*iC)[1] + kSpring*(*iS)[1];
        (*iVc)[2] += f * kCentroid*(*iC)[2] + kSpring*(*iS)[2];
      }
    }    
  }

  std::cout << "WMesh size : " << wmesh.size() << std::endl;
  // write result
  writer.write(wmesh);
  wt.write( tex );
}

void testPush()
{
  Mesh1 mesh;
  {
    std::cout << "Starting here" << std::endl;
    AimsSurfaceTriangle *s = makeSphere(Point3df(1.0, 2.0, 3.0), 0.7, 7);
    convert(*s, mesh);
  }

  MeshTraits<Mesh1>::Vertex centr;
  centroid(getVertices(mesh), centr);
  
}



/// A simple C array
template <int N>
class A
{
public:
  A() { for (int i = 0; i < N; ++i) m_data[i] = std::sqrt(double(i)); }
  double operator[](int i) const { return m_data[i]; }

private:

  double m_data[N];
};

struct B
{
  B(int n) : vA(n) {}
  std::vector<A<5> > vA;
};

const std::vector<A<5> > &
getVA(const B & b)
{
  return b.vA;
}

void testSpeed()
{
  const int n = 1000000;
  B *b = new B(n);
  {
    double d;
    PRINT_TIME3("first", ({
      for (int i = 0; i < n; ++i)
      {
        d = b->vA[i][3];
      }
    }));
  }

  {
    double d;
    PRINT_TIME3("second", ({
      for (int i = 0; i < n; ++i)
      {
        d = getVA(*b)[i][3];
      }
    }));
  }
  {
    double d;
    PRINT_TIME3("first", ({
      for (int i = 0; i < n; ++i)
      {
        d = b->vA[i][3];
      }
    }));
  }
};


void testConvert(int argc, char* argv[])
{
  AimsApplication app( argc, aims_const_hack(argv), "myinflate" );
  aims::Reader<AimsSurfaceTriangle> reader;
  aims::Writer<AimsSurfaceTriangle> writer;
  aims::Writer<Texture1d> wt;
  Mesh2_N mesh;
  //Mesh_N mesh;

  AimsSurfaceTriangle wmesh;
  // Read and convert input mesh
  {
    //Mesh1 mesh0;
    Mesh2 mesh0;
    app.addOption(reader, "-i", "input mesh" );
    app.addOption(writer, "-o", "output mesh" );
    app.initialize();
    
    //AimsSurfaceTriangle *s = makeSphere(Point3df(0,0,0), 50, 4);
    //AimsSurfaceTriangle aimsmesh = *s;
    //delete s;
    AimsSurfaceTriangle aimsmesh;
    std::cout << "Reading mesh..."<< std::flush;
    reader.read(aimsmesh);
    std::cout << "Done" << std::endl;
    std::cout << "Read mesh with" << size(getVertices(aimsmesh)) << " vertices" << std::endl;    
    
    std::cout << "First face" << std::endl;
    for (int i = 0; i < 3; ++i)
    {
      std::cout << getFaceIndices(aimsmesh)[0][i] << " " << getVertices(aimsmesh)[getFaceIndices(aimsmesh)[0][i]] << " " << &(getVertices(aimsmesh)[getFaceIndices(aimsmesh)[0][i]]) << std::endl;
    }
    std::cout << "size " << sizeof(MeshTraits<AimsSurfaceTriangle>::Vertex) << std::endl;
    convert(aimsmesh, mesh0);
    
    std::cout << "First face" << std::endl;
    for (int i = 0; i < 3; ++i)
    {
      std::cout << *(getFaceIndices(mesh0)[0][i]) << " " << (getFaceIndices(mesh0)[0][i]) << std::endl;
    }    
    mesh = addNeighborsToMesh(mesh0);
    std::cout << "size " << sizeof(MeshTraits<Mesh2_N>::Vertex) << std::endl;

    std::cout <<(getFaceIndices(mesh)[0][1] - &(getVertices(mesh)[0])) << std::endl;
  }
  std::cout << "WMesh size : " << wmesh.size() << std::endl;
  // write result

  AimsSurface<3,Void> tmp;
  convert(mesh, tmp);
  wmesh[0] = tmp;

  std::cout << "First face" << std::endl;
  for (int i = 0; i < 3; ++i)
  {
    std::cout << getFaceIndices(wmesh)[0][i] << " " << getVertices(wmesh)[getFaceIndices(wmesh)[0][i]] << " " << &(getVertices(wmesh)[getFaceIndices(wmesh)[0][i]]) << std::endl;
  }
  std::cout << "size " << sizeof(MeshTraits<AimsSurface<3,Void> >::Vertex) << std::endl;
  
  writer.write(wmesh);  
}

void testMyMesh(int argc, char* argv[])
{
  //std::cout << Is_BoostArray_N<boost::array<std::size_t, 3>, 3>::value << std::endl;


  Reader<AimsSurfaceTriangle> r;
  Writer<AimsSurfaceTriangle> w;
  AimsApplication app( argc, aims_const_hack(argv), "testMyMesh" );
  app.addOption(r, "-i", "input mesh" );
  app.addOption(w, "-o", "output mesh");
  app.initialize();

  Mesh1 mesh;
  std::cout << "Aims -> mesh" << std::endl;
  {
    AimsSurfaceTriangle s;
    r.read( s );
    convert(s, mesh);
  }
  std::cout << "Mesh -> aims" << std::endl;
  {
    AimsSurfaceTriangle s;
    convert(mesh, s);
    w.write(s);
  }  
}


void testKDTreeSerious(int argc, char * argv[])
{
  Reader<AimsSurfaceTriangle> r;
  AimsApplication app( argc, aims_const_hack(argv), "testMyMesh" );
  app.addOption(r, "-i", "input mesh" );
  app.initialize();

  Mesh2 mesh;
  std::cout << "Aims -> mesh" << std::endl;
  {
    AimsSurfaceTriangle s;
    r.read( s );
    convert(s, mesh);
  }

  //boost::uniform_real<float> random(-150,150);    
  Point3D<double> center;
  centroid(getVertices(mesh), center);
  std::cout << "Centroid: " << center << std::endl;
  Point3D<double> std;
  stdev(getVertices(mesh), std);
  std::cout << "Stdev: " << std << std::endl;
  //std::vector<Point3D<float>*> res;
  typedef KDTree<Point3D<float>*, MeshTraits<Mesh2>::VertexCollection> MyKDTree;
  MyKDTree res;
  makeKDTree(getVertices(mesh), res);
  std::cout << "kdtree.size " << size(res) << std::endl;

  struct timeval start, finish;
  struct timezone tz;
  //double timetake;
  const int N = 50;
  std::vector<double> timeTaken_kd(N);
  std::vector<double> niter_kd(N);
  std::vector<double> timeTaken_std(N);
  const int imsize = 100;
  AimsData<float> im(imsize, imsize, imsize);
  //for (int i = 0; i < N; ++i)
  for (int i = 0; i < imsize; ++i)
  {
    std::cout << i << " / " << imsize << std::endl;
  for (int j = 0; j < imsize; ++j)
  for (int k = 0; k < imsize; ++k)
  {
    Point3D<float> target = Point3D<float>
    (
/*     ((std::rand() / double(RAND_MAX)) * 2 - 1) * std[0] + center[0],
     ((std::rand() / double(RAND_MAX)) * 2 - 1) * std[1] + center[1],
     ((std::rand() / double(RAND_MAX)) * 2 - 1) * std[2] + center[2]*/
     3*(2.0*i-imsize) / imsize * std[0] + center[0],
     3*(2.0*j-imsize) / imsize * std[1] + center[1],
     3*(2.0*k-imsize) / imsize * std[2] + center[2]
    );
    //std::cout << "Target: " << target << std::endl;

    gettimeofday(&start, &tz);
    FindClosest<double, MyKDTree> fc(res);
    fc(target);
    gettimeofday(&finish, &tz);
    //im(i,j,k) = (double)((finish.tv_sec-start.tv_sec) * 1000000L + (finish.tv_usec-start.tv_usec));
    im(i,j,k) = fc.niter();
    //timeTaken_kd[i] = (double)((finish.tv_sec-start.tv_sec) * 1000000L + (finish.tv_usec-start.tv_usec));
    //niter_kd[i] = fc.niter();
  
    /*
    Point3D<float> const * p;
    gettimeofday(&start, &tz);
    {
      std::vector<Point3D<float> >::const_iterator iVertex = getVertices(mesh).begin();
      double minDist = std::numeric_limits<double>::max();
      for (; iVertex != getVertices(mesh).end(); ++iVertex)
      {
        double tmp = dist2<double>(*iVertex, target);
        if (minDist > tmp)
        {
          minDist = tmp;
          p = &*iVertex;
          //std::cout << *iVertex << " " << target << std::endl;
          //std::cout << minDist << " " << *p << std::endl;
        }
      }
    }
    gettimeofday(&finish, &tz);
    //timeTaken_std[i] = (double)((finish.tv_sec-start.tv_sec) * 1000000L + (finish.tv_usec-start.tv_usec));
    */
    
    //std::cout << im(i,j,k) << std::endl;
        
    //std::cout << "kdtree: " << timeTaken_kd[i] << " (" << niter_kd[i] << ")   std: " << timeTaken_std[i] << std::endl;
  }
  }
  /*
  double mean_kd = std::accumulate(timeTaken_kd.begin(), timeTaken_kd.end(), 0.0)/N;
  double mean_std = std::accumulate(timeTaken_std.begin(), timeTaken_std.end(), 0.0)/N;
  double var_kd = std::inner_product(timeTaken_kd.begin(), timeTaken_kd.end(), timeTaken_kd.begin(), 0)/N - square(mean_kd);
  double var_std = std::inner_product(timeTaken_std.begin(), timeTaken_std.end(), timeTaken_std.begin(), 0)/N - square(mean_std);
  */
  
  //std::cout << "Means: kdtree: " << mean_kd << " (" << std::sqrt(var_kd) << ")  std: " << mean_std << " (" << std::sqrt(var_std) << ")" << std::endl;
  
  Writer<AimsData<float> > wim("/home/cathier/tmp/out");
  wim.write(im);
}

template < typename T >
T findClosest(const std::vector<T> & c, T v)
{
  T res;
  double dist = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < size(c); ++i)
  {
    if (dist2<double>(c[i], v) < dist)
    {
      dist = dist2<double>(c[i], v);
      res = c[i];
    }
  }
  return res;
}

void testKDTreeSerious3(int argc, char * argv[])
{
  
  Reader<AimsSurfaceTriangle> r;
  Writer<AimsSurfaceTriangle> w;
  Mesh2_N mesh;
  {
    AimsApplication app( argc, aims_const_hack(argv), "testKDTreeSerious3" );
    app.addOption(r, "-i", "input mesh" );
    app.addOption(w, "-o", "output mesh" );
    app.initialize();
    AimsSurfaceTriangle s;
    r.read( s );
    Mesh2 mesh0;
    convert(s, mesh0);
    mesh = addNeighborsToMesh(mesh0);
  }
  Mesh2_N mesh2;
  //til::deepCopy(mesh, mesh2);
  mesh2.deepCopy(mesh);
  std::cout << size(getVertices(mesh)) << " " << size(getVertices(mesh2)) << std::endl;
  std::cout << size(getFaceIndices(mesh)) << " " << size(getFaceIndices(mesh2)) << std::endl;
  
  TriangleMeshGeodesicMap< Mesh2_N, double > geomap(mesh);
  std::vector<MeshTraits<Mesh2_N>::FaceIndex::value_type> startPoints;
  startPoints.push_back(getFaceIndices(mesh)[0][0]);
  std::vector<double> dist(1, 0.0);
  geomap.init(startPoints, dist);
  geomap.compute();
  std::vector<double> res = geomap.getDistanceMap();
  {
    Texture1d t;
    convert(res, t);
    Writer<Texture1d> w("/home/cathier/tmp/texgeo");
    w.write(t);
  }  
  
  for (std::size_t i = 0; i < size(getVertices(mesh)); ++i)
  {
    getVertices(mesh2)[i] += Vector<float,3>(20.0, 20.0, 20.0); 
  }  
  {
    AimsSurfaceTriangle s;
    convert(mesh2, s);
    w.write(s);
  }

  typedef KDTree<std::size_t, MeshTraits<Mesh2_N>::VertexCollection> MyKDTree;
  MyKDTree kdt(getVertices(mesh));
  makeKDTree(getVertices(mesh), kdt);
  
  std::vector<double> res2(size(res), 0);
  Mesh2_N::VertexCollection::iterator iVertex = getVertices(mesh2).begin();
  int count = 0;
  for (; iVertex != getVertices(mesh2).end(); ++iVertex)
  {
    FindClosest< double, MyKDTree > fc(kdt);
    //std::cout << ++count << " " << std::flush;
    std::size_t i = fc(*iVertex);
    /*
    Point3D<float> p = findClosest(getVertices(mesh), *iVertex);
    if (p != getVertices(mesh)[i])
    {
      std::cout << "ERROR" << std::endl;
      std::cout << p << std::endl;
      std::cout << getVertices(mesh)[i] << std::endl;
      exit(1);
    }
    */
    //res2[getVertexNumber(mesh2, &*iVertex)] = res[i];
    res2[count] = res[i];
    ++count;
  }
  std::cout << "toto2" << std::endl;
  {
    Texture1d t;
    convert(res2, t);
    Writer<Texture1d> w("/home/cathier/tmp/texgeo2");
    w.write(t);
  }  
}


void testKDTreeSerious2(int argc, char * argv[])
{
  Reader<AimsSurfaceTriangle> r;
  AimsApplication app( argc, aims_const_hack(argv), "testKDTreeSerious2" );
  app.addOption(r, "-i", "input mesh" );
  app.initialize();

  Mesh2 mesh;
  std::cout << "Aims -> mesh" << std::endl;
  {
    AimsSurfaceTriangle s;
    r.read( s );
    convert(s, mesh);
  }

  /*
  std::vector<Point3D<float>*> res;
  PRINT_TIME3("Building kdtree, simple", ({
  makeKDTree(getVertices(mesh), res);
  }));  
  std::cout << "kdtree.size " << res.size() << std::endl;
  */
  
  KDTree<Point3D<float>*, MeshTraits<Mesh2>::VertexCollection> kdt;
  PRINT_TIME3("Building kdtree, full", ({
  makeKDTree(getVertices(mesh), kdt);
  }));

  /*
  PRINT_TIME3("Searching, simple", ({
  Mesh2::VertexCollection::iterator iVertex = getVertices(mesh).begin();
  for (; iVertex != getVertices(mesh).end(); ++iVertex)
  {
    Point3D<float> target;// = *iVertex + Vector<float,3>(1.0,1.0,1.0);
    add(*iVertex, Vector<float,3>(1.0,1.0,1.0), target);
    //std::cout << target << std::endl;

    FindClosest<double, Point3D<float> > fc(res);
    fc(target);    
  }
  }));
  */
  
  
  PRINT_TIME3("Searching, full", ({
  Mesh2::VertexCollection::iterator iVertex = getVertices(mesh).begin();
  for (; iVertex != getVertices(mesh).end(); ++iVertex)
  {
    Point3D<float> target;// = *iVertex + Vector<float,3>(1.0,1.0,1.0);
    add(*iVertex, Vector<float,3>(1.0,1.0,1.0), target);
    //std::cout << target << std::endl;

    FindClosest<double, KDTree<Point3D<float>*, MeshTraits<Mesh2>::VertexCollection> > fc(kdt);
    fc(target);
  }
  }));
  
  //std::cout << "out" << std::endl;
}


/*
void testKDTree()
{
  std::vector<Vector<float,3> > v;
  
  v.push_back(Vector<float,3>(8,7,9));
  v.push_back(Vector<float,3>(28,7,7));
  v.push_back(Vector<float,3>(18,7,1));
  v.push_back(Vector<float,3>(-8,7,0));
  v.push_back(Vector<float,3>(9,27,3));
  v.push_back(Vector<float,3>(8,37,5));
  v.push_back(Vector<float,3>(2,17,2));
  v.push_back(Vector<float,3>(4,-7,2));
  v.push_back(Vector<float,3>(8,5,22));
  v.push_back(Vector<float,3>(4,3,23));
  
  Vector<float,3> v0(120, -2, 1);
  //std::vector<Vector<float,3>*> res = makeKDTree(v);
  {
    std::vector<std::size_t> res;
    makeKDTree(v, res);  
    printTree(res);
    std::cout << "ClosertPoint: ";
    //std::cout << findClosest<double>(res, v0) << std::endl;
  }
  {
    std::vector<Vector<float,3>*> res;
    makeKDTree(v, res);  
    printTree(res);
    std::cout << "ClosertPoint: ";
    std::cout << findClosest<double>(res, v0) << std::endl;
  }
  
  {
    KDTree<Vector<float,3>*, std::vector<Vector<float,3> > > k;
    makeKDTree(v, k);
    //makeKDTree<Vector<float,3> >(v, k);
    print(k);
    std::cout << "ClosertPoint: " << std::flush;
    std::cout << findClosest<double>(k, v0) << std::endl;
  }
}

*/

void testPrintStats2(int argc, char * argv[])
{
  AimsApplication app( argc, aims_const_hack(argv), "testPrintStats2" );
  aims::Reader<AimsSurfaceTriangle> reader;
  typedef Mesh2_N MeshType;
  MeshType mesh;

  AimsSurfaceTriangle wmesh;
  // Read and convert input mesh
  {
    //Mesh1 mesh0;
    Mesh2 mesh0;
    app.addOption(reader, "-i", "input mesh" );
    app.initialize();
    
    AimsSurfaceTriangle aimsmesh;
    std::cout << "Reading mesh..."<< std::flush;
    reader.read(aimsmesh);
    std::cout << "Done" << std::endl;
    std::cout << "Read mesh with " << size(getVertices(aimsmesh)) << " vertices" << std::endl;
    convert(aimsmesh, mesh0);
    std::cout << "Adding neighbors" << std::endl;
    mesh = addNeighborsToMesh(mesh0);
  }
 
 std::cout << "stats" << std::endl;
  
  meshStats(mesh);
}


void doo(int argc, char* argv[])
{
  Reader<AimsSurfaceTriangle> r;
  Writer<AimsSurfaceTriangle> w;
  AimsApplication app( argc, aims_const_hack(argv), "testMyMesh" );
  app.addOption(r, "-i", "input mesh" );
  app.addOption(w, "-o", "output mesh");
  app.initialize();

  r.read();
  AimsSurfaceTriangle s;
  r.read( s );
  AimsSurfaceTriangle s2;
  s2[0] = s[49];
  w.write(s2);
}


void testBinaryTree()
{
  NaryTree<double, 2> bt;
  NaryTree<double, 2>::iterator i;

  i = bt.addChild(0, 0, 2.34);
  bt.addChild(i, 1, 8.99);
  i = bt.addChild(i, 0, 4.56);
  i = bt.addChild(i, 1, 7.89);
  print(bt);
    
  /*
  std::cout << "b" << std::endl;
  bt.insert(i, 2.34);
  std::cout << "c" << std::endl;
  i = i.child(0);
  std::cout << "d" << std::endl;
  bt.insert(i, 4.56);
  std::cout << "e" << std::endl;
  i = i.child(1);
  std::cout << "f" << std::endl;
  bt.insert(i, 7.89);
  std::cout << "g" << std::endl;
  print(bt);
  */
}

void testAdd1()
{
  Vector<double,3> v1, v2(23.45, 287.9823, 374397924.2938), v3(12098409834.23, 323.283, 0.000000001);
  for (int i = 0; i < 100000; ++i)
  {
    v1 = v2 + v3;
  }
}

void testAdd2()
{
  Vector<double,3> v1, v2(23.45, 287.9823, 374397924.2938), v3(12098409834.23, 323.283, 0.000000001);
  for (int i = 0; i < 100000; ++i)
  {
    add(v2, v3, v1);
  }
}

void testAdd3()
{
  Vector<double,3> v1, v2(23.45, 287.9823, 374397924.2938), v3(12098409834.23, 323.283, 0.000000001);
  for (int i = 0; i < 100000; ++i)
  {
    Vector<double,3> v(v2);
  }
}



void testAdds()
{
  PRINT_TIME3("add1", testAdd1());
  PRINT_TIME3("add2", testAdd2());
  PRINT_TIME3("add3", testAdd2());
  PRINT_TIME3("add1", testAdd1());
  PRINT_TIME3("add2", testAdd2());
  PRINT_TIME3("add3", testAdd2());
  PRINT_TIME3("add1", testAdd1());
  PRINT_TIME3("add2", testAdd2());
  PRINT_TIME3("add3", testAdd2());
}


void testTILFunctors()
{
  using namespace til;
  std::vector<double> a(10,0);
  std::list<double> b(10,3);
  
  loop1(*_1 = *_2 **_2, a, b);
  
  std::cout << "*" << a[0] << std::endl;
}

void imageJoin(int argc, char* argv[])
{
  Reader<AimsData<short> > r1;
  Reader<AimsData<short> > r2;
  Writer<AimsData<unsigned char> > w;
  AimsApplication app( argc, aims_const_hack(argv), "imageJoin" );
  app.addOption(r1, "-i1", "first image" );
  app.addOption(r2, "-i2", "second image" );
  app.addOption(w, "-o", "output image");
  app.initialize();
  
  // Read images and convert them in U8
  AimsData<unsigned char> *pim1;
  AimsData<unsigned char> *pim2;
  {
    AimsData<short> _im1;
    std::cout << "Reading first file..." << std::flush;
    r1.read(_im1);
    std::cout << "OK" << std::endl;
    AimsData<short> _im2;
    std::cout << "Reading second file..." << std::flush;
    r2.read(_im2);
    std::cout << "OK" << std::endl;
    carto::Converter<AimsData<short>, AimsData<unsigned char> > c;
    pim1 = c(_im1);
    pim2 = c(_im2);
  }
  AimsData<unsigned char> &im1 = *pim1;
  AimsData<unsigned char> &im2 = *pim2;
  
  
  std::cout << "Merging..." << std::flush;
  for (int i = 0; i < im1.dimX()*im1.dimY()*im1.dimZ(); ++i)
  {
    if (im2[i] > im1[i]) im1[i] = im2[i];
  }
  std::cout << "OK" << std::endl;

  std::cout << "Writing output..." << std::flush;
  w.write(im1);
  std::cout << "OK" << std::endl;
}


void meshTrace(int argc, char* argv[])
{
  Reader<AimsData<short> > r;
  Reader<AimsSurfaceTriangle> rm;
  Writer<AimsData<short> > w;
  AimsApplication app( argc, aims_const_hack(argv), "meshTrace" );
  app.addOption(r, "-i", "anatomic image" );
  app.addOption(rm, "-m", "surface mesh" );
  app.addOption(w, "-o", "output image");
  app.initialize();
  
  AimsData<short> im;
  r.read(im);
  AimsSurfaceTriangle mesh;
  rm.read(mesh);

  // clear image
  for (int i = 0; i < im.dimX()*im.dimY()*im.dimZ(); ++i)
  {
    im[i] = 0;
  }
  
  // Set mesh trace
  MeshTraits<AimsSurfaceTriangle>::VertexCollection vertices = getVertices(mesh);
  MeshTraits<AimsSurfaceTriangle>::VertexCollection::const_iterator iVertex = vertices.begin();
  for (; iVertex != vertices.end(); ++iVertex)
  {
    im(int((*iVertex)[0]/im.sizeX()),
       int((*iVertex)[1]/im.sizeY()),
       int((*iVertex)[2]/im.sizeZ())) = 1;
/*    im(int(rint((*iVertex)[0]/im.sizeX())),
       int(rint((*iVertex)[1]/im.sizeY())),
       int(rint((*iVertex)[2]/im.sizeZ()))) = 1;
       */
  }
  w.write(im);
}

void testSparseVectors()
{
  til::SparseVector<double> v1(10);
  til::SparseVector<double> v2(10);
  
  v1.set(2, 30);
  v1.set(7, 15);
  v1.set(3, 12);
  v1.set(1, 39);
  v1.set(1, -29);
  v1.set(9, 95.2);
  v1.set(0, 1);

  v2.set(3, 330);
  v2.set(5, 315);
  v2.set(4, 212);
  v2.set(2, 239);
  v2.set(8, -429);
  v2.set(5, 295.2);
  v2.set(4, 11);
  
  std::cout << v1 << std::endl;
  std::cout << v2 << std::endl;
  
  std::cout << v1*v2 << std::endl;
  
  std::cout << v1+v2 << std::endl;

  std::cout << til::dot<double>(v1,v2) << std::endl;

//  std::cout << (v1 == v2) << std::endl;
//  std::cout << (v1 == v1) << std::endl;

  
  std::cout << "Starting now" << std::endl;
  const int N = 1000000;
  
  PRINT_TIME3("Rand", ({
  for (int i = 0; i < 10*N; ++i)
  {
    std::rand() % N;
    std::rand() % N;
  }
  }));

  PRINT_TIME3("Rand write2", ({
  std::vector<til::SparseVector<double> > v(N, til::SparseVector<double>(N));
  for (int i = 0; i < 10*N; ++i)
  {
    v[std::rand() % N][std::rand() % N] = 1.0;
  }
  }));

  PRINT_TIME3("Rand write", ({
  std::vector<til::SparseVector<double> > v(N, til::SparseVector<double>(N));
  for (int i = 0; i < 10*N; ++i)
  {
    v[std::rand() % N].set(std::rand()%N, 1.0);
  }
  }));

  double d;
  PRINT_TIME3("Rand read", ({
  std::vector<til::SparseVector<double> > v(N, til::SparseVector<double>(N));
  for (int i = 0; i < 10*N; ++i)
  {
    d = v[std::rand() % N].get(std::rand()%N);
  }
  }));

  PRINT_TIME3("Rand write2", ({
  std::vector<til::SparseVector<double> > v(N, til::SparseVector<double>(N));
  for (int i = 0; i < 10*N; ++i)
  {
    v[std::rand() % N][std::rand() % N] = 1.0;
  }
  }));
}

void mask(int argc, char* argv[])
{
  Reader<AimsData<unsigned char> > rim;
  Reader<AimsData<short> > rmask;
  Writer<AimsData<unsigned char> > wim;
  AimsApplication app( argc, aims_const_hack(argv), "mask" );
  app.addOption(rim, "-i", "input image" );
  app.addOption(rmask, "-m", "mask" );
  app.addOption(wim, "-o", "output image");
  app.initialize();
  
  std::cout << "Reading images..." << std::flush;
  AimsData<unsigned char> im;
  rim.read(im);
  AimsData<short> mask;
  rmask.read(mask);
  std::cout << "OK" << std::endl;
  
  // Get mask value.
  // We suppose here, and don't check, that the mask is a binary image, and that the background
  // is zero.
  short maskValue = 0;
  {
    for (int i = 0; i < mask.dimX()*mask.dimY()*mask.dimZ(); ++i)
    {
      if (mask[i])
      {
        maskValue = mask[i];
        break;
      }
    }
  }
  std::cout << "Mask value found : " << maskValue << std::endl;
  
  //carto::rc_ptr<aims::Interpolator> interp = getLinearInterpolator(mask);
  
  //MaskIteratorOf<AimsData<short> > maskIter(mask);
  
  rc_ptr<MaskIterator> maskIter;
  maskIter = getMaskIterator(mask);
  
  std::cout << "Masking..." << std::flush;
  int i, j, k;
  int mi, mj, mk;
  for ( k = 0; k < im.dimZ(); k++ )
  {
    if (k % int(im.dimZ()/10.0) == 0)
    {
      std::cout << 10 * int((k / int(im.dimZ()/10.0))) << "%..." << std::flush;
    }
    for ( j = 0; j < im.dimY(); j++ )
    for ( i = 0; i < im.dimX(); i++ )
    {
      /*
      if ((*interp)(i*im.sizeX(), j*im.sizeY(), k*im.sizeZ()) < 0.9 * maskValue)
      {
        im(i,j,k) = 0;
      }
      */


      bool test1 = maskIter->contains(Point3df(i*im.sizeX(),
                                       j*im.sizeY(),
                                       k*im.sizeZ()));


      /*
      if (!maskIter->contains(Point3df(i*im.sizeX(),
                                       j*im.sizeY(),
                                       k*im.sizeZ())))
      / *
      if (!maskIter.contains(Point3df(i*im.sizeX(),
                                      j*im.sizeY(),
                                      k*im.sizeZ())))
      * /
      {
        im(i,j,k) = 0;
      }
      */

      mi = int(std::floor(i*im.sizeX() / mask.sizeX() + 0.5));
      mj = int(std::floor(j*im.sizeY() / mask.sizeY() + 0.5));
      mk = int(std::floor(k*im.sizeZ() / mask.sizeZ() + 0.5));

      bool test2 = (mi < 0 || mi >= mask.dimX() ||
          mj < 0 || mj >= mask.dimY() ||
          mk < 0 || mk >= mask.dimZ() ||
          mask(mi,mj,mk) == 0);

      if (test1 == test2)
      {
        std::cout << "arg" << test1 << " " << test2 << std::endl;
        std::cout << i*im.sizeX() / mask.sizeX() << " " << j*im.sizeY() / mask.sizeY() << " " << k*im.sizeZ() / mask.sizeZ() << std::endl;
        std::cout << mi << " " << mj << " " << mk << std::endl;
        std::cout << rint(i*im.sizeX() / mask.sizeX()) << " " << rint(j*im.sizeY() / mask.sizeY()) << " " << rint(k*im.sizeZ() / mask.sizeZ()) << std::endl;
      }

      /*
      if (mi < 0 || mi >= mask.dimX() ||
          mj < 0 || mj >= mask.dimY() ||
          mk < 0 || mk >= mask.dimZ() ||
          mask(mi,mj,mk) == 0)
      {
        im(i,j,k) = 0;
      }
      //else im(i,j,k) = 1;
       */

    }
  }
  std::cout << "OK" << std::endl;

  std::cout << "Writing image..." << std::flush;
  wim.write(im);  
  std::cout << "OK" << std::endl;
}

void testCircularNeighbors()
{
  AimsSurfaceTriangle *s1 = makeSphere(Point3df(1.0, 2.0, 3.0), 0.7, 7);
  //AimsSurfaceTriangle &s = *s1;
  Mesh2_N mesh;
  {
    Mesh2 mesh0;
    convert(*s1, mesh0);
    std::cout << "Adding neighbors" << std::endl;
    mesh = addNeighborsToMesh(mesh0);
  }

  std::cout << "A" << std::endl;
  std::vector<std::vector<Point3D<float>*> > neigh = getNeighborIndices(mesh);
  std::cout << "B" << std::endl;
  shared_ptr<std::vector<std::vector<Point3D<float>*> > > neighc
    = getCircularNeighborIndices(mesh);
  std::cout << "C" << std::endl;
    
  std::cout << size(*neighc) << std::endl;
  std::cout << size(neigh) << std::endl;
  
  for (std::size_t i = 0; i < size(*neighc); ++i)
  {
    if (size((*neighc)[i])!=size(neigh[i]))
    {
      std::cout << "ERROR " << i << std::endl;
    }
  }
}

long double imnorm(const AimsData<float> & im)
{
  long double res = 0.0;
  int i, j, k;
  for ( k = 1; k < im.dimZ(); k++ )
  for ( j = 1; j < im.dimY(); j++ )
  for ( i = 1; i < im.dimX(); i++ )
  {
    res +=
      til::square<long double>(im(i,j,k) - im(i-1,j,k)) +
      til::square<long double>(im(i,j,k) - im(i,j-1,k)) +
      til::square<long double>(im(i,j,k) - im(i,j,k-1));
  }
  return res;
}

void vfnorm(int argc, char * argv[])
{
  std::string name;
  AimsApplication app( argc, aims_const_hack(argv), "vfnorm" );
  app.addOption(name, "-i", "input image" );
  app.initialize();
  
  long double res = 0.0;
  {
    Reader<AimsData<float> > r(name + ".x.ima");
    AimsData<float> im;
    r.read(im);
    res += imnorm(im);
  }
  {
    Reader<AimsData<float> > r(name + ".y.ima");
    AimsData<float> im;
    r.read(im);
    res += imnorm(im);
  }
  {
    Reader<AimsData<float> > r(name + ".z.ima");
    AimsData<float> im;
    r.read(im);
    res += imnorm(im);
  }
  std::cout << res << std::endl;
}

void testsizeof()
{
  std::cout << "double " << sizeof(double) << std::endl;
  std::cout << "long double " << sizeof(long double) << std::endl;
}

void testGeodesicDistance(int argc, char * argv[])
{
  //AimsSurfaceTriangle &s = *s1;
  Mesh2_N mesh;
  {
    Reader<AimsSurfaceTriangle> r;
    AimsApplication app( argc, aims_const_hack(argv), "testMyMesh" );
    app.addOption(r, "-i", "input mesh" );
    app.initialize();
    AimsSurfaceTriangle s;
    r.read( s );
    Mesh2 mesh0;
    convert(s, mesh0);

    //AimsSurfaceTriangle *s1 = makeSphere(Point3df(1.0, 2.0, 3.0), 0.7, 7);
    //Mesh2 mesh0;
    //convert(*s1, mesh0);
    
    std::cout << "Adding neighbors" << std::endl;
    mesh = addNeighborsToMesh(mesh0);
  }

  std::cout << "Computing geomap" << std::endl;
  ghost::GMapStop_AboveThreshold<double> stopGhost(10.0);
  std::cout << "Finishing" << std::endl;

  TriangleMeshGeodesicMap<Mesh2_N, double, ghost::GMapStop_AboveThreshold<double>, policy::GMap_DefaultStorage<til::sparse_vector, Mesh2_N, double > > geomap(mesh, stopGhost);
  std::vector<MeshTraits<Mesh2_N>::FaceIndex::value_type> startPoints;
  startPoints.push_back(getFaceIndices(mesh)[0][0]);
  std::vector<double> dist(1, 0.0);
  geomap.init(startPoints, dist);
  geomap.compute();
  til::sparse_vector<double> sres = geomap.getDistanceMap();

  {
    til::sparse_vector<double>::sparse_iterator iRes = sres.sparse_begin();
    for (; iRes != sres.sparse_end(); ++iRes)
    {
      std::cout << iRes->second << " ";
    }
    std::cout << std::endl;
  }
  
  
  std::vector<double> res(size(sres));
  {
    std::vector<double>::iterator iRes = res.begin();
    til::sparse_vector<double>::const_iterator iRes2 = sres.begin();
    for (; iRes2 != sres.end(); ++iRes, ++iRes2)
    {
      *iRes = *iRes2;
    }
  }

  /*  
  TriangleMeshGeodesicMap<Mesh2_N, double, ghost::GMapStop_AboveThreshold<double> > geomap(mesh, stopGhost);
  std::vector<MeshTraits<Mesh2_N>::FaceIndex::value_type> startPoints;
  startPoints.push_back(getFaceIndices(mesh)[0][0]);
  std::vector<double> dist(1, 0.0);
  geomap.init(startPoints, dist);
  geomap.compute();
  std::vector<double> res = geomap.getDistanceMap();
  */
  
  
  for (std::size_t i = 0; i < size(res); ++i) res[i] = til::min(res[i], 100.0);
  
  {
    std::vector<double> res2 = res;
    std::sort  (res2.begin(), res2.end());
    std::vector<double>::iterator newend = std::unique(res2.begin(), res2.end());
    std::cout << "Unique : " << std::distance(res2.begin(), newend) << std::endl;
    for (int i = 0; i < 10; ++i) std::cout << res[i] << " ";
    std::cout << std::endl;
  }  
  {
    AimsSurfaceTriangle s;
    convert(mesh, s);
    Writer<AimsSurfaceTriangle> w("sphere");
    w.write(s);
  }
  {
    Texture1d t; //(1, size(res));
    t.reserve(size(res));
    for (std::size_t i=0; i<size(res); ++i) //t[i] = res[i];
      t.push_back(res[i]);
    Writer<Texture1d> w("spheredist");
    w.write(t);
  }
  
  
}

void testGeodesicDistance2(int argc, char * argv[])
{
  //AimsSurfaceTriangle &s = *s1;
  Mesh2_N mesh;
  {
    Reader<AimsSurfaceTriangle> r;
    AimsApplication app( argc, aims_const_hack(argv), "testMyMesh" );
    app.addOption(r, "-i", "input mesh" );
    app.initialize();
    AimsSurfaceTriangle s;
    r.read( s );
    Mesh2 mesh0;
    convert(s, mesh0);

    //AimsSurfaceTriangle *s1 = makeSphere(Point3df(1.0, 2.0, 3.0), 0.7, 7);
    //Mesh2 mesh0;
    //convert(*s1, mesh0);
    
    std::cout << "Adding neighbors" << std::endl;
    mesh = addNeighborsToMesh(mesh0);
  }


  double distance = 10.0;
  std::cout << "Computing geomap" << std::endl;
  ghost::GMapStop_AboveThreshold<double> stopGhost(distance);
  TriangleMeshGeodesicMap<Mesh2_N, double, ghost::GMapStop_AboveThreshold<double>, policy::GMap_DefaultStorage<til::sparse_vector, Mesh2_N, double > > geomap(mesh, stopGhost);
  std::vector<const MeshTraits<Mesh2_N>::Vertex *> startPoints(1);
  std::vector<double> dist(1, 0.0);

  std::vector<double> res(size(getVertices(mesh)), 0.0);

  int count = 0;

  MeshTraits<Mesh2_N>::VertexCollection::const_iterator iVertex = getVertices(mesh).begin();
  for(; iVertex != getVertices(mesh).end(); ++iVertex)
  {
    ++count;
    //if (++count > 10) break;
    startPoints[0] = &*iVertex;
    geomap.init(startPoints, dist);
    geomap.compute();
    til::sparse_vector<double> sres = geomap.getDistanceMap();
    //sres.setDefaultValue(0);
    //std::transform(sres.begin(), sres.end(), res.begin(), res.begin(), std::plus<double>());
  }
  std::cout << count << std::endl;

  std::cout << "Finishing" << std::endl;

  /*
  {
    til::sparse_vector<double>::sparse_iterator iRes = sres.sparse_begin();
    for (; iRes != sres.sparse_end(); ++iRes)
    {
      std::cout << iRes->second << " ";
    }
    std::cout << std::endl;
  }
  */
  /*
  std::vector<double> res(size(sres));
  {
    std::vector<double>::iterator iRes = res.begin();
    til::sparse_vector<double>::const_iterator iRes2 = sres.begin();
    for (; iRes2 != sres.end(); ++iRes, ++iRes2)
    {
      *iRes = *iRes2;
    }
  }
  */
  
  /*  
  TriangleMeshGeodesicMap<Mesh2_N, double, ghost::GMapStop_AboveThreshold<double> > geomap(mesh, stopGhost);
  std::vector<MeshTraits<Mesh2_N>::FaceIndex::value_type> startPoints;
  startPoints.push_back(getFaceIndices(mesh)[0][0]);
  std::vector<double> dist(1, 0.0);
  geomap.init(startPoints, dist);
  geomap.compute();
  std::vector<double> res = geomap.getDistanceMap();
  */
  
  
  //for (std::size_t i = 0; i < size(res); ++i) res[i] = til::min(res[i], 100.0);
  
  {
    std::vector<double> res2 = res;
    std::sort  (res2.begin(), res2.end());
    std::vector<double>::iterator newend = std::unique(res2.begin(), res2.end());
    std::cout << "Unique : " << std::distance(res2.begin(), newend) << std::endl;
    for (int i = 0; i < 10; ++i) std::cout << res[i] << " ";
    std::cout << std::endl;
  }  
  {
    AimsSurfaceTriangle s;
    convert(mesh, s);
    Writer<AimsSurfaceTriangle> w("sphere");
    w.write(s);
  }
  {
    Texture1d t; //(1, size(res));
    t.reserve(size(res));
    for (std::size_t i=0; i<size(res); ++i) //t[i] = res[i];
      t.push_back(res[i]);
    Writer<Texture1d> w("spheredist");
    w.write(t);
  }  
}


void testSolver()
{
  std::cout << norm<double>(Vector<float,3>(1, 2, 3) - Vector<float,3>(4, 5, 6)) << std::endl;
  til::math::PolySolver_real<double, til::math::policy::InfinitySolutions_None> solver;
  solver.solve(1, 1, 1);
  for (int i = 0; i < solver.nsols(); ++i) std::cout << solver.sols()[i] << " "; std::cout << std::endl;
  solver.solve(1, 1, -1);
  for (int i = 0; i < solver.nsols(); ++i) std::cout << solver.sols()[i] << " "; std::cout << std::endl;
  solver.solve(0, 1, -1);
  for (int i = 0; i < solver.nsols(); ++i) std::cout << solver.sols()[i] << " "; std::cout << std::endl;
  solver.solve(1, 0, -2);
  for (int i = 0; i < solver.nsols(); ++i) std::cout << solver.sols()[i] << " "; std::cout << std::endl;
}



template < typename TMesh, typename TPrec, typename TConn, typename Finder >
class FiberDistance
{
public: // typedefs

  typedef typename Finder::index_type index_type;
  
public: // constructors & destructor

  FiberDistance(const TConn & conn1, const TConn & conn2, Finder finder) : m_conn1(conn1), m_conn2(conn2), m_finder(finder) {}

  TPrec operator()(const TMesh & mesh1, const TMesh & mesh2)
  {
    TPrec d = 0.0;
    
    for (typename MeshTraits<TMesh>::VertexCollection iVertex1 = getVertices(mesh1).begin();
         iVertex1 != getVertices(mesh1).end(); ++iVertex1)
    {
      index_type i2 = finder(*iVertex1);

      d += norm2<TPrec>(*iVertex1, getVertices(mesh2)[i2]);      
    }
    
    return d;
  }
  
private: // data

  const TConn & m_conn1;
  const TConn & m_conn2;
  Finder m_finder;
};


void testBundles(int argc, char * argv[])
{
  Mesh2_N mesh;
  {
    Reader<AimsSurfaceTriangle> r("/volatile/cathier/data/icbm/icbm100T_Lwhite.mesh");
    AimsSurfaceTriangle s;
    r.read( s );
    Mesh2 mesh0;
    convert(s, mesh0);
    mesh = addNeighborsToMesh(mesh0);
  }
  
  const double DIST = 5.0;
  std::cout << "Computing geomap" << std::endl;
  ghost::GMapStop_AboveThreshold<double> stopGhost(DIST);
  TriangleMeshGeodesicMap<Mesh2_N, double, ghost::GMapStop_AboveThreshold<double>, policy::GMap_DefaultStorage<til::sparse_vector, Mesh2_N, double > > geomap(mesh, stopGhost);
  std::vector<const MeshTraits<Mesh2_N>::Vertex *> startPoints(1);
  std::vector<double> dist(1, 0.0);
  typedef std::vector< std::pair<std::size_t, double> > QuickMap;
  std::vector< QuickMap > res(size(getVertices(mesh)));
  //MeshTraits<Mesh2_N>::VertexCollection::const_iterator iVertex = getVertices(mesh).begin();
  //for(; iVertex != getVertices(mesh).end(); ++iVertex)
  std::cout << "wow" << std::endl;
  {
    til::Histogram<std::size_t> h;
    std::vector<std::size_t> nneigh(size(getVertices(mesh)));
    for (std::size_t i = 0; i < size(getVertices(mesh)); ++i)
    {
      startPoints[0] = &(getVertices(mesh)[i]);
      geomap.init(startPoints, dist);
      geomap.compute();
      //std::cout << count << std::endl;
      til::sparse_vector<double> tmp = geomap.getDistanceMap();
      res[i].resize(size(tmp.getMap()));
      convert(tmp.getMap(), res[i]);
      //std::cout << count << std::endl;
      h.accumulate(size(res[i]));
      nneigh[i] = size(res[i]);
      //res[i] = geomap.getDistanceMap();
      //sres.setDefaultValue(0);
      //std::transform(sres.begin(), sres.end(), res.begin(), res.begin(), std::plus<double>());
    }
    std::cout << "Finishing" << std::endl;

    // writing number of neighbors in the DIST-cell
    {
      Texture1d t(1, size(nneigh));
      convert(nneigh, t);
      Writer<Texture1d> w("/home/cathier/tmp/nneigh");
      w.write(t);
    }
  }

  // writing distances in the DIST-cell of the first vertex
  {
    Texture1d t(1, size(getVertices(mesh)));
    for (std::size_t i = 0; i < size(res[0]); ++i)
    {
      t.item(res[0][i].first) = res[0][i].second;
    }
    Writer<Texture1d> w("/home/cathier/tmp/toubou");
    w.write(t);
  }

  // loading fibers
  boost::shared_ptr<std::vector<std::vector<Point3D<float> > > > pfibers;
  {
    aims::BundleReader r("/home/Panabase/pascal/fibers/res/old/res.bundles");
    BundleLoader loader;
    r.addBundleListener(loader);
    r.read();
    pfibers = loader.getFibers();
  }
  typedef std::vector<Point3D<float> > Fiber;
  typedef std::vector<Fiber> Fibers;
  Fibers & fibers = *pfibers;

  std::cout << "Number of fibers: " << size(fibers) << std::endl;
  /*
  for (std::size_t i = 0; i < size(fibers[0]); ++i)
  {
    std::cout << fibers[0][i] << std::endl;
  }
  */

  std::cout << "Generating kdtree" << std::endl;
  typedef KDTree<std::size_t, MeshTraits<Mesh2_N>::VertexCollection> MyKDTree;
  MyKDTree kdt(getVertices(mesh));
  makeKDTree(getVertices(mesh), kdt);

  std::cout << "Looking for closest points" << std::endl;
  typedef til::sparse_vector<double> Connectivity;
  typedef std::vector< Connectivity > Connectivities;
  Connectivities conn(size(getVertices(mesh)), til::sparse_vector<double>(size(getVertices(mesh))));
  for (Fibers::const_iterator iFiber = fibers.begin(); iFiber != fibers.end(); ++iFiber)
  {
    FindClosest< double, MyKDTree > fc(kdt);
    std::size_t A = fc(iFiber->front());
    std::size_t B = fc(iFiber->back());
    for (QuickMap::const_iterator i = res[A].begin(); i != res[A].end(); ++i)
    for (QuickMap::const_iterator j = res[B].begin(); j != res[B].end(); ++j)
    {
      double e1 = til::square(i->second);
      double e2 = til::square(j->second);
      double w = std::exp(- ( e1 + e2 ) / ( 2*DIST*DIST ));
      conn[i->first][j->first] += w;
      conn[j->first][i->first] += w;
    }
  }
  
  // To compress a little bit, remove points that are below some threshold
  std::cout << "Removing weak connections" << std::endl;
  {
    const double THRESH = 1.0;
    Connectivities::iterator i = conn.begin();
    int count1 = 0;
    int count2 = 0;
    for (; i != conn.end(); ++i)
    {
      Connectivity::sparse_iterator j = i->sparse_begin();
      for (; j != i->sparse_end(); ++j)
      {
        ++count1;
        if (j->second <= THRESH)
        {
          ++count2;
          i->erase(j);
        }
      }
    }
    std::cout << "Removed " << count2 << " elements out of " << count1 << std::endl;
  }
  
  /*
  std::cout << "Writing textures" << std::endl;
  for (int V = 0; V < 10; ++V)
  {
    Texture1d t(1, size(getVertices(mesh)));
    til::sparse_vector<double>::Map::const_iterator i = conn[V].getMap().begin();
    for (; i != conn[V].getMap().end(); ++i)
    {
      t.item(i->first) = i->second;
    }
    char s[128];
    std::sprintf(s, "/home/cathier/tmp/wi%i", V);
    Writer<Texture1d> w(s);
    w.write(t);
  }
  */  

  // Registration between brains
  {
    typedef FindClosest< double, MyKDTree > Finder;
    Finder finder(kdt);
    FiberDistance< Mesh2_N, double, Connectivities, Finder > dist(conn, conn, finder);
  }

  
}

void testsqrt()
{
  {
    typedef double prec;
    std::cout << std::numeric_limits<prec>::min() << std::endl;
    std::cout << std::sqrt(std::numeric_limits<prec>::min()) << std::endl;
  }  
  {
    typedef float prec;
    std::cout << std::numeric_limits<prec>::min() << std::endl;
    std::cout << std::sqrt(std::numeric_limits<prec>::min()) << std::endl;
  }
  
  /*
  std::cout << til::max(std::sqrt(std::numeric_limits<float>::min()), 
                   std::sqrt(std::numeric_limits<double>::min())) << std::endl;


  std::cout << std::max<long double>(std::sqrt(std::numeric_limits<float>::min()), 
                                     std::sqrt(std::numeric_limits<double>::min())) << std::endl;
    */
}


struct Toto : public std::unary_function<Vector<double, 5>, float>
{
  float operator()(const Vector<double, 5> & v)
  {
    return v[0]*v[0] + 100*v[1]*v[1] + 1000*v[2]*v[2] + 10*v[3]*v[3] + 50*v[4]*v[4] + v[3]*v[4];
  }
};

void testPowell()
{
  til::Powell<Toto> p;
  Vector<double, 5> min;
  min[0] = 1000.0;
  min[1] = -1000.0;
  min[2] = 2000.0;
  min[3] = -10000.0;
  min[4] = 2;
  
  min = p(min);
  std::cout << min << std::endl;
}



int main( int argc, char* argv[] )
{
  try
  {
    std::cout << "main starts now!" << std::endl;
    
    //foo();
    //testMTL();
    //testConnectivityMatrix(argc, argv);
    //testMul();
    //testCentroid();
    //testCentroid3();
    //testCentroid2();
    //testSpeed();
    //testMyMesh(argc, argv);
    //myinflate(argc, argv);
    //testConvert(argc, argv);
    //testKDTree();
    //doo(argc, argv);
	  //testKDTreeSerious(argc, argv);
    //testKDTreeSerious3(argc, argv);
    //testKDTreeSerious2(argc, argv);
    //testBinaryTree();
    //testPrintStats2(argc, argv);
    //testAdds();
    //testTILFunctors();
    //imageJoin(argc, argv);
    //meshTrace(argc, argv);
    //testSparseVectors();
    //mask(argc, argv);
    //testCircularNeighbors();
    //vfnorm(argc, argv);
    //testsizeof();
    //testGeodesicDistance(argc, argv);
    //testGeodesicDistance2(argc, argv);
    testBundles(argc, argv);
    //testPowell();
    //testsqrt();
  	return 0;

		std::cout << "Fresh start" << std::endl;
		
		Reader<AimsSurfaceTriangle> r;
		Writer<Texture1d> w;
		Writer<AimsSurfaceTriangle> r2;
		AimsApplication app( argc, aims_const_hack(argv), "voila voila" );
		app.addOption( r, "-i", "input mesh" );
		app.addOption( w, "-o", "output texture" );
		app.addOption(r2, "-u", "output mesh");
		app.initialize();

		AimsSurfaceTriangle s;
		r.read( s );
		Texture1d tex;
		std::vector<Point3df> & vert = s.vertex();
		unsigned i, n = vert.size();
		tex.reserve( n );
		for( i=0; i<n; ++i )
			tex.push_back( vert[i][0] );
		
		w.write( tex );
		r2.write(s);
	}
	catch( user_interruption & )
	{
		// Exceptions thrown by command line parser (already handled, simply exit)
	}
	catch( exception & e )
	{
		cerr << e.what() << endl;      
	}
}

