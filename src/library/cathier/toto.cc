// includes from STL
//#include <math.h>
#include <cmath>
//#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include <sstream>
#include <sys/times.h>
#include <sys/time.h>
#include <vector>

// includes from BOOST
#include <boost/logic/tribool_io.hpp>
//#include <boost/utility.hpp>
//#include <boost/utility/result_of.hpp>
#include <boost/random/uniform_real.hpp>

// includes from MTL
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
// do not include MTL for GCC 3.4.4 as it does not compile
#if GCC_VERSION < 030404
#include <mtl/matrix.h>
#include <mtl/mtl.h>
#include <mtl/utils.h>
#endif


// includes from CARTO
//#include "carto/Converter.h"

// includes from AIMS
#include "connectomist/fibertracking/bundles.h"
#include "aims/getopt/getopt2.h"
#include "aims/io/motionR.h"
#include "aims/io/reader.h"
#include "aims/io/writer.h"
#include "aims/mesh/inflate.h"
#include "aims/mesh/surface.h"
#include "aims/mesh/surfacegen.h"
#include "aims/mesh/texture.h"
#include "aims/resampling/linearInterpolator.h"
#include "aims/resampling/motion.h"
#include "aims/roi/maskIterator.h"
#include "aims/signalfilter/g3dsmooth.h"
#include "aims/signalfilter/gjacobian.h"
#include "aims/utility/anytype_reader.h"
#include "aims/utility/converter_volume.h"
#include "aims/vector/vector.h"

// includes from TIL
#include "til/proba_distributions.h"
#include "til/AffineMap.h"
#include "til/numeric_array.h"
#include "til/TExpr.h"

// includes from TIL
#include "accessors.h"
#include "aims_wrap.h"
#include "binary_tree.h"
#if GCC_VERSION < 030404
//#include "connectivityMatrix.h"
#endif
//#include "containers.h"
//#include "domino.h"
#include "dwt.h"
#include "Forces.h"
#include "func_iterator.h"
#include "fuzzy_logic.h"
#include "geometrics.h"
#include "globalTraits.h"
#include "histogram.h"
#include "io.h"
#include "kdtree.h"
#include "math_functions.h"
#include "mesh_conversion.h"
#include "mesh_decimation.h"
#include "MeshTraits.h"
#include "meshUtils.h"
#include "minTools.h"
#include "miscUtils.h"
#include "progress_indicator.h"
#include "triangle_mesh_geodesic_map.h"
#include "SparseVector.h"

#include "totodebug.h"
#include "cathier/toto.h"

extern "C"
{
#include "lapack.h"
}



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
using namespace comist;


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
/// Columns of mat are summed and added to vec.181
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
    
  mtl::mult(mat, vec, res);
}

/*
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
  std::cout << "Read mesh with" << til::size(getVertices(mesh)) << " vertices" << std::endl;
  
  std::cout << "Building connectivity matrix..." << std::flush;
  //BuildConnectivityMatrix::SuggestedType mat(til::size(getVertices(mesh)), til::size(getVertices(mesh)));
  //buildConnectivityMatrix(mesh, mat);
  TIL_FOR_ALL_CONST_FACES(mesh)
  {
    std::size_t i;
    // NB: actually, container should be STL-compliant...
    for (i=0; i<size(*iFace)-1; ++i)
    {
      mat((*iFace)[i], (*iFace)[i+1]) = 1;
    }
    mat((*iFace)[i], (*iFace)[0]) = 1;
  }
  
  std::cout << "OK" << std::endl;
  
  int n = (til::size(getVertices(mesh))-1) - 5;
  std::cout << "Printing neighbors of " << n << std::endl;
  for (int i = 0; i < mat.ncols(); ++i)
  {
    if (mat(n,i)) std::cout << i << "(" << mat(n,i) << ")" <<" "<< std::flush;
  }
  std::cout << std::endl;
  std::cout << "Printing neighbors (2) of " << n << std::endl;
  for (int i = 0; i < mat.ncols(); ++i)
  {
    if (mat(i,n)) std::cout << i << "(" << mat(i,n) << ")" <<" "<< std::flush;
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
*/
#endif

template < typename TMesh >
std::vector<til::sparse_vector<double> >
//std::vector<std::vector<double> >
getDistToNeighbors(const TMesh & mesh, double distance)
{
  std::cout << "Computing geomap (1)..." << std::flush;
  std::size_t n = size(getVertices(mesh));
  til::ghost::GMapStop_AboveThreshold<double> stopGhost(distance);
  //til::Triangle_mesh_geodesic_map<TMesh, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage<til::sparse_vector, TMesh, double > > 
  //til::Graph_distance_map<TMesh, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage<til::sparse_vector, TMesh, double > > 
  typedef std::vector<std::vector<std::size_t> > CNeighborhoods;
  shared_ptr<CNeighborhoods> pneighc = til::circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  til::Graph_distance_map<typename TMesh::VertexCollection, CNeighborhoods, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage<til::sparse_vector, double > >
  //til::Graph_distance_map<TMesh, double, til::ghost::GMapStop_AboveThreshold<double> > 
    geomap(getVertices(mesh), *pneighc, stopGhost);
  std::vector<std::size_t> startPoints(1);
  std::vector<double> dist(1, 0.0);
  //std::vector<double> res(n, 0.0);
  std::vector<til::sparse_vector<double> > res(n);
  shared_ptr<til::sparse_vector<unsigned char> > labels;
  //std::vector<std::vector<double> > res(n);
  //std::vector<unsigned char> labels;
  for (std::size_t i=0; i < n; ++i)
  {
    startPoints[0] = i;
    geomap.init(startPoints, dist);
    geomap.process();
    res[i] = *(geomap.distanceMap());
    labels = geomap.labels();

    //for (til::sparse_vector<double>::sparse_iterator j = res[i].sparse_begin(); j != res[i].sparse_end();)
    /*
    for (std::size_t j = 0; j < til::size(res[i]); ++j)
    {
      if (labels[j] != (unsigned char)(1)) res[i][j] = 0.0;
    }
    */
    for (til::sparse_vector<double>::Map::iterator j = res[i].getMap().begin(); j != res[i].getMap().end();)
    {
      if ((*labels)[j->first] != (unsigned char)(1))
      {
        res[i].getMap().erase(j++);
      }
      else
      {
        ++j;
      }
    }
  }
  std::cout << "OK" << std::endl;
  /*
  {
    til::sparse_vector<double>::Map::const_iterator i = res[0]->getMap().begin();
    for (; i != res[0]->getMap().end(); ++i)
    {
      std::cout << i->first << " " << i->second << std::endl;
    }
  }
  */
  return res;
}

template < typename TMesh >
std::vector<shared_ptr<til::sparse_vector<double> > >
getDistToNeighbors(const TMesh & mesh, const std::vector<std::vector<std::size_t> > & nindex)
{
  std::cout << "Computing geomap..." << std::flush;
  std::size_t n = size(getVertices(mesh));
  typedef til::ghost::GMapStop_MyPointsDone<std::vector<std::size_t>::const_iterator> Plugin;
  //Plugin stopGhost(nindex.begin(), nindex.end());
  //til::Triangle_mesh_geodesic_map<TMesh, double, Plugin, til::policy::GMap_DefaultStorage<til::sparse_vector, TMesh, double > > 
  //til::Graph_distance_map<TMesh, double, Plugin, til::policy::GMap_DefaultStorage<til::sparse_vector, TMesh, double > > 
  //  geomap(mesh);
  typedef std::vector<std::vector<std::size_t> > CNeighborhoods;
  shared_ptr<CNeighborhoods> pneighc = til::circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  til::Graph_distance_map<typename TMesh::VertexCollection, CNeighborhoods, double, Plugin, til::policy::GMap_DefaultStorage<til::sparse_vector, double > > 
    geomap(getVertices(mesh), *pneighc);
  std::vector<std::size_t> startPoints(1);
  std::vector<double> dist(1, 0.0);
  //std::vector<double> res(n, 0.0);
  std::vector<shared_ptr<til::sparse_vector<double> > > res(n);
  for (std::size_t i=0; i < n; ++i)
  {
    startPoints[0] = i;
    geomap.init(startPoints, dist);
    geomap.stopGhost().init(nindex[i].begin(), nindex[i].end());
    geomap.process();
    res[i] = geomap.distanceMap();
  }
  std::cout << "OK" << std::endl;
  return res;
}

template < typename TMesh >
shared_ptr<til::sparse_vector<double> >
getDistToNeighbor
(
 const TMesh & mesh,
 const shared_ptr<std::vector<std::vector<std::size_t> > > & neighc,
 const std::vector<std::vector<std::size_t> > & nindex,
 std::size_t i)
{
  //std::cout << "Computing geomap..." << std::flush;
  std::size_t n = size(getVertices(mesh));
  typedef til::ghost::GMapStop_MyPointsDone<std::vector<std::size_t>::const_iterator> Plugin;
  //Plugin stopGhost(nindex.begin(), nindex.end());
  //til::Triangle_mesh_geodesic_map<TMesh, double, Plugin, til::policy::GMap_DefaultStorage<til::sparse_vector, TMesh, double > > 
  til::Graph_distance_map<typename TMesh::VertexCollection, std::vector<std::vector<std::size_t> >, double, Plugin, til::policy::GMap_DefaultStorage<til::sparse_vector, double > >
    geomap(getVertices(mesh), *neighc);
  std::vector<std::size_t> startPoints(1);
  std::vector<double> dist(1, 0.0);
  //std::vector<double> res(n, 0.0);
  shared_ptr<til::sparse_vector<double> > res(n);
  //for (std::size_t i=0; i < n; ++i)
  {
    startPoints[0] = i;
    geomap.init(startPoints, dist);
    geomap.stopGhost().init(nindex[i].begin(), nindex[i].end());
    geomap.process();
    res = geomap.distanceMap();
  }
  //std::cout << "OK" << std::endl;
  return res;
}

template < typename TMesh >
shared_ptr<std::vector<double> >
getDistToNeighbor2
(
 const TMesh & mesh,
 const shared_ptr<std::vector<std::vector<std::size_t> > > & neighc,
 const std::vector<std::vector<std::size_t> > & nindex,
 std::size_t i)
{
  typedef til::ghost::GMapStop_MyPointsDone<std::vector<std::size_t>::const_iterator> Plugin;
  //Plugin stopGhost(nindex.begin(), nindex.end());
  //til::Triangle_mesh_geodesic_map<TMesh, double, Plugin, til::policy::GMap_DefaultStorage<til::sparse_vector, TMesh, double > > 
  //til::Graph_distance_map<TMesh, double, Plugin, til::policy::GMap_DefaultStorage<til::sparse_vector, TMesh, double > > 
  til::Graph_distance_map<typename TMesh::VertexCollection, std::vector<std::vector<std::size_t> >, double, Plugin  >
    geomap(getVertices(mesh), *neighc);
  std::vector<std::size_t> startPoints(1);
  std::vector<double> dist(1, 0.0);
  //std::vector<double> res(n, 0.0);
  shared_ptr<std::vector<double> > res;
  //for (std::size_t i=0; i < n; ++i)
  {
    startPoints[0] = i;
    geomap.init(startPoints, dist);
    geomap.stopGhost().init(nindex[i].begin(), nindex[i].end());
    geomap.process();
    res = geomap.distanceMap();
  }
  //std::cout << "OK" << std::endl;
  return res;
}

void indicator(std::size_t i, std::size_t n, double parts = 10.0)
{
  if (i % int(n/parts) == 0)
  {
    std::cout << int(parts * int(i / int(n/parts))) << "%..." << std::flush;
  }
}


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
    typename til::MeshTraits<TMesh>::Vertex tmp;
    typename til::MeshTraits<TMesh>::FaceIndexCollection::const_iterator iFaceIndex;
    for (iFaceIndex = til::getFaceIndices(mesh).begin(); iFaceIndex != til::getFaceIndices(mesh).end(); ++iFaceIndex)
    {
      tmp = til::getFaceVertex(mesh, iFaceIndex, 0);
    }
    }));
  }
  {
    PRINT_TIME3("standard4", ({
    typename til::MeshTraits<TMesh>::Vertex tmp;
    typename til::MeshTraits<TMesh>::FaceIndexCollection::const_iterator iFaceIndex;
    for (iFaceIndex = til::getFaceIndices(mesh).begin(); iFaceIndex != til::getFaceIndices(mesh).end(); ++iFaceIndex)
    {
      tmp = til::getFaceVertex(mesh, iFaceIndex, 0);
    }
    }));
  }
  {
    PRINT_TIME3("superfast?", ({
    typename til::MeshTraits<TMesh>::Vertex tmp;
    typename til::MeshTraits<TMesh>::FaceIndexCollection::const_iterator iFaceIndex = til::getFaceIndices(mesh).begin();
    typename til::MeshTraits<TMesh>::FaceIndexCollection::const_iterator iFaceIndexEnd = til::getFaceIndices(mesh).end();
    for (; iFaceIndex != iFaceIndexEnd; ++iFaceIndex)
    {
      tmp = til::getFaceVertex(mesh, iFaceIndex, 0);
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
    til::MeshTraits<AimsSurfaceTriangle>::FaceIndexCollection::const_iterator iFaceIndex;
    //for (iFaceIndex = getFaceIndices(*s).begin(); iFaceIndex != getFaceIndices(*s).end(); ++iFaceIndex)
    for (iFaceIndex = til::getFaceIndices(*s).begin(); iFaceIndex != s->polygon().end(); ++iFaceIndex)
    {
      tmp = til::getVertices(s2)[(*iFaceIndex)[0]];
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
      tmp = til::getFaceVertex(s2, iFaceIndex, 0);
      //std::cout << tmp << std::endl;
    }
    }));
  }

  /*
  {
    til::MeshTraits<AimsSurface<3,Void> >::Vertex tmp;
    PRINT_TIME3("traits", ({    
    til::MeshTraits<AimsSurface<3,Void> >::FaceCollection::const_iterator iFaces = til::getFaces(*s).begin();
    for ( ; iFaces != til::getFaces(*s).end(); ++iFaces)
    {
      tmp = (*iFaces)[0];      
      //std::cout << tmp << std::endl;
    }
    }));
  }
  */
  
/////////////////////////////////////////////////////////////////////////
// MESH1
/////////////////////////////////////////////////////////////////////////


  std::cout << "til::Mesh1" << std::endl;

  {
    til::Mesh1 mesh;
    //const AimsSurface<3,Void> &ss = *s;
    //til::getFaceIndices(ss);
    til::convert(mesh, *s);
  
    {
      PRINT_TIME3("standard", ({
      til::MeshTraits<til::Mesh1>::Vertex tmp;
      std::vector<til::numeric_array<std::size_t, 3> >::iterator iFaceIndex;
      for (iFaceIndex = mesh.getFaceIndices().begin(); iFaceIndex != mesh.getFaceIndices().end(); ++iFaceIndex)
      {
        tmp = mesh.getVertices()[(*iFaceIndex)[0]];
        //std::cout << tmp << std::endl;
      }
      }));
    }  
  
    {
      PRINT_TIME3("standard2", ({
      til::MeshTraits<til::Mesh1>::Vertex tmp;
      til::MeshTraits<til::Mesh1>::FaceIndexCollection::const_iterator iFaceIndex;
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
      til::MeshTraits<til::Mesh1>::Vertex tmp;
      std::vector<til::numeric_array<std::size_t, 3> >::iterator iFaceIndex;
      for (iFaceIndex = mesh.getFaceIndices().begin(); iFaceIndex != mesh.getFaceIndices().end(); ++iFaceIndex)
      {
        tmp = mesh.getVertices()[(*iFaceIndex)[0]];
        //std::cout << tmp << std::endl;
      }
      }));
    }  
  
    {
      PRINT_TIME3("superfast?", ({
      til::MeshTraits<til::Mesh1>::Vertex tmp;
      til::MeshTraits<til::Mesh1>::FaceIndexCollection::iterator iFaceIndex;
      til::MeshTraits<til::Mesh1>::FaceIndexCollection & face = getFaceIndices(mesh);
      til::MeshTraits<til::Mesh1>::VertexCollection & vertex = getVertices(mesh);
      for (iFaceIndex = face.begin(); iFaceIndex != face.end(); ++iFaceIndex)
      {
        tmp = vertex[(*iFaceIndex)[0]];
        //std::cout << tmp << std::endl;
      }
      }));
    }  
    {
      PRINT_TIME3("standard4", ({
      til::MeshTraits<til::Mesh1>::Vertex tmp;
      std::vector<til::numeric_array<std::size_t, 3> >::iterator iFaceIndex;
      for (iFaceIndex = mesh.getFaceIndices().begin(); iFaceIndex != mesh.getFaceIndices().end(); ++iFaceIndex)
      {
        tmp = getFaceVertex(mesh, iFaceIndex, 0);
        //std::cout << tmp << std::endl;
      }
      }));
    }
  
    /*
    {
      til::MeshTraits<til::Mesh1>::Vertex tmp;
      PRINT_TIME3("traits", {    
      til::MeshTraits<til::Mesh1>::FaceCollection::const_iterator iFaces = getFaces(mesh).begin();
      for ( ; iFaces != getFaces(mesh).end(); ++iFaces)
      {
        tmp = (*iFaces)[0];
        //std::cout << tmp << std::endl;
      }
      });
    }
    */
    
    {
      const int N = 10000;
      til::numeric_array<double, N> p;
      til::numeric_array<double, N> v;
      std::fill(p.begin(), p.end(), 1.23);
      std::fill(v.begin(), v.end(), 3.45);
      PRINT_TIME3("addto", p += v);
      PRINT_TIME3("direct", ({
      til::numeric_array<double,N>::iterator iP = p.begin();
      til::numeric_array<double,N>::const_iterator iV = v.begin();
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

  std::cout << "til::Mesh2" << std::endl;

  {
    til::Mesh2 mesh;
    til::convert(mesh, *s);
  
    {
      PRINT_TIME3("standard", ({
      til::MeshTraits<til::Mesh2>::Vertex tmp;
      til::MeshTraits<til::Mesh2>::FaceIndexCollection::iterator iFaceIndex;
      for (iFaceIndex = mesh.getFaceIndices().begin(); iFaceIndex != mesh.getFaceIndices().end(); ++iFaceIndex)
      {
        tmp = *((*iFaceIndex)[0]);
        //std::cout << tmp << std::endl;
      }
      }));
    }  
  
    {
      PRINT_TIME3("standard2", ({
      til::MeshTraits<til::Mesh2>::Vertex tmp;
      til::MeshTraits<til::Mesh2>::FaceIndexCollection::const_iterator iFaceIndex;
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
      til::MeshTraits<til::Mesh2>::Vertex tmp;
      til::MeshTraits<til::Mesh2>::FaceIndexCollection::iterator iFaceIndex;
      for (iFaceIndex = mesh.getFaceIndices().begin(); iFaceIndex != mesh.getFaceIndices().end(); ++iFaceIndex)
      {
        tmp = *((*iFaceIndex)[0]);
        //std::cout << tmp << std::endl;
      }
      }));
    }  
  
    {
      PRINT_TIME3("superfast?", ({
      til::MeshTraits<til::Mesh2>::Vertex tmp;
      til::MeshTraits<til::Mesh2>::FaceIndexCollection::iterator iFaceIndex;
      til::MeshTraits<til::Mesh2>::FaceIndexCollection & face = getFaceIndices(mesh);
      //til::MeshTraits<til::Mesh2>::VertexCollection & vertex = getVertices(mesh);
      for (iFaceIndex = face.begin(); iFaceIndex != face.end(); ++iFaceIndex)
      {
        tmp = *((*iFaceIndex)[0]);
        //std::cout << tmp << std::endl;
      }
      }));
    }  

    {
      PRINT_TIME3("standard4", ({
      til::MeshTraits<til::Mesh2>::Vertex tmp;
      til::MeshTraits<til::Mesh2>::FaceIndexCollection::iterator iFaceIndex;
      for (iFaceIndex = mesh.getFaceIndices().begin(); iFaceIndex != mesh.getFaceIndices().end(); ++iFaceIndex)
      {
        tmp = getFaceVertex(mesh, iFaceIndex, 0);
        //std::cout << tmp << std::endl;
      }
      }));
    }
  
    /*
    {
      til::MeshTraits<til::Mesh2>::Vertex tmp;
      PRINT_TIME3("traits", {    
      til::MeshTraits<til::Mesh2>::FaceCollection::const_iterator iFaces = getFaces(mesh).begin();
      for ( ; iFaces != getFaces(mesh).end(); ++iFaces)
      {
        tmp = (*iFaces)[0];
        //std::cout << tmp << std::endl;
      }
      });
    }
    */
    
    {
      const int N = 10000;
      til::numeric_array<double, N> p;
      til::numeric_array<double, N> v;
      std::fill(p.begin(), p.end(), 1.23);
      std::fill(v.begin(), v.end(), 3.45);
      PRINT_TIME3("addto", p += v);
      PRINT_TIME3("direct", ({
      til::numeric_array<double,N>::iterator iP = p.begin();
      til::numeric_array<double,N>::const_iterator iV = v.begin();
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
  til::Mesh1 mesh;
  til::convert(mesh, s);
  til::Mesh2 mesh2;
  til::convert(mesh2, s);
  _gra(s,mesh,mesh2);
}


void testCentroid2()
{
  til::Mesh1 mesh;
  {
    std::cout << "Starting here" << std::endl;
    AimsSurfaceTriangle *s = makeSphere(Point3df(1.0, 2.0, 3.0), 0.7, 0);
    std::cout << "zero" << std::endl;
    til::convert(mesh, *s);
  }

  std::cout << "one" << std::endl;

  til::MeshTraits<til::Mesh1>::Vertex centr;
  centroid(getVertices(mesh), centr);

  std::cout << "two" << std::endl;

  std::cout << "Centroid : " << centr << " (should be 1.0, 2.0, 3.0 )" << std::endl;
  
  til::Mesh_N mesh2 = addNeighborsToMesh(mesh);

  std::cout << "three" << std::endl;
  
  std::cout << getVertices(mesh2)[0] << std::endl;
  std::cout << getVertices(mesh)[0] << std::endl;
  std::cout << &getVertices(mesh2) << std::endl;
  std::cout << &getVertices(mesh) << std::endl;

  til::MeshTraits<til::Mesh_N>::NeighborIndexCollection::const_iterator iNic = getNeighborIndices(mesh2).begin();
  til::MeshTraits<til::Mesh_N>::NeighborIndex::const_iterator iNi = iNic->begin();
  std::cout << "Voici" << std::endl;
  for (; iNi != iNic->end(); ++iNi)
  {
    std::cout << *iNi << std::endl;
  }
  til::SpringForce<til::Mesh_N, til::numeric_array<float,3> > sf;
  sf.initializeLengths(mesh2);
  std::vector<til::numeric_array<float,3> > f;
  sf.getForces(mesh2, f);
  std::cout << f[0] << std::endl;
}

template < typename T >
bool foo(T)
{
  return boost::is_same<typename til::value_type_of<T>::type, std::size_t>::value;
}

template < typename TMesh >
//void meshStats(til::Mesh2_N & mesh)
void meshStats(TMesh & mesh)
{
  int nVertices = getVertices(mesh).size();
  std::cout << "Number of vertices : " << nVertices << std::endl;
  int nFaces = til::size(getFaceIndices(mesh));
  std::cout << "Number of faces: " << nFaces << std::endl;

  int nEdges;
  
  {
    shared_ptr<std::vector<std::pair<std::size_t, std::size_t> > > pedges
      = til::faces2edges(til::getFaceIndices(mesh));
    nEdges = pedges->size();
    //typedef til::numeric_array<typename TMesh::FaceIndex::value_type, 2> Edge;
    //std::set<Edge, til::Lexicographical_compare<Edge> > e = getEdges(mesh);
    //nEdges = til::size(e);
    std::cout << "Number of edges: " << nEdges << std::endl;
  }
  
  std::cout << "Euler number: " << nVertices + nFaces - nEdges << std::endl;
  
  {
    std::cout << "Histogram: number of neighbors" << std::endl;
    //std::map<std::size_t, int> h = histogram(getNeighborIndices(mesh), functor::Size<til::MeshTraits<til::Mesh2_N>::NeighborIndexCollection>());
    til::MapHistogram<std::size_t> h;
    foo(til::func_it<til::functor::Size>(getNeighborIndices(mesh).begin()));
    h.accumulate(til::func_it<til::functor::Size>(getNeighborIndices(mesh).begin()), 
                 til::func_it<til::functor::Size>(getNeighborIndices(mesh).end()));
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


/*
namespace til { namespace functor {

    template <  >
    class CastTo<AimsSurface<3,::Void>, til::Mesh2_N >
//     : public std::binary_function<AimsSurface<3, ::Void> &, const Mesh<TParam> &, void>
    {
    public:
      void
      //operator()(AimsSurface<3,Void> & aimsMesh, const Mesh<TParam> & mesh) const
      operator()(AimsSurface<3,Void> & aimsMesh, const Mesh2_N & mesh) const
      {
        //if (boost::is_pointer<typename TParam::FaceIndex::value_type>::value)
          detail::convert_mesh_3(mesh, aimsMesh);
        //else
        //  detail::convert_mesh_1(mesh, aimsMesh);
      }
    };        
}}
*/

void myinflate(int argc, char* argv[])
{
  AimsApplication app( argc, aims_const_hack(argv), "myinflate" );
  aims::Reader<AimsSurfaceTriangle> reader;
  aims::Writer<AimsSurfaceTriangle> writer;
  aims::Writer<Texture1d> wt;
  typedef til::Mesh2_N MeshType;
  MeshType mesh;
  float kSpring = 0.2;
  float kCentroid = 0.05;
  //float kLapl = 0.01;
  int niter = 10000;
  int dump = 50;

  AimsSurfaceTriangle wmesh;
  // Read and convert input mesh
  {
    //til::Mesh1 mesh0;
    til::Mesh2 mesh0;
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
    std::cout << "Read mesh with" << til::size(til::getVertices(aimsmesh)) << " vertices" << std::endl;
    til::convert(mesh0, aimsmesh);
    mesh = addNeighborsToMesh(mesh0);
    std::cout << til::getVertices(aimsmesh)[0] << std::endl;
  }
  
  //meshStats(mesh);
  
  Texture1d tex;
  
  // write texture
  /*
  {
    tex.reserve(til::size(getVertices(mesh)));
    for( unsigned i=0; i<til::size(getVertices(mesh)); ++i )
      tex.push_back( i==0 );
    wt.write( tex );
  }
  */
  
  // initialize spring force functor with initial lengths
  til::SpringForce<MeshType, til::numeric_array<float,3> > sf;
  //til::SpringForce<til::Mesh_N, float> sf;
  sf.initializeLengths(mesh);
  
  // Main loop
  til::MeshTraits<MeshType>::Vertex center;
  //til::MeshTraits<til::Mesh_N>::Vertex center;
  typedef til::numeric_array<float,3> Vec3D;
  std::vector<Vec3D> centroidForces(til::size(getVertices(mesh)));
  std::vector<Vec3D> totalForces(til::size(getVertices(mesh)));
  std::vector<Vec3D> springForces(til::size(getVertices(mesh)));
  std::vector<Vec3D> lForces(til::size(getVertices(mesh)));
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
        til::convert(tmp, mesh);
        wmesh[c] = tmp;
        ++c;
      }
      // write texture
      {
        Texture<float> tmp(til::size(getVertices(mesh)));
        for_each_neighbors(mesh, sf.getLengths(), tmp.data(), til::functor::SpringEnergy());
        tex[itex] = tmp;
        ++itex;
      }
    }
    std::cout << "Iteration " << i << std::endl;
    
    // compute centroid forces
    centroid(getVertices(mesh), center);
    {    
      til::MeshTraits<MeshType>::VertexCollection::const_iterator iVc = getVertices(mesh).begin();
      std::vector<til::numeric_array<float,3> >::iterator iC = centroidForces.begin();
      for (; iVc != getVertices(mesh).end(); ++iVc, ++iC)
      {
        *iC = (*iVc - center) * (1.0f/til::norm<float>(*iVc - center));
      }
    }

    // compute spring forces
    {
      Vec3D z; z[0] = z[1] = z[2] = 0;
      std::fill(springForces.begin(), springForces.end(), z);
    }
    sf.getForces(mesh, springForces);
    //std::cout << "sf " << springForces[0] << std::endl;

    // compute laplacian forces
    //laplacianForces(mesh, lForces);

    /*
    std::cout << til::size(getNeighborIndices(mesh)) << "*" << std::flush;
    std::cout << til::size(getVertices(mesh)) << "*" << std::flush;
    for (int j = 0; j < til::size(getNeighborIndices(mesh)); ++j)
    {
      std::cout << til::size(getNeighborIndices(mesh)[j]) << " " << std::flush;
    }
    */
    
    /*
    
    // hand-compute spring force for first voxel
    std::vector<til::Vector<float,3> > f(til::size(getVertices(mesh)));
    til::MeshTraits<til::Mesh_N>::VertexCollection v = getVertices(mesh);
    for (int j = 0; j < til::size(getVertices(mesh)); ++j)
    {
      for (std::size_t k = 0; k < til::size(getNeighborIndices(mesh)[j]); ++k)
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
      til::MeshTraits<MeshType>::VertexCollection::iterator iVc = getVertices(mesh).begin();
      std::vector<til::numeric_array<float,3> >::const_iterator iC = centroidForces.begin();
      std::vector<til::numeric_array<float,3> >::const_iterator iS = springForces.begin();
      //std::vector<til::Vector<float,3> >::const_iterator iL = lForces.begin();
      //for (; iVc != getVertices(mesh).end(); ++iVc, ++iC, ++iS, ++iL)
      float f = exp(-i/1000.0);
      for (; iVc != getVertices(mesh).end(); ++iVc, ++iC, ++iS)
      //std::vector<til::Vector<float,3> >::const_iterator iF = f.begin();
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

  std::cout << "WMesh til::size : " << wmesh.size() << std::endl;
  // write result
  writer.write(wmesh);
  wt.write( tex );
}



void testPush()
{
  til::Mesh1 mesh;
  {
    std::cout << "Starting here" << std::endl;
    AimsSurfaceTriangle *s = makeSphere(Point3df(1.0, 2.0, 3.0), 0.7, 7);
    til::convert(mesh, *s);
  }

  til::MeshTraits<til::Mesh1>::Vertex centr;
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
  AimsApplication app( argc, aims_const_hack(argv), "testConvert" );
  aims::Reader<AimsSurfaceTriangle> reader;
  aims::Writer<AimsSurfaceTriangle> writer;
  aims::Writer<Texture1d> wt;
  til::Mesh2_N mesh;
  //til::Mesh_N mesh;

  AimsSurfaceTriangle wmesh;
  // Read and convert input mesh
  {
    //til::Mesh1 mesh0;
    til::Mesh2 mesh0;
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
    std::cout << "Read mesh with" << til::size(til::getVertices(aimsmesh)) << " vertices" << std::endl;    
    
    std::cout << "First face" << std::endl;
    for (int i = 0; i < 3; ++i)
    {
      std::cout << til::getFaceIndices(aimsmesh)[0][i] << " " << til::getVertices(aimsmesh)[til::getFaceIndices(aimsmesh)[0][i]] << " " << &(til::getVertices(aimsmesh)[til::getFaceIndices(aimsmesh)[0][i]]) << std::endl;
    }
    std::cout << "size " << sizeof(til::MeshTraits<AimsSurfaceTriangle>::Vertex) << std::endl;
    til::convert(mesh0, aimsmesh);
    
    std::cout << "First face" << std::endl;
    for (int i = 0; i < 3; ++i)
    {
      std::cout << *(getFaceIndices(mesh0)[0][i]) << " " << (getFaceIndices(mesh0)[0][i]) << std::endl;
    }    
    mesh = addNeighborsToMesh(mesh0);
    std::cout << "size " << sizeof(til::MeshTraits<til::Mesh2_N>::Vertex) << std::endl;

    std::cout <<(getFaceIndices(mesh)[0][1] - &(getVertices(mesh)[0])) << std::endl;
  }
  std::cout << "WMesh size : " << wmesh.size() << std::endl;
  // write result

  AimsSurface<3,Void> tmp;
  til::convert(tmp, mesh);
  wmesh[0] = tmp;

  std::cout << "First face" << std::endl;
  for (int i = 0; i < 3; ++i)
  {
    std::cout << til::getFaceIndices(wmesh)[0][i] << " " << til::getVertices(wmesh)[til::getFaceIndices(wmesh)[0][i]] << " " << &(til::getVertices(wmesh)[til::getFaceIndices(wmesh)[0][i]]) << std::endl;
  }
  std::cout << "size " << sizeof(til::MeshTraits<AimsSurface<3,Void> >::Vertex) << std::endl;
  
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

  til::Mesh1 mesh;
  std::cout << "Aims -> mesh" << std::endl;
  {
    AimsSurfaceTriangle s;
    r.read( s );
    til::convert(mesh, s);
  }
  std::cout << "Mesh -> aims" << std::endl;
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh);
    w.write(s);
  }  
}


void testKDTreeSerious(int argc, char * argv[])
{
  Reader<AimsSurfaceTriangle> r;
  AimsApplication app( argc, aims_const_hack(argv), "testMyMesh" );
  app.addOption(r, "-i", "input mesh" );
  app.initialize();

  til::Mesh2 mesh;
  std::cout << "Aims -> mesh" << std::endl;
  {
    AimsSurfaceTriangle s;
    r.read( s );
    til::convert(mesh, s);
  }

  //boost::uniform_real<float> random(-150,150);    
  til::numeric_array<double,3> center;
  centroid(getVertices(mesh), center);
  std::cout << "Centroid: " << center << std::endl;
  til::numeric_array<double,3> std;
  stdev(getVertices(mesh), std);
  std::cout << "Stdev: " << std << std::endl;
  //std::vector<til::Point<float,3>*> res;
  typedef til::KDTree<til::numeric_array<float,3>*, til::MeshTraits<til::Mesh2>::VertexCollection> MyKDTree;
  MyKDTree res;
  makeKDTree(getVertices(mesh), res);
  std::cout << "kdtree.size " << til::size(res) << std::endl;

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
    til::numeric_array<float,3> target = til::numeric_array<float,3>
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
    til::Find_closest<double, MyKDTree> fc(res);
    fc(target);
    gettimeofday(&finish, &tz);
    //im(i,j,k) = (double)((finish.tv_sec-start.tv_sec) * 1000000L + (finish.tv_usec-start.tv_usec));
    im(i,j,k) = fc.niter();
    //timeTaken_kd[i] = (double)((finish.tv_sec-start.tv_sec) * 1000000L + (finish.tv_usec-start.tv_usec));
    //niter_kd[i] = fc.niter();
  
    /*
    til::Point<float,3> const * p;
    gettimeofday(&start, &tz);
    {
      std::vector<til::Point<float,3> >::const_iterator iVertex = getVertices(mesh).begin();
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
  
  til::aimswrite(im, "/home/cathier/tmp/out");
}

template < typename T >
T find_closest(const std::vector<T> & c, T v)
{
  T res;
  double dist = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < til::size(c); ++i)
  {
    if (til::dist2(c[i], v, til::prec<double>()) < dist)
    {
      dist = til::dist2(c[i], v, til::prec<double>());
      res = c[i];
    }
  }
  return res;
}


/*
 * OOps, code is now broken because of deepCopy
void testKDTreeSerious3(int argc, char * argv[])
{
  
  Reader<AimsSurfaceTriangle> r;
  Writer<AimsSurfaceTriangle> w;
  til::Mesh2_N mesh;
  {
    AimsApplication app( argc, aims_const_hack(argv), "testKDTreeSerious3" );
    app.addOption(r, "-i", "input mesh" );
    app.addOption(w, "-o", "output mesh" );
    app.initialize();
    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh2 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  til::Mesh2_N mesh2;
  //til::deepCopy(mesh, mesh2);
  mesh2.deepCopy(mesh);
  std::cout << til::size(getVertices(mesh)) << " " << til::size(getVertices(mesh2)) << std::endl;
  std::cout << til::size(getFaceIndices(mesh)) << " " << til::size(getFaceIndices(mesh2)) << std::endl;
  
  til::Triangle_mesh_geodesic_map<til::Mesh2_N, double > geomap(mesh);
  std::vector<til::MeshTraits<til::Mesh2_N>::FaceIndex::value_type> startPoints;
  startPoints.push_back(getFaceIndices(mesh)[0][0]);
  std::vector<double> dist(1, 0.0);
  geomap.init(startPoints, dist);
  geomap.process();
  shared_ptr<std::vector<double> > res = geomap.distanceMap();
  {
    Texture1d t;
    til::convert(t, *res);
    Writer<Texture1d> w("/home/cathier/tmp/texgeo");
    w.write(t);
  }  
  
  for (std::size_t i = 0; i < til::size(getVertices(mesh)); ++i)
  {
    getVertices(mesh2)[i] += til::Vector<float,3>(20.0, 20.0, 20.0); 
  }  
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh2);
    w.write(s);
  }

  typedef til::KDTree<std::size_t, til::MeshTraits<til::Mesh2_N>::VertexCollection> MyKDTree;
  MyKDTree kdt(getVertices(mesh));
  makeKDTree(getVertices(mesh), kdt);
  
  std::vector<double> res2(til::size(*res), 0);
  til::Mesh2_N::VertexCollection::iterator iVertex = getVertices(mesh2).begin();
  int count = 0;
  for (; iVertex != getVertices(mesh2).end(); ++iVertex)
  {
    til::Find_closest< double, MyKDTree > fc(kdt);
    //std::cout << ++count << " " << std::flush;
    std::size_t i = fc(*iVertex);
    / *
    Point<float,3> p = find_closest(getVertices(mesh), *iVertex);
    if (p != getVertices(mesh)[i])
    {
      std::cout << "ERROR" << std::endl;
      std::cout << p << std::endl;
      std::cout << getVertices(mesh)[i] << std::endl;
      exit(1);
    }
    * /
    //res2[getVertexNumber(mesh2, &*iVertex)] = res[i];
    res2[count] = (*res)[i];
    ++count;
  }
  std::cout << "toto2" << std::endl;
  {
    Texture1d t;
    til::convert(t, res2);
    Writer<Texture1d> w("/home/cathier/tmp/texgeo2");
    w.write(t);
  }  
}
*/

void testKDTreeSerious2(int argc, char * argv[])
{
  Reader<AimsSurfaceTriangle> r;
  AimsApplication app( argc, aims_const_hack(argv), "testKDTreeSerious2" );
  app.addOption(r, "-i", "input mesh" );
  app.initialize();

  til::Mesh2 mesh;
  std::cout << "Aims -> mesh" << std::endl;
  {
    AimsSurfaceTriangle s;
    r.read( s );
    til::convert(mesh, s);
  }

  /*
  std::vector<Point<float,3>*> res;
  PRINT_TIME3("Building kdtree, simple", ({
  makeKDTree(getVertices(mesh), res);
  }));  
  std::cout << "kdtree.size " << res.size() << std::endl;
  */
  
  til::KDTree<til::numeric_array<float,3>*, til::MeshTraits<til::Mesh2>::VertexCollection> kdt;
  PRINT_TIME3("Building kdtree, full", ({
  makeKDTree(getVertices(mesh), kdt);
  }));

  /*
  PRINT_TIME3("Searching, simple", ({
  til::Mesh2::VertexCollection::iterator iVertex = getVertices(mesh).begin();
  for (; iVertex != getVertices(mesh).end(); ++iVertex)
  {
    til::Point<float,3> target;// = *iVertex + til::Vector<float,3>(1.0,1.0,1.0);
    add(*iVertex, til::Vector<float,3>(1.0,1.0,1.0), target);
    //std::cout << target << std::endl;

    til::Find_closest<double, til::Point<float,3> > fc(res);
    fc(target);    
  }
  }));
  */
  
  
  PRINT_TIME3("Searching, full", ({
  til::Mesh2::VertexCollection::iterator iVertex = getVertices(mesh).begin();
  for (; iVertex != getVertices(mesh).end(); ++iVertex)
  {
    til::numeric_array<float,3> target;// = *iVertex + til::Vector<float,3>(1.0,1.0,1.0);
    //add(*iVertex, til::Vector<float,3>(1.0,1.0,1.0), target);
    target = *iVertex + til::numeric_array<float,3>(1.0,1.0,1.0);
    //std::cout << target << std::endl;

    til::Find_closest<double, til::KDTree<til::numeric_array<float,3>*, til::MeshTraits<til::Mesh2>::VertexCollection> > fc(kdt);
    fc(target);
  }
  }));
  
  //std::cout << "out" << std::endl;
}


/*
void testKDTree()
{
  std::vector<til::Vector<float,3> > v;
  
  v.push_back(til::Vector<float,3>(8,7,9));
  v.push_back(til::Vector<float,3>(28,7,7));
  v.push_back(til::Vector<float,3>(18,7,1));
  v.push_back(til::Vector<float,3>(-8,7,0));
  v.push_back(til::Vector<float,3>(9,27,3));
  v.push_back(til::Vector<float,3>(8,37,5));
  v.push_back(til::Vector<float,3>(2,17,2));
  v.push_back(til::Vector<float,3>(4,-7,2));
  v.push_back(til::Vector<float,3>(8,5,22));
  v.push_back(til::Vector<float,3>(4,3,23));
  
  til::Vector<float,3> v0(120, -2, 1);
  //std::vector<til::Vector<float,3>*> res = makeKDTree(v);
  {
    std::vector<std::size_t> res;
    makeKDTree(v, res);  
    printTree(res);
    std::cout << "Closerttil::Point: ";
    //std::cout << find_closest<double>(res, v0) << std::endl;
  }
  {
    std::vector<til::Vector<float,3>*> res;
    makeKDTree(v, res);  
    printTree(res);
    std::cout << "Closerttil::Point: ";
    std::cout << find_closest<double>(res, v0) << std::endl;
  }
  
  {
    til::KDTree<til::Vector<float,3>*, std::vector<til::Vector<float,3> > > k;
    makeKDTree(v, k);
    //makeKDTree<til::Vector<float,3> >(v, k);
    print(k);
    std::cout << "ClosertPoint: " << std::flush;
    std::cout << find_closest<double>(k, v0) << std::endl;
  }
}

*/

void testPrintStats2(int argc, char * argv[])
{
  AimsApplication app( argc, aims_const_hack(argv), "testPrintStats2" );
  aims::Reader<AimsSurfaceTriangle> reader;
  typedef til::Mesh2_N MeshType;
  MeshType mesh;

  AimsSurfaceTriangle wmesh;
  // Read and convert input mesh
  {
    //til::Mesh1 mesh0;
    til::Mesh2 mesh0;
    app.addOption(reader, "-i", "input mesh" );
    app.initialize();
    
    AimsSurfaceTriangle aimsmesh;
    std::cout << "Reading mesh..."<< std::flush;
    reader.read(aimsmesh);
    std::cout << "Done" << std::endl;
    std::cout << "Read mesh with " << til::size(til::getVertices(aimsmesh)) << " vertices" << std::endl;
    til::convert(mesh0, aimsmesh);
    std::cout << "Adding neighbors" << std::endl;
    mesh = addNeighborsToMesh(mesh0);
  }
 
 std::cout << "stats" << std::endl;
  
  //meshStats(mesh);
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
  til::NaryTree<double, 2> bt;
  til::NaryTree<double, 2>::iterator i;

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
  til::numeric_array<double,3> v1, v2(23.45, 287.9823, 374397924.2938), v3(12098409834.23, 323.283, 0.000000001);
  for (int i = 0; i < 100000; ++i)
  {
    v1 = v2 + v3;
  }
}

void testAdd2()
{
  til::numeric_array<double,3> v1, v2(23.45, 287.9823, 374397924.2938), v3(12098409834.23, 323.283, 0.000000001);
  for (int i = 0; i < 100000; ++i)
  {
    //v1 = v2 + v3;
    // formerly add(v2,v3,v1)
  }
}

void testAdd3()
{
  til::numeric_array<double,3> v1, v2(23.45, 287.9823, 374397924.2938), v3(12098409834.23, 323.283, 0.000000001);
  for (int i = 0; i < 100000; ++i)
  {
    til::numeric_array<double,3> v(v2);
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

/*
void testTILFunctors()
{
  using namespace til;
  std::vector<double> a(10,0);
  std::list<double> b(10,3);
  
  loop1(*_1 = *_2 **_2, a, b);
  
  std::cout << "*" << a[0] << std::endl;
}
*/

/*
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
*/

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
  
  std::cout << "Reading image..." << std::flush;
  AimsData<short> im;
  r.read(im);
  std::cout << "OK" << std::endl;
  
  std::cout << "Reading mesh..." << std::flush;
  AimsSurfaceTriangle mesh;
  rm.read(mesh);
  std::cout << "OK" << std::endl;

  std::cout << "Clear image..." << std::flush;
  // clear image
  for (int i = 0; i < im.dimX()*im.dimY()*im.dimZ(); ++i)
  {
    im[i] = 0;
  }
  std::cout << "OK" << std::endl;
  
  
  // Set mesh trace
  std::cout << "Computing trace..." << std::flush;
  til::MeshTraits<AimsSurfaceTriangle>::VertexCollection vertices = til::getVertices(mesh);
  til::MeshTraits<AimsSurfaceTriangle>::VertexCollection::const_iterator iVertex = vertices.begin();
  std::size_t count = 0;
  for (; iVertex != vertices.end(); ++iVertex, ++count)
  {
    indicator(count, vertices.size());
    
    im(int((*iVertex)[0]/im.sizeX()),
       int((*iVertex)[1]/im.sizeY()),
       int((*iVertex)[2]/im.sizeZ())) = 1;
/*    im(int(rint((*iVertex)[0]/im.sizeX())),
       int(rint((*iVertex)[1]/im.sizeY())),
       int(rint((*iVertex)[2]/im.sizeZ()))) = 1;
       */
  }
  std::cout << "OK" << std::endl;

  std::cout << "Writing image..." << std::flush;
  w.write(im);
  std::cout << "OK" << std::endl;  
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

  // Commented out to avoid a warning
/*  
  PRINT_TIME3("Rand", ({
  for (int i = 0; i < 10*N; ++i)
  {
    std::rand() % N;
    std::rand() % N;
  }
  }));
*/

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

/*
void testCircularNeighbors()
{
  AimsSurfaceTriangle *s1 = makeSphere(Point3df(1.0, 2.0, 3.0), 0.7, 7);
  //AimsSurfaceTriangle &s = *s1;
  til::Mesh2_N mesh;
  {
    til::Mesh2 mesh0;
    til::convert(mesh0, *s1);
    std::cout << "Adding neighbors" << std::endl;
    mesh = addNeighborsToMesh(mesh0);
  }

  std::cout << "A" << std::endl;
  std::vector<std::vector<til::numeric_array<float,3>*> > neigh = getNeighborIndices(mesh);
  std::cout << "B" << std::endl;
  shared_ptr<std::vector<std::vector<til::numeric_array<float,3>*> > > neighc
    = circular_neighborhoods(getVertices(mesh), getFaces(mesh));
  std::cout << "C" << std::endl;
    
  std::cout << til::size(*neighc) << std::endl;
  std::cout << til::size(neigh) << std::endl;
  
  for (std::size_t i = 0; i < til::size(*neighc); ++i)
  {
    if (til::size((*neighc)[i])!=til::size(neigh[i]))
    {
      std::cout << "ERROR " << i << std::endl;
    }
  }
}
*/

/*
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
*/

void testsizeof()
{
  std::cout << "double " << sizeof(double) << std::endl;
  std::cout << "long double " << sizeof(long double) << std::endl;
}

void testGeodesicDistance(int argc, char * argv[])
{
  //AimsSurfaceTriangle &s = *s1;
  typedef til::Mesh_N MyMesh;
  MyMesh mesh;
  {
    Reader<AimsSurfaceTriangle> r;
    AimsApplication app( argc, aims_const_hack(argv), "testMyMesh" );
    app.addOption(r, "-i", "input mesh" );
    app.initialize();
    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);

    //AimsSurfaceTriangle *s1 = makeSphere(Point3df(1.0, 2.0, 3.0), 0.7, 7);
    //til::Mesh2 mesh0;
    //convert(*s1, mesh0);
    
    std::cout << "Adding neighbors" << std::endl;
    mesh = addNeighborsToMesh(mesh0);
  }

  std::cout << "Computing geomap" << std::endl;
  til::ghost::GMapStop_AboveThreshold<double> stopGhost(10.0);
  std::cout << "Finishing" << std::endl;

  //til::Triangle_mesh_geodesic_map<til::Mesh2_N, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage<til::sparse_vector, til::Mesh2_N, double > >
    //geomap(mesh, stopGhost);
  typedef std::vector<std::vector<std::size_t> > CNeighborhoods;
  shared_ptr<CNeighborhoods> pneighc = til::circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  til::Triangle_mesh_geodesic_map<MyMesh::VertexCollection, CNeighborhoods, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage<til::sparse_vector, double > >
    geomap(getVertices(mesh), *pneighc, stopGhost);
  std::vector<til::MeshTraits<MyMesh>::FaceIndex::value_type> startPoints;
  startPoints.push_back(getFaceIndices(mesh)[0][0]);
  std::vector<double> dist(1, 0.0);
  geomap.init(startPoints, dist);
  geomap.process();
  shared_ptr<til::sparse_vector<double> > sres = geomap.distanceMap();

  {
    til::sparse_vector<double>::sparse_iterator iRes = sres->sparse_begin();
    for (; iRes != sres->sparse_end(); ++iRes)
    {
      std::cout << iRes->second << " ";
    }
    std::cout << std::endl;
  }
  
  
  std::vector<double> res(til::size(*sres));
  {
    std::vector<double>::iterator iRes = res.begin();
    til::sparse_vector<double>::const_iterator iRes2 = sres->begin();
    for (; iRes2 != sres->end(); ++iRes, ++iRes2)
    {
      *iRes = *iRes2;
    }
  }

  /*  
  til::Triangle_mesh_geodesic_map<til::Mesh2_N, double, ghost::GMapStop_AboveThreshold<double> > geomap(mesh, stopGhost);
  std::vector<til::MeshTraits<til::Mesh2_N>::FaceIndex::value_type> startPoints;
  startPoints.push_back(getFaceIndices(mesh)[0][0]);
  std::vector<double> dist(1, 0.0);
  geomap.init(startPoints, dist);
  geomap.process();
  std::vector<double> res = geomap.distanceMap();
  */
  
  
  for (std::size_t i = 0; i < til::size(res); ++i) res[i] = til::min(res[i], 100.0);
  
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
    til::convert(s, mesh);
    til::aimswrite(s, "sphere");
  }
  {
    Texture1d t; //(1, til::size(res));
    t.reserve(til::size(res));
    for (std::size_t i=0; i<til::size(res); ++i) //t[i] = res[i];
      t.push_back(res[i]);
    til::aimswrite(t, "spheredist");
  }  
}

void testGeodesicDistance2(int argc, char * argv[])
{
  //AimsSurfaceTriangle &s = *s1;
  typedef til::Mesh_N MyMesh;
  MyMesh mesh;
  {
    Reader<AimsSurfaceTriangle> r;
    AimsApplication app( argc, aims_const_hack(argv), "testMyMesh" );
    app.addOption(r, "-i", "input mesh" );
    app.initialize();
    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);

    //AimsSurfaceTriangle *s1 = makeSphere(Point3df(1.0, 2.0, 3.0), 0.7, 7);
    //til::Mesh2 mesh0;
    //til::convert(*s1, mesh0);
    
    std::cout << "Adding neighbors" << std::endl;
    mesh = addNeighborsToMesh(mesh0);
  }


  double distance = 10.0;
  std::cout << "Computing geomap" << std::endl;
  til::ghost::GMapStop_AboveThreshold<double> stopGhost(distance);
  //til::Triangle_mesh_geodesic_map<til::Mesh2_N, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage<til::sparse_vector, til::Mesh2_N, double > > geomap(mesh, stopGhost);
  typedef std::vector<std::vector<std::size_t> > CNeighborhoods;
  shared_ptr<CNeighborhoods> pneighc = til::circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  til::Triangle_mesh_geodesic_map<MyMesh::VertexCollection, CNeighborhoods, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage<til::sparse_vector, double > >
    geomap(getVertices(mesh), *pneighc, stopGhost);
  //std::vector<const til::MeshTraits<MyMesh>::Vertex *> startPoints(1);
  std::vector<std::size_t> startPoints(1);
  std::vector<double> dist(1, 0.0);

  std::vector<double> res(getVertices(mesh).size(), 0.0);

  int count = 0;

  til::MeshTraits<MyMesh>::VertexCollection::const_iterator iVertex = getVertices(mesh).begin();
  for(; iVertex != getVertices(mesh).end(); ++iVertex)
  {
    //if (++count > 10) break;
    //startPoints[0] = &*iVertex;
    startPoints[0] = count;
    geomap.init(startPoints, dist);
    geomap.process();
    shared_ptr<til::sparse_vector<double> > sres = geomap.distanceMap();
    //sres.setDefaultValue(0);
    //std::transform(sres.begin(), sres.end(), res.begin(), res.begin(), std::plus<double>());
    ++count;
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
  std::vector<double> res(til::size(sres));
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
  til::Triangle_mesh_geodesic_map<til::Mesh2_N, double, ghost::GMapStop_AboveThreshold<double> > geomap(mesh, stopGhost);
  std::vector<til::MeshTraits<til::Mesh2_N>::FaceIndex::value_type> startPoints;
  startPoints.push_back(getFaceIndices(mesh)[0][0]);
  std::vector<double> dist(1, 0.0);
  geomap.init(startPoints, dist);
  geomap.process();
  std::vector<double> res = geomap.distanceMap();
  */
  
  
  //for (std::size_t i = 0; i < til::size(res); ++i) res[i] = til::min(res[i], 100.0);
  
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
    til::convert(s, mesh);
    til::aimswrite(s, "sphere");
  }
  {
    Texture1d t; //(1, til::size(res));
    t.reserve(til::size(res));
    for (std::size_t i=0; i<til::size(res); ++i) //t[i] = res[i];
      t.push_back(res[i]);
    til::aimswrite(t, "spheredist");
  }  
}


void testSolver()
{
  std::cout << til::norm(til::numeric_array<float,3>(1, 2, 3) - til::numeric_array<float,3>(4, 5, 6), til::prec<double>()) << std::endl;
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


void param2affine(til::numeric_array<float, 12> const & params, til::AffineMap<float> & a)
{
  a.transfo().setTransl(til::numeric_array<float,3>(params[0],params[1],params[2]));
  a.transfo().setMatrix(
    til::Matrix3<float>(
      params[3],params[6],params[8],
      params[7],params[4],params[10],
      params[9],params[11],params[5]
      ));
    /*
    a.transfo().setMatrix(
      til::Matrix3<float>(
        params[0],params[1],params[2],
        params[3],params[4],params[5],
        params[6],params[7],params[8]
        ));
    a.transfo().setTransl(til::Vector<float,3>(params[9],params[10],params[11]));
    */
}



template < typename TMesh, typename TPrec, typename TConn, typename Finder >
class FiberDistance : public std::unary_function< til::numeric_array<float,12>, float >
{
private: // classes

  struct Dist2
  {
    Dist2() { res = 0; }
    void operator()(std::pair<std::size_t, double> p, std::pair<std::size_t, double> b) { return this->operator()(p.second, b.second); }
    void operator()(std::pair<std::size_t, double> p, double b) { return this->operator()(p.second, b); }
    void operator()(double a, std::pair<std::size_t, double> p) { return this->operator()(a, p.second); }
    void operator()(double a, double b) { res += square(a - b); }
    double res;
  };
 
public: // typedefs

  typedef typename Finder::index_type index_type;
  typedef typename til::MeshTraits<TMesh>::VertexCollection::const_iterator const_vertex_iterator;

public: // constructors & destructor

  FiberDistance(const TMesh & mesh1, const TMesh & mesh2, const TConn & conn1, const TConn & conn2, double alpha, Finder finder)
   : m_mesh1(mesh1)
   , m_mesh2(mesh2)
   , m_conn1(conn1)
   , m_conn2(conn2)
   , m_finder(finder)
   , m_alpha(alpha)
   , m_nn(til::size(getVertices(mesh2)))
   {
    //std::cout << "bbst2 " << getVertices(mesh2)[0] << " " << getVertices(mesh1)[0] << std::endl;
    //std::cout << "bbst " << getVertices(m_mesh2)[0] << " " << getVertices(m_mesh1)[0] << std::endl;
   }

  //float operator()(const til::Vector<float,12> & params)
  float operator()(til::numeric_array<float,12> const & params)
  {
    std::cout << "." << std::flush;
    typedef typename TMesh::Vertex Vertex;
    til::AffineMap<float> a;
    param2affine(params, a);

    double resdist = 0.0;
    double resconn = 0.0;

    //for (const_vertex_iterator iVertex2 = getVertices(mesh2).begin(); iVertex2 != getVertices(mesh2).end(); ++iVertex2)
    for (std::size_t i2 = 0; i2 < getVertices(m_mesh2).size(); ++i2)
    {
      Vertex av = a(getVertices(m_mesh2)[i2]);
      //std::size_t i1 = m_finder(a(getVertices(m_mesh2)[i2]));
      //resdist += til::dist2<double>(getVertices(m_mesh2)[i2], getVertices(m_mesh1)[i1]);
      std::size_t i1 = m_finder(av);
      resdist += til::dist2(av, getVertices(m_mesh1)[i1], til::prec<double>());
      m_nn[i2] = i1;
      /*
      if (i2 == 0)
      {
        std::cout << "@ " << getVertices(m_mesh2)[i2] << " " << a(getVertices(m_mesh2)[i2]) << " " << getVertices(m_mesh1)[i1] << " " << getVertices(m_mesh1)[0] << std::endl;
      }
      */
    }
    
    //resdist /= std::sqrt(max(0.0f,til::det(a.transfo().getMatrix()))) + 128 * std::numeric_limits<double>::epsilon();

    //unsigned int countMiss = 0;
    //unsigned int countHit = 0;
    for (std::size_t i2 = 0; i2 < getVertices(m_mesh2).size(); ++i2)
    {
      //std::cout << "i2=" << i2 << std::endl;
      til::sparse_vector<double> tmp;
      for (til::sparse_vector<double>::Map::const_iterator c1 = m_conn1[m_nn[i2]].getMap().begin();
           c1 != m_conn1[m_nn[i2]].getMap().end(); ++c1)
      {
        tmp[m_nn[c1->first]] += c1->second;
      }
      //std::cout << "A" << std::endl;
      Dist2 d;
      til::loop_mapEach(tmp.getMap(), m_conn2[i2].getMap(), d);
      resconn += d.res;
      //std::cout << "B" << std::endl;
      
      /*
      for (til::sparse_vector<double>::Map::const_iterator c2 = m_conn2[i2].getMap().begin();
        c2 != m_conn2[i2].getMap().end(); ++c2)
      {
        resconn += c2->second * m_conn1[m_nn[i2]][m_nn[c2->first]];
        if (m_conn1[m_nn[i2]][m_nn[c2->first]]) ++countHit;
        else ++countMiss;
        / *
        if (m_nn[c2->first] < 10)
        {
          std::cout << getVertices(m_mesh2)[c2->first] << " " << getVertices(m_mesh1)[m_nn[c2->first]] << " " << m_conn1[m_nn[i2]][m_nn[c2->first]] << " " << i2 << " " << m_nn[c2->first] << std::endl;
          for (int i = 0; i < 10; ++i) std::cout << m_conn1[m_nn[i2]][i] << std::endl;
        }
        * /
      }
      */
    }
    //std::cout << "Hit: " << countHit << "  Miss: " << countMiss << std::endl;

    double denom = std::sqrt(max(0.0f,til::det(a.transfo().getMatrix()))) + 128 * std::numeric_limits<double>::epsilon();
    resdist /= denom;
    
    //double res = -m_alpha * resconn + (1-m_alpha) * resdist;
    double res = m_alpha * resconn + (1-m_alpha) * resdist;
    //res /= denom;
    std::cout << "Func called at " << params << " : " << resdist/denom << " " << -resconn/denom << " " << res << std::endl;
    return res;
  }
  
private: // data

  const TMesh & m_mesh1;
  const TMesh & m_mesh2;
  const TConn & m_conn1;
  const TConn & m_conn2;
  Finder m_finder;
  double m_alpha;
  //std::vector<const_vertex_iterator> m_nn;
  std::vector<std::size_t> m_nn;
};



/*
void checkBundleTransfo(int argc, char * argv[])
{
  std::string bundleFilename;
  Motion motion;
  {
    std::string mname;
    AimsApplication app( argc, aims_const_hack(argv), "checkBundleTransfo" );
    app.addOption(bundleFilename, "-bundles", "input bundles" );
    app.addOption(mname, "-trs", "transform from anat to t2");
    app.initialize();
    MotionReader mreader(mname);
    mreader.read(motion);
    motion = motion.inverse();
  }
  til::AffineMap<float> a;
  {
    AimsData<float> rot = motion.rotation();
    Point3df trs = motion.translation();
    a.transfo().setMatrix(til::Matrix3<float>(rot[0], rot[1], rot[2], rot[3], rot[4], rot[5], rot[6], rot[7], rot[8]));
    a.transfo().setTransl(til::numeric_array<float,3>(trs[0], trs[1], trs[2]));
  }
  
  // loading fibers
  std::cout << "Loading fibers..." << std::flush;
  typedef std::vector<til::numeric_array<float,3> > Fiber;
  typedef std::vector<Fiber> Fibers;
  boost::shared_ptr<Fibers> pfibers;
  {
    aims::BundleReader bundleReader(bundleFilename);
    //aims::BundleReader r("/home/Panabase/pascal/fibers/res/old/res.bundles");
    //aims::BundleReader r("/home/Panabase/pascal/fibers/res/res.bundles");
    til::BundleLoader loader;
    bundleReader.addBundleListener(loader);
    bundleReader.read();
    //r.addBundleListener(loader);
    //r.read();
    pfibers = loader.getFibers();
  }
  std::cout << "OK" << std::endl;
  Fibers & fibers = *pfibers;

  std::cout << "Number of fibers: " << fibers.size() << std::endl;

  {
    typedef AimsData<unsigned char> Image;
    Image im(256, 256, 256);
    for (std::size_t i = 0; i < fibers.size(); ++i)
    {
      //Point3df p2 = motion.transform(fibers[i].front()[0], fibers[i].front()[1], fibers[i].front()[2]);
      //til::Point<float, 3> p(p2[0], p2[1], p2[2]);
      til::numeric_array<float, 3> p = a(fibers[i].front());
      if (int(p[0]) < 0 || int(p[1]) < 0 || int(p[2]) < 0)
      {
        std::cout << "w!";
        continue;
      }
      im(int(p[0]), int(p[1]), int(p[2])) = 1;
    }
    Writer<Image> w("/home/cathier/tmp/front");
    w.write(im);
  }
  {
    typedef AimsData<unsigned char> Image;
    Image im(256, 256, 256);
    for (std::size_t i = 0; i < fibers.size(); ++i)
    {
      //Point3df p2 = motion.transform(fibers[i].front()[0], fibers[i].front()[1], fibers[i].front()[2]);
      //til::Point<float, 3> p(p2[0], p2[1], p2[2]);
      til::numeric_array<float, 3> p = a(fibers[i].back());
      if (int(p[0]) < 0 || int(p[1]) < 0 || int(p[2]) < 0)
      {
        std::cout << "w!";
        continue;
      }
      im(int(p[0]), int(p[1]), int(p[2])) = 1;
    }
    Writer<Image> w("/home/cathier/tmp/back");
    w.write(im);
  }
}
*/



void testBundles(int argc, char * argv[])
{
  typedef til::Mesh_N MyMesh;
  MyMesh mesh;
  TimeTexture<short> tex;
  //double alpha = 0.5;
  double distthresh = 5.0;
  double wthresh = 1.0;
  bool writecm = false;
  std::size_t gyrus = 0;
  Motion motion;
  std::string bundleFilename;
  std::string outputDirname;
  const double G_THRESH = 0.001;
  const std::size_t NUM_LABELS = 36;
  {
    Reader<AimsSurfaceTriangle> r;
    Reader< TimeTexture<short> > texR;
    char mname[512];
    //Reader<Motion> rmotion;
    AimsApplication app( argc, aims_const_hack(argv), "testBundles" );
    app.addOption(r, "-i", "input mesh" );
    app.addOption(texR, "-tex", "input gyri tex" );
    app.addOption(gyrus, "-gyrus", "input gyrus label, 0 to calculate the connectivity for all the gyri" );
    app.addOption(writecm, "-cmatrix", "write vertex connectivity matrix", true );
    app.addOption(bundleFilename, "-bundles", "input bundles" );
    //app.addOption(alpha, "-alpha", "alpha");
    app.addOption(distthresh, "-dist", "dist");
    app.addOption(wthresh, "-wthresh", "weight threshold");
    app.addOption(mname, "-trs", "transform from anat to t2");
    app.addOption(outputDirname, "-outdir", "output dir" );
    app.initialize();

    std::cout << "input bundles: " << bundleFilename << std::endl;
    std::cout << "dist: " << distthresh << std::endl;
    std::cout << "output dir: " << outputDirname << std::endl;
    std::cout << "gyrus : " << gyrus << endl;

    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
    
    MotionReader mreader(mname);
    mreader.read(motion);
    
    cout << "reading texture..." << flush;
    texR.read( tex );
    cout << "done" << endl;
  }

  std::cout << "First vertice: " << getVertices(mesh)[0] << std::endl;
  std::cout << "# vertices : " << getVertices(mesh).size() << std::endl;
  std::cout << "# faces : " << getFaceIndices(mesh).size() << std::endl;
  std::cout << "texture dim : " << tex[0].nItem() << endl;

  std::cout << "Computing geomap..." << std::flush;
  til::ghost::GMapStop_AboveThreshold<double> stopGhost(distthresh);
  //til::Triangle_mesh_geodesic_map<MyMesh, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage<til::sparse_vector, MyMesh, double > >
  //  geomap(mesh, stopGhost);
  typedef std::vector<std::vector<std::size_t> > CNeighborhoods;
  shared_ptr<CNeighborhoods> pneighc = til::circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));

/*
  std::cout << "pneighc size : " << pneighc->size() << std::endl;
  std::cout << "Writing Circular Neighborhoods..." << std::endl;
  {
    std::fstream f;
    std::string tmpString = outputDirname + "CNeigh.txt";
    f.open(tmpString.c_str(), std::fstream::out);

    for (std::size_t i = 0; i < pneighc->size(); ++i)
    {
      f << i << ": ";
      for (std::size_t j = 0; j < ((*pneighc)[i]).size(); ++j)
      {
        f << (*pneighc)[i][j] << " ";
      }
      f << std::endl;
    }
  }
*/

  til::Triangle_mesh_geodesic_map<MyMesh::VertexCollection, CNeighborhoods, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage<til::sparse_vector, double > >
    geomap(getVertices(mesh), *pneighc, stopGhost);
  std::vector<std::size_t> startPoints(1);
  std::vector<double> dist(1, 0.0);
  typedef std::vector< std::pair<std::size_t, double> > QuickMap;
  std::vector<QuickMap> res(getVertices(mesh).size());
  std::vector<std::size_t> nneigh(til::size(getVertices(mesh)));
  til::MapHistogram<std::size_t> h;

//std::cout << "res size : " << res.size() << std::endl;
std::cout << "Writing GeoMap..." << std::endl;

  {
    std::fstream f;
    std::string tmpString = outputDirname + "GeoMap.txt";
//    f.open(tmpString.c_str(), std::fstream::out);
  
    for (std::size_t i = 0; i < til::size(getVertices(mesh)); ++i)
    {
      startPoints[0] = i;
      geomap.init(startPoints, dist);
      geomap.process();
      shared_ptr<til::sparse_vector<double> > tmp = geomap.distanceMap();
      res[i].resize(tmp->getMap().size());

      {
        using namespace til::expr;
        til::detail::loop_xx(castTo(*_1, *_2), res[i], tmp->getMap());
      }

      nneigh[i] = res[i].size();
//PAM//std::cout << "nneigh[" << i <<  "] : " << nneigh[i] << std::endl;
      h.accumulate(res[i].size());

/*      f << i << ": ";
      for (std::size_t j = 0; j < (res[i]).size(); ++j)
      {
        f << "(" << res[i][j].first << "," << res[i][j].second << ") ";
      }
      f << std::endl;
*/
    }
  }
  std::cout << "OK" << std::endl;

  /*
  std::cout << "HISTOGRAM:" << std::endl;
  std::map<std::size_t, unsigned int>::const_iterator iH = h.get().begin();
  for (; iH != h.get().end(); ++iH)
  {
    std::cout << iH->first << " : " << iH->second << std::endl;
  }
  */
  
  // writing number of neighbors in the DIST-cell
  std::string tmpString;
  {
    Texture1d t(1, til::size(nneigh));
    til::convert(t, nneigh);

    tmpString = outputDirname + "nneigh";
    til::aimswrite(t, tmpString);
  }

  // writing distances in the DIST-cell of the vertex N
  {
    int N = 1000;
    Texture1d t(1, til::size(getVertices(mesh)));
    for (std::size_t i = 0; i < til::size(res[N]); ++i)
    {
      t.item(res[N][i].first) = res[N][i].second;
    }
    tmpString = outputDirname + "toubou";
    til::aimswrite(t, tmpString);
  }

  // loading fibers
  std::cout << "Loading fibers..." << std::flush;
  typedef std::vector<til::numeric_array<float,3> > Fiber;
  typedef std::vector<Fiber> Fibers;
  boost::shared_ptr<Fibers> pfibers;
  {
    comist::BundleReader bundleReader(bundleFilename);
    //aims::BundleReader r("/home/Panabase/pascal/fibers/res/old/res.bundles");
    //aims::BundleReader r("/home/Panabase/pascal/fibers/res/res.bundles");
    til::BundleLoader loader;
    bundleReader.addBundleListener(loader);
    bundleReader.read();
    //r.addBundleListener(loader);
    //r.read();
    pfibers = loader.getFibers();
  }
  std::cout << "OK" << std::endl;
  Fibers & fibers = *pfibers;

  std::cout << "Number of fibers: " << fibers.size() << std::endl;
//0.939423 -2.49658 -16.9287
//1 -5.68286e-05 -0.000116213
//5.75763e-05 0.999979 0.00644452
//0.000115845 -0.00644452 0.999979
//    1.0000    0.0001    0.0001   -0.9373
//   -0.0001    1.0000   -0.0064    2.3875
//   -0.0001    0.0064    1.0000   16.9446
//         0         0         0    1.0000
  til::AffineMap<float> a;
  a.transfo().setMatrix(til::Matrix3<float>(1.0000, 0.0001, 0.0001, -0.0001, 1.0000, -0.0064, -0.0001, 0.0064, 1.0000));
  a.transfo().setTransl(til::numeric_array<float,3>(-0.9373, 2.3875, 16.9446));

/*////
  {
    std::string tmpString = outputDirname + "vertices";
    typedef AimsData<unsigned char> Image;
    Image im(256, 256, 256);
    for (std::size_t i = 0; i < getVertices(mesh).size(); ++i)
    {
      im(int(getVertices(mesh)[i][0]), int(getVertices(mesh)[i][1]), int(getVertices(mesh)[i][2])) = 1;
    }
    Writer<Image> w(tmpString.c_str());
    w.write(im);
  }
*/////

  {
    typedef AimsData<unsigned char> Image;
    Image im(256, 256, 256);
    for (std::size_t i = 0; i < fibers.size(); ++i)
    {
      til::numeric_array<float, 3> p = a(fibers[i].front());
      Point3df p2 = motion.transform(fibers[i].front()[0], fibers[i].front()[1], fibers[i].front()[2]);
      if (p2[0] != p[0])
      {
        std::cout << p2[0] << " " << p[0] << std::endl;
      }
      
      if (int(p[0]) < 0 || int(p[1]) < 0 || int(p[2]) < 0)
      {
        std::cout << "w!";
        continue;
      }
      im(int(p[0]), int(p[1]), int(p[2])) = 1;
    }
    tmpString = outputDirname + "front";
    til::aimswrite(im, tmpString);
  }
  /*
  {
    typedef AimsData<unsigned char> Image;
    Image im(256, 256, 256);
    for (int i = 0; i < fibers.size(); ++i)
    {
      til::Point<float, 3> p = a(fibers[i].back());
      if (int(p[0]) < 0 || int(p[1]) < 0 || int(p[2]) < 0)
      {
        std::cout << "w!";
        continue;
      }
      im(int(p[0]), int(p[1]), int(p[2])) = 1;
    }
    
    tmpString = outputDirname + "back";
    Writer<Image> w(tmpString.c_str());
    w.write(im);
  }
  */
  
  std::cout << "Generating kdtree" << std::endl;
  typedef til::KDTree<std::size_t, til::MeshTraits<MyMesh>::VertexCollection> MyKDTree;
  MyKDTree kdt(getVertices(mesh));
  makeKDTree(getVertices(mesh), kdt);

  std::cout << "Looking for closest points" << std::endl;
  typedef til::sparse_vector<double> Connectivity;
  typedef std::vector< Connectivity > Connectivities;
  Connectivities conn(getVertices(mesh).size(), til::sparse_vector<double>(getVertices(mesh).size()));
  int countRemoved = 0;
  std::size_t fiberCount = 0, nFibers = fibers.size();
  
  for (Fibers::const_iterator iFiber = fibers.begin(); iFiber != fibers.end(); ++iFiber, ++fiberCount)
  {
    if (fiberCount % int(nFibers/10.0) == 0)
    {
      std::cout << 10 * int((fiberCount / int(nFibers/10.0))) << "%..." << std::flush;
    }
    
    til::Find_closest< double, MyKDTree > fc(kdt);
    std::size_t A = fc(a(iFiber->front()));
    std::size_t B = fc(a(iFiber->back()));
    
    //std::cout << "front " << iFiber->front() << " " << getVertices(mesh)[A] << std::endl;
    //std::cout << "back " << iFiber->back() << " " << getVertices(mesh)[B] << std::endl;

    if (til::dist2(a(iFiber->front()), getVertices(mesh)[A], til::prec<float>()) > 25.0 ||
        til::dist2(a(iFiber->back()), getVertices(mesh)[B], til::prec<float>()) > 25.0)
    {
      //std::cout << "!";
      ++countRemoved;
      continue;
    }
    
    for (QuickMap::const_iterator j = res[B].begin(); j != res[B].end(); ++j)
    {
      double e2 = til::square(j->second);
      double w = std::exp( -e2  / ( 2*distthresh*distthresh ));
      if (w > G_THRESH)
      {
        conn[A][j->first] += w;
        conn[j->first][A] += w;
      }
    }
    //*/
  }
  std::cout << std::endl << "NB: " << countRemoved << " fibers out of " << fibers.size() << " have been discarded" << std::endl;
  
  // To compress a little bit, remove points that are below some threshold
  std::cout << "Removing weak connections" << std::endl;
  {
    int count1 = 0;
    int count2 = 0;
    for (Connectivities::iterator i = conn.begin(); i != conn.end(); ++i)
    {
      for (Connectivity::Map::iterator j = i->getMap().begin(); j != i->getMap().end(); )
      //for (Connectivity::sparse_iterator j = i->sparse_begin(); j != i->sparse_end(); ++j)
      {
        ++count1;
        if (j->second <= wthresh)
        {
          ++count2;
          //i->erase(j);
          i->getMap().erase(j++);
        }
        else
        {
          ++j;
        }
      }
    }
    std::cout << "Removed " << count2 << " elements out of " << count1 << std::endl;
  }

  std::vector< std::size_t > labels;
  labels.resize(NUM_LABELS + 1);
  for (std::size_t i = 0; i < labels.size(); ++i)
  {
     labels[i] = 0;
  }

  std::cout << "Reading Gyrus Labels: " << std::endl;
  {
    for (std::size_t i = 0; i < tex[0].nItem(); ++i)
    {
       labels[tex[0].item(i)]++;
    }

    for (std::size_t i = 0; i < labels.size(); ++i)
    {
       std::cout << i << ": " << labels[i] << std::endl;
    }
  }
  std::cout << "OK" << std::endl;

/*
  std::cout << "Reading Gyrus Labels: " << std::endl;
  {
    std::size_t countLabel = 0;
    std::size_t label;
    
    for (std::size_t i = 0; i < tex[0].nItem(); ++i)
    {
       label = tex[0].item(i);
       std::cout << "item " << i << ": " << label << std::endl;
       if ((label >= 0) && (label < NUM_LABELS)) 
         labels[label]++;
       else labels[NUM_LABELS]++;
       labels[tex[0].item(i)]++;
    }

    for (std::size_t i = 0; i < labels.size(); ++i)
    {
       std::cout << i << ": " << labels[i] << std::endl;
       if (labels[i] != 0) countLabel++;
    }
    std::cout << "Number of non-empty labels:" << countLabel << std::endl;
  }
  std::cout << "OK" << std::endl;
*/

/*
  std::cout << "Converting Gyrus Labels: " << std::endl;
  {
    std::vector< std::size_t > labelsMap;
    labelsMap.resize(NUM_LABELS + 1);
    labelsMap[0] = 0;
    std::size_t label;
    std::size_t countMap = 1;

    for (std::size_t i = 0; i < labels.size(); ++i)
    {
      if (labels[i] != 0)
      {
        labelsMap[i] = countMap;
        countMap++;
      }
      else labelsMap[i] = 0;

    }

    std::cout << "--Labels Map: " << std::endl;
    for (std::size_t i = 0; i < labelsMap.size(); ++i)
    {
      std::cout << i << ": " << labelsMap[i] << std::endl;
    }

    std::cout << "--Writing Converted Texture: " << std::endl;
    TimeTexture<short> texOut(1,getVertices(mesh).size());

    for (std::size_t i = 0; i < tex[0].nItem(); ++i)
    {
       label = tex[0].item(i);
       texOut[0].item(i) = labelsMap[label];
    }   
    tmpString = outputDirname + "outputTexture";
    Writer< TimeTexture<short> > texW(tmpString.c_str());
    texW.write( texOut );
  }
  std::cout << "OK" << std::endl;
*/
  std::vector< std::size_t > gyri;
  if (gyrus > NUM_LABELS) gyrus = 0;
  if (gyrus == 0)
  {
    for (std::size_t i = 1; i <= NUM_LABELS; ++i)
    {
      gyri.push_back(i);
    }
  }
  else
  {
    gyri.push_back(gyrus);
  }

  if (writecm)
  {
    for (std::size_t k = 0; k < gyri.size(); ++k)
    {
      std::cout << "Writing vertex connectivity matrix for gyrus: " << gyri[k] << ", #vertices: " << labels[gyri[k]] << std::endl;
      std::fstream f;
      char s[128];
      std::sprintf(s, "%sVertexConnMatrixGyrus%i.txt", outputDirname.c_str(), gyri[k]);
      f.open(s, std::fstream::out);

      std::size_t countLabel = 0;
      //cMatrix.reserve(labels[gyri[k]]); 

      std::size_t label;
      //for (std::size_t i = 0; i < getVertices(mesh).size(); ++i)
      for (std::size_t i = 0; i < 100; ++i)
      {
        label = tex[0].item(i);
        if (label == gyri[k])  // vertice appartient au gyrus
        {
           std::vector<float> line;
           //line.reserve(getVertices(mesh).size());
           for (std::size_t j = 0; j < getVertices(mesh).size(); ++j)
           {
             //line.push_back(conn[i][j]);
             f << conn[i][j] << " ";
           } 
           //cMatrix[countLabel] = line;
           countLabel++;
           f << std::endl;
        }
      }
    }
    std::cout << "OK" << std::endl;  
  }

  std::vector< std::vector< float > > cMatrix; // cMatrix[gyrus_vertices][gyri]
  
  for (std::size_t k = 0; k < gyri.size(); ++k)
  {
    cMatrix.resize(labels[gyri[k]]);

    std::cout << "Writing connectivity matrix for gyrus " << gyri[k];
    std::cout << ", size (#vertices_gyrus x #gyri): " << labels[gyri[k]] << " x "<< NUM_LABELS << std::endl;
    std::fstream f;
    char s[128];
    std::sprintf(s, "%sConnMatrixGyrus%i.txt", outputDirname.c_str(), gyri[k]);
    f.open(s, std::fstream::out);

    // Create and Initialize connectivity matrix for gyrus gyri[k]
    std::size_t label;
    for (std::size_t j = 0; j < labels[gyri[k]]; ++j)
    {
      std::vector<float> line;
      line.resize(NUM_LABELS);
      for (std::size_t i = 0; i < NUM_LABELS; ++i)
      {
        line[i] = 0;
      }
      cMatrix[j] = line;
    }

    // Calculate connectivity matrix for gyrus gyri[k]
    std::size_t label2;
    std::size_t countLabel = 0;
    for (std::size_t i = 0; i < getVertices(mesh).size(); ++i)
    {
      label = tex[0].item(i);
      if (label == gyri[k])  // vertice appartient au gyrus
      {
         for (std::size_t j = 0; j < getVertices(mesh).size(); ++j)
         {
            label2 = tex[0].item(j);  // gyrus auquel gyri[k] est connecté
            cMatrix[countLabel][label2-1] += conn[i][j];
         }
         countLabel++;
      }
    }

    // Write connectivity matrix for gyrus gyri[k]
    
    for (std::size_t i = 0; i < cMatrix.size(); ++i)
    {
      for (std::size_t j = 0; j < cMatrix[i].size(); ++j)
      {
          f << cMatrix[i][j] << " ";
      }
      f << std::endl;
    }
    std::cout << "Done. Matrix size : " << cMatrix.size() << " x "<< cMatrix[0].size() << std::endl;
  }

exit(0); 

  
  std::cout << "Writing connectivity matrix..." << std::endl;
  {
    std::fstream f;
    tmpString = outputDirname + "ConnMatrix.txt";
    f.open(tmpString.c_str(), std::fstream::out);
    for (std::size_t i = 0; i < 5000; ++i)
    {
      for (std::size_t j = 0; j < 5000; ++j)
      {
        f << conn[i][j] << " ";
      }
      f << std::endl;
    }
    // not needed!
    //f.close();
  }
  std::cout << "OK" << std::endl;

  std::cout << "Writing textures" << std::endl;
  for (int V = 0; V < 1000; V+=100)
  {
    Texture1d t(1, getVertices(mesh).size());
    for (til::sparse_vector<double>::Map::const_iterator i = conn[V].getMap().begin(); i != conn[V].getMap().end(); ++i)
    {
      //std::cout << i->first << " " <<  i->second << "*" << std::flush;
      t.item(i->first) = i->second;
    }
    
    char s[128];
    std::sprintf(s, "%swi%i", outputDirname.c_str(), V);
    til::aimswrite(t, s);
  }

}




/*
template < typename KDT, typename Conn, typename QuickMap, typename Mesh >
void smoothConnectivity(const Mesh & mesh, const KDT & kdt, const std::vector<QuickMap> & distmap, Conn & conn, double distthresh, double wthresh, double maxdist)
{
  int countRemoved = 0;
  std::size_t fiberCount = 0, nFibers = fibers.size();
  for (Fibers::const_iterator iFiber = fibers.begin(); iFiber != fibers.end(); ++iFiber, ++fiberCount)
  {
    indicator(fiberCount, nFibers);
    
    til::Find_closest< double, KDT > fc(kdt);
    //std::size_t A = fc(a(iFiber->front()));
    //std::size_t B = fc(a(iFiber->back()));
    std::size_t A = fc(iFiber->front());
    std::size_t B = fc(iFiber->back());
    
    //if (til::dist2(a(iFiber->front()), getVertices(mesh)[A], til::prec<float>()) > maxdist ||
    //    til::dist2(a(iFiber->back()), getVertices(mesh)[B], til::prec<float>()) > maxdist)
    if (til::dist2(iFiber->front(), getVertices(mesh)[A], til::prec<float>()) > maxdist ||
        til::dist2(iFiber->back(), getVertices(mesh)[B], til::prec<float>()) > maxdist)
    {
      ++countRemoved;
      continue;
    }
    
    / *
    for (QuickMap::const_iterator i = distmap[A].begin(); i != distmap[A].end(); ++i)
    for (QuickMap::const_iterator j = distmap[B].begin(); j != distmap[B].end(); ++j)
    {
      double e1 = til::square(i->second);
      double e2 = til::square(j->second);
      double w = std::exp(- ( e1 + e2 ) / ( 2*distthresh*distthresh ));
      if (w > G_THRESH)
      {
        conn[i->first][j->first] += w;
        conn[j->first][i->first] += w;
      }
    }
    * /
    /// *
    for (typename QuickMap::const_iterator j = distmap[B].begin(); j != distmap[B].end(); ++j)
    {
      double e2 = til::square(j->second);
      double w = std::exp( -e2  / ( 2*distthresh*distthresh ));
      if (w > G_THRESH)
      {
        conn[A][j->first] += w;
        conn[j->first][A] += w;
      }
    }
    // * /
  }
  std::cout << "OK" << std::endl;
  std::cout <<  countRemoved << " fibers out of " << fibers.size() << " have been discarded" << std::endl;
  
  // To compress a little bit, remove points that are below some threshold
  {
    int count1 = 0;
    int count2 = 0;
    for (typename Conn::iterator i = conn.begin(); i != conn.end(); ++i)
    {
      for (typename Conn::value_type::Map::iterator j = i->getMap().begin(); j != i->getMap().end(); )
      {
        ++count1;
        if (j->second <= wthresh)
        {
          ++count2;
          i->getMap().erase(j++);
        }
        else
        {
          ++j;
        }
      }
    }
    std::cout << "Removed " << count2 << " elements out of " << count1 << std::endl;
  }
}
*/


// Starting from a list of points, computes the associated normalized RBFs.
template < typename VertexCollection, typename Quant, typename RBF >
void computingNormalizedRBFWithSupport(const VertexCollection & vertices, const Quant & quant, const RBF & gaussian, std::vector<std::vector<std::pair<std::size_t, double> > > & defoBasis)
{
  std::size_t n = vertices.size();
  std::vector<std::vector<double> > rbf(quant->size(), std::vector<double>(n));
  for (std::size_t j = 0; j < rbf.size(); ++j)
  {
    //std::transform(meshes[level].first.begin(), meshes[level].first.end(), rbf[j].begin(), bind2nd(gaussian, *(*quant)[j]));
    for (std::size_t i = 0; i < n; ++i)
    {
      rbf[j][i] = gaussian(*(*quant)[j], vertices[i]);
    }
  }

  // It's probably better to threshold the RBF before normalization; first because thresholding does not
  // strictly preserve normalization (even though for mild thresholds it should not impact much); second
  // because before normalization we know the shape of the kernel and in particular its maximum of 1, meaning that
  // the threshold can be taken more confidently; last, because normalization does not undo thresholding -- zeroed
  // coefficients remain null.
  for (std::size_t j = 0; j < rbf.size(); ++j)
  {
    for (std::size_t i = 0; i < n; ++i)
    {
      if (rbf[j][i] < 10e-4) rbf[j][i] = 0.0;
    }
  }

  // normalization: unity sum
  {
    std::vector<double> sum(n, 0.0);
    for (std::size_t j = 0; j < rbf.size(); ++j)
    {
      std::transform(rbf[j].begin(), rbf[j].end(), sum.begin(), sum.begin(), std::plus<double>());
    }
    for (std::size_t j = 0; j < rbf.size(); ++j)
    {
      std::transform(rbf[j].begin(), rbf[j].end(), sum.begin(), rbf[j].begin(), std::divides<double>());
    }
  }

  // Conversion into compact storage
  for (std::size_t j = 0; j < rbf.size(); ++j)
  {
    for (std::size_t i = 0; i < n; ++i)
    {
      if (rbf[j][i] > 0.0)
      {
        defoBasis[j].push_back(std::make_pair(i, rbf[j][i]));
      }
    }
  }
}




template < typename TVertexCollection, typename TFaceCollection >
void meshMultiresStep2
(
  const TVertexCollection & inputvertices,
  const TFaceCollection & inputfaces,
  TVertexCollection & vertices,
  TFaceCollection & faces,
  int niter, float lambda, float mu)
{
  // copy input
  vertices = inputvertices;
  faces = inputfaces;
  
  std::size_t n = vertices.size();
  
  // computing neighborhoods
  //std::vector<std::vector<std::size_t> > neigh;
  //getNeighborIndices(vertices, faces, neigh);

  // smoothing mesh
  std::cout << "Smoothing mesh..." << std::flush;
  {
    typedef typename TVertexCollection::iterator                                        vertex_iterator;
    typedef til::xsr::graph::Position<vertex_iterator>                                  VertexIndexing;
    //typedef til::xsr::Integer_access<std::vector<std::vector<std::size_t> > >           NeighborIndexing;
    typedef til::xsr::graph::Neighbors<vertex_iterator>                                 NeighborIndexing;
    //typedef til::xsr::Iterator_access<typename TVertexCollection::iterator>             VertexIndexing2;
    typedef til::LaplacianSmoothing<VertexIndexing, VertexIndexing, NeighborIndexing>   Smoother;

    TVertexCollection verticesTmp(n);
    //VertexIndexing tmp1;
    //NeighborIndexing indNeigh;
    Smoother lambdaSmoothing((VertexIndexing()), (VertexIndexing()), (NeighborIndexing()), lambda);
    //VertexIndexing tmp4(verticesTmp);
    Smoother muSmoothing((VertexIndexing()), (VertexIndexing()), (NeighborIndexing()), -mu);

    for (int i = 0; i < niter; ++i)
    {
      //std::cout << "." << std::endl;
      //lambdaSmoothing(0, n, verticesTmp.begin());
      //muSmoothing(0, n, vertices.begin());      
      lambdaSmoothing(vertices.begin(), vertices.end(), verticesTmp.begin());
      muSmoothing(vertices.begin(), vertices.end(), vertices.begin());      
    }
  }
  std::cout << "OK" << std::endl;
  
  // quantizing mesh
  std::cout << "quantizing mesh..." << std::flush;
  std::map<std::size_t, std::size_t> q;
  {
    std::list<std::size_t> q2 = til::quantizer(vertices, faces, std::size_t(vertices.size() / 4.0));
    std::size_t i = 0;
    for (std::list<std::size_t>::iterator iQ = q2.begin(); iQ != q2.end(); ++iQ, ++i)
    {
      q[*iQ] = i;
    }
  }
  std::cout << "OK" << std::endl;  
}


/*
template < typename TMesh >
void meshMultires(const TMesh & mesh, int nlevel, int niter, double lambda, double mu,
std::vector<std::vector<til::numeric_array<float,3> > > & vertices,
std::vector<std::vector<til::numeric_array<std::size_t, 3> > > & faces
)
{
  // building multires pyr  
  typedef std::vector<std::vector<til::numeric_array<float,3> > > VertexCollections;
  typedef std::vector<std::vector<til::numeric_array<std::size_t, 3> > > FaceCollections;
  vertices.resize(nlevel);
  faces.resize(nlevel);
  vertices[0] = getVertices(mesh);
  faces[0] = getFaceIndices(mesh);
  for (int i = 1; i < nlevel; ++i)
  {
    std::cout << "[level-" << i << "]" << std::endl;
    meshMultiresStep(vertices[i-1], faces[i-1], vertices[i], faces[i], niter, lambda, mu);
  }
}
*/


template < typename TVertices, typename TFaces, typename TVerticesOut, typename TFacesOut >
void meshMultires3
(
  TVertices             const &   verticesIn,
  TFaces                const &   facesIn,
  int nlevel, int niter, double lambda, double mu,
  std::vector<TVerticesOut>   &   verticesOut,
  std::vector<TFacesOut>      &   facesOut
)
{
  // building multires pyr  
  verticesOut.resize(nlevel);
  facesOut.resize(nlevel);

  til::list2graph_mesh_conversion(verticesIn, facesIn, verticesOut[0], facesOut[0]);
  
  for (int i = 1; i < nlevel; ++i)
  {
    std::cout << "[level-" << i << "]" << std::endl;
    meshMultiresStep2(verticesOut[i-1], facesOut[i-1], verticesOut[i], facesOut[i], niter, lambda, mu);
  }
}




// bin/toto -i1 ~/data/cleaned_meshes/jeff_Lwhite_clean.mesh -i2 ~/data/cleaned_meshes/annieclaude_9398_Lwhite.clean.mesh -bundles1 ~Panabase/pascal/fibers/jeff/res/res.bundles -bundles2 /home/Panabase/vincent/data/annieclaude_9398/diffusion/annieclaude_9398_trackingALL_ASCII.bundles -trs /home/Panabase/vincent/data/annieclaude_9398/anniclaude_9398_T1_TO_T2.trm -alpha 0.2 -dist 10.0 -wthresh 0.5 -maxdist 25
// /home/cathier/code/c++/myproject-main-linux-release/bin/toto -i1 /home/cathier/data/cleaned_meshes/jeanpierre_9430_Lwhite.clean.transfo.mesh -i2 /home/cathier/data/cleaned_meshes/annieclaude_9398_Lwhite.clean.transfo.mesh -bundles1 /home/Panabase/vincent/data/jeanpierre_9430/diffusion/jeanpierre_9430_trackingALL_ASCII.bundles -bundles2 /home/Panabase/vincent/data/annieclaude_9398/diffusion/annieclaude_9398_trackingALL_ASCII.bundles -alpha 0.0001 -dist 10.0 -wthresh 0.5 -maxdist 25 -initstd 35.0 -sigma 20.0 -beta 0.8 -niter 100 -lambda 0.01 -mu 0.01 -nlevel 4
// regBundles


void testApplyAffineTransform()
{
  til::Mesh2_N mesh;
  {
    Reader<AimsSurfaceTriangle> r("/volatile/cathier/data/icbm/icbm100T_Lwhite.mesh");
    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh2 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  
  til::AffineMap<float> a;
  a.transfo().setMatrix(til::Matrix3<float>(1.1, -0.1, 0.1, 0.1, 0.9, 0.1, -0.1, 0, 1));
  a.transfo().setTransl(til::numeric_array<float,3>(20, 20, 20));
  
  
  PRINT_TIME3("loop",(
  {
    using namespace til::expr;
    til::detail::loop_x(*_1 = func(a)(*_1), getVertices(mesh));
  }
  ));
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh);
    til::aimswrite(s, "/home/cathier/tmp/aff");
  } 

  PRINT_TIME3("direct",(
  {
    til::MeshTraits<til::Mesh2_N>::VertexCollection::iterator iVertex = getVertices(mesh).begin();
    for (; iVertex != getVertices(mesh).end(); ++iVertex)
    {
      *iVertex = a(*iVertex);
    }
  }
  ));
  PRINT_TIME3("loop",(
  {
    using namespace til::expr;
    //til::detail::loop_x(*_1 = bind(a,*_1), getVertices(mesh));
    til::detail::loop_x(*_1 = func(a)(*_1), getVertices(mesh));
  }
  ));

  PRINT_TIME3("direct",(
  {
    til::MeshTraits<til::Mesh2_N>::VertexCollection::iterator iVertex = getVertices(mesh).begin();
    for (; iVertex != getVertices(mesh).end(); ++iVertex)
    {
      *iVertex = a(*iVertex);
    }
  }
  ));

  
  typedef til::KDTree<til::MeshTraits<til::Mesh2_N>::Vertex*, til::MeshTraits<til::Mesh2_N>::VertexCollection> MyKDTree;
  //typedef til::KDTree<std::size_t, til::MeshTraits<til::Mesh2_N>::VertexCollection> MyKDTree;
  MyKDTree kdtree(getVertices(mesh));
  makeKDTree(getVertices(mesh), kdtree);
  typedef til::Find_closest< double, MyKDTree > Finder;
  Finder finder(kdtree);
  {    
    double res = 0.0;
    //std::size_t i;
    til::MeshTraits<til::Mesh2_N>::Vertex *i;
    til::MeshTraits<til::Mesh2_N>::VertexCollection::iterator iVertex = getVertices(mesh).begin();
    for (; iVertex != getVertices(mesh).end(); ++iVertex)
    {
      i = finder(a(*iVertex));
      //res += til::dist2<double>(*iVertex, getVertices(mesh)[i])
      res += til::dist2(*iVertex, *i, til::prec<double>());
    }
  }

  /*
  {    
    double res = 0.0;
    til::detail::loop_x(var(res) +=  (_1 - func(finder)(func(a)(*_1))));
    til::detail::loop_x(
    (
      var(i) = func(finder)(func(a)(*_1)),
      var(res) += func(Norm2<double>)(*_1, var(getVertices(mesh))[var(i)])
    )
    , getVertices(mesh));
  }
  */

}


template < class TMesh, class TFinder >
struct Flabadel : public std::unary_function< til::numeric_array<float,12>, float>
{
  Flabadel
  (
    TMesh const & mesh,
    TMesh const & mesh2,
    TFinder const & finder
  )
    : m_mesh(mesh)
    , m_mesh2(mesh2)
    , m_finder(finder)
  {}
  
  //float operator()(const boost::array<float,12> & params)
  //float operator()(const til::Vector<float,12> & params)
  float operator()(til::numeric_array<float,12> const & params)
  {
    // TODO: When the Grand Task is done, there should be an Affine (mathematical) object
    // that accepts a container of size 12 as the underlying storage. And since that's obviously
    // extensible (say, who knows what crazy parameters people will come up with), that's
    // probably where a factory would be usefull.
    typedef typename til::MeshTraits<TMesh>::VertexCollection   VertexCollection;
    typedef typename til::MeshTraits<TMesh>::Vertex             Vertex;

    til::AffineMap<float> a;
    param2affine(params, a);
    double res = 0.0;
    typename VertexCollection::const_iterator iVertex = til::getVertices(m_mesh).begin();
    //std::cout << "Function called..." << std::flush;
    for (; iVertex != getVertices(m_mesh).end(); ++iVertex)
    {
      //std::cout << *iVertex << " " << a(*iVertex) << " " << *m_finder(a(*iVertex)) << std::endl;
      //res += til::dist2<double>(*iVertex, *m_finder(a(*iVertex)));
      //res += til::dist2<double>(*iVertex, getVertices(m_mesh)[m_finder(a(*iVertex))]);
      Vertex av = a(*iVertex);
      res += til::dist2(av, getVertices(m_mesh2)[m_finder(av)], til::prec<double>());
      /*
      if (iVertex == til::getVertices(m_mesh).begin())
      {
        std::cout << "@ " << *iVertex << " " << a(*iVertex) << " " << getVertices(m_mesh2)[m_finder(a(*iVertex))] << " " << getVertices(m_mesh2)[0] << std::endl;
      }
      if (params[0] == 10) std::cout << res << " ";
      */
    }
    //std::cout << "OK" << std::endl;
    //std::cout << a.transfo().getMatrix() << std::endl;
    //std::cout << til::det(a.transfo().getMatrix()) << std::endl;
    //std::cout << size(til::getVertices(m_mesh)) << " " << res << std::endl;
    res /= std::sqrt(max(0.0f,til::det(a.transfo().getMatrix()))) + 128 * std::numeric_limits<double>::epsilon();

    //std::cout << "Func called at " << params << " : " << res << std::endl;
    
    return float(res);
  }
  const TMesh & m_mesh;
  const TMesh & m_mesh2;
  TFinder m_finder;
};

template < class TMesh, class TFinder >
struct Flabadel2
  : public std::unary_function<til::numeric_array<float,3>, float>
{
  Flabadel2(const TMesh & mesh, const TMesh & mesh2, const TFinder & finder) : m_mesh(mesh), m_mesh2(mesh2), m_finder(finder) {}
  
  //float operator()(const boost::array<float,12> & params)
  //float operator()(const til::Vector<float,3> & params)
  float operator()(til::numeric_array<float,3> const & params)
  {
    // TODO: When the Grand Task is done, there should be an Affine (mathematical) object
    // that accepts a container of size 12 as the underlying storage. And since that's obviously
    // extensible (say, who knows what crazy parameters people will come up with), that's
    // probably where a factory would be usefull.
    til::AffineMap<float> a;
    a.transfo().setMatrix(
      til::Matrix3<float>(
        1,0,0,
        0,1,0,
        0,0,1
        ));
    a.transfo().setTransl(til::numeric_array<float,3>(params[0],params[1],params[2]));
    double res = 0.0;
    typename til::MeshTraits<TMesh>::VertexCollection::const_iterator iVertex;
    iVertex = til::getVertices(m_mesh).begin();
    //std::cout << "Function called..." << std::flush;
    for (; iVertex != getVertices(m_mesh).end(); ++iVertex)
    {
      //std::cout << *iVertex << " " << a(*iVertex) << " " << *m_finder(a(*iVertex)) << std::endl;
      //res += til::dist2<double>(*iVertex, *m_finder(a(*iVertex)));
      Vertex av = a(*iVertex);
      res += til::dist2(av, getVertices(m_mesh2)[m_finder(av)], til::prec<double>());
      //res += til::dist2<double>(*iVertex, getVertices(m_mesh2)[m_finder(a(*iVertex))]);
    }
    res /= std::sqrt(max(0.0f,til::det(a.transfo().getMatrix()))) + 128 * std::numeric_limits<double>::epsilon();
    //std::cout << "Func called at " << params << " : " << res << std::endl;
    return res;
  }
  const TMesh & m_mesh;
  const TMesh & m_mesh2;
  TFinder m_finder;
};

namespace boost
{
  std::ostream &
  operator<<(std::ostream & os, const array<float,12> & a)
  {
    os << "[ ";
    for (int i = 0; i < 12; ++i)
    {
      os << a[i] << " ";
    }
    os << "]";
    return os;
  }
}

//    0.8911    0.0990   -0.0990   -8.9109
//   -0.1089    1.0990   -0.0990   -8.9109
//    0.0891    0.0099    0.9901  -10.8911
//         0         0         0    1.0000
void testSimpleRegistration()
{
  std::cout << "Loading mesh..." << std::flush;
  til::Mesh_N mesh;
  
  {
    Reader<AimsSurfaceTriangle> r("/volatile/cathier/data/icbm/icbm100T_Lwhite.mesh");
    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  
  /*
  {
    AimsSurfaceTriangle *s = makeSphere(Point3df(1.0, 2.0, 3.0), 20, 7);
    til::Mesh1 mesh0;
    til::convert2(*s).into(mesh0);
    mesh = addNeighborsToMesh(mesh0);
  }
  */
  
  std::cout << "OK" << std::endl;
  
  til::AffineMap<float> a;
  //a.transfo().setMatrix(til::Matrix3<float>(1,0,0,0,1,0,0,0,1));
  //a.transfo().setTransl(til::Vector<float,3>(-10,0,0));

  a.transfo().setMatrix(til::Matrix3<float>(1.1, -0.1, 0.1, 0.1, 0.9, 0.1, -0.1, 0, 1));
  a.transfo().setTransl(til::numeric_array<float,3>(10, 10, 10));
  
  //a.transfo().setMatrix(til::Matrix3<float>(1.1, -0.1, 0.1, 0.1, 0.9, 0.1, -0.1, 0, 1));
  //a.transfo().setTransl(til::Vector<float,3>(1, 1, 1));
  //a.transfo().setMatrix(til::Matrix3<float>(1, 0, 0, 0, 1, 0, 0, 0, 1));
  //a.transfo().setTransl(til::Vector<float,3>(5, 5, 5));
  
  /*
  std::cout << "Smoothing mesh..." << std::flush;
  {
    std::cout << getVertices(mesh)[0] << std::endl;
    std::vector<til::Vector<float,3> > laplF;
    float k = 0.2;
    for (int i = 0; i < 1000; ++i)
    {
      laplacianForces(mesh, laplF);
      {
        using namespace til::expr;
        til::detail::loop_xx(*_1 += *_2 * k, getVertices(mesh), laplF);
      }
    }
    std::cout << getVertices(mesh)[0] << std::endl;
  }
  std::cout << "OK" << std::endl;
  */
  
  std::cout << "Generating new mesh..." << std::flush;
  // shallow copy
  til::Mesh_N mesh2(mesh);
  {
    using namespace til::expr;
    mesh2.vertices() = 
    boost::shared_ptr<til::MeshTraits<til::Mesh_N>::VertexCollection>(new til::MeshTraits<til::Mesh_N>::VertexCollection(size(getVertices(mesh))));
    //mesh2.vertices() = shared_ptr(new std::vector<til::Point<float,3> >(size(getVertices(mesh))));
    //til::detail::loop_xx(*_2 = bind(a,*_1), getVertices(mesh), getVertices(mesh2));
    til::detail::loop_xx(*_2 = func(a)(*_1), getVertices(mesh), getVertices(mesh2));
  }
  std::cout << "OK" << std::endl;
  std::cout << getVertices(mesh)[0] << " " << getVertices(mesh2)[0] << std::endl;

  std::cout << "Writing mesh..." << std::flush;
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh2);
    til::aimswrite(s, "/home/cathier/tmp/reginit");
  } 
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh);
    til::aimswrite(s, "/home/cathier/tmp/reginit0");
  } 
  std::cout << "OK" << std::endl;
  
  std::cout << "Generating KDTree..." << std::flush;
  //typedef til::KDTree<til::MeshTraits<til::Mesh_N>::Vertex*, til::MeshTraits<til::Mesh_N>::VertexCollection> MyKDTree;
  typedef til::KDTree<std::size_t, til::MeshTraits<til::Mesh_N>::VertexCollection> MyKDTree;
  MyKDTree kdtree(getVertices(mesh));
  makeKDTree(getVertices(mesh), kdtree);
  std::cout << "OK" << std::endl;

  typedef til::Find_closest< double, MyKDTree > Finder;
  Finder finder(kdtree);
  //typedef Flabadel2<til::Mesh_N, Finder> Energy;
  typedef Flabadel<til::Mesh_N, Finder> Energy;
  Energy energy(mesh2, mesh, finder);

  std::cout << "Starting minimization" << std::endl;
  til::Powell<Energy> minalgo(energy);
  // initial scaling estimates of the parameters
  {
    std::vector<float> initStd(12, 0.1f);
    //initStd[9] = initStd[10] = initStd[11] = 50;
    initStd[0] = initStd[1] = initStd[2] = 10;
    minalgo.initStd() = initStd;
  }

  //til::Vector<float,12> params;
  til::numeric_array<float,12> params;
  params[3] = 1;
  params[4] = 1;
  params[5] = 1;
  /*
  params[0] = 1;
  params[4] = 1;
  params[8] = 1;
  */
  //til::Vector<float,3> params;
  params = minalgo(params);
  /*
  params[0] = 0.894439;
  params[1] = 0.102035;
  params[2] = -0.100592;
  params[3] = -0.101044;
  params[4] = 1.09844;
  params[5] = -0.0982294;
  params[6] = 0.0802507;
  params[7] = 0.00803045;
  params[8] = 0.992575;
  params[9] = -9.6078;
  params[10] = -9.99204;
  params[11] = -9.65953;
  */
  std::cout << "Minimization finished: " << params << std::endl;
  {
    using namespace til::expr;
    til::AffineMap<float> a;
    param2affine(params, a);
    for (std::size_t i = 0; i < size(getVertices(mesh2)); ++i)
    {
      getVertices(mesh2)[i] = a(getVertices(mesh2)[i]);
      //a(getVertices(mesh2)[i]);
    }
    //til::detail::loop_x(*_1 = func(a)(*_1), getVertices(mesh2));
  }
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh2);
    til::aimswrite(s, "/home/cathier/tmp/regres");
  } 
  std::cout << "b" << std::endl;
}

struct Count
{
  Count() : m_c(0) {}
  template < typename T >
  void operator()(T) { ++m_c; }
  std::size_t get() { return m_c; }
  std::size_t m_c;
};

void testKDTree0()
{
  std::cout << "Loading mesh..." << std::flush;
  til::Mesh_N mesh;
  
  {
    Reader<AimsSurfaceTriangle> r("/volatile/cathier/data/icbm/icbm100T_Lwhite.mesh");
    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  
  std::cout << "OK" << std::endl;

  std::cout << til::size(getVertices(mesh)) << std::endl;
  
  til::AffineMap<float> a;
  a.transfo().setMatrix(til::Matrix3<float>(1,0,0,0,1,0,0,0,1));
  a.transfo().setTransl(til::numeric_array<float,3>(0,0,0));
  
  std::cout << "Generating new mesh..." << std::flush;
  // shallow copy
  til::Mesh_N mesh2(mesh);
  /*
  {
    using namespace til::expr;
    mesh2.vertices() = 
    boost::shared_ptr<til::MeshTraits<til::Mesh_N>::VertexCollection>(new til::MeshTraits<til::Mesh_N>::VertexCollection(size(getVertices(mesh))));
    //mesh2.vertices() = shared_ptr(new std::vector<til::Point<float,3> >(size(getVertices(mesh))));
    //til::detail::loop_xx(*_2 = bind(a,*_1), getVertices(mesh), getVertices(mesh2));
    til::detail::loop_xx(*_2 = func(a)(*_1), getVertices(mesh), getVertices(mesh2));
  }
  */
  std::cout << "OK" << std::endl;
  
  std::cout << "Generating KDTree..." << std::flush;
  //typedef til::KDTree<til::MeshTraits<til::Mesh_N>::Vertex*, til::MeshTraits<til::Mesh_N>::VertexCollection> MyKDTree;
  typedef til::KDTree<std::size_t, til::MeshTraits<til::Mesh_N>::VertexCollection> MyKDTree;
  MyKDTree kdtree(getVertices(mesh));
  makeKDTree(getVertices(mesh), kdtree);
  std::cout << "OK" << std::endl;

  std::cout << til::size(kdtree) << std::endl;

  Count c;  
  til::pre_order_scan(kdtree.root(), c);
  std::cout << c.get() << std::endl;

  typedef til::Find_closest< double, MyKDTree > Finder;
  Finder finder(kdtree);

  
  til::numeric_array<float,3> p;
  p[0] = 106.748;
  p[1] = 121.989;
  p[2] = 29.9544;


  
  til::MeshTraits<til::Mesh2_N>::VertexCollection::const_iterator iVertex;
  iVertex = til::getVertices(mesh2).begin();
  for (; iVertex != getVertices(mesh2).end(); ++iVertex)
  {
    if (til::dist2(*iVertex, getVertices(mesh2)[finder(a(*iVertex))], til::prec<double>()))
    {
      std::cout << *iVertex << " " << getVertices(mesh2)[finder(a(*iVertex))] << " " << find_closest(getVertices(mesh2), *iVertex) << std::endl;
    }
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


struct Toto : public std::unary_function<til::numeric_array<double, 5>, double>
{
  double operator()(const til::numeric_array<double, 5> & v)
  {
    float res = v[0]*v[0] + 100*v[1]*v[1] + 1000*v[2]*v[2] + 10*v[3]*v[3] + 50*v[4]*v[4] + v[3]*v[4];
    //std::cout << "Toto: " << res << std::endl;
    return res;
  }
};

struct DToto : public std::binary_function<const til::numeric_array<double, 5> &, til::numeric_array<double, 5> &, void>
{
  void operator()(const til::numeric_array<double, 5> & v, til::numeric_array<double, 5> & d)
  {
    d[0] = 2 * v[0];
    d[1] = 2 * 100 * v[1];
    d[2] = 2 * 1000 * v[2];
    d[3] = 2 * 10 * v[3] + v[4];
    d[4] = 2 * 50 * v[4] + v[3];
  }
};


void testPowell()
{
  til::Powell<Toto> p((Toto()));
  til::numeric_array<double, 5> min;
  min[0] = 1000.0;
  min[1] = -1000.0;
  min[2] = 2000.0;
  min[3] = -10000.0;
  min[4] = 2;
  
  std::cout << "Powell" << std::endl;
  min = p(min);
  std::cout << "Powell result: " << min << std::endl;
}

void testCG()
{
  til::PRConjugateGradient<Toto, til::GradientEstimator<Toto> > p((Toto()), (til::GradientEstimator<Toto>(Toto(), 0.001)), til::LineMin<Toto>((Toto())));
  til::numeric_array<double, 5> min;
  min[0] = 1000.0;
  min[1] = -1000.0;
  min[2] = 2000.0;
  min[3] = -10000.0;
  min[4] = 2;
  
  std::cout << "CG" << std::endl;
  min = p(min);
  std::cout << "CG result: " << min << std::endl;
}

void testCG2()
{
  til::PRConjugateGradient<Toto, DToto> p((Toto()), (DToto()), til::LineMin<Toto>((Toto())));
  til::numeric_array<double, 5> min;
  min[0] = 1000.0;
  min[1] = -1000.0;
  min[2] = 2000.0;
  min[3] = -10000.0;
  min[4] = 2;
  
  std::cout << "CG" << std::endl;
  min = p(min);
  std::cout << "CG result: " << min << std::endl;
}

void testGD()
{
  til::GradientDescent<Toto, DToto> p((Toto()), (DToto()), til::LineMin<Toto>((Toto())));
  til::numeric_array<double, 5> min;
  min[0] = 1000.0;
  min[1] = -1000.0;
  min[2] = 2000.0;
  min[3] = -10000.0;
  min[4] = 2;
  
  std::cout << "GD" << std::endl;
  min = p(min);
  std::cout << "GD result: " << min << std::endl;
}


void testMatrixDet()
{
  {
    til::Matrix3<float> m(1,0,0,0,1,0,0,0,1);
    std::cout << til::det<float>(m) << " (1)" << std::endl;
  }
  {
    til::Matrix3<float> m(1,0,0,0,1,0,0,0,-1);
    std::cout << til::det<float>(m) << " (-1)"  << std::endl;
  }
  {
    til::Matrix3<float> m(1,0,0,0,-1,0,0,0,1);
    std::cout << til::det<float>(m) << " (-1)"  << std::endl;
  }
  {
    til::Matrix3<float> m(-1,0,0,0,-1,0,0,0,1);
    std::cout << til::det<float>(m) << " (1)"  << std::endl;
  }
  {
    til::Matrix3<float> m(0,1,0,0,0,1,1,0,0);
    std::cout << til::det<float>(m) << " (-1)"  << std::endl;
  }
  {
    til::Matrix3<float> m(2,0,0,0,1,0,0,0,1);
    std::cout << til::det<float>(m) << " (2)" << std::endl;
  }
  {
    til::Matrix3<float> m(2,3,4,0,3,0,0,0,4);
    std::cout << til::det<float>(m) << " (24)" << std::endl;
  }
}



template < typename TVertexCollection >
void printBoundingBox(const TVertexCollection & vertices)
{
  typedef typename TVertexCollection::value_type::value_type prec_type;
  prec_type minx, miny, minz;
  prec_type maxx, maxy, maxz;
  minx = miny = minz = std::numeric_limits<prec_type>::max();
  maxx = maxy = maxz = -minx;
  
  for (typename TVertexCollection::const_iterator i = vertices.begin(); i != vertices.end(); ++i)
  {
    minx = min(minx, (*i)[0]);
    miny = min(miny, (*i)[1]);
    minz = min(minz, (*i)[2]);

    maxx = max(maxx, (*i)[0]);
    maxy = max(maxy, (*i)[1]);
    maxz = max(maxz, (*i)[2]);
  }

  std::cout << "Bounding box = " << minx << " " << miny << " " << minz << " ; " << maxx << " " << maxy << " " << maxz << " " << std::endl;
}

void geoinflate(int argc, char * argv[])
{
  til::Mesh_N mesh;
  std::cout << "Reading mesh..." << std::flush;
  float kcentroid = 0.02f;
  float kgeo = 0.02f;
  {
    Reader<AimsSurfaceTriangle> r;

    AimsApplication app( argc, aims_const_hack(argv), "geoinflate" );
    app.addOption(r, "-i", "input mesh" );
    app.addOption(kcentroid, "-kcentroid", "centroid forces" );
    app.addOption(kgeo, "-kgeo", "geodesic forces" );
    app.initialize();

    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  /*
  {
    AimsSurfaceTriangle * s = aims::SurfaceGenerator::sphere(Point3df(0,0,0), 100.0f, 100);
    til::Mesh1 mesh0;
    til::convert(mesh0, *s);
    mesh = addNeighborsToMesh(mesh0);
    delete s;    
  }
  */
  /*
  {
    AimsSurfaceTriangle *s = makeSphere(Point3df(1.0, 2.0, 3.0), 100, 5);
    til::Mesh1 mesh0;
    til::convert2(*s).into(mesh0);
    mesh = addNeighborsToMesh(mesh0);
  }
  */
  
  std::size_t n = size(getVertices(mesh));
  
  printBoundingBox(getVertices(mesh));
  
  /*
  std::cout << "Writing mesh..." << std::flush;
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh);
    Writer<AimsSurfaceTriangle> w("/home/cathier/tmp/sphere");
    w.write(s);
  }
  std::cout << "OK" << std::endl;
  */
  
  std::cout << "Get circular neighbor..." << std::flush;
  shared_ptr<std::vector<std::vector<std::size_t> > > neighc = circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  std::cout << "OK" << std::endl;
  
  // computing initial geodesic distance
  double distthresh = 5.0;
  std::vector<til::sparse_vector<double> > dist = getDistToNeighbors(mesh, distthresh);  
  //std::vector<std::vector<double> > dist = getDistToNeighbors(mesh, 10.0);
  //std::vector<shared_ptr<til::sparse_vector<double> > > dist = getDistToNeighbors(mesh, 10.0);
   
  std::cout << "Write dist..." << std::flush;
  /*
  {
    
    Texture1d t(1, til::size(getVertices(mesh)));
    for (std::size_t i = 0; i < til::size(dist[0]); ++i)
    {
      t.item(i) = dist[0][i];
    }
    Writer<Texture1d> w("/home/cathier/tmp/toubou");
    w.write(t);
  }
  */
  
  /*
  {
    Texture1d t(1, til::size(getVertices(mesh)));
    til::sparse_vector<double>::Map::const_iterator i = dist[0].getMap().begin();
    std::cout << getVertices(mesh)[0] << std::endl;
    for (; i != dist[0].getMap().end(); ++i)
    {
      t.item(i->first) = i->second;
      std::cout << getVertices(mesh)[i->first] << " " << i->second << std::endl;
    }
    Writer<Texture1d> w("/home/cathier/tmp/toubou");
    w.write(t);
  }
  std::cout << "OK" << std::endl;
  */
  
  std::vector<std::vector<std::size_t> >  geoneighbors(n);
  std::vector<std::vector<float> >        geodist0(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    til::sparse_vector<double>::Map::const_iterator j = dist[i].getMap().begin();
    for (; j != dist[i].getMap().end(); ++j)
    {
      // don't push the point itself
      if (j->first == i)
      {
        assert(j->second == 0.0);
        continue;
      }
      if (j->second == 0.0)
      {
        // weird!!!
        std::cerr << "!";
        continue;
      }
      geoneighbors[i].push_back(j->first);
      geodist0[i].push_back(j->second);
    }
    /*
    for (std::size_t j = 0; j < til::size(dist[i]); ++j)
    {
      if (j == i) continue;
      if (dist[i][j] == 0.0) continue;
      geoneighbors[i].push_back(j);
      geodist0[i].push_back(dist[i][j]);    
    }
    */
  }
  
  /*
  for (std::size_t i = 0; i < til::size(geoneighbors[0]); ++i) std::cout << geoneighbors[0][i] << " "; std::cout << std::endl;
  for (std::size_t i = 0; i < til::size(geoneighbors[n-1]); ++i) std::cout << geoneighbors[n-1][i] << " "; std::cout << std::endl;
  for (std::size_t i = 0; i < til::size(geodist0[0]); ++i) std::cout << geodist0[0][i] << " "; std::cout << std::endl;
  for (std::size_t i = 0; i < til::size(geodist0[n-1]); ++i) std::cout << geodist0[n-1][i] << " "; std::cout << std::endl;
  */
  
  til::MeshTraits<til::Mesh_N>::Vertex center;
  std::vector<til::numeric_array<float,3> > fcenter(n);
  std::vector<til::numeric_array<float,3> > fdist(n);
  AimsSurfaceTriangle wmesh;
  //int niter = 100;
  int niter = 1000;
  int dump = 100;
  int c = 0;
  for (int iter = 0; iter < niter; ++iter)
  {
    std::cout << "iter " << iter << std::endl;
    
    if (iter % (niter/dump) == 0)
    //if (0)
    //if (1)
    {
      // write mesh
      {       
        AimsSurfaceTriangle tmp;
        til::convert(tmp, mesh);
        char s[256];
        sprintf(s, "/home/cathier/tmp/geoinf-%i", c);
        til::aimswrite(tmp, s);
        
        //std::cout << "copying" << std::endl;
        //AimsSurface<3,Void> tmp;
        //til::convert(tmp, mesh);
        //wmesh[c] = tmp;
        ++c;
      }
    }

    struct timeval start, finish;
    struct timezone tz;
    gettimeofday(&start, &tz);
    
    
    // Compute centroidal forces
    std::cout << "centroidal forces..." << std::flush;
    centroid(getVertices(mesh), center);
    {
      til::MeshTraits<til::Mesh_N>::VertexCollection::const_iterator iVc = getVertices(mesh).begin();
      std::vector<til::numeric_array<float,3> >::iterator iC = fcenter.begin();
      for (; iVc != getVertices(mesh).end(); ++iVc, ++iC)
      {
        //*iC = (*iVc - center) * (1.0f/til::norm<float>(*iVc - center));
        *iC = (*iVc - center);
        *iC *= (1.0f/til::norm<float>(*iC));
      }
    }
    std::cout << "OK" << std::endl;
    std::cout << "center " << center << std::endl;

    
    // Compute geodesic forces
    std::cout << "geodesic forces..." << std::flush;
    {
      const float delta = 0.05;
      //til::sparse_vector<double> geodist;
      shared_ptr<std::vector<double> > geodist;
      for (std::size_t i = 0; i < n; ++i)
      {
        //std::cout << i << std::endl;

        // Get current distances to neighbors
        geodist = getDistToNeighbor2(mesh, neighc, geoneighbors, i);
        // Sum all these distances to get current energy
        float e = 0;
        for (std::size_t j = 0; j < til::size(geoneighbors[i]); ++j)
        {
          e += square(geodist0[i][j] - (*geodist)[geoneighbors[i][j]]);
        }

        // save position of vertex
        til::MeshTraits<til::Mesh_N>::Vertex v0 = getVertices(mesh)[i];
        
        
        // Do the same thing as before, but shifting current point by a small amount in the
        // x, y, then z direction, to get the gradient of the energy
        for (int k = 0; k < 3; ++ k)
        {
          // Shift current vertex by delta along the current axis
          getVertices(mesh)[i] = v0;
          getVertices(mesh)[i][k] += delta;
          // get distances to neighbors
          geodist = getDistToNeighbor2(mesh, neighc, geoneighbors, i);
          // Sum distances to get energy
          float tmp = 0;
          for (std::size_t j = 0; j < til::size(geoneighbors[i]); ++j)
          {
            tmp += square(geodist0[i][j] - (*geodist)[geoneighbors[i][j]]);
          }
          fdist[i][k] = e - tmp;
        }
        fdist[i] *= 1/delta;
        // Restore original vertex position
        getVertices(mesh)[i] = v0;
      }
    }
    std::cout << "OK" << std::endl;

    centroid(getVertices(mesh), center);
    std::cout << "center2 " << center << std::endl;

    // apply forces
    double tmp1, tmp2;
    tmp1 = 0.0;
    tmp2 = 0.0;
    til::numeric_array<float,3> sumcf(0,0,0);
    std::cout << "Applying forces..." << std::flush;
    for (std::size_t i = 0; i < size(getVertices(mesh)); ++i)
    {
      sumcf += fcenter[i];
      tmp1 += til::norm(fcenter[i], til::prec<double>());
      tmp2 += til::norm(fdist[i], til::prec<double>());
      getVertices(mesh)[i] += kcentroid * fcenter[i];
      getVertices(mesh)[i] += kgeo * fdist[i];
    }
    std::cout << "OK" << std::endl;

    std::cout << "fcenter " << tmp1 << std::endl;
    std::cout << "fdist " << tmp2 << std::endl;
    sumcf *= 1.0f / size(getVertices(mesh));
    std::cout << "sumcf " << sumcf << std::endl;

    gettimeofday(&finish, &tz);
    std::cout << "iter time : " << (double)((finish.tv_sec-start.tv_sec) * 1000000L + (finish.tv_usec-start.tv_usec)) << std::endl;
  }
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh);
    til::aimswrite(s, "/home/cathier/tmp/finalevel");
  } 
  til::aimswrite(wmesh, "/home/cathier/tmp/mevel");
}

struct MeshVertexNode;

struct MeshFaceNode
{
  typedef boost::array<std::list<MeshVertexNode>::iterator,3>   Face;
  boost::array<std::list<MeshVertexNode>::iterator,3>           face;
};

class MeshVertexNode
{
public: // typedefs

  typedef MeshVertexNode                        Self;
  typedef til::numeric_array<float,3>                   Vertex;
  typedef std::list<Self>                       VertexCollection;
  typedef VertexCollection::iterator            VertexIndex;
  typedef std::list<VertexIndex>                VertexIndexCollection;
  
  typedef MeshFaceNode                          Face;
  typedef std::list<Face>                       FaceCollection;
  typedef FaceCollection::iterator              FaceIndex;
  typedef std::list<FaceIndex>                  FaceIndexCollection;
 
public: // set & get
 
  Vertex       & pos()       { return m_pos; }
  Vertex const & pos() const { return m_pos; }

  VertexIndexCollection       & neighbors()       { return m_neighbors; }
  VertexIndexCollection const & neighbors() const { return m_neighbors; }

  FaceIndexCollection       & faces()       { return m_faces; }
  FaceIndexCollection const & faces() const { return m_faces; }
  
private: // data

  Vertex                  m_pos;
  VertexIndexCollection   m_neighbors;
  FaceIndexCollection     m_faces;
};


void testcurv()
{ 
  til::Mesh_N mesh;
  
  std::cout << "Reading mesh..." << std::flush;
  {
    Reader<AimsSurfaceTriangle> r("/volatile/cathier/data/icbm/icbm100T_Lwhite.mesh");
    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }

  /*
  {
    AimsSurfaceTriangle *s = makeSphere(Point3df(1.0, 2.0, 3.0), 20, 7);
    til::Mesh1 mesh0;
    til::convert2(*s).into(mesh0);
    mesh = addNeighborsToMesh(mesh0);
  }
  */

  std::size_t n = size(getVertices(mesh));

  /*
  std::cout << "Smoothing mesh..." << std::flush;
  {
    std::cout << getVertices(mesh)[0] << std::endl;
    std::vector<til::Vector<float,3> > laplF;
    float k = 0.2;
    for (int i = 0; i < 1000; ++i)
    {
      laplacianForces(mesh, laplF);
      {
        using namespace til::expr;
        til::detail::loop_xx(*_1 += *_2 * k, getVertices(mesh), laplF);
      }
    }
    std::cout << getVertices(mesh)[0] << std::endl;
  }
  std::cout << "OK" << std::endl;
  */
  
  
  std::cout << "Get circular neighbor..." << std::flush;
  shared_ptr<std::vector<std::vector<std::size_t> > > neighc = circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  std::cout << "OK" << std::endl;
  
  
  std::vector<std::list<std::size_t> > faceIndices = til::invertIndices(getFaceIndices(mesh));
  {
    std::list<std::size_t>::const_iterator i = faceIndices[0].begin();
    for (; i != faceIndices[0].end(); ++i)
    {
      std::cout << *i << std::endl;
    }
  }
  {
    std::list<std::size_t>::const_iterator i = faceIndices.back().begin();
    for (; i != faceIndices.back().end(); ++i)
    {
      std::cout << *i << std::endl;
    }
  }
  

  /*
  std::cout << "Converting into graph..." << std::flush;
  std::list<MeshVertexNode> graph_vertices;
  std::list<MeshFaceNode> graph_faces;
  std::vector<MeshVertexNode*> index_vertex(size(getVertices(mesh)));
  std::vector<MeshFaceNode*> index_face(size(getFaceIndices(mesh)));
  {
    // First pass: push all points and get an index2iterator translation
    for (std::size_t i = 0; i < n; ++i)
    {
      MeshVertexNode m;
      m.pos() = getVertices(mesh)[i];
      index_vertex[i] = graph_vertices.insert(graph_vertices.end(), m);
    }
    // second pass: finish filling vertex structure
    for (std::size_t i = 0; i < n; ++i)
    {
      for (std::size_t j = 0; j < size((*neighc)[i]); ++j)
      {
        index_vertex[i].pos().append(index_vertex[(*neighc)[i][j]]);
      }
    }
    // faces
    for (std::size_t i = 0; i < size(getFaceIndices(mesh)); ++i)
    {
      graph_vertices();
    }    
  }
  std::cout << "OK" << std::endl;
  */
  
  {
    const unsigned char KEEP = 0;
    const unsigned char DELETE = 1;
    const unsigned char FREEZE = 2;
    std::cout << "Suppressing 3-connected vertices..." << std::endl;
    std::vector<unsigned char> suppressVertex(size(getVertices(mesh)), 0);
    std::vector<unsigned char> suppressFace(size(getFaceIndices(mesh)), 0);
    int nDelVertex = 0;
    int nDelFace = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
      if ((suppressVertex[i] == KEEP) && (til::size((*neighc)[i]) == 3))
      {
        ++nDelVertex;        
        // mark vertex as deleted
        suppressVertex[i] = DELETE;
        // mark its neighbors as undeletable
        for (std::size_t j = 0; j < til::size((*neighc)[i]); ++j)
        {
          suppressVertex[(*neighc)[i][j]] = FREEZE;
        }
        
        std::list<std::size_t>::iterator j = faceIndices[i].begin();
        for (;j != faceIndices[i].end(); ++j)
        {
          if (suppressFace[*j] == KEEP)
          {
            suppressFace[*j] = DELETE;
            ++nDelFace;
          }
        }
      }
    }
    std::cout << "Deleted " << nDelVertex << " vertices and " << nDelFace << " faces" << std::endl;    
  }

  std::cout << "Writing mesh..." << std::flush;
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh);
    til::aimswrite(s, "/home/cathier/tmp/smoothed1");
  } 

  //for (std::size_t i = 0; i < til::size((*neighc)[0]); ++i) std::cout << (*neighc)[0][i] << " "; std::cout << std::endl;
  //for (std::size_t i = 0; i < til::size((*neighc)[n-1]); ++i) std::cout << (*neighc)[n-1][i] << " "; std::cout << std::endl;
  
  std::cout << "Get curvature..." << std::flush;
  typedef float prec_type;
  til::Mesh_curvature<til::MeshTraits<til::Mesh_N>::VertexCollection, std::vector<std::vector<std::size_t> >, prec_type>
    mc(getVertices(mesh), *neighc);
  std::vector<float> curv(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    mc.process(i);
    curv[i] = min(prec_type(3000.0), 1000*mc.meanCurvature());
    curv[i] = max(float(-3000.0), curv[i]);
  }
  std::cout << "OK" << std::endl;

  std::cout << "Write curvature..." << std::flush;
  {
    Texture1d tex(1,n);
    for (std::size_t i = 0; i < n; ++i)
    {
      tex.item(i) = curv[i];
      //std::cout << curv[i] << " ";
    }
    //til::convert(tex, curv);
    til::aimswrite(tex, "/home/cathier/tmp/curv");
  }
  std::cout << "OK" << std::endl;

}





void printPointNeighborhood(std::list<MeshVertexNode>::iterator i)
{
  std::cout << &*i << std::endl;
  for (std::list<std::list<MeshVertexNode>::iterator>::iterator g = i->neighbors().begin(); g != i->neighbors().end(); ++g)
  {
    //std::cout << (*g)->pos() << std::endl;
    std::cout << &**g << std::endl;
    for (std::list<std::list<MeshVertexNode>::iterator>::iterator j = (*g)->neighbors().begin(); j != (*g)->neighbors().end(); ++j)
    {
      std::cout << &**j << " ";
    }
    std::cout << std::endl;
  }
}



/// Returns the number of points that have been removed
template < typename TMesh, typename TCNeighborhoods, typename TIndexCollection >
std::size_t remove_vertices(TMesh & mesh, const TCNeighborhoods & neighc, const TIndexCollection & removed)
{
  // Convert mesh into a graph
  typedef til::MeshVertexNodeX<std::size_t> VertexNode;
  typedef til::MeshFaceNodeX<VertexNode> FaceNode;
  std::cout << "Converting into graph..." << std::flush;
  std::list<VertexNode> graph_vertices;
  std::list<FaceNode> graph_faces;
  til::list2graph_mesh_conversion(getVertices(mesh), getFaceIndices(mesh), neighc, graph_vertices, graph_faces);
  /*
  std::vector<std::list<std::size_t> > faceIndices = til::invertIndices(getFaceIndices(mesh));
  til::List2GraphMeshConvertor<VertexNode, FaceNode> l2g;
  l2g(getVertices(mesh), getFaceIndices(mesh), faceIndices, neighc);
  std::list<VertexNode> & graph_vertices = l2g.graph_vertices();
  std::list<FaceNode> & graph_faces = l2g.graph_faces();
  */
  std::cout << "OK" << std::endl;
  
  std::cout << "Adding attribute..." << std::flush;
  {
    std::list<VertexNode>::iterator it = graph_vertices.begin();
    std::size_t gvsize = graph_vertices.size();
    for (std::size_t i = 0; i < gvsize; ++i, ++it)
    {
      it->attribute() = i;
    }
  }
  std::cout << "OK" << std::endl;;
  
  std::cout << "Removing vertices..." << std::flush;
  std::size_t count = 0;
  bool complete;
  {
    int iter = 0;
    bool done;
    til::Vertex_remover<VertexNode, FaceNode> vertexRemover(graph_vertices, graph_faces);
    do
    {
      std::cout << "£" << std::flush;
      complete = true;
      done = true;
      std::list<VertexNode>::iterator i = graph_vertices.begin();
      std::list<FaceNode>::iterator itmp;
      std::cout << "Q" << std::flush;
      while (i != graph_vertices.end())
      {
        std::cout << "@" << std::flush;
        if (til::size(i->neighbors()) != til::size(i->faces()))
        {
          std::cout << "wO!" << til::size(i->neighbors()) << "-" << til::size(i->faces()) << "-";
        }
        //std::cout << i->attribute() << ":" << int(removed[i->attribute()]) << " ";
        if (removed[i->attribute()] == 1)
        {
          std::cout << "Y" << std::flush;
          //std::cout << "." << std::flush;
          // Here, the 'done' boolean should be here.
          // Note that it means that to be robust, I should add a check that the number of points
          // has strictly decreased during two iterations.
          complete = false;
          if (vertexRemover(i))
          {
            done = false;
            ++count;
            continue;
          }
        }
        ++i;
      }
      ++iter;
      std::cout << "pass " << iter << " : " << count << std::endl;
    } while (!done);
  }
  std::cout << "OK" << std::endl;
  std::cout << "Removed " << count << " points" << std::endl;
  if (!complete)
  {
    std::cout << "Warning: not all desired points could be removed" << std::endl;
  }
  
  std::cout << "Converting into mesh..." << std::flush;
  til::Graph2ListMeshConvertor<TMesh> g2l;
  g2l(graph_vertices, graph_faces);
  mesh = g2l.mesh();
  std::cout << "OK" << std::endl;

  /*
  {
    til::Mesh_N mesh1n = addNeighborsToMesh(mesh1);
    meshStats(mesh1n);
  }
  */
  return count;
}



void testGraphConvertion()
{
  til::Mesh_N mesh;
  
  std::cout << "Reading mesh..." << std::flush;
  {
    Reader<AimsSurfaceTriangle> r("/volatile/cathier/data/icbm/icbm100T_Lwhite.mesh");
    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  std::cout << "OK" << std::endl;
  std::size_t n = til::size(getVertices(mesh));
  std::cout << "mesh size " << n << std::endl;

  std::cout << "Get circular neighbor..." << std::flush;
  shared_ptr<std::vector<std::vector<std::size_t> > > neighc = circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  std::cout << "OK" << std::endl;
  
  
  std::cout << "Inverting face indices..." << std::flush;
  std::vector<std::list<std::size_t> > faceIndices = til::invertIndices(getFaceIndices(mesh));
  std::cout << "OK" << std::endl;
  {
    std::list<std::size_t>::const_iterator i = faceIndices[0].begin();
    for (; i != faceIndices[0].end(); ++i)
    {
      std::cout << getFaceIndices(mesh)[*i][0] << " " << getFaceIndices(mesh)[*i][1] << " " << getFaceIndices(mesh)[*i][2] << std::endl;
    }
  }
  {
    std::list<std::size_t>::const_iterator i = faceIndices.back().begin();
    for (; i != faceIndices.back().end(); ++i)
    {
      std::cout << getFaceIndices(mesh)[*i][0] << " " << getFaceIndices(mesh)[*i][1] << " " << getFaceIndices(mesh)[*i][2] << std::endl;
    }
  }
  
  std::cout << "Get curvature..." << std::flush;
  typedef float prec_type;
  til::Mesh_curvature<til::MeshTraits<til::Mesh_N>::VertexCollection, std::vector<std::vector<std::size_t> >, prec_type>
    mc(getVertices(mesh), *neighc);
  std::vector<float> curv(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    mc.process(i);
    curv[i] = min(prec_type(3000.0), 1000*mc.meanCurvature());
    curv[i] = max(float(-3000.0), curv[i]);
  }
  std::cout << "OK" << std::endl;
    
  std::cout << "Checking graph..." << std::flush;
  {
    for (std::size_t i = 0; i < til::size(getVertices(mesh)); ++i)
    {
      if (til::size(faceIndices[i]) != til::size((*neighc)[i]))
      {
        std::cout << "wR!" << til::size(faceIndices[i]) << "-" << til::size((*neighc)[i]) << "-";
      }
    }
  }
  std::cout << "OK" << std::endl;
  
  
  std::cout << "Converting into graph..." << std::flush;
  std::list<MeshVertexNode> graph_vertices;
  std::list<MeshFaceNode> graph_faces;
  til::list2graph_mesh_conversion(getVertices(mesh), getFaceIndices(mesh), *neighc, graph_vertices, graph_faces);
  /*
  til::List2GraphMeshConvertor<MeshVertexNode,MeshFaceNode> l2g;
  l2g(getVertices(mesh), getFaceIndices(mesh), faceIndices, *neighc);
  std::list<MeshVertexNode> & graph_vertices = l2g.graph_vertices();
  std::list<MeshFaceNode> & graph_faces = l2g.graph_faces();
  */
  
  /*
  {
    std::vector<std::list<MeshVertexNode>::iterator> index_vertex(til::size(getVertices(mesh)));
    std::vector<std::list<MeshFaceNode>::iterator> index_face(til::size(getFaceIndices(mesh)));
    {
      // First pass: push all points and get an index2iterator translation
      for (std::size_t i = 0; i < n; ++i)
      {
        MeshVertexNode m;
        m.pos() = getVertices(mesh)[i];
        index_vertex[i] = graph_vertices.insert(graph_vertices.end(), m);
      }
      // faces
      for (std::size_t i = 0; i < til::size(getFaceIndices(mesh)); ++i)
      {
        MeshFaceNode f;
        for (int j = 0; j < 3; ++j)
          f.face[j] = index_vertex[getFaceIndices(mesh)[i][j]];
        index_face[i] = graph_faces.insert(graph_faces.end(), f);
      }    
      // second pass: finish filling vertex structure
      for (std::size_t i = 0; i < n; ++i)
      {
        for (std::list<std::size_t>::const_iterator j = faceIndices[i].begin(); j != faceIndices[i].end(); ++j)      
        {
          index_vertex[i]->faces().push_back(index_face[*j]);
        }
        for (std::size_t j = 0; j < til::size((*neighc)[i]); ++j)
        {
          index_vertex[i]->neighbors().push_back(index_vertex[(*neighc)[i][j]]);
        }
      }
    }
  }
  */
  std::cout << "OK" << std::endl;
  
  /*
  {
    std::cout << "Removing three-connected vertices..." << std::flush;
    int count = 0;
    std::list<MeshVertexNode>::iterator i = graph_vertices.begin();
    int ind = 0;
    while (i != graph_vertices.end())
    {
      if (til::size(i->neighbors()) == 3)
      {
        ++count;

        // save list of neighbors
        std::list<std::list<MeshVertexNode>::iterator> neighbors = i->neighbors();

        // remove vertex
        i = remove_vertex(i, graph_vertices, graph_faces);

        // Add new face
        {
          MeshFaceNode f;
          {
            std::list<std::list<MeshVertexNode>::iterator>::const_iterator n = neighbors.begin();
            f.face[0] = *(n++);
            f.face[1] = *(n++);
            f.face[2] = *n;
          }
          std::list<MeshFaceNode>::iterator newf = graph_faces.insert(graph_faces.end(), f);
          {
            for (std::size_t n = 0; n < 3; ++n)
            {
              f.face[n]->faces().push_back(newf);
            }
          }
          // Check normal          
          std::cout << til::dot(
            til::cross(f.face[0]->pos()-f.face[1]->pos(),f.face[0]->pos()-f.face[2]->pos()),
            til::cross(i->faces().front()->face[0]->pos()-i->faces().front()->face[1]->pos(),i->faces().front()->face[0]->pos()-i->faces().front()->face[2]->pos())) << "*";
        }
        
        //std::cout << "Removed point of index " << ind << std::endl;
      }
      else
      {
        ++i;
      }
      ++ind;
    }
    std::cout << "OK" << std::endl;
    std::cout << "Removed " << count << " 3-connected points" << std::endl;
  }
  */
  std::cout << "Checking graph..." << std::flush;
  {
    std::list<MeshVertexNode>::iterator i = graph_vertices.begin();
    for (; i != graph_vertices.end(); ++i)
    {
      if (til::size(i->neighbors()) != til::size(i->faces()))
      {
        std::cout << "wM!" << til::size(i->neighbors()) << "-" << til::size(i->faces()) << "-";
      }
    }
  }
  std::cout << "OK" << std::endl;

  std::list<MeshVertexNode>::iterator i0 = graph_vertices.begin();
  {
    std::cout << "Removing 4-connected vertices..." << std::flush;
    int count = 0;
    bool done;
    int iter = 0;
    til::Vertex_remover<MeshVertexNode,MeshFaceNode> vertexRemover(graph_vertices, graph_faces);
    do
    {
      done = true;
      std::list<MeshVertexNode>::iterator i = graph_vertices.begin();
      std::list<MeshVertexNode>::iterator itmp;
      while (i != graph_vertices.end())
      {
        if (til::size(i->neighbors()) != til::size(i->faces()))
        {
          std::cout << "wO!" << til::size(i->neighbors()) << "-" << til::size(i->faces()) << "-";
        }
        
        if (til::size(i->neighbors()) <= 4)
        {
          if (vertexRemover(i))
          {
            done = false;
            ++count;
            continue;
          }
        }
        ++i;
      }
      ++iter;
      std::cout << "pass: " << count << std::endl;
    } while (!done);
    //} while (0);
    //} while (iter <= 8);
    {
      Texture1d t;
      t.reserve(graph_vertices.size());
      for (std::list<MeshVertexNode>::iterator j = graph_vertices.begin(); j != graph_vertices.end(); ++j)
      {
        if (i0 == j) t.push_back(100);
        else t.push_back(0);
      }
      til::aimswrite(t, "mods");
    }


    std::cout << "OK" << std::endl;
    std::cout << "Removed " << count << " 3-connected points" << std::endl;
  }
 
  
  std::cout << "Converting into mesh..." << std::flush;
  til::Graph2ListMeshConvertor<til::Mesh1> g2l;
  g2l(graph_vertices, graph_faces);
  til::Mesh1 & mesh1 = g2l.mesh();
  /*
  til::Mesh1 mesh1;
  {
    std::map<std::list<MeshVertexNode>::const_iterator, std::size_t, ItComp<std::list<MeshVertexNode>::const_iterator> > index_vertex;
    {
      getVertices(mesh1).resize(til::size(graph_vertices));
      std::size_t c = 0;
      for (std::list<MeshVertexNode>::const_iterator i = graph_vertices.begin(); i != graph_vertices.end(); ++i)
      {
        getVertices(mesh1)[c] = i->pos();
        //index_vertex[i] = c;
        index_vertex.insert(std::make_pair(i,c));
        ++c;
      }
    }
    {
      int c = 0;
      getFaceIndices(mesh1).resize(til::size(graph_faces));
      for (std::list<MeshFaceNode>::const_iterator i = graph_faces.begin(); i != graph_faces.end(); ++i)
      {
        for (int j = 0; j < 3; ++j)
          getFaceIndices(mesh1)[c][j] = index_vertex[i->face[j]];
        ++c;
      }
    }
  }
  */
  std::cout << "OK" << std::endl;

  {
    til::Mesh_N mesh1n = addNeighborsToMesh(mesh1);
    meshStats(mesh1n);
  }

  std::cout << "Writing mesh..." << std::flush;
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh1);
    til::aimswrite(s, "/home/cathier/tmp/biconverted");
  }
}

template < typename T2 >
double quant_dist_energy(const til::sparse_vector<double>::Map & distmap, const T2 & points, double distthresh)
{
  double e = 0.0;
  for (til::sparse_vector<double>::Map::const_iterator j = distmap.begin(); j != distmap.end(); ++j)
  {
    if (points[j->first])
    {
      double tmp = til::square(distthresh - j->second) / (0.001 + j->second);
      //double tmp = 1.0 / (0.001 + j->second) - 1.0 / (0.001 + distthresh);
      //double tmp = 1.0 / (0.001 + til::square(j->second)) - 1.0 / (0.001 + til::square(distthresh));
      //double tmp = til::square(j->second) - til::square(distthresh);
      e += tmp;
      /*
      e += til::square(distthresh) - til::square(j->second);
      if (tmp < 0)
      {
        std::cout << "wH!";
      }
      */
    }
  }
  return e;
}

void testSimplification(int argc, char * argv[])
{
  til::Mesh_N mesh;
  std::cout << "Reading mesh..." << std::flush;
  double kcurv;
  double distthresh = 30.0;
  std::size_t N = 3000;
  {
    Reader<AimsSurfaceTriangle> r;

    AimsApplication app( argc, aims_const_hack(argv), "testSimplification" );
    app.addOption(r, "-i", "input mesh" );
    app.addOption(kcurv, "-kcurv", "kcurv" );
    app.addOption(distthresh, "-distthresh", "distthresh" );
    app.addOption(N, "-N", "N" );
    app.initialize();

    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  std::size_t n = size(getVertices(mesh));

  std::cout << "Get circular neighbor..." << std::flush;
  shared_ptr<std::vector<std::vector<std::size_t> > > neighc = circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  std::cout << "OK" << std::endl;
  
  std::cout << "Get curvature..." << std::flush;
  til::Mesh_curvature<til::MeshTraits<til::Mesh_N>::VertexCollection, std::vector<std::vector<std::size_t> >, float>
    mc(getVertices(mesh), *neighc);
  std::vector<float> curv(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    mc.process(i);
    curv[i] = min(3.0f, mc.meanCurvature());
    curv[i] = max(-3.0f, curv[i]);
  }
  std::cout << "OK" << std::endl;

  //double distthresh = 30.0;
  std::vector<til::sparse_vector<double> > dist = getDistToNeighbors(mesh, distthresh);
  //std::vector<std::vector<double> > dist = getDistToNeighbors(mesh, 10.0);
  //std::vector<shared_ptr<til::sparse_vector<double> > > dist = getDistToNeighbors(mesh, 10.0);
  
  // Random initialization
  //std::vector<std::size_t> points(N);
  std::vector<std::size_t> center_labels(getVertices(mesh).size(), 0);
  std::vector<std::size_t> center_indices(N);
  {
    std::set<std::size_t> tmp;
    while (tmp.size() != N)
    {
      tmp.insert(til::rand<std::size_t>(0, getVertices(mesh).size()-1));
    }
    std::set<std::size_t>::const_iterator i = tmp.begin();
    std::vector<std::size_t>::iterator iLabel = center_indices.begin();
    for (; i != tmp.end(); ++i, ++iLabel)
    {
      *iLabel = *i;
    }
    for (std::size_t i = 0; i < N; ++i)
    {
      center_labels[center_indices[i]] = i;
    }
  }
  //std::vector<unsigned char> points2 = center_labels;
  
  std::cout << "Write random initialization..." << std::flush;
  {
    Texture1d tex(1,n);
    for (std::size_t i = 0; i < n; ++i)
    {
      tex.item(i) = center_labels[i];
    }
    til::aimswrite(tex, "/home/cathier/tmp/random_init");
  }
  std::cout << "OK" << std::endl;  


  {
    bool flagNotFinished;
    std::vector<std::size_t> labels(n);
    std::vector<unsigned char> points2(n);
    std::vector<double> minDist(N);
    std::vector<std::size_t> minIndex(N);
    int pass = 0;
    do
    {
      std::cout << "Pass " << ++pass << std::endl;

      /*
      std::cout << "Writing centers..." << std::flush;
      {
        Texture1d tex(1,n);
        for (std::size_t i = 0; i < n; ++i)
        {
          tex.item(i) = center_labels[i];
        }
        char s[256];
        sprintf(s, "/home/cathier/tmp/center_labels-%i", pass);
        Writer<Texture1d> w(s);
        w.write(tex);
      }
      std::cout << "OK" << std::endl;  
      */
      {
        int tmp = 0;
        for (std::size_t i = 0; i < n; ++i)
        {
          if (center_labels[i] > 0) ++tmp;
        }
        std::cout << "npoitns " << tmp << std::endl;
      }
      
      // Compute voronoi cells
      for (std::size_t i = 0; i < n; ++i)
      {        
        double mind = std::numeric_limits<double>::max();
        std::size_t mylabel = std::numeric_limits<std::size_t>::max();
        for (til::sparse_vector<double>::Map::const_iterator k = dist[i].getMap().begin(); k != dist[i].getMap().end(); ++k)
        {
          if (center_labels[k->first] != 0 && k->second < mind)
          {
            mylabel = center_labels[k->first];
            mind = k->second;
          }
        }
        labels[i] = mylabel;
      }

      /*
      std::cout << "Writing parcels..." << std::flush;
      {
        Texture1d tex(1,n);
        for (std::size_t i = 0; i < n; ++i)
        {
          tex.item(i) = labels[i];
        }
        char s[256];
        sprintf(s, "/home/cathier/tmp/parcels-%i", pass);
        Writer<Texture1d> w(s);
        w.write(tex);
      }
      std::cout << "OK" << std::endl;  
      */
      
      // Compute intra-cell distances
      std::fill(minDist.begin(), minDist.end(), std::numeric_limits<double>::max());
      for (std::size_t i = 0; i < n; ++i)
      {
        // Skip unlabeled points
        if (labels[i] >= N) continue;
        int count = 0;
        double d = 0.0;
        for (til::sparse_vector<double>::Map::const_iterator k = dist[i].getMap().begin(); k != dist[i].getMap().end(); ++k)
        {
          if (labels[k->first] == labels[i])
          {
            ++count;
            d += til::square(k->second);
          }
        }
        if (count == 0)
        {
          std::cout << "wZ!";
          d = std::numeric_limits<double>::max();
        }
        else
        {
          d /= count;
        }
        d -= kcurv * curv[i];
        if (d < minDist[labels[i]])
        {
          minDist[labels[i]] = d;
          minIndex[labels[i]] = i;
        }
      }
      
      flagNotFinished = false;
      for (std::size_t i = 0; i < N; ++i)
      {
        if (center_indices[i] != minIndex[i])
        {
          flagNotFinished = true;
          break;
        }
      }
      
      if (flagNotFinished)
      {
        center_indices = minIndex;
        std::fill(center_labels.begin(), center_labels.end(), 0);
        for (std::size_t i = 0; i < N; ++i)
        {
          center_labels[center_indices[i]] = i;
        }        
      }
      
      /*
      std::fill(points2.begin(), points2.end(), 0);
      for (std::size_t i = 0; i < N; ++i)
      {
        points2[minIndex[i]] = 1;
      }
      
      flagNotFinished = false;
      for (std::size_t i = 0; i < n; ++i)
      {
        if (center_labels[i] != points2[i])
        {
          flagNotFinished = true;
          break;
        }
      }
      */
      
    } while (flagNotFinished);
  }

  std::cout << "Writing centers..." << std::flush;
  {
    Texture1d tex(1,n);
    for (std::size_t i = 0; i < n; ++i)
    {
      tex.item(i) = (center_labels[i] > 0);
    }
    til::aimswrite(tex, "/home/cathier/tmp/quant");
  }
  std::cout << "OK" << std::endl;  

  /*  
  std::cout << "Writing parcels..." << std::flush;
  {
    Texture1d tex(1,n);
    for (std::size_t i = 0; i < n; ++i)
    {
      tex.item(i) = parcels[i];
    }
    Writer<Texture1d> w("/home/cathier/tmp/parcels");
    w.write(tex);
  }
  std::cout << "OK" << std::endl;  
  */


  /*  
  {
    bool flagNotFinished;
    int count = 0;
    do
    {
      std::cout << "Pass " << ++count << std::endl;
      
      if (count > 40) break;
      
      flagNotFinished = false;
      
      double eTot = 0.0;
      
      for (std::size_t i = 0; i < points.size(); ++i)
      {
        if (points[i] == 0) continue;
        
        double e = quant_dist_energy(dist[i].getMap(), points, distthresh) - kcurv * til::square(curv[i]);

        double emin = std::numeric_limits<double>::max();
        std::size_t kmin = 0;
        for (std::vector<std::size_t>::const_iterator k = (*neighc)[i].begin(); k != (*neighc)[i].end(); ++k)
        {
          if (points[*k] == 1) continue;
          double tmp = quant_dist_energy(dist[*k].getMap(), points, distthresh) - kcurv * til::square(curv[*k]);
          if (tmp < emin)
          {
            emin = tmp;
            kmin = *k;
          }
        }
        
        if (emin < e)
        {
          e = emin;
          points[i] = 0;
          assert(points[kmin] == 0);
          points[kmin] = 1;
          flagNotFinished = true;
        }
       eTot += e;
      }
      std::cout << "Energy: " << eTot << std::endl;
    } while (flagNotFinished);
  }
  std::cout << "Remaining points : " << std::accumulate(points.begin(), points.end(), 0) << std::endl;

  std::cout << "Write final quantization..." << std::flush;
  {
    Texture1d tex(1,n);
    for (std::size_t i = 0; i < n; ++i)
    {
      tex.item(i) = points[i];
    }
    Writer<Texture1d> w("/home/cathier/tmp/final_quant");
    w.write(tex);
  }
  std::cout << "OK" << std::endl;  

  {
    bool flagNotFinished;
    int count = 0;
    do
    {
      std::cout << "Pass " << ++count << std::endl;

      if (count > 40) break;
      
      flagNotFinished = false;
      
      for (std::size_t i = 0; i < points.size(); ++i)
      {
        if (points[i] == 0) continue;
        
        double e = quant_dist_energy(dist[i].getMap(), points, distthresh) - kcurv * til::square(curv[i]);

        double emin = std::numeric_limits<double>::max();
        std::size_t kmin = 0;
        //for (std::vector<std::size_t>::const_iterator k = (*neighc)[i].begin(); k != (*neighc)[i].end(); ++k)
        for (til::sparse_vector<double>::Map::const_iterator k = dist[i].getMap().begin(); k != dist[i].getMap().end(); ++k)
        {
          if (points[k->first] == 1) continue;
          double tmp = quant_dist_energy(dist[k->first].getMap(), points, distthresh) - kcurv * til::square(curv[k->first]);
          if (tmp < emin)
          {
            emin = tmp;
            kmin = k->first;
          }
        }
        if (emin < e)
        {
          points[i] = 0;
          assert(points[kmin] == 0);
          points[kmin] = 1;
          flagNotFinished = true;
        }
      }
    } while (flagNotFinished);
  }
  std::cout << "Remaining points : " << std::accumulate(points.begin(), points.end(), 0) << std::endl;

  std::cout << "Write final quantization (2)..." << std::flush;
  {
    Texture1d tex(1,n);
    for (std::size_t i = 0; i < n; ++i)
    {
      tex.item(i) = points[i];
    }
    Writer<Texture1d> w("/home/cathier/tmp/final_quant2");
    w.write(tex);
  }
  std::cout << "OK" << std::endl;  
  */ 
}


void testSimplification2(int argc, char * argv[])
{
  til::Mesh_N mesh;
  std::cout << "Reading mesh..." << std::flush;
  double kcurv;
  double distthresh = 30.0;
  std::size_t N = 3000;
  {
    Reader<AimsSurfaceTriangle> r;

    AimsApplication app( argc, aims_const_hack(argv), "testSimplification2" );
    app.addOption(r, "-i", "input mesh" );
    app.addOption(kcurv, "-kcurv", "kcurv" );
    app.addOption(distthresh, "-distthresh", "distthresh" );
    app.addOption(N, "-N", "N" );
    app.initialize();

    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  std::size_t n = size(getVertices(mesh));

  std::cout << "Get circular neighbor..." << std::flush;
  shared_ptr<std::vector<std::vector<std::size_t> > > neighc = circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  std::cout << "OK" << std::endl;
  
  std::cout << "Get curvature..." << std::flush;
  til::Mesh_curvature<til::MeshTraits<til::Mesh_N>::VertexCollection, std::vector<std::vector<std::size_t> >, float>
    mc(getVertices(mesh), *neighc);
  std::vector<float> curv(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    mc.process(i);
    curv[i] = min(3.0f, mc.meanCurvature());
    curv[i] = max(-3.0f, curv[i]);
  }
  std::cout << "OK" << std::endl;

  //double distthresh = 30.0;
  std::vector<til::sparse_vector<double> > dist = getDistToNeighbors(mesh, distthresh);
  //std::vector<std::vector<double> > dist = getDistToNeighbors(mesh, 10.0);
  //std::vector<shared_ptr<til::sparse_vector<double> > > dist = getDistToNeighbors(mesh, 10.0);
  
  /*
  std::cout << "Extracting edges..." << std::flush;
  typedef til::numeric_array<std::size_t,2> Edge;
  std::vector<Edge> edges;
  getEdges(mesh, edges);
  std::cout << "OK" << std::endl;
  
  std::cout << "Number of edges: " << edges.size() << std::endl;
  
  typedef std::priority_queue<std::pair<Edge, double>, std::vector<std::pair<Edge, double> >, til::Lesser_Pair2<Edge, double> >
    EdgeQueue;
  EdgeQueue q;
  for (std::vector<Edge>::const_iterator iEdge = edges.begin(); iEdge != edges.end(); ++iEdge)
  {
    q.push(std::make_pair(*iEdge, til::dist2(getVertices(mesh)[(*iEdge)[0]], getVertices(mesh)[(*iEdge)[1]], til::prec<double>())));
  }
  */
}

void testQuantization(int argc, char * argv[])
{
  til::Mesh_N mesh;
  std::cout << "Reading mesh..." << std::flush;
  {
    Reader<AimsSurfaceTriangle> r;

    AimsApplication app( argc, aims_const_hack(argv), "testQuantization" );
    app.addOption(r, "-i", "input mesh" );
    app.initialize();

    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  /*
  {
    AimsSurfaceTriangle * s = aims::SurfaceGenerator::sphere(Point3df(0,0,0), 100.0f, 100);
    til::Mesh1 mesh0;
    til::convert(mesh0, *s);
    mesh = addNeighborsToMesh(mesh0);
    delete s;    
  }
  */
  /*
  {
    AimsSurfaceTriangle *s = makeSphere(Point3df(1.0, 2.0, 3.0), 100, 5);
    til::Mesh1 mesh0;
    til::convert2(*s).into(mesh0);
    mesh = addNeighborsToMesh(mesh0);
  }
  */
  
  std::size_t n = size(getVertices(mesh));
  
  printBoundingBox(getVertices(mesh));
  
  /*
  std::cout << "Writing mesh..." << std::flush;
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh);
    Writer<AimsSurfaceTriangle> w("/home/cathier/tmp/sphere");
    w.write(s);
  }
  std::cout << "OK" << std::endl;
  */
  
  std::cout << "Get circular neighbor..." << std::flush;
  shared_ptr<std::vector<std::vector<std::size_t> > > neighc = circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  std::cout << "OK" << std::endl;
  
  //for (std::size_t i = 0; i < til::size((*neighc)[0]); ++i) std::cout << (*neighc)[0][i] << " "; std::cout << std::endl;
  //for (std::size_t i = 0; i < til::size((*neighc)[n-1]); ++i) std::cout << (*neighc)[n-1][i] << " "; std::cout << std::endl;
  
  
  std::cout << "Get curvature..." << std::flush;
  til::Mesh_curvature<til::MeshTraits<til::Mesh_N>::VertexCollection, std::vector<std::vector<std::size_t> >, float>
    mc(getVertices(mesh), *neighc);
  std::vector<float> curv(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    mc.process(i);
    curv[i] = min(3.0f, mc.meanCurvature());
    curv[i] = max(-3.0f, curv[i]);
  }
  std::cout << "OK" << std::endl;

  std::cout << "Write curvature..." << std::flush;
  {
    Texture1d tex(1,n);
    for (std::size_t i = 0; i < n; ++i)
    {
      tex.item(i) = curv[i];
    }
    //til::convert(tex, curv);
    til::aimswrite(tex, "/home/cathier/tmp/curv");
  }
  std::cout << "OK" << std::endl;
  
  
  // computing initial geodesic distance
  double distthresh = 5.0;
  std::vector<til::sparse_vector<double> > dist = getDistToNeighbors(mesh, distthresh);  
  //std::vector<std::vector<double> > dist = getDistToNeighbors(mesh, 10.0);
  //std::vector<shared_ptr<til::sparse_vector<double> > > dist = getDistToNeighbors(mesh, 10.0);
   
  std::cout << "Write dist..." << std::flush;
  /*
  {    
    Texture1d t(1, til::size(getVertices(mesh)));
    for (std::size_t i = 0; i < til::size(dist[0]); ++i)
    {
      t.item(i) = dist[0][i];
    }
    Writer<Texture1d> w("/home/cathier/tmp/toubou");
    w.write(t);
  }
  */  
  {
    Texture1d t(1, til::size(getVertices(mesh)));
    til::sparse_vector<double>::Map::const_iterator i = dist[0].getMap().begin();
    std::cout << getVertices(mesh)[0] << std::endl;
    for (; i != dist[0].getMap().end(); ++i)
    {
      t.item(i->first) = i->second;
      std::cout << getVertices(mesh)[i->first] << " " << i->second << std::endl;
    }
    til::aimswrite(t, "/home/cathier/tmp/toubou");
  }
  std::cout << "OK" << std::endl;

  int count = 0;
  std::cout << "Symmetrizing distance-based neighborhood..." << std::flush;
  {
    for (std::size_t i = 0; i < til::size(dist); ++i)
    {
      til::sparse_vector<double>::Map::iterator j = dist[i].getMap().begin();
      while (j != dist[i].getMap().end())
      {
        if (dist[j->first].getMap().find(i) == dist[j->first].getMap().end())
        {
          dist[i].erase(j++);
          ++count;
        }
        else
        {
          ++j;
        }
      }
    }
  }
  std::cout << "OK" << std::endl;

  std::cout << "Removed " << count << " unsymmetric neighbors" << std::endl;
  {
    for (std::size_t i = 0; i < til::size(dist); ++i)
    {
      til::sparse_vector<double>::Map::iterator j = dist[i].getMap().begin();
      for (; j != dist[i].getMap().end(); ++j)
      {
        if (dist[j->first].getMap().find(i) == dist[j->first].getMap().end())
        {
          std::cout << "w1!" << std::flush;
        }
      }
    }
  }
  
  // this quantization removes a point and locks its neighbors.
  // obviously, the criterion needs to hold on a neighborhood only, otherwise, strictly speaking, a recomputation
  // of the criterion would be needed everytime a point is removed.
  std::cout << "Vector quantization..." << std::endl;
  std::vector<unsigned char> removed(n, 0);
  count = 0;
  {
    // loop through all vertices
    for (std::size_t i = 0; i < n; ++i)
    {
      // skip points that are not labeled as 'unprocessed'
      if (removed[i]) continue;
      //double maxdist = 0.0;
      double mindist = std::numeric_limits<double>::max();
      std::size_t i2 = i;
      bool flag;
      do
      {
        // remove current point
        removed[i2] = 2;
        // check that none of its neighbors are labeled as being removed -- this would be a bug
        {
          for (std::size_t t = 0; t < (*neighc)[i2].size(); ++t)
          {
            if (removed[(*neighc)[i2][t]] == 2)
            {
              std::cout << "w8! " << i2 << " " << (*neighc)[i2][t] << " " << dist[i2][(*neighc)[i2][t]] << " " << dist[(*neighc)[i2][t]][i2] << std::endl;
            }
          }
        }
        flag = false;
        // Look for a point in the dist-neighborhood that minimize a criteria for being the next point to be removed
        til::sparse_vector<double>::Map::iterator j = dist[i2].getMap().begin();
        til::sparse_vector<double>::Map::iterator maxPos;
        for (; j != dist[i2].getMap().end(); ++j)
        {
          if (removed[j->first]) continue;
          
          float myc = (til::size((*neighc)[j->first])>3 ? abs(curv[j->first]) : 0);
          //if (maxdist < (j->second + 5*myc) && j->second > 0.7 * distthresh)
          float e = max(0.0, distthresh - j->second) * max(0.0f, 1 - myc);
          if (mindist > e && j->second > 0.7 * distthresh)
          {
            if (flag)
            {
              if (removed[maxPos->first])
              {
                std::cerr << "w6!" << std::flush;
              }
              else
              {
                ++count;
                removed[maxPos->first] = 1;
              }
            }
            mindist = e;
            //maxdist = e;
            maxPos = j;
            flag = true;
          }
          else
          {
            ++count;
            if (removed[j->first])
            {
              std::cerr << "w7!" << std::endl;
            }
            removed[j->first] = 1;
          }
        }
        if (flag)
        {
          i2 = maxPos->first;
          if (removed[i2]) std::cerr << "w5!" << std::flush;
        }
      } while (flag);
    }
  }
  std::cout << "OK" << std::endl;

  /*
  // this quantization removes a point and lock its neighbors, and then goes to a point in the dist-neighbor
  // of the current point, if possible.
  std::cout << "Vector quantization..." << std::endl;
  std::vector<unsigned char> removed(n, 0);
  count = 0;
  {
    // loop through all vertices
    for (std::size_t i = 0; i < n; ++i)
    {
      // skip points that are not labeled as 'unprocessed'
      if (removed[i]) continue;
      //double maxdist = 0.0;
      double mindist = std::numeric_limits<double>::max();
      std::size_t i2 = i;
      bool flag;
      do
      {
        // remove current point
        removed[i2] = 2;
        // check that none of its neighbors are labeled as being removed -- this would be a bug
        {
          for (std::size_t t = 0; t < (*neighc)[i2].size(); ++t)
          {
            if (removed[(*neighc)[i2][t]] == 2)
            {
              std::cout << "w8! " << i2 << " " << (*neighc)[i2][t] << " " << dist[i2][(*neighc)[i2][t]] << " " << dist[(*neighc)[i2][t]][i2] << std::endl;
            }
          }
        }
        flag = false;
        // Look for a point in the dist-neighborhood that minimize a criteria for being the next point to be removed
        til::sparse_vector<double>::Map::iterator j = dist[i2].getMap().begin();
        til::sparse_vector<double>::Map::iterator maxPos;
        for (; j != dist[i2].getMap().end(); ++j)
        {
          if (removed[j->first]) continue;
          
          float myc = (til::size((*neighc)[j->first])>3 ? abs(curv[j->first]) : 0);
          //if (maxdist < (j->second + 5*myc) && j->second > 0.7 * distthresh)
          float e = max(0.0, distthresh - j->second) * max(0.0f, 1 - myc);
          if (mindist > e && j->second > 0.7 * distthresh)
          {
            if (flag)
            {
              if (removed[maxPos->first])
              {
                std::cerr << "w6!" << std::flush;
              }
              else
              {
                ++count;
                removed[maxPos->first] = 1;
              }
            }
            mindist = e;
            //maxdist = e;
            maxPos = j;
            flag = true;
          }
          else
          {
            ++count;
            if (removed[j->first])
            {
              std::cerr << "w7!" << std::endl;
            }
            removed[j->first] = 1;
          }
        }
        if (flag)
        {
          i2 = maxPos->first;
          if (removed[i2]) std::cerr << "w5!" << std::flush;
        }
      } while (flag);
    }
  }
  std::cout << "OK" << std::endl;
  */
  
  /*
  // This quantization removed a point and locks its neighbors
  std::cout << "Vector quantization..." << std::endl;
  std::vector<unsigned char> removed(n, 0);
  count = 0;
  {
    std::size_t k = 3;
    std::size_t s;
    bool active;
    std::size_t maxK;
    std::size_t minK;
    do
    {
      std::cout << "k = " << k << " " << count << std::endl;
      active = false;
      maxK = 0;
      minK = std::numeric_limits<int>::max();
      for (std::size_t i = 0; i < til::size(dist); ++i)
      {
        if (removed[i]) continue;

        s = til::size(dist[i].getMap());
        maxK = max(maxK, s);
        minK = min(minK, s);
        if (s > k) continue;

        bool del = true;
        til::sparse_vector<double>::Map::iterator j = dist[i].getMap().begin();
        for (; j != dist[i].getMap().end(); ++j)
        {
          //if (j->first == i) continue;
          int c = 0;
          til::sparse_vector<double>::Map::iterator x = dist[j->first].getMap().begin();
          for (; x != dist[j->first].getMap().end(); ++x)
          {
            //if (x->first == j->first) continue;
            if (!removed[x->first])
            {
              if (++c == 2) break;
            }
          }
          if (c == 1)
          {
            del = false;
            break;
          }
        }
        if (!del) continue;
        active = true;
        removed[i] = 1;
        ++count;
      }
      if (!active) ++k;
    } while (active || k <= maxK);
  }
  std::cout << "OK" << std::endl;
  //*/
  
  std::cout << "Before: " << til::size(dist) << std::endl;
  std::cout << "Removed: " << count << std::endl;
  std::cout << "Now: " << til::size(dist) - count << std::endl;
  count = 0;
  for (std::size_t i = 0; i < til::size(dist); ++i)
  {
    if (removed[i] == 1) ++count;
  }
  std::cout << "count: " << count << std::endl;
  
  std::cout << "Write remaining points..." << std::flush;
  {
    Texture1d tex(1,n);
    for (std::size_t i = 0; i < n; ++i)
    {
      tex.item(i) = removed[i];
    }
    til::aimswrite(tex, "/home/cathier/tmp/removed");
  }
  std::cout << "OK" << std::endl;  

  std::cout << "Converting into graph..." << std::flush;

  typedef til::MeshVertexNodeX<std::size_t> VertexNode;
  typedef til::MeshFaceNodeX<VertexNode> FaceNode;
  std::list<VertexNode> graph_vertices;
  std::list<FaceNode> graph_faces;
  til::list2graph_mesh_conversion(getVertices(mesh), getFaceIndices(mesh), *neighc, graph_vertices, graph_faces);

  /*    
  std::cout << "Inverting face indices..." << std::flush;
  std::vector<std::list<std::size_t> > faceIndices = til::invertIndices(getFaceIndices(mesh));
  std::cout << "OK" << std::endl;
  
  til::List2GraphMeshConvertor<VertexNode, FaceNode> l2g;
  l2g(getVertices(mesh), getFaceIndices(mesh), faceIndices, *neighc);
  std::list<VertexNode> & graph_vertices = l2g.graph_vertices();
  std::list<FaceNode> & graph_faces = l2g.graph_faces();
  */
  
  /*
  std::list<MeshVertexNode> graph_vertices;
  std::list<MeshFaceNode> graph_faces;
  listmesh2graphmesh(getVertices(mesh), getFaceIndices(mesh), *neighc, faceIndices, graph_vertices, graph_faces);
  */
  std::cout << "OK" << std::endl;
  
  std::cout << graph_vertices.size() << std::endl;
  
  std::cout << "Adding attribute..." << std::flush;
  std::list<VertexNode>::iterator it = graph_vertices.begin();
  {
    std::size_t gvsize = graph_vertices.size();
    for (std::size_t i = 0; i < gvsize; ++i, ++it)
    {
      it->attribute() = i;
    }
  }
  std::cout << "OK" << std::endl;;
  
    
  std::cout << "Removing vertices..." << std::flush;
  count = 0;
  bool complete;
  {
    int iter = 0;
    bool done;
    til::Vertex_remover<VertexNode, FaceNode> vertexRemover(graph_vertices, graph_faces);
    do
    {
      complete = true;
      done = true;
      std::list<VertexNode>::iterator i = graph_vertices.begin();
      std::list<FaceNode>::iterator itmp;
      while (i != graph_vertices.end())
      {
        if (i->neighbors().size() != i->faces().size())
        {
          std::cout << "wO!" << til::size(i->neighbors()) << "-" << til::size(i->faces()) << "-";
        }
        //std::cout << i->attribute() << ":" << int(removed[i->attribute()]) << " ";
        if (removed[i->attribute()] == 1)
        {
          //std::cout << "." << std::flush;
          // Here, the 'done' boolean should be here.
          // Note that it means that to be robust, I should add a check that the number of points
          // has strictly decreased during two iterations.
          complete = false;
          if (vertexRemover(i))
          {
            done = false;
            ++count;
            continue;
          }
        }
        ++i;
      }
      ++iter;
      std::cout << "pass " << iter << " : " << count << std::endl;
    } while (!done);
  }
  std::cout << "OK" << std::endl;
  std::cout << "Removed " << count << " points" << std::endl;
  if (!complete)
  {
    std::cout << "Warning: not all desired points could be removed" << std::endl;
  }
  
  std::cout << "Converting into mesh..." << std::flush;
  til::Graph2ListMeshConvertor<til::Mesh1> g2l;
  g2l(graph_vertices, graph_faces);
  til::Mesh1 & mesh1 = g2l.mesh();
  std::cout << "OK" << std::endl;

  {
    til::Mesh_N mesh1n = addNeighborsToMesh(mesh1);
    meshStats(mesh1n);
  }

  std::cout << "Writing mesh..." << std::flush;
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh1);
    til::aimswrite(s, "/home/cathier/tmp/biconverted");
  }
}


// void removingPointsWithHighCurvature(int argc, char * argv[])
// moved into commands/meshCleaner/meshCleaner.cc


void removingPointsWithHighCurvature2(int argc, char * argv[])
{
  til::Mesh_N mesh;
  std::cout << "Reading mesh..." << std::flush;
  Writer<AimsSurfaceTriangle> w;
  int N;
  {
    Reader<AimsSurfaceTriangle> r;
    AimsApplication app( argc, aims_const_hack(argv), "removingPointsWithHighCurvature");
    app.addOption(r, "-i", "input mesh" );
    app.addOption(w, "-o", "output mesh" );
    app.addOption(N, "-n", "");
    app.initialize();

    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  std::size_t n =  getVertices(mesh).size();

  std::cout << "Get circular neighbor..." << std::flush;
  shared_ptr<std::vector<std::vector<std::size_t> > > neighc = circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  std::cout << "OK" << std::endl;

  std::cout << "Converting into graph..." << std::flush;
  typedef til::MeshVertexNodeX<unsigned int>  VertexNode;
  typedef til::MeshFaceNodeX<VertexNode>      FaceNode;
  typedef std::list<VertexNode>::iterator     iterator;
  std::list<VertexNode> graph_vertices;
  std::list<FaceNode> graph_faces;
  til::list2graph_mesh_conversion(getVertices(mesh), getFaceIndices(mesh), *neighc, graph_vertices, graph_faces);
  /*
  std::vector<std::list<std::size_t> > faceIndices = til::invertIndices(getFaceIndices(mesh));
  til::List2GraphMeshConvertor<VertexNode, FaceNode> l2g;
  l2g(getVertices(mesh), getFaceIndices(mesh), faceIndices, *neighc);
  std::list<VertexNode> & graph_vertices = l2g.graph_vertices();
  std::list<FaceNode> & graph_faces = l2g.graph_faces();
  */
  std::cout << "OK" << std::endl;

  std::cout << "Adding attribute..." << std::flush;
  {
    for (iterator it = graph_vertices.begin(); it != graph_vertices.end(); ++it)
    {
      it->attribute() = 0;
    }
  }
  /*
  unsigned long top = 0;
  {
    std::list<VertexNode>::iterator it = graph_vertices.begin();
    std::size_t gvsize = graph_vertices.size();
    for (; top < gvsize; ++top, ++it)
    {
      it->attribute() = top;
    }
  }
  */
  std::cout << "OK" << std::endl;

  // Declare curvature calculator
  typedef float prec_type;
  typedef til::xsr::graph::Position<VertexNode::VertexIndex> A;
  typedef til::xsr::graph::Neighbors<VertexNode::VertexIndex> B;
  A a; B b; til::MeshCurvature2<A,B,float> mc(a,b);
  // NB: this hack doesn't work in gcc3.3 :(
  //til::MeshCurvature2<A,B,prec_type> mc((A()),(B()));
  
  // Declare curvature queue
  typedef std::pair<iterator, unsigned int> CurvIndex;
  typedef std::pair<CurvIndex, prec_type> IndexedCurv;
  std::priority_queue<IndexedCurv, std::vector<IndexedCurv>, til::Lesser_Pair2<CurvIndex, prec_type> >
    curvqueue;

  std::cout << "Initializing queue..." << std::flush;
  for (iterator i = graph_vertices.begin(); i != graph_vertices.end(); ++i)
  {
    mc.process(i);
    curvqueue.push(std::make_pair(std::make_pair(i, i->attribute()), mc.unorientedMeanCurvature()));
  }
  std::cout << "OK" << std::endl;

  unsigned int count = 0;
  {
    std::set<iterator, til::Lesser_PointeeAddress<iterator> > removedVertices;
    til::Vertex_remover<VertexNode, FaceNode> vertexRemover(graph_vertices, graph_faces);
    while (count < n-N)
    {
      //std::cout << "dodi" << std::flush;
      IndexedCurv curv = curvqueue.top();
      curvqueue.pop();
      //std::cout << "doda" << std::flush;
      // Check whether curvature is obsolete
      if (removedVertices.find(curv.first.first) != removedVertices.end() ||
          curv.first.second != curv.first.first->attribute()
      )
      {
        //std::cout << "Removed curv " << curv.second << std::endl;
        continue;
      }
      //std::cout << "dodu " << &*curv.first.first << std::endl;
      //std::cout << "curv " << curv.second << std::endl;
      // Try to remove point
      if (vertexRemover.isRemovable(curv.first.first))
      {
        // Add fraction of position to neighbors
        for (std::vector<iterator>::iterator iNeigh = vertexRemover.neighbors().begin(); iNeigh != vertexRemover.neighbors().end(); ++iNeigh)
        {
          (*iNeigh)->pos() *= 1.0f - 1.0f/vertexRemover.neighbors().size();
          (*iNeigh)->pos() += curv.first.first->pos() * (1.0f / vertexRemover.neighbors().size());
        }
        //std::cout << "dodon" << std::flush;
        // Signal point as removed
        removedVertices.insert(curv.first.first);
        // Remove point
        vertexRemover.remove(curv.first.first);

        
        // Recompute the curvature of neighbors
        for (std::vector<iterator>::iterator iNeigh = vertexRemover.neighbors().begin(); iNeigh != vertexRemover.neighbors().end(); ++iNeigh)
        {
          //std::cout << "e" << std::endl;
          mc.process(*iNeigh);
          //std::cout << "adding curv " << mc.unorientedMeanCurvature() << std::endl;
          curvqueue.push(std::make_pair(std::make_pair(*iNeigh, ++(*iNeigh)->attribute()), mc.unorientedMeanCurvature()));
        }
        // if successfull, increase counter
        //std::cout << "f" << std::endl;
        ++count;
      }
      /*
      else
      {
        std::cout << "skip" << std::flush;
      }
      */
    }
  }
  std::cout << "Remove " << count << " points" << std::endl;
  {
    Texture1d t(1, graph_vertices.size());
    iterator ivertex = graph_vertices.begin();
    std::size_t i = 0;
    for (; ivertex != graph_vertices.end(); ++i, ++ivertex)
    {
      t.item(i) = (ivertex->attribute() == 1234 ? 1 : 0);
    }
    til::aimswrite(t, "/home/cathier/tmp/toubou");
  }

  std::cout << "Converting into mesh..." << std::flush;
  til::Graph2ListMeshConvertor<til::Mesh1> g2l;
  g2l(graph_vertices, graph_faces);
  til::Mesh1 & mesh1 = g2l.mesh();
  std::cout << "OK" << std::endl;

  std::cout << "Writing mesh..." << std::flush;
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh1);
    w.write(s);
  }
}


void testQuantization3(int argc, char * argv[])
{
  til::Mesh_N mesh;
  std::cout << "Reading mesh..." << std::flush;
  std::size_t N;
  Writer<AimsSurfaceTriangle> w;
  {
    Reader<AimsSurfaceTriangle> r;
    AimsApplication app( argc, aims_const_hack(argv), "testQuantization" );
    app.addOption(r, "-i", "input mesh" );
    app.addOption(w, "-o", "output mesh" );
    app.addOption(N, "-n", "desired number of vertices");
    app.initialize();

    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  std::cout << "OK" << std::endl;

  meshStats(mesh);
    
  // this quantization removes a point and locks its neighbors.
  // obviously, the criterion needs to hold on a neighborhood only, otherwise, strictly speaking, a recomputation
  // of the criterion would be needed everytime a point is removed.
  std::cout << "Vector quantization..." << std::endl;
  std::size_t bigCount = 0;
  do
  {
    std::size_t n = getVertices(mesh).size();
    std::vector<unsigned char> removed(n, 0);
    // computing circular neighborhoods
    shared_ptr<std::vector<std::vector<std::size_t> > > neighc = circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
    // compute mesh curvature
    /*
    til::Mesh_curvature<til::MeshTraits<til::Mesh_N>::VertexCollection, std::vector<std::vector<std::size_t> >, float>
      mc(getVertices(mesh), *neighc);
    */
    til::MeshWaveletEnergy<til::MeshTraits<til::Mesh_N>::VertexCollection, std::vector<std::vector<std::size_t> >, float>
      mwe(getVertices(mesh), *neighc);
    std::vector<std::pair<std::size_t, float> > sortedCurv(n);
    for (std::size_t i = 0; i < n; ++i)
    {
      /*
      mc.process(i);
      sortedCurv[i] = std::make_pair(i, min(1.0f, std::abs(mc.gaussianCurvature())));
      */
      sortedCurv[i] = std::make_pair(i, mwe(i));
    }
    /*
    {
      Texture1d t;
      t.reserve(n);
      for (std::size_t i=0; i<n; ++i)
      {
        if (sortedCurv[i].second < 0)
          t.push_back(-std::sqrt(-sortedCurv[i].second));
        else
          t.push_back(std::sqrt(sortedCurv[i].second));
      }
      Writer<Texture1d> w("mwe");
      w.write(t);
    }    
    exit(0);
    */
    //*/
    // order points according to their curvatures, higher curvatures first.
    std::sort(sortedCurv.begin(), sortedCurv.end(), til::Greater_Pair2<std::size_t, float>());
    //for (std::size_t i = 0; i < n; ++i) cout << sortedCurv[i].second << " ";
    //std::cout << std::endl;
    // loop through all vertices, higher curvature first
    std::size_t count = 0;
    for (std::size_t i0 = 0; i0 < n; ++i0)
    {
      // get real index of current vertex
      std::size_t i = sortedCurv[i0].first;
      // skip point if not labeled as 'unprocessed'
      if (removed[i]) continue;
      // label point as removed
      removed[i] = 2;
      ++count;
      if (n - count < N)
      {
        std::cout << "Warning: reached max number" << std::endl;
        break;
      }
      // mark all of its neighbors as kept
      for (std::size_t j = 0; j < (*neighc)[i].size(); ++j)
      {
        if (removed[(*neighc)[i][j]] == 2)
        {
          std::cout << "w8! ";
        }
        removed[(*neighc)[i][j]] = 1;
      }
      // add a fraction of the displacement to neighbors
      {
        for (std::size_t j = 0; j < (*neighc)[i].size(); ++j)
        {
          getVertices(mesh)[(*neighc)[i][j]] += (getVertices(mesh)[i] - getVertices(mesh)[(*neighc)[i][j]]) * (1.0f / (*neighc)[i].size());
        }
      }
    }
    // Put 'removed' in standard binary form
    for (std::size_t i = 0; i < n; ++i) removed[i] = (removed[i] == 2 ? 1 : 0) ;
    std::cout << "Marked " << count << " points out of " << n << " for deletion" << std::endl;
    count = remove_vertices(mesh, *neighc, removed);
    std::cout << count << "points could be effectively removed" << std::endl;
    bigCount += count;
    mesh.getNeighborIndices() = getNeighborIndices(mesh);
  //} while (0);
  //} while (bigCount < N);
  } while (getVertices(mesh).size() > N);
  std::cout << "OK" << std::endl;

  std::cout << "Writing mesh..." << std::flush;
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh);
    w.write(s);
  }
}


void testQuantization2(int argc, char * argv[])
{
  til::Mesh_N mesh;
  std::cout << "Reading mesh..." << std::flush;
  {
    Reader<AimsSurfaceTriangle> r;

    AimsApplication app( argc, aims_const_hack(argv), "testQuantization2" );
    app.addOption(r, "-i", "input mesh" );
    app.initialize();

    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  std::size_t n = size(getVertices(mesh));
    
  std::cout << "Get circular neighbor..." << std::flush;
  shared_ptr<std::vector<std::vector<std::size_t> > > neighc = circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  std::cout << "OK" << std::endl;
  
  std::cout << "Get curvature..." << std::flush;
  til::Mesh_curvature<til::MeshTraits<til::Mesh_N>::VertexCollection, std::vector<std::vector<std::size_t> >, float>
    mc(getVertices(mesh), *neighc);
  std::vector<float> curv(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    mc.process(i);
    curv[i] = min(3.0f, mc.meanCurvature());
    curv[i] = max(-3.0f, curv[i]);
  }
  std::cout << "OK" << std::endl;

  std::cout << "Sort curvature" << std::endl;
  std::vector<std::pair<std::size_t, float> > indexed_curv(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    indexed_curv[i] = std::make_pair(i, std::abs(curv[i]));
  }
  std::sort(indexed_curv.begin(), indexed_curv.end(), til::Lesser_Pair2<std::size_t, float>());
  std::cout << "OK" << std::endl;
   
  std::cout << "Vector quantization..." << std::endl;
  const unsigned char UNPROCESSED = 0;
  const unsigned char KEEP = 1;
  const unsigned char REMOVE = 2;
  std::vector<unsigned char> label(n, UNPROCESSED);
  std::size_t count = 0;
  {
    // loop from lowest to highest curvature
    for (std::size_t i0 = 0; i0 < n; ++i0)
    {
      // get index of point
      std::size_t i = indexed_curv[i0].first;      
      // skip point if it has already been labeled
      if (label[i] != UNPROCESSED) continue;
      // otherwise, label point as removable...
      label[i] = REMOVE;
      ++count;
      // .. and label its neighbors as unremovable      
      for (std::size_t j = 0; j < (*neighc)[i].size(); ++j)
      {
        assert(label[(*neighc)[i][j]] != REMOVE);
        label[(*neighc)[i][j]] = KEEP;
      }
    }
  }
  std::cout << "OK" << std::endl;

  std::cout << "Before: " << n << std::endl;
  std::cout << "Removed: " << count << std::endl;
  std::cout << "Now: " << n - count << std::endl;
  count = 0;
  for (std::size_t i = 0; i < n; ++i)
  {
    if (label[i] == KEEP) ++count;
  }
  std::cout << "count: " << count << std::endl;
  
  std::cout << "Write remaining points..." << std::flush;
  {
    Texture1d tex(1,n);
    for (std::size_t i = 0; i < n; ++i)
    {
      tex.item(i) = label[i];
    }
    til::aimswrite(tex, "/home/cathier/tmp/removed");
  }
  std::cout << "OK" << std::endl;  
  

  std::cout << "Converting into graph..." << std::flush;
  typedef til::MeshVertexNodeX<std::size_t> VertexNode;
  typedef til::MeshFaceNodeX<VertexNode> FaceNode;
  std::list<VertexNode> graph_vertices;
  std::list<FaceNode> graph_faces;
  til::list2graph_mesh_conversion(getVertices(mesh), getFaceIndices(mesh), *neighc, graph_vertices, graph_faces);
  std::cout << "OK" << std::endl;

  /*
  std::cout << "Inverting face indices..." << std::flush;
  std::vector<std::list<std::size_t> > faceIndices = til::invertIndices(getFaceIndices(mesh));
  std::cout << "OK" << std::endl;
  
  til::List2GraphMeshConvertor<VertexNode, FaceNode> l2g;
  l2g(getVertices(mesh), getFaceIndices(mesh), faceIndices, *neighc);
  std::list<VertexNode> & graph_vertices = l2g.graph_vertices();
  std::list<FaceNode> & graph_faces = l2g.graph_faces();
  */
  
  std::cout << "Adding attribute..." << std::flush;
  std::list<VertexNode>::iterator it = graph_vertices.begin();
  {
    std::size_t gvsize = graph_vertices.size();
    for (std::size_t i = 0; i < gvsize; ++i, ++it)
    {
      it->attribute() = i;
    }
  }
  std::cout << "OK" << std::endl;;
  
    
  std::cout << "Removing vertices..." << std::flush;
  count = 0;
  bool complete;
  {
    int iter = 0;
    bool done;
    til::Vertex_remover<VertexNode, FaceNode> vertexRemover(graph_vertices, graph_faces);
    do
    {
      complete = true;
      done = true;
      std::list<VertexNode>::iterator i = graph_vertices.begin();
      std::list<FaceNode>::iterator itmp;
      while (i != graph_vertices.end())
      {
        if (til::size(i->neighbors()) != til::size(i->faces()))
        {
          std::cout << "wO!" << til::size(i->neighbors()) << "-" << til::size(i->faces()) << "-";
        }
        //std::cout << i->attribute() << ":" << int(removed[i->attribute()]) << " ";
        if (label[i->attribute()] == REMOVE)
        {
          //std::cout << "." << std::flush;
          // Here, the 'done' boolean should be here.
          // Note that it means that to be robust, I should add a check that the number of points
          // has strictly decreased during two iterations.
          complete = false;
          if (vertexRemover(i))
          {
            done = false;
            ++count;
            continue;
          }
        }
        ++i;
      }
      ++iter;
      std::cout << "pass " << iter << " : " << count << std::endl;
    } while (!done);
  }
  std::cout << "OK" << std::endl;
  std::cout << "Removed " << count << " points" << std::endl;
  if (!complete)
  {
    std::cout << "Warning: not all desired points could be removed" << std::endl;
  }
  
  std::cout << "Converting into mesh..." << std::flush;
  til::Graph2ListMeshConvertor<til::Mesh1> g2l;
  g2l(graph_vertices, graph_faces);
  til::Mesh1 & mesh1 = g2l.mesh();
  std::cout << "OK" << std::endl;

  {
    til::Mesh_N mesh1n = addNeighborsToMesh(mesh1);
    meshStats(mesh1n);
  }

  std::cout << "Writing mesh..." << std::flush;
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh1);
    til::aimswrite(s, "/home/cathier/tmp/biconverted");
  }
  std::cout << "OK" << std::endl;
}

class Pol : public std::unary_function<boost::array<double,1>, double>
{
  double operator()(double x)
  {
    double p = x * ( (x + 9/8)*x + 201/1280) + 11/640;
    return p*p;
  }
};


void writeEmptyFloatImage()
{
  AimsData<float> im(256,256,256);
  for (unsigned int i = 0; i < 256*256*256; ++i)
  {
    im[i] = 0.0f;
  }
  til::aimswrite(im, "/home/cathier/tmp/empty");
}


/*
void dwtImage(int argc, char * argv[])
{
  Writer<AimsData<float> > w;
//  Writer<AimsData<float> > w2;
  AimsData<float> im;
  double alpha = 0.001;
  {
    Reader<AimsData<float> > r;
    AimsApplication app( argc, aims_const_hack(argv), "dwtImage" );
    app.addOption( r, "-i", "input image" );
    app.addOption( w, "-o", "output image" );
    //app.addOption( w2, "-o2", "output image2" );
    app.addOption( alpha, "-alpha", "alpha" );
    app.initialize();
    std::cout << "Loading image..." << std::flush;
    r.read(im);
    std::cout << "OK" << std::endl;
  }
  til::numeric_array<std::size_t, 3> dim;
  dim[0] = im.dimX();
  dim[1] = im.dimY();
  dim[2] = im.dimZ();
  
  std::cout << "DWT..." << std::flush;
  til::multi_dwtND_cubic(&im[0], dim);
  std::cout << "OK" << std::endl;
  
  std::cout << "Abs..." << std::flush;
  AimsData<float> im2 = im.clone();
  std::cout << std::endl;
  for (int i = 0; i < im2.dimX()*im2.dimY()*im2.dimZ(); ++i)
  {
    im2[i] = std::abs(im2[i]);
  }
  std::cout << "OK" << std::endl;

  std::cout << "DWT power..." << std::flush;
  til::multi_dwtND_power(&im2[0], dim);
  std::cout << "OK" << std::endl;

  std::cout << "Computing power sum..." << std::flush;
  double sum = std::accumulate(&im2[0], &im2[im2.dimX()*im2.dimY()*im2.dimZ()], 0.0, std::plus<double>());
  std::cout << "OK" << std::endl;

  std::cout << "Sorting initialization..." << std::flush;  
  std::vector<std::pair<std::size_t, float> > im3(im2.dimX()*im2.dimY()*im2.dimZ());
  for (std::size_t i = 0; i < im3.size(); ++i)
  {
    im3[i].first = i;
    im3[i].second = im2[i];
  }
  std::cout << "OK" << std::endl;
  
  std::cout << "Sorting..." << std::flush;  
  std::sort(im3.begin(), im3.end(), til::Lesser_Pair2<std::size_t, float>());
  std::cout << "OK" << std::endl;
  
  std::cout << "Flattening spectrum..." << std::flush;
  double psum = 0.0;
  int count = 0;
  for (std::size_t i = 0; i < im3.size(); ++i)
  {
    psum += im3[i].second;
    if (psum < alpha * sum)
    {
      im[im3[i].first] = 0.0f;
      ++count;
    }
    else break;
  }
  std::cout << "OK" << std::endl;
  std::cout << "Removed " << count << " points out of " << im3.size() << std::endl;
  
  std::cout << "DWT shuffle..." << std::flush;
  til::multi_dwtND_shuffle(&im[0], dim);
  std::cout << "OK" << std::endl;
  
  std::cout << "Writing image..." << std::flush;
  w.write(im);
  std::cout << "OK" << std::endl;

  / *  
  std::cout << "DWT unshuffle..." << std::flush;
  til::multi_dwtND_unshuffle(&im[0], dim);
  std::cout << "OK" << std::endl;

  std::cout << "IDWT..." << std::flush;
  til::multi_idwtND_cubic(&im[0], dim);
  std::cout << "OK" << std::endl;

  std::cout << "Writing image..." << std::flush;
  w2.write(im);
  std::cout << "OK" << std::endl;
  * /  
}
*/

void testdwt2D()
{
  til::numeric_array<std::size_t, 2> N;
  N[0] = 11;
  N[1] = 13;
  til::numeric_array<std::size_t, 2> step;
  step[0] = 1;
  step[1] = 1;
  std::vector<double> im(N[0]*N[1]);
  for (std::size_t i = 0; i < N[0]*N[1]; ++i)
  {
    im[i] = i;
  }
  //til::dwt2D_cubic(&im.front(), N, 1, step);
  //til::idwt2D_cubic(&im.front(), N, 1, step);
  til::multi_dwtND_cubic(&im.front(), N);
  std::copy(im.begin(), im.end(), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;

  til::multi_dwtND_shuffle(&im.front(), N);
  std::copy(im.begin(), im.end(), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;

  til::multi_dwtND_unshuffle(&im.front(), N);
  std::copy(im.begin(), im.end(), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;

  til::multi_idwtND_cubic(&im.front(), N);
  std::copy(im.begin(), im.end(), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;  
}

void testdwt()
{
  for (int N = 4; N < 47; ++N)
  {
    std::vector<double> v(N);
    for (int i = 0; i < N; ++i) v[i] = i;
    std::copy(v.begin(), v.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

    til::multi_dwt_cubic(&v.front(), N);
    std::copy(v.begin(), v.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

    til::multi_dwt_shuffle(&v.front(), N);
    std::copy(v.begin(), v.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

    til::multi_dwt_unshuffle(&v.front(), N);
    std::copy(v.begin(), v.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

    til::multi_idwt_cubic(&v.front(), N);
    std::copy(v.begin(), v.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
  }
  /*
  int N = 13;
  std::vector<double> v(N);
  for (int i = 0; i < N; ++i) v[i] = i;
  std::copy(v.begin(), v.end(), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;
  
  til::dwt_cubic(&v.front(), N, 1);
  std::copy(v.begin(), v.end(), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;
  
  til::idwt_cubic(&v.front(), N, 1);
  std::copy(v.begin(), v.end(), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;
  
  til::multi_dwt_cubic(&v.front(), N, 1);
  std::copy(v.begin(), v.end(), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;
  
  til::multi_idwt_cubic(&v.front(), N, 1);
  std::copy(v.begin(), v.end(), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;
  */
}

void testUnsignedDiff()
{
  unsigned int i = 10;
  unsigned int j = 11;
  int k = i - j;
  std::cout << "k=" << k << std::endl;
}


void testQuant2(int argc, char* argv[])
{
  typedef til::Mesh_N MyMesh;
  MyMesh mesh;
  std::size_t N = 10;
  std::cout << "Reading mesh..." << std::flush;
  {
    Reader<AimsSurfaceTriangle> r;
    AimsApplication app( argc, aims_const_hack(argv), "testQuant2" );
    app.addOption(r, "-i", "input mesh" );
    app.addOption(N, "-n", "number of vertices" );
    app.initialize();

    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  std::cout << "OK" << std::endl;
  if (N >= getVertices(mesh).size())
  {
    std::cerr << "Number of final vertices cannot exceed the initial number" << std::endl;
    return;
  }
  if (N >= getVertices(mesh).size() / 2)
  {
    std::cout << "Warning: the number of vertices is relatively high -- this method is not appropriate for low decimation rates" << std::endl;
  }
  
  std::cout << "Maxdist quantization..." << std::flush;
  std::vector<std::size_t> vIndex(N);
  vIndex[0] = 0;
  //til::Triangle_mesh_geodesic_map<til::Mesh_N, double > geomap(mesh);
  typedef std::vector<std::vector<std::size_t> > CNeighborhoods;
  shared_ptr<CNeighborhoods> pneighc = til::circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  til::Triangle_mesh_geodesic_map<til::Mesh_N::VertexCollection, CNeighborhoods, double >
    geomap(getVertices(mesh), *pneighc);
  for (std::size_t i = 1; i < N; ++i)
  {
    std::vector<std::size_t> startPoints(i);
    std::copy(vIndex.begin(), vIndex.begin() + i, startPoints.begin());
    std::vector<double> dist(i, 0.0);
    geomap.init(startPoints, dist);
    geomap.process();
    shared_ptr<std::vector<double> > res = geomap.distanceMap();
    double dmax = *std::max_element(res->begin(), res->end());
    std::size_t iAdd = std::distance(res->begin(), std::max_element(res->begin(), res->end()));
    std::cout << "Adding point " << iAdd << " at a distance " << dmax << std::endl;
    vIndex[i] = iAdd;
  }
  std::cout << "OK" << std::endl;

  {
    std::vector<unsigned char> tmp(getVertices(mesh).size(), 0);
    for (std::size_t i = 0; i < N; ++i)
      tmp[vIndex[i]] = 1;
    Texture1d t; //(1, til::size(res));
    t.reserve(getVertices(mesh).size());
    for (std::size_t i = 0; i < getVertices(mesh).size(); ++i)
    {
      //std::cout << tmp[i] << " ";
      t.push_back(tmp[i]);
    }
    til::aimswrite(t, "quant2");
  }

  std::cout << "Computing Voronoi cells..." << std::flush;  
  std::vector<unsigned int> vcells;
  {
    std::vector<double> dist(N,0.0);
    //til::Voronoi_triangle_mesh_geodesic_map<til::Mesh_N, double> geomap(mesh);
    typedef std::vector<std::vector<std::size_t> > CNeighborhoods;
    shared_ptr<CNeighborhoods> pneighc = til::circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
    til::Voronoi_triangle_mesh_geodesic_map<til::Mesh_N::VertexCollection, CNeighborhoods, double>
      geomap(getVertices(mesh), *pneighc);
    geomap.init(vIndex, dist);
    geomap.process();
    vcells = *geomap.clusterLabels();
  }
  std::cout << "OK" << std::endl;
  
  std::cout << "Writing cells..." << std::flush;
  {
    Texture1d t;
    t.reserve(getVertices(mesh).size());
    for (std::size_t i = 0; i < getVertices(mesh).size(); ++i)
    {
      t.push_back(vcells[i]);
    }
    Writer<Texture1d> w("clust");
    w.write(t);
  }
  std::cout << "OK" << std::endl;

}

/*
void addBorders( int argc, char* argv[] )
{
  typedef AimsData<short> Image; 
  Reader<Image> r;
  Writer<Image> w;
  AimsApplication app( argc, aims_const_hack(argv), "clearBorders" );
  app.addOption(r, "-i", "input image" );
  app.addOption(margin, "-margin", "margin");
  app.addOption(w, "-o", "output image");
  app.initialize();

  Image im;
  r.read(im);

  Image imout(im.dimX()+2*margin, im.dimY()+2*margin, im.dimZ()+2*margin);

  for (int k = 0; k < im.dimZ(); ++k)
  for (int j = 0; j < im.dimY(); ++j)
  for (int i = 0; i < im.dimX(); ++i)
  {
    imout(i+margin, j+margin, k+margin) = im(i, j, k);
  }
}
*/


void clearBorders( int argc, char* argv[] )
{
  typedef AimsData<short> Image; 
  Reader<Image> r;
  Writer<Image> w;
  AimsApplication app( argc, aims_const_hack(argv), "clearBorders" );
  app.addOption(r, "-i", "input image" );
  app.addOption(w, "-o", "output image");
  app.initialize();
  
  Image im;
  r.read(im);
  
  for (int j = 0; j < im.dimY(); ++j)
  for (int i = 0; i < im.dimX(); ++i)
  {
    im(i, j, 0) = 0;
    im(i, j, im.dimZ()-1) = 0;
  }

  for (int k = 0; k < im.dimZ(); ++k)
  for (int i = 0; i < im.dimX(); ++i)
  {
    im(i, 0, k) = 0;
    im(i, im.dimY()-1, k) = 0;
  }

  for (int k = 0; k < im.dimZ(); ++k)
  for (int j = 0; j < im.dimY(); ++j)
  {
    im(0, j, k) = 0;
    im(im.dimX()-1, j, k) = 0;
  }

  w.write(im);
}

/*
/// To deform a data set with RBFs.
template < typename TDataIterator, typename TRBFCenterIterator, typename TRBFSupportIterator, typename TRBF, typename TCoeffIterator, typename TOutputIterator >
class VectQuantRBFTransfo
{
public: // constructors

  VectQuantRBFTransfo(
  TDataIterator dataBegin, TDataIterator dataEnd,
  TRBFCenterIterator centerBegin, TRBFCenterIterator centerEnd, 
  TRBFSupportIterator supportBegin, TRBFSupportIterator supportEnd,
  TRBF rbf)
    : m_dataBegin(dataBegin)
    , m_dataEnd(dataEnd)
    , m_centerBegin(centerBegin)
    , m_centerEnd(centerEnd)
    , m_supportBegin(supportBegin)
    , m_supportEnd(supportEnd)
    , m_rbf(rbf)
  {}

public: //

  void operator()(TCoeffIterator coeffBegin, TOutputIterator outBegin)
  {
    for (TDataIterator iData = m_dataBegin; iData != m_dataEnd; ++iData)
    {
      for (TCenterIterator iCenter = m_centerBegin, TCoeffIterator iCoeffs = coeffBegin; iCenter != centerEnd; ++iCenter, ++iCoeffs)
      {
        
      }
    }
  }

private: // data, input

  // data
  TDataIterator m_dataBegin;
  TDataIterator m_dataEnd;
  
  // RBF centers
  TRBFCenterIterator m_centerBegin;
  TRBFCenterIterator m_centerEnd;
  
  // RBF support
  TRBFSupportIterator m_supportBegin;
  TRBFSupportIterator m_supportEnd;
  
  // RBF function
  TRBF m_rbf;
};
*/


template < class TVertexCollection, class TNeighborhoodCollection, class TSupport, class TFinder, class TQuant, class TKernel >
class RBFReg : public std::unary_function< std::vector<float>, double >
{
public: // typedefs

  //typedef typename til::MeshTraits<TMesh>::VertexCollection   VertexCollection;
  //typedef typename til::MeshTraits<TMesh>::Vertex             Vertex;
  typedef TVertexCollection                       VertexCollection;
  typedef typename TVertexCollection::value_type  Vertex;

public: // constructors

  RBFReg
  (
    const TVertexCollection & vertices,
    const TVertexCollection & vertices2,
    const TNeighborhoodCollection & neighc2,
    const TFinder & finder,
    const TSupport & support,
    const TQuant & quant,
    TKernel kernel
  )
//    : m_mesh(mesh)
//    , m_mesh2(mesh2)
    : m_vertices(vertices)
    , m_vertices2(vertices2)
    , m_neighc2(neighc2)
    , m_finder(finder)
    , m_support(support)
    , m_quant(quant)
    , m_kernel(kernel)
  {}

public: // operators

  double operator()(std::vector<float> const & params)
  {
    double res = 0.0;    
    std::size_t n = m_vertices.size();
    for (std::size_t i = 0; i < n; ++i)
    //for (typename VertexCollection::const_iterator iVertex = begin; iVertex != end; ++iVertex)
    {
      m_tmp = m_vertices[i];
      std::size_t supportSize = m_support[i].size();
      for (std::size_t j = 0; j < supportSize; ++j)
      {
        std::size_t sij = m_support[i][j];
        double coeff = m_kernel(*m_quant[sij], m_vertices[i]);
        //double dexp = std::exp( - til::dist2(*m_quant[sij], m_vertices[i]) / (2 * sigma * sigma) );
        sij *= 3;
        m_tmp[0] += coeff * params[sij++];
        m_tmp[1] += coeff * params[sij++];
        m_tmp[2] += coeff * params[sij];
      }
      std::size_t cpi = m_finder(m_tmp);
      res += til::dist2_surface<double>(m_tmp, cpi, m_vertices2, m_neighc2[cpi]);
      //res += til::dist2<double>(m_tmp, m_vertices2[]);

    }

    //res /= std::sqrt(max(0.0f,til::det(a.transfo().getMatrix()))) + 128 * std::numeric_limits<double>::epsilon();

    std::cout << "Func called : " << res << std::endl;
    
    return res;
  }

private: // data, input

//  const TMesh & m_mesh;
//  const TMesh & m_mesh2;
  const TVertexCollection & m_vertices;
  const TVertexCollection & m_vertices2;
  const TNeighborhoodCollection & m_neighc2;
  TFinder m_finder;
  const TSupport & m_support;
  const TQuant & m_quant;
  TKernel m_kernel;

private: // data, internals
  Vertex m_tmp;  
};

template < class TVertexCollection, class TNeighborhoodCollection, class TFinder >
class RBFReg2 : public std::unary_function< std::vector<float>, double >
{
public: // typedefs

  //typedef typename til::MeshTraits<TMesh>::VertexCollection   VertexCollection;
  //typedef typename til::MeshTraits<TMesh>::Vertex             Vertex;
  typedef TVertexCollection                       VertexCollection;
  typedef typename TVertexCollection::value_type  Vertex;
  typedef std::vector<std::vector<std::pair<std::size_t, double> > >       Kernel;

public: // constructors

  RBFReg2
  (
    const TVertexCollection & vertices,
    const TVertexCollection & vertices2,
    const TNeighborhoodCollection & neighc2,
    const TFinder & finder,
    const Kernel & kernel
  )
    : m_vertices(vertices)
    , m_vertices2(vertices2)
    , m_neighc2(neighc2)
    , m_finder(finder)
    , m_kernel(kernel)
    , m_movedVertices(vertices.size())
  {}

/*
public: // set & get

  void setReferenceBrain(const TVertexCollection & vertices, const TFinder & finder)
  {
    m_vertices = vertices
  }
  */

public: // operators

  double operator()(std::vector<float> const & params)
  {
    // check that the number of parameters matches the number of deformation functions.
    assert(params.size() == 3 * m_kernel.size());
    // apply the transformation
    this->applyTransform(params);    
    return this->computeEnergy();
  }

private: // function

  double computeEnergy()
  {
    double res = 0.0;
    // loop through all moved vertices
    for (typename std::vector<Vertex>::iterator iMovedVertex = m_movedVertices.begin(); iMovedVertex != m_movedVertices.end(); ++iMovedVertex)
    {
      std::size_t i = m_finder(*iMovedVertex);
      res += til::dist2_surface<double>(*iMovedVertex, i, m_vertices2, m_neighc2[i]);
      //res += til::dist2<double>(m_tmp, m_vertices2[]);
    }
    //std::cout << "Func called : " << res << std::endl;
    return res;
  }

  /// Apply the transform with given coefficients and set result in m_movedVertices.
  void applyTransform(std::vector<float> const & params)
  {
    // initialization
    //std::fill(m_movedVertices.begin(), m_movedVertices.end(), til::Point<float,3>(0,0,0));
    std::copy(m_vertices.begin(), m_vertices.end(), m_movedVertices.begin());

    // loop through all kernels
    std::size_t iParam = 0;
    for (typename Kernel::iterator iKernel = m_kernel.begin(); iKernel != m_kernel.end(); ++iKernel, iParam += 3)
    {
      // for each kernel, loop on its definition domain
      for (typename Kernel::value_type::iterator iiKernel = iKernel->begin(); iiKernel != iKernel->end(); ++iiKernel)
      {
        for (std::size_t j = 0; j < 3; ++j)
        {
          m_movedVertices[iiKernel->first][j] += iiKernel->second * params[iParam + j];
        }
      }
    }
  }

private: // data, input

  /*
  TVertexIterator m_iVertexBegin;
  TVertexIterator m_iVertexEnd;

  TVertexIterator m_iVertex2Begin;
  TVertexIterator m_iVertex2End;
  */
  
  const TVertexCollection & m_vertices;
  const TVertexCollection & m_vertices2;
  const TNeighborhoodCollection & m_neighc2;
  TFinder m_finder;
  Kernel m_kernel;

private: // data, internals

  //Vertex m_tmp;
  std::vector<Vertex> m_movedVertices;
};


/*
template < typename TKernel, typename TCenterIterator, typename TCoeffIterator, typename TRes = typename TKernel::result_type >
class KernelMixture
  : std::unary_function<typename TKernel::second_argument_type, TRes>
{
public: // typedef
  
  // NB: by convention, the first argument of the kernel is used for centers, and the second for mixture arguments.
  
  typedef std::binary_function<typename TKernel::first_argument_type, typename TKernel::second_argument_type, TRes> Base;
  //typedef typename Base::first_argument_type first_argument_type;
  //typedef typename Base::second_argument_type second_argument_type;
  typedef typename Base::second_argument_type argument_type;
  typedef typename Base::result_type result_type;

public: // constructors

  KernelMixture(TKernel kernel)
    : m_kernel(kernel)
  {}

public: // set & get

  void setCenters(TCenterIterator begin, TCenterIterator end) { m_begin = begin; m_end = end; }
  void setCoeffs(TCoeffIterator cbegin, TCoeffIterator cend) { m_cbegin = cbegin; m_cend = cend; }

public: // functions

  result_type operator()(argument_type x)
  {
    result_type res = result_type();
    for (TCenterIterator iCenter = m_begin, TCoeffIterator iCoeff = m_cbegin; iCenter != m_end; ++iCenter, ++iCoeff)
    {
      res += *iCoeff * m_kernel(*iCenter, x);
    }
  }
  
private: // data, input

  TCenterIterator m_begin;
  TCenterIterator m_end;
  TCoeffIterator m_cbegin;
  TCoeffIterator m_cend;  
  TKernel m_kernel;
};

template < typename TKernel, typename TCenterIterator, typename TCoeffIterator, typename TRes = typename TKernel::result_type >
class NormalizedKernelMixture
  : std::unary_function<typename TKernel::second_argument_type, TRes>
{
public: // typedef
  typedef KernelMixture<TKernel, TCenterIterator, TCoeffIterator, TRes> KernelMix;
public: // constructors

  NormalizedKernelMixture(TKernel kernel) : m_kernelmix(kernel) {}

public: // operator
  
  TRes operator()(typename KernelMix::argument_type)
  {
    return 
  }

private: // data
  KernelMixt m_kernelmix;
};
*/

void testGonzCluster( int argc, char* argv[] )
{
  typedef til::Mesh_N MyMesh;
  MyMesh mesh;
  float sigma;
  
  std::cout << "Reading mesh..." << std::flush;
  {
    Reader<AimsSurfaceTriangle> r;
    AimsApplication app( argc, aims_const_hack(argv), "testGonzCluster" );
    app.addOption(r, "-i", "input mesh" );
    //app.addOption(maxDist, "-maxDist", "maxdist" );
    app.addOption(sigma, "-sigma", "maxdist" );
    app.initialize();

    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  std::size_t n = getVertices(mesh).size();
  std::cout << "OK" << std::endl;

  std::cout << "Duplicating mesh..." << std::flush;
  til::Mesh_N mesh2;
  getVertices(mesh2) = getVertices(mesh);
  getFaceIndices(mesh2) = getFaceIndices(mesh);
  std::cout << "OK" << std::endl;
    
  std::cout << "Computing neighborhoods..." << std::flush;
  shared_ptr<std::vector<std::vector<std::size_t> > > neighc = circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  std::cout << "OK" << std::endl;
        
  std::cout << "Computing Gonzalez clustering..." << std::flush;
  typedef std::vector<til::numeric_array<float,3> >::const_iterator iterator;
  //til::GonzalezClustering<std::vector<til::Point<float,3> >, float > gc(getVertices(mesh));
  //std::cout << til::square(sigma / 0.6) << std::endl;
  //gc.clusterize_maxDiam(til::square(sigma / 0.6));
  //shared_ptr<std::vector<iterator> > quant = gc.quantization();
  shared_ptr<std::vector<iterator> > quant = til::gonzalez_clustering(getVertices(mesh), til::square(sigma / 0.6));
  std::cout << "OK" << std::endl;
  
    
  {
    Texture1d t(1, n);
    for (std::size_t i = 0; i < quant->size(); ++i)
    {
      t.item(std::distance(iterator(getVertices(mesh).begin()),(*quant)[i])) = 1;
    }
    Writer<Texture1d> w("/home/cathier/tmp/quantgonz");
    w.write(t);
  }
  exit(0);
  /*
  shared_ptr<std::vector<std::size_t> > labels = gc.labels();
  {
    Texture1d t(1, n);
    for (std::size_t i = 0; i < n; ++i)
    {
      t.item(i) = (*labels)[i];
    }
    Writer<Texture1d> w("/home/cathier/tmp/gonzclust");
    w.write(t);
  }
  */
  
  std::cout << "Computing RBF support..." << std::flush;
  std::vector<std::vector<std::size_t> > support(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t j = 0; j < quant->size(); ++j)
    {
      if (til::dist2(getVertices(mesh)[i], *(*quant)[j], til::prec<float>()) < til::square(4*sigma))
      {
        support[i].push_back(j);
      }
    }
  }
  std::cout << "OK" << std::endl;

  {
    std::vector<std::size_t>::iterator iS;
    Texture1d t(1, n);
    for (std::size_t i = 0; i < n; ++i)
    {
      if ((iS = std::find(support[i].begin(), support[i].end(), std::size_t(5))) != support[i].end())
      //if ((std::find(support[i].begin(), support[i].end(), 5)) != support[i].end())
      {
        t.item(i) = std::exp( - til::dist2(*(*quant)[*iS], getVertices(mesh)[i], til::prec<float>()) / (2 * sigma * sigma) );
        //t.item(i) = 1;
      }
    }
    Writer<Texture1d> w("/home/cathier/tmp/support0");
    w.write(t);
  }
  
  std::vector<til::numeric_array<float, 3> > coeffs(quant->size());
  
  for (std::size_t i = 0; i < coeffs.size(); ++i)
  {
    for (std::size_t j = 0; j < 3; ++j)
    {
      coeffs[i][j] = 0;
    }
  }
  coeffs[0][2] = -10.0f;
  

  /*
  std::cout << "Taking random coefficients..." << std::flush;
  til::UniformRandomDistribution<float> urand(-10.0f, 10.0f);
  for (std::size_t i = 0; i < coeffs.size(); ++i)
  {
    for (std::size_t j = 0; j < 3; ++j)
    {
      coeffs[i][j] = urand.draw();
    }
  }
  std::cout << "OK" << std::endl;
  */
  
  std::cout << "Applying transformation..." << std::flush;
  std::vector<til::numeric_array<float,3> > newVertices(n);
  std::copy(getVertices(mesh).begin(), getVertices(mesh).end(), newVertices.begin());
  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t j = 0; j < support[i].size(); ++j)
    {
      newVertices[i] += std::exp( - til::dist2(*(*quant)[support[i][j]], getVertices(mesh)[i], til::prec<float>()) / (2 * sigma * sigma) )
        * coeffs[support[i][j]];
    }
  }
  std::copy(newVertices.begin(), newVertices.end(), getVertices(mesh).begin());
  std::cout << "OK" << std::endl;
  
  std::cout << "Writing mesh..." << std::flush;
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh);
    Writer<AimsSurfaceTriangle> w("/home/cathier/tmp/deformed");
    w.write(s);
  }
  std::cout << "OK" << std::endl;
  
  std::cout << "Generating KDTree..." << std::flush;
  //typedef til::KDTree<til::MeshTraits<til::Mesh_N>::Vertex*, til::MeshTraits<til::Mesh_N>::VertexCollection> MyKDTree;
  typedef til::KDTree<std::size_t, til::MeshTraits<til::Mesh_N>::VertexCollection> MyKDTree;
  MyKDTree kdtree(getVertices(mesh));
  makeKDTree(getVertices(mesh), kdtree);
  std::cout << "OK" << std::endl;
  
  std::cout << "Starting minimization (" << 3*quant->size() << " degrees of freedom)" << std::endl;
  typedef til::Find_closest< double, MyKDTree > Finder;
  Finder finder(kdtree);
  typedef RBFReg <
    std::vector<til::numeric_array<float,3> >,
    std::vector<std::vector<std::size_t> >,
    std::vector<std::vector<std::size_t> >,
    //std::vector<std::vector<til::Point<float,3>*> >,
    Finder,
    std::vector<iterator>,
    til::math::IsotropicGaussianKernel<til::numeric_array<float,3>, float>
  > Energy;
  til::math::IsotropicGaussianKernel<til::numeric_array<float,3>, float> gaussianK(sigma);
  Energy energy(getVertices(mesh2), getVertices(mesh), *neighc, finder, support, *quant, gaussianK);
  til::Powell<Energy> minalgo(energy);
  // initial scaling estimates of the parameters
  minalgo.initStd() = std::vector<float>(3*quant->size(), 2.0f);
  std::vector<float> params(3*quant->size(), 0.0f);
  params = minalgo(params);
  std::cout << "Minimization finished" << std::endl;

  std::cout << "Applying transformation..." << std::flush;
  {
    std::size_t n = getVertices(mesh2).size();
    std::vector<til::numeric_array<float,3> > newVertices(n);
    for (std::size_t i = 0; i < n; ++i)
    {
      newVertices[i] = getVertices(mesh2)[i];
      for (std::size_t j = 0; j < support[i].size(); ++j)
      {
        std::size_t sij = support[i][j];
        double dexp = gaussianK(*(*quant)[sij], getVertices(mesh2)[i]);
        //double dexp = std::exp( - til::dist2(*m_quant[sij], m_vertices[i]) / (2 * sigma * sigma) );
        sij *= 3;
        newVertices[i][0] += dexp * params[sij++];
        newVertices[i][1] += dexp * params[sij++];
        newVertices[i][2] += dexp * params[sij];
      }
    }
    std::copy(newVertices.begin(), newVertices.end(), getVertices(mesh2).begin());
  }
  std::cout << "OK" << std::endl;

  std::cout << "Writing result..." << std::flush;
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh2);
    Writer<AimsSurfaceTriangle> w("/home/cathier/tmp/regres");
    w.write(s);
  }
  std::cout << "OK" << std::endl;
}


void testTriangleIntersection( int argc, char* argv[] )
{
  double a0, a1, a2;
  double b0, b1, b2;
  double c0, c1, c2;
  double d0, d1, d2;
  double e0, e1, e2;
  double f0, f1, f2;
  {
    AimsApplication app( argc, aims_const_hack(argv), "testTriangleIntersection" );
    app.addOption(a0, "-a0", "" );
    app.addOption(a1, "-a1", "" );
    app.addOption(a2, "-a2", "" );
    app.addOption(b0, "-b0", "" );
    app.addOption(b1, "-b1", "" );
    app.addOption(b2, "-b2", "" );
    app.addOption(c0, "-c0", "" );
    app.addOption(c1, "-c1", "" );
    app.addOption(c2, "-c2", "" );
    app.addOption(d0, "-d0", "" );
    app.addOption(d1, "-d1", "" );
    app.addOption(d2, "-d2", "" );
    app.addOption(e0, "-e0", "" );
    app.addOption(e1, "-e1", "" );
    app.addOption(e2, "-e2", "" );
    app.addOption(f0, "-f0", "" );
    app.addOption(f1, "-f1", "" );
    app.addOption(f2, "-f2", "" );
    app.initialize();
  }
  
  til::geo::AreIntersecting::Triangles3D<til::numeric_array<double,3> > ai;
  til::numeric_array<double,3> A;
  A[0] = a0;
  A[1] = a1;
  A[2] = a2;
  til::numeric_array<double,3> B;
  B[0] = b0;
  B[1] = b1;
  B[2] = b2;
  til::numeric_array<double,3> C;
  C[0] = c0;
  C[1] = c1;
  C[2] = c2;
  til::numeric_array<double,3> D;
  D[0] = d0;
  D[1] = d1;
  D[2] = d2;
  til::numeric_array<double,3> E;
  E[0] = e0;
  E[1] = e1;
  E[2] = e2;
  til::numeric_array<double,3> F;
  F[0] = f0;
  F[1] = f1;
  F[2] = f2;
  
  std::cout << A << " " << B << " " << C << " " << D << " " << E << " " << F << std::endl;
  
  std::cout << ai(A,B,C,D,E,F) << std::endl;
  
}

void testSelfIntersection( int argc, char* argv[] )
{
  til::Mesh_N mesh;
  std::cout << "Reading mesh..." << std::flush;
  {
    Reader<AimsSurfaceTriangle> r;
    AimsApplication app( argc, aims_const_hack(argv), "testSelfIntersection" );
    app.addOption(r, "-i", "input mesh" );
    app.initialize();

    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  std::size_t n = getVertices(mesh).size();
  std::cout << "OK" << std::endl;
  
  til::geo::AreIntersecting::Triangles3D<til::numeric_array<double,3> > ai;
  //til::AreIntersecting<til::numeric_array<float,3> > ai;
  til::numeric_array<double,3> v1, v2, v3, w1, w2, w3;
  for (std::size_t i = 0; i < getFaceIndices(mesh).size(); ++i)
  {
    std::cout << "." << std::flush;
    for (std::size_t j = 0; j < getFaceIndices(mesh).size(); ++j)
    {
      {
        til::numeric_array<float,3> & tmp = getVertices(mesh)[getFaceIndices(mesh)[i][0]];
        std::copy(tmp.begin(), tmp.end(), v1.begin());
      }
      {
        til::numeric_array<float,3> & tmp = getVertices(mesh)[getFaceIndices(mesh)[i][1]];
        std::copy(tmp.begin(), tmp.end(), v2.begin());
      }
      {
        til::numeric_array<float,3> & tmp = getVertices(mesh)[getFaceIndices(mesh)[i][2]];
        std::copy(tmp.begin(), tmp.end(), v3.begin());
      }
      {
        til::numeric_array<float,3> & tmp = getVertices(mesh)[getFaceIndices(mesh)[j][0]];
        std::copy(tmp.begin(), tmp.end(), w1.begin());
      }
      {
        til::numeric_array<float,3> & tmp = getVertices(mesh)[getFaceIndices(mesh)[j][1]];
        std::copy(tmp.begin(), tmp.end(), w2.begin());
      }
      {
        til::numeric_array<float,3> & tmp = getVertices(mesh)[getFaceIndices(mesh)[j][2]];
        std::copy(tmp.begin(), tmp.end(), w3.begin());
      }
      
      /*
      if (ai(
      getVertices(mesh)[getFaceIndices(mesh)[i][0]],
      getVertices(mesh)[getFaceIndices(mesh)[i][1]],
      getVertices(mesh)[getFaceIndices(mesh)[i][2]],
      getVertices(mesh)[getFaceIndices(mesh)[j][0]],
      getVertices(mesh)[getFaceIndices(mesh)[j][1]],
      getVertices(mesh)[getFaceIndices(mesh)[j][2]]))
      */

      if (ai(v1,v2,v3,w1,w2,w3))
      {        
        if (
        v1 != w1 &&
        v1 != w2 &&
        v1 != w3 &&
        v2 != w1 &&
        v2 != w2 &&
        v2 != w3 &&
        v3 != w1 &&
        v3 != w2 &&
        v3 != w3)
        {
          std::cout << "*" << std::endl;
          
          Texture1d t(1, n);
  
          t.item(getFaceIndices(mesh)[i][0]) = 1;
          t.item(getFaceIndices(mesh)[i][1]) = 1;
          t.item(getFaceIndices(mesh)[i][2]) = 1;
          
          t.item(getFaceIndices(mesh)[j][0]) = 1;
          t.item(getFaceIndices(mesh)[j][1]) = 1;
          t.item(getFaceIndices(mesh)[j][2]) = 1;
  
          Writer<Texture1d> w("/home/cathier/tmp/inters");
          w.write(t);
          exit(0);

          /*        
          std::cout << "Intersecting triangles " << i << " and " << j << std::endl;
          std::cout << 
            getVertices(mesh)[getFaceIndices(mesh)[i][0]] << " " << 
            getVertices(mesh)[getFaceIndices(mesh)[i][1]] << " " << 
            getVertices(mesh)[getFaceIndices(mesh)[i][2]] << " " << 
            getVertices(mesh)[getFaceIndices(mesh)[j][0]] << " " << 
            getVertices(mesh)[getFaceIndices(mesh)[j][1]] << " " << 
            getVertices(mesh)[getFaceIndices(mesh)[j][2]] << std::endl;
          */
        }
      }
    }
  }
}

void testCurvatureSmoothing( int argc, char* argv[] )
{
  // NB: it seems that curvature smooting cannot work. When surface is smoothing, unfortunately, vertices
  // lying on top of high curvature ridges tend to collapse. Curvature estimation becomes unstable in these areas.
  til::Mesh_N mesh;
  int niter;
  float alpha;
  Writer<AimsSurfaceTriangle> w;
  std::cout << "Reading mesh..." << std::flush;
  {
    Reader<AimsSurfaceTriangle> r;
    AimsApplication app( argc, aims_const_hack(argv), "testCurvatureSmoothing" );
    app.addOption(r, "-i", "input mesh" );
    app.addOption(w, "-o", "output mesh" );
    app.addOption(niter, "-niter", "niter");
    app.addOption(alpha, "-alpha", "alpha");
    app.initialize();

    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  std::cout << "OK" << std::endl;
  std::size_t n = getVertices(mesh).size();

  std::cout << "Get circular neighbor..." << std::flush;
  shared_ptr<std::vector<std::vector<std::size_t> > > neighc = circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  std::cout << "OK" << std::endl;

  std::cout << "Surface evolution" << std::flush;
  typedef float prec_type;
  til::Mesh_curvature<til::MeshTraits<til::Mesh_N>::VertexCollection, std::vector<std::vector<std::size_t> >, prec_type>
    mc(getVertices(mesh), *neighc);
  std::vector<prec_type> curv(n);
  std::vector<til::numeric_array<float,3> > normals(n);
  std::vector<til::numeric_array<float,3> > vertices(n);
  for (int iter = 0; iter < niter; ++iter)
  {
    std::copy(getVertices(mesh).begin(), getVertices(mesh).end(), vertices.begin());
    std::cout << "." << std::flush;
    for (std::size_t i = 0; i < n; ++i)
    {
      mc.process(i);
      curv[i] = max(prec_type(-3), min(prec_type(3), mc.meanCurvature()));
      normals[i] = mc.normal();
    }
    prec_type curvmean = std::accumulate(curv.begin(), curv.end(), prec_type(0));
    curvmean /= curv.size();
    for (std::size_t i = 0; i < n; ++i)
    {
      vertices[i] -= alpha * normals[i] * (curv[i] - curvmean);
    }
    std::copy(vertices.begin(), vertices.end(), getVertices(mesh).begin());
  }
  std::cout << "OK" << std::endl;

  
  std::cout << "Writing mesh..." << std::flush;
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh);
    w.write(s);
  }   
  std::cout << "OK" << std::endl;
}


void testCurvatureSmoothing2( int argc, char* argv[] )
{
  til::Mesh_N mesh;
  int niter;
  float alpha;
  double beta;
  Writer<AimsSurfaceTriangle> w;
  std::cout << "Reading mesh..." << std::flush;
  {
    Reader<AimsSurfaceTriangle> r;
    AimsApplication app( argc, aims_const_hack(argv), "testCurvatureSmoothing2" );
    app.addOption(r, "-i", "input mesh" );
    app.addOption(w, "-o", "output mesh" );
    app.addOption(niter, "-niter", "niter");
    app.addOption(alpha, "-alpha", "alpha");
    app.addOption(beta, "-beta", "alpha");    
    app.initialize();

    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  std::cout << "OK" << std::endl;
  std::size_t n = getVertices(mesh).size();

  std::cout << "Get circular neighbor..." << std::flush;
  shared_ptr<std::vector<std::vector<std::size_t> > > neighc = circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  std::cout << "OK" << std::endl;


  typedef til::MeshVertexNodeX<std::size_t> VertexNode;
  typedef til::MeshFaceNodeX<VertexNode> FaceNode;
  std::cout << "Converting into graph..." << std::flush;
  std::list<VertexNode> graph_vertices;
  std::list<FaceNode> graph_faces;
  til::list2graph_mesh_conversion(getVertices(mesh), getFaceIndices(mesh), *neighc, graph_vertices, graph_faces);
  /*
  std::vector<std::list<std::size_t> > faceIndices = til::invertIndices(getFaceIndices(mesh));
  til::List2GraphMeshConvertor<VertexNode, FaceNode> l2g;
  l2g(getVertices(mesh), getFaceIndices(mesh), faceIndices, *neighc);
  std::list<VertexNode> & graph_vertices = l2g.graph_vertices();
  std::list<FaceNode> & graph_faces = l2g.graph_faces();
  */
  std::cout << "OK" << std::endl;
  
  std::cout << "Adding attribute..." << std::flush;
  {
    std::list<VertexNode>::iterator it = graph_vertices.begin();
    std::size_t gvsize = graph_vertices.size();
    for (std::size_t i = 0; i < gvsize; ++i, ++it)
    {
      it->attribute() = i;
    }
  }
  std::cout << "OK" << std::endl;;

  std::cout << "Surface evolution" << std::flush;
  typedef float prec_type;
  typedef til::xsr::graph::Position<VertexNode::VertexIndex> A;
  typedef til::xsr::graph::Neighbors<VertexNode::VertexIndex> B;
  A a; B b; til::MeshCurvature2<A,B,float> mc(a,b);
  //til::MeshCurvature2<A,B,float> mc((A()),(B()));
  //til::Mesh_curvature<til::MeshTraits<til::Mesh_N>::VertexCollection, std::vector<std::vector<std::size_t> >, prec_type> mc(getVertices(mesh), *neighc);
  std::vector<prec_type> curv(n);
  std::vector<til::numeric_array<float,3> > normals(n);
  std::vector<til::numeric_array<float,3> > vertices(n);
  til::Vertex_remover<VertexNode, FaceNode> vertexRemover(graph_vertices, graph_faces);
  for (int iter = 0; iter < niter; ++iter)
  {
    std::cout << "." << std::flush;
    std::size_t n = graph_vertices.size();
    {
      std::vector<til::numeric_array<float,3> >::iterator iVertex = vertices.begin();
      std::list<VertexNode>::iterator iGraphVertex = graph_vertices.begin();
      for (; iGraphVertex != graph_vertices.end(); ++iVertex, ++iGraphVertex)
      {
        *iVertex = iGraphVertex->pos();
      }
    }
    {
      std::list<VertexNode>::iterator i = graph_vertices.begin();
      std::size_t index = 0;
      for (; i != graph_vertices.end(); ++i, ++index)
      {
        mc.process(i);
        curv[index] = max(prec_type(-3), min(prec_type(3), mc.meanCurvature()));
        normals[index] = mc.normal();
      }
    }
    prec_type curvmean = std::accumulate(curv.begin(), curv.end(), prec_type(0));
    curvmean /= curv.size();
    for (std::size_t i = 0; i < n; ++i)
    {
      vertices[i] -= alpha * normals[i] * (curv[i] - curvmean);
    }
    {
      std::vector<til::numeric_array<float,3> >::iterator iVertex = vertices.begin();
      std::list<VertexNode>::iterator iGraphVertex = graph_vertices.begin();
      for (; iGraphVertex != graph_vertices.end(); ++iVertex, ++iGraphVertex)
      {
        iGraphVertex->pos() = *iVertex;
      }
    }
    ///*
    {
      std::list<VertexNode>::iterator i = graph_vertices.begin();
      std::size_t index = 0;
      int count = 0;
      bool notdone;
      do
      {
        notdone = false;
        for (; i != graph_vertices.end(); ++index)
        {
          double perimeter = 0.0;
          double area = 0.0;
          VertexNode::VertexIndexCollection::const_iterator j = i->neighbors().begin();
          til::const_cyclic_iterator<VertexNode::VertexIndexCollection> j2(i->neighbors(), i->neighbors().begin());
          ++j2;
          for (; j != i->neighbors().end(); ++j, ++j2)
          {
            perimeter += til::dist((*j)->pos(), (*j2)->pos(), til::prec<double>());
            area = 0.5 * til::norm(til::cross((*j)->pos(), (*j2)->pos(), til::prec<double>()));
          }
          if (4 * M_PI * area  < beta * (perimeter*perimeter) )
          {
            if (vertexRemover(i))
            {
              notdone = true;
              ++count;
              continue;
            }
          }
          ++i;
        }
      } while (notdone);
      if (count) std::cout << "(removed " << count << " points)";
    }
    //*/
  }
  std::cout << "OK" << std::endl;

  std::cout << "Converting into mesh..." << std::flush;
  til::Graph2ListMeshConvertor<til::Mesh1> g2l;
  g2l(graph_vertices, graph_faces);
  til::Mesh1 & mesh1 = g2l.mesh();
  std::cout << "OK" << std::endl;
  
  std::cout << "Writing mesh..." << std::flush;
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh1);
    w.write(s);
  }   
  std::cout << "OK" << std::endl;
}

bool absLess(double x, double y)
{
  return std::abs(x) < std::abs(y);
}

void testGaussPatchApprox( int argc, char* argv[] )
{
  til::Mesh_N mesh;
  Writer<AimsSurfaceTriangle> w;
  std::cout << "Reading mesh..." << std::flush;
  float sigma;
  float beta;
  {
    Reader<AimsSurfaceTriangle> r;
    AimsApplication app( argc, aims_const_hack(argv), "testGaussPatchApprox" );
    app.addOption(r, "-i", "input mesh" );
    app.addOption(w, "-o", "output mesh" );
    app.addOption(sigma, "-sigma", "sigma" );
    app.addOption(beta, "-beta", "beta" );
    app.initialize();

    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  std::cout << "OK" << std::endl;
  std::size_t n = getVertices(mesh).size();

  std::vector<til::numeric_array<float, 3> > vertices(n);
  std::copy(getVertices(mesh).begin(), getVertices(mesh).end(), vertices.begin());

  /*
  std::cout << "Computing surface pseudo-coordinates..." << std::flush;
  std::vector<double> cx(n);
  std::vector<double> cy(n);
  std::vector<double> cz(n);
  {
    til::Triangle_mesh_geodesic_map<til::Mesh_N, double > geomap(mesh);
    std::vector<double> dist(1, 0.0);
    shared_ptr<std::vector<double> > res;
    double mean;
    
    std::vector<std::size_t> startPoints(1, 0);
    geomap.init(startPoints, dist);
    geomap.process();
    res = geomap.distanceMap();
    mean = std::accumulate(res->begin(), res->end(), 0.0);
    mean /= res->size();
    std::transform(res->begin(), res->end(), res->begin(), std::bind2nd(std::minus<double>(), mean));
    std::copy(res->begin(), res->end(), cx.begin());
    
    startPoints[0] = std::distance(res->begin(), std::min_element(res->begin(), res->end(), absLess));
    geomap.init(startPoints, dist);
    geomap.process();
    res = geomap.distanceMap();
    mean = std::accumulate(res->begin(), res->end(), 0.0);
    mean /= res->size();
    std::transform(res->begin(), res->end(), res->begin(), std::bind2nd(std::minus<double>(), mean));
    std::copy(res->begin(), res->end(), cy.begin());
    
    double mincrit = std::numeric_limits<double>::max();
    int imin;
    for (int i = 0; i < n; ++i)
    {
      double tmp = std::abs(cx[i]) + std::abs(cy[i]);
      if (tmp < mincrit)
      {
        mincrit = tmp;
        imin = i;
      }
    }
    startPoints[0] = imin;
    geomap.init(startPoints, dist);
    geomap.process();
    res = geomap.distanceMap();
    mean = std::accumulate(res->begin(), res->end(), 0.0);
    mean /= res->size();
    std::transform(res->begin(), res->end(), res->begin(), std::bind2nd(std::minus<double>(), mean));
    std::copy(res->begin(), res->end(), cz.begin());
  }
  std::cout << "OK" << std::endl;

  std::cout << "Writing distance maps..." << std::flush;
  {
    Texture1d t;
    t.reserve(getVertices(mesh).size());
    for (std::size_t i = 0; i < getVertices(mesh).size(); ++i)
    {
      t.push_back(cx[i]);
    }
    Writer<Texture1d> w("cx");
    w.write(t);
  }
  {
    Texture1d t;
    t.reserve(getVertices(mesh).size());
    for (std::size_t i = 0; i < getVertices(mesh).size(); ++i)
    {
      t.push_back(cy[i]);
    }
    Writer<Texture1d> w("cy");
    w.write(t);
  }
  {
    Texture1d t;
    t.reserve(getVertices(mesh).size());
    for (std::size_t i = 0; i < getVertices(mesh).size(); ++i)
    {
      t.push_back(cz[i]);
    }
    Writer<Texture1d> w("cz");
    w.write(t);
  }
  std::cout << "OK" << std::endl;
  */

  
  std::cout << "Computing Gonzalez clustering..." << std::flush;
  typedef std::vector<til::numeric_array<float,3> >::const_iterator iterator;
  shared_ptr<std::vector<iterator> > quant;
  {
    til::GonzalezClustering<std::vector<til::numeric_array<float,3> >, float > gc(getVertices(mesh));
    gc.clusterize_maxDiam(til::square(sigma / beta));
    //shared_ptr<std::vector<iterator> > quant = gc.quantization();
    quant = gc.quantization();
  }
  std::cout << "OK" << std::endl;  
  std::cout << "Found " << quant->size() << " clusters" << std::endl;
  
  std::cout << "Computing distance maps" << std::flush;
  //til::Triangle_mesh_geodesic_map<til::Mesh_N, double > geomap(mesh);
  typedef std::vector<std::vector<std::size_t> > CNeighborhoods;
  shared_ptr<CNeighborhoods> pneighc = til::circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  til::Triangle_mesh_geodesic_map<til::Mesh_N::VertexCollection, CNeighborhoods, double >
    geomap(getVertices(mesh), *pneighc);
  std::vector<std::size_t> startPoints(1);
  std::vector<double> dist(1, 0.0);
  std::vector<std::vector<double> > A(quant->size(), std::vector<double>(n));
  til::math::Gaussian<double> gaussian(sigma);
  for (std::size_t i = 0; i < quant->size(); ++i)
  {
    std::cout << "." << std::flush;
    startPoints[0] = std::distance(iterator(getVertices(mesh).begin()), (*quant)[i]);
    geomap.init(startPoints, dist);
    geomap.process();
    shared_ptr<std::vector<double> > tmp = geomap.distanceMap();
    std::transform(tmp->begin(), tmp->end(), A[i].begin(), gaussian);
  }
  std::cout << "OK" << std::endl;

  std::cout << "Computing sum of bases..." << std::flush;
  std::vector<double> sum(n, 0.0);
  {
    for (std::size_t i = 0; i < quant->size(); ++i)
    {
      for (std::size_t j = 0; j < n; ++j)
      {
        sum[j] += A[i][j];
      }
    }
  }
  std::cout << "OK" << std::endl;
  std::cout << "Sum:  " << *std::min_element(sum.begin(), sum.end()) << " " << *std::max_element(sum.begin(), sum.end()) << std::endl;

  std::cout << "Normalizing..." << std::flush;
  {
    for (std::size_t i = 0; i < quant->size(); ++i)
    {
      for (std::size_t j = 0; j < n; ++j)
      {
        A[i][j] /= sum[j];
      }
    }
  }
  std::cout << "OK" << std::endl;
  

  /*
  std::cout << "Writing sum..." << std::flush;
  {
    Texture1d t;
    t.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
    {
      t.push_back(sum[i]);
    }
    Writer<Texture1d> w("sum");
    w.write(t);
  }
  std::cout << "OK" << std::endl;
  */
  
  // TODO: could also add "x-y-z"-like sphere coordinates
  //A.push_back(std::vector<double>(n, 1.0));
  
  /*
  std::vector<std::vector<double> > A(4, std::vector<double>(n));
  A[0] = std::vector<double>(n, 10.0);
  A[1] = cx;
  A[2] = cy;
  A[3] = cz;
  */
  
  /*
  std::cout << "Writing distance maps..." << std::flush;
  {
    Texture1d t;
    t.reserve(getVertices(mesh).size());
    for (std::size_t i = 0; i < getVertices(mesh).size(); ++i)
    {
      t.push_back(A[1][i]);
    }
    Writer<Texture1d> w("surfgauss");
    w.write(t);
  }
  std::cout << "OK" << std::endl;
  */
  
  /*
  std::cout << "Constructing linear system..." << std::flush;
  std::size_t m = A.size();
  double * AA = new double[m*m];
  for (std::size_t i = 0; i < m; ++i)
  {
    for (std::size_t j = i; j < m; ++j)
    {
      AA[i + j*m] = AA[j + i*m] = std::inner_product(A[i].begin(), A[i].end(), A[j].begin(), 0.0);
    }
  }
  double * AB = new double[3*m];
  for (std::size_t i = 0; i < m; ++i)
  {
    AB[i]     = 0.0;
    AB[i+m]   = 0.0;
    AB[i+2*m] = 0.0;
    for (std::size_t j = 0; j < n; ++j)
    {
      AB[i]     += A[i][j] * getVertices(mesh)[j][0];
      AB[i+m]   += A[i][j] * getVertices(mesh)[j][1];
      AB[i+2*m] += A[i][j] * getVertices(mesh)[j][2];
    }
  }
  std::cout << "OK" << std::endl;
  */
  
  /*
  std::cout << "A = [" << std::endl;
  for (std::size_t i = 0; i < m; ++i)
  {
    for (std::size_t j = 0; j < m; ++j)
    {
      std::cout << AA[i + j*m] << " ";
    }
  std::cout << std::endl;
  }
  std::cout << "]" << std::endl;
  */
  
  /*
  std::cout << "B = [" << std::endl;
  for (std::size_t i = 0; i < m; ++i)
  {
    std::cout << AB[i] << std::endl;
  }
  std::cout << "]" << std::endl;
  */
  
  /*
  std::cout << "Computing eigenvectors..." << std::flush;
  {
    char JOBZ = 'V';
    char UPLO = 'U';
    int N = m;
    double * W = new double[m];
    int LWORK = 5*m;
    double * WORK = new double[LWORK];
    int INFO;
    
    dsyev_( &JOBZ, &UPLO, &N, AA, &N, W, WORK, & LWORK, &INFO );
    delete [] WORK;
    
    if (INFO < 0)
      std::cout << "Invalid argument number " << -INFO << std::endl;
    else if (INFO > 0)
      std::cout << "Failure, minor " << INFO << std::endl;

    for (std::size_t i = 0; i < m; ++i) std::cout << W[i] << std::endl;
        
    delete [] W;
  }
  std::cout << "OK" << std::endl;
  */

  /*
  for (std::size_t i = 0; i < m; ++i)
  {
    for (std::size_t j = 0; j < m; ++j)
    {
      if (i == j)
        AA[i + j*m] = 1.0;
      else 
        AA[i + j*m] = 0.0;
    }
  }
  AA[m] = 1.0;
  */

  /*
  std::cout << "A = [" << std::endl;
  for (std::size_t i = 0; i < m; ++i)
  {
    for (std::size_t j = 0; j < m; ++j)
    {
      std::cout << AA[i + j*m] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "]" << std::endl;
  */

  /*
  std::cout << "Constructing eigenmaps..." << std::flush;
  std::vector<std::vector<double> > eigs(m, std::vector<double>(n, 0.0));
  for (std::size_t j = 0; j < m; ++j)
  {
    for (std::size_t i = 0; i < n; ++i)
    {
      for (std::size_t k = 0; k < m; ++k)
      {
        //eigs[j][i] +=  A[k][i] * AA[j+m*k];
        eigs[m-1-j][i] +=  A[k][i] * AA[k+m*j];
      }
    }
  }
  std::cout << "OK" << std::endl;
  */

  /*
  std::cout << "checking orthogonality..." << std::flush;
  {
    for (std::size_t i = 0; i < m; ++i)
    {
      for (std::size_t j = i; j < m; ++j)
      {
        std::cout << i << " " << j << " " << std::inner_product(eigs[i].begin(), eigs[i].end(), eigs[j].begin(), 0.0) << std::endl;
      }
    }
  }
  std::cout << "OK" << std::endl;
  */
  
  /*
  std::cout << "Writing maps..." << std::flush;
  {
    Texture1d t(1, getVertices(mesh).size());
    char s[500];
    for (int j = 0; j < m; ++j)
    {
      for (std::size_t i = 0; i < getVertices(mesh).size(); ++i)
      {
        t.item(i) = eigs[j][i];
      }
      sprintf(s, "eigenvec%i", j);
      Writer<Texture1d> w(s);
      w.write(t);      
    }
  }
  std::cout << "OK" << std::endl;
  */
  /*
  for (int d = 4; d < 10; ++d)
  {
    std::cout << "Constructing linear system..." << std::flush;
    std::size_t m = d;
    double * AA = new double[m*m];
    for (std::size_t i = 0; i < m; ++i)
    {
      for (std::size_t j = i; j < m; ++j)
      {
        AA[i + j*m] = AA[j + i*m] = std::inner_product(eigs[i].begin(), eigs[i].end(), eigs[j].begin(), 0.0);
      }
    }
    double * AB = new double[3*m];
    for (std::size_t i = 0; i < m; ++i)
    {
      AB[i]     = 0.0;
      AB[i+m]   = 0.0;
      AB[i+2*m] = 0.0;
      for (std::size_t j = 0; j < n; ++j)
      {
        AB[i]     += eigs[i][j] * vertices[j][0];
        AB[i+m]   += eigs[i][j] * vertices[j][1];
        AB[i+2*m] += eigs[i][j] * vertices[j][2];
      }
    }
    std::cout << "OK" << std::endl;
    std::cout << "Solving linear system..." << std::flush;
    {
      char UPLO = 'U';
      int N = m;
      int NRHS = 3;
      int INFO;
      dposv_(&UPLO, &N, &NRHS, AA, &N, AB, &N, &INFO);
      if (INFO < 0)
        std::cout << "Invalid argument number " << -INFO << std::endl;
      else if (INFO > 0)
        std::cout << "System is not definite positive, minor " << INFO << std::endl;
    }
    std::cout << "OK" << std::endl;
    std::cout << "Reconstructing mesh..." << std::flush;
    for (std::size_t i = 0; i < n; ++i)
    {
      getVertices(mesh)[i][0] = getVertices(mesh)[i][1] = getVertices(mesh)[i][2] = 0.0f;
      for (std::size_t j = 0; j < m; ++j)
      {
        getVertices(mesh)[i][0] += A[j][i] * AB[j];
        getVertices(mesh)[i][1] += A[j][i] * AB[j+m];
        getVertices(mesh)[i][2] += A[j][i] * AB[j+2*m];
      }
    }
    std::cout << "OK" << std::endl;
    
    delete [] AB;
    delete [] AA;
  
    
    std::cout << "Writing mesh..." << std::flush;
    {
      char name[512];
      sprintf(name, "rec%i", d);
      Writer<AimsSurfaceTriangle> w(name);
      AimsSurfaceTriangle s;
      til::convert(s, mesh);
      w.write(s);
    }   
    std::cout << "OK" << std::endl;
  }
  */
  
  {
    std::cout << "Constructing linear system..." << std::flush;
    std::size_t m = A.size();
    double * AA = new double[m*m];
    for (std::size_t i = 0; i < m; ++i)
    {
      for (std::size_t j = i; j < m; ++j)
      {
        AA[i + j*m] = AA[j + i*m] = std::inner_product(A[i].begin(), A[i].end(), A[j].begin(), 0.0);
      }
    }
    double * AB = new double[3*m];
    for (std::size_t i = 0; i < m; ++i)
    {
      AB[i]     = 0.0;
      AB[i+m]   = 0.0;
      AB[i+2*m] = 0.0;
      for (std::size_t j = 0; j < n; ++j)
      {
        AB[i]     += A[i][j] * getVertices(mesh)[j][0];
        AB[i+m]   += A[i][j] * getVertices(mesh)[j][1];
        AB[i+2*m] += A[i][j] * getVertices(mesh)[j][2];
      }
    }
    std::cout << "OK" << std::endl;
  
    std::cout << "Solving linear system..." << std::flush;
    {
      char UPLO = 'U';
      int N = m;
      int NRHS = 3;
      int INFO;
      dposv_(&UPLO, &N, &NRHS, AA, &N, AB, &N, &INFO);
      if (INFO < 0)
        std::cout << "Invalid argument number " << -INFO << std::endl;
      else if (INFO > 0)
        std::cout << "System is not definite positive, minor " << INFO << std::endl;
    }
    
    std::cout << "Reconstructing mesh..." << std::flush;
    for (std::size_t i = 0; i < n; ++i)
    {
      getVertices(mesh)[i][0] = getVertices(mesh)[i][1] = getVertices(mesh)[i][2] = 0.0f;
      for (std::size_t j = 0; j < m; ++j)
      {
        getVertices(mesh)[i][0] += A[j][i] * AB[j];
        getVertices(mesh)[i][1] += A[j][i] * AB[j+m];
        getVertices(mesh)[i][2] += A[j][i] * AB[j+2*m];
      }
    }
    std::cout << "OK" << std::endl;
    
    delete [] AB;
    delete [] AA;
  
    
    std::cout << "Writing mesh..." << std::flush;
    {
      AimsSurfaceTriangle s;
      til::convert(s, mesh);
      w.write(s);
    }   
    std::cout << "OK" << std::endl;
  }
}


void testGaussPatchApprox2( int argc, char* argv[] )
{
  til::Mesh_N mesh;
  Writer<AimsSurfaceTriangle> w;
  std::cout << "Reading mesh..." << std::flush;
  float sigma;
  float beta;
  float alpha;
  {
    Reader<AimsSurfaceTriangle> r;
    AimsApplication app( argc, aims_const_hack(argv), "testGaussPatchApprox2" );
    app.addOption(r, "-i", "input mesh" );
    app.addOption(w, "-o", "output mesh" );
    app.addOption(sigma, "-sigma", "sigma" );
    app.addOption(beta, "-beta", "beta" );
    app.addOption(alpha, "-alpha", "alpha" );
    
    app.initialize();

    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  std::cout << "OK" << std::endl;
  std::size_t n = getVertices(mesh).size();

  std::vector<til::numeric_array<float, 3> > vertices(n);
  std::copy(getVertices(mesh).begin(), getVertices(mesh).end(), vertices.begin());

  
  std::cout << "Computing Gonzalez clustering..." << std::flush;
  typedef std::vector<til::numeric_array<float,3> >::const_iterator iterator;
  shared_ptr<std::vector<iterator> > quant;
  {
    til::GonzalezClustering<std::vector<til::numeric_array<float,3> >, float > gc(getVertices(mesh));
    gc.clusterize_maxDiam(til::square(sigma / beta));
    quant = gc.quantization();
  }
  std::cout << "OK" << std::endl;  
  std::cout << "Found " << quant->size() << " clusters" << std::endl;
  
  std::cout << "Computing distance maps" << std::flush;
  typedef std::vector<std::vector<std::size_t> > CNeighborhoods;
  shared_ptr<CNeighborhoods> pneighc = til::circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  //til::Triangle_mesh_geodesic_map<til::Mesh_N, double > geomap(mesh);
  til::Triangle_mesh_geodesic_map<til::Mesh_N::VertexCollection, CNeighborhoods, double >
    geomap(getVertices(mesh), *pneighc);
  std::vector<std::size_t> startPoints(1);
  std::vector<double> dist(1, 0.0);
  std::vector<std::vector<double> > A(quant->size(), std::vector<double>(n));
  til::math::Gaussian<double> gaussian(sigma);
  for (std::size_t i = 0; i < quant->size(); ++i)
  {
    std::cout << "." << std::flush;
    startPoints[0] = std::distance(iterator(getVertices(mesh).begin()), (*quant)[i]);
    geomap.init(startPoints, dist);
    geomap.process();
    shared_ptr<std::vector<double> > tmp = geomap.distanceMap();
    std::transform(tmp->begin(), tmp->end(), A[i].begin(), gaussian);
  }
  std::cout << "OK" << std::endl;

  A.push_back(std::vector<double>(n, 1.0/5.0));
  
  std::cout << "Constructing linear system..." << std::flush;
  std::size_t m = A.size();
  double * AA = new double[m*m];
  for (std::size_t i = 0; i < m; ++i)
  {
    for (std::size_t j = i; j < m; ++j)
    {
      AA[i + j*m] = AA[j + i*m] = std::inner_product(A[i].begin(), A[i].end(), A[j].begin(), 0.0);
    }
  }
  std::cout << "OK" << std::endl;

  std::cout << "Computing eigenvectors..." << std::flush;
  std::size_t K;
  double * W = new double[m];
  {
    char JOBZ = 'V';
    char UPLO = 'U';
    int N = m;
    int LWORK = 5*m;
    double * WORK = new double[LWORK];
    int INFO;
    
    dsyev_( &JOBZ, &UPLO, &N, AA, &N, W, WORK, & LWORK, &INFO );
    delete [] WORK;
    
    if (INFO < 0)
      std::cout << "Invalid argument number " << -INFO << std::endl;
    else if (INFO > 0)
      std::cout << "Failure, minor " << INFO << std::endl;

    std::cout << "Eigenvalues: ";
    for (std::size_t i = 0; i < m; ++i) std::cout << W[i] << std::endl;

    double sum = std::accumulate(W, W+m, 0.0);
    double partialSum = 0;
    for (K = 0; K < m; ++K)
    {
      partialSum += W[m-1-K];
      if (partialSum > alpha * sum)
      {
        ++K;
        break;
      }
    }
  }
  std::cout << "OK" << std::endl;
  std::cout << "Cropping at " << K << " eigenvectors" << std::endl;

  std::cout << "Constructing eigenmaps..." << std::flush;
  std::vector<std::vector<double> > eigs(K, std::vector<double>(n, 0.0));
  for (std::size_t j = 0; j < K; ++j)
  {
    for (std::size_t i = 0; i < n; ++i)
    {
      for (std::size_t k = 0; k < m; ++k)
      {
        eigs[j][i] +=  A[k][i] * AA[k+m*(m-1-j)] / std::sqrt(W[m-1-j]);
      }
    }
  }
  std::cout << "OK" << std::endl;

  std::cout << "Writing maps..." << std::flush;
  {
    Texture1d t(1, getVertices(mesh).size());
    char s[500];
    for (std::size_t j = 0; j < K; ++j)
    {
      for (std::size_t i = 0; i < getVertices(mesh).size(); ++i)
      {
        t.item(i) = eigs[j][i];
      }
      sprintf(s, "eigenvec%i", j);
      Writer<Texture1d> w(s);
      w.write(t);      
    }
  }
  std::cout << "OK" << std::endl;
  
  std::cout << "checking orthonormality..." << std::flush;
  {
    for (std::size_t i = 0; i < K; ++i)
    {
      for (std::size_t j = i; j < K; ++j)
      {
        std::cout << i << " " << j << " " << std::inner_product(eigs[i].begin(), eigs[i].end(), eigs[j].begin(), 0.0) << std::endl;
      }
    }
  }
  std::cout << "OK" << std::endl;
  

  std::vector<double> coeffx(K, 0.0);
  std::vector<double> coeffy(K, 0.0);
  std::vector<double> coeffz(K, 0.0);
  for (std::size_t i = 0; i < K; ++i)
  {
    for (std::size_t j = 0; j < n; ++j)
    {
      coeffx[i] += eigs[i][j] * getVertices(mesh)[j][0];
      coeffy[i] += eigs[i][j] * getVertices(mesh)[j][1];
      coeffz[i] += eigs[i][j] * getVertices(mesh)[j][2];
    }
  }
  
  for (std::size_t d = 4; d < K; ++d)
  {
    std::cout << "Reconstructing mesh..." << std::flush;
    for (std::size_t i = 0; i < n; ++i)
    {
      getVertices(mesh)[i][0] = getVertices(mesh)[i][1] = getVertices(mesh)[i][2] = 0.0f;
      for (std::size_t j = 0; j < d; ++j)
      {
        getVertices(mesh)[i][0] += eigs[j][i] * coeffx[j];
        getVertices(mesh)[i][1] += eigs[j][i] * coeffy[j];
        getVertices(mesh)[i][2] += eigs[j][i] * coeffz[j];
      }
    }
    std::cout << "OK" << std::endl;
    
    std::cout << "Writing mesh..." << std::flush;
    {
      char name[512];
      sprintf(name, "rec%i", d);
      Writer<AimsSurfaceTriangle> w(name);
      AimsSurfaceTriangle s;
      til::convert(s, mesh);
      w.write(s);
    }   
    std::cout << "OK" << std::endl;
  }
  
  delete [] AA;
  delete [] W;  
}



void testLambdaMu( int argc, char* argv[] )
{
  til::Mesh_N mesh;
  Writer<AimsSurfaceTriangle> w;
  std::cout << "Reading mesh..." << std::flush;
  float lambda;
  float mu;
  int niter;
  {
    Reader<AimsSurfaceTriangle> r;
    AimsApplication app( argc, aims_const_hack(argv), "testLambdaMu" );
    app.addOption(r, "-i", "input mesh" );
    app.addOption(w, "-o", "output mesh" );
    app.addOption(lambda, "-lambda", "" );
    app.addOption(mu, "-mu", "" );
    app.addOption(niter, "-niter", "");
    app.initialize();
    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  std::cout << "OK" << std::endl;
  std::size_t n = getVertices(mesh).size();

  /*
   shared_ptr<std::vector<std::vector<til::Point<float,3>*> > > neighc
    = circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
 */

  
  {
    std::vector<til::numeric_array<float, 3> > vertices(n);
  
    typedef til::xsr::Integer_access<std::vector<til::numeric_array<float,3> > >              VertexIndexing;
    typedef til::xsr::Integer_access<std::vector<std::vector<std::size_t> > >                 NeighborIndexing;
    typedef til::xsr::Iterator_access<std::vector<til::numeric_array<float,3> >::iterator>    VertexIndexing2;
    
    //typedef til::LaplacianSmoothing<VertexIndexing, NeighborIndexing, VertexIndexing2> Smoother;
    typedef til::ScaleIndependantLaplacianSmoothing<VertexIndexing, NeighborIndexing, VertexIndexing2> Smoother;
    VertexIndexing tmp1(getVertices(mesh));
    NeighborIndexing tmp2(getNeighborIndices(mesh));
    VertexIndexing2 tmp3;
    Smoother lambdaSmoothing(tmp1, tmp2, tmp3, lambda);
    VertexIndexing tmp4(vertices);
    NeighborIndexing tmp5(getNeighborIndices(mesh));
    Smoother muSmoothing(tmp4, tmp5, tmp3, -mu);

    //til::LambdaMuSmoothing<VertexIndexing, NeighborIndexing> lambdaMuSmoothing(VertexIndexing(getVertices(mesh)), NeighborIndexing(getNeighborIndices(mesh)), lambda, mu, getVertices(mesh).size());
    
    for (int i = 0; i < niter; ++i)
    {
      std::cout << "iter " << i << std::endl;
      lambdaSmoothing(0, getVertices(mesh).size(), vertices.begin());
      muSmoothing(0, vertices.size(), getVertices(mesh).begin());
      //lambdaMuSmoothing(0, getVertices(mesh).size());
    }
  }
    
  /* 
  float coeff;
  for (int i = 0; i < 2*niter; ++i)
  {
    coeff = ( i%2 ? -mu : lambda );
    std::cout << "iter " << i << std::endl;
    std::copy(getVertices(mesh).begin(), getVertices(mesh).end(), vertices.begin());
    for (std::size_t j = 0; j < n; ++j)
    {
      vertices[j] *= (1-coeff);
      std::size_t nneigh = getNeighborIndices(mesh)[j].size();
      //std::size_t nneigh = (*neighc)[i].size();
      for (std::size_t k = 0; k < nneigh; ++k)
      {
        vertices[j] += (coeff / nneigh) * getVertices(mesh)[getNeighborIndices(mesh)[j][k]];
      }
    }
    std::copy(vertices.begin(), vertices.end(), getVertices(mesh).begin());
  }
  */
  
  std::cout << "Writing mesh..." << std::flush;
  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh);
    w.write(s);
  }   
  std::cout << "OK" << std::endl;
}

// bin/toto -i 

void testMultiscaleRegistration( int argc, char* argv[] )
{
  til::Mesh1 mesh;
  Writer<AimsSurfaceTriangle> w;
  std::cout << "Reading mesh..." << std::flush;
  float lambda;
  float mu;
  int niter;
  int nlevel;
  double freq;
  double ampl;
  double sigma;
  double beta;
  {
    Reader<AimsSurfaceTriangle> r;
    AimsApplication app( argc, aims_const_hack(argv), "testMultiscaleRegistration" );
    app.addOption(r, "-i", "input mesh" );
    app.addOption(w, "-o", "output mesh" );
    app.addOption(lambda, "-lambda", "" );
    app.addOption(mu, "-mu", "" );
    app.addOption(niter, "-niter", "");
    app.addOption(nlevel, "-nlevel", "");
    app.addOption(freq, "-freq", "");
    app.addOption(ampl, "-ampl", "");
    app.addOption(sigma, "-sigma", "");
    app.addOption(beta, "-beta", "");
    
    app.initialize();
    AimsSurfaceTriangle s;
    r.read( s );
    til::convert(mesh, s);
  }
  std::cout << "OK" << std::endl;
  
  til::Mesh1 mesh2;
  mesh2.vertices() = shared_ptr<std::vector<til::numeric_array<float, 3> > >(new std::vector<til::numeric_array<float, 3> >(getVertices(mesh).size()));
  for (std::size_t i = 0; i < getVertices(mesh).size(); ++i)
  {
    getVertices(mesh2)[i] = getVertices(mesh)[i];
  }
  mesh2.faces() = shared_ptr<std::vector<til::numeric_array<std::size_t, 3> > >(new std::vector<til::numeric_array<std::size_t, 3> >(getFaceIndices(mesh).size()));
  for (std::size_t i = 0; i < getFaceIndices(mesh).size(); ++i)
  {
    getFaceIndices(mesh2)[i] = getFaceIndices(mesh)[i];
  }
  
  for (std::size_t i = 0; i < getVertices(mesh).size(); ++i)
  {
    getVertices(mesh)[i][0] += ampl * std::sin( getVertices(mesh)[i][1] / freq ) * std::sin( getVertices(mesh)[i][2] / freq );
    getVertices(mesh)[i][1] += ampl * std::sin( getVertices(mesh)[i][0] / freq ) * std::sin( getVertices(mesh)[i][2] / freq );
    getVertices(mesh)[i][2] += ampl * std::sin( getVertices(mesh)[i][0] / freq ) * std::sin( getVertices(mesh)[i][1] / freq );
  }
  
  typedef std::vector<std::pair<std::vector<til::numeric_array<float,3> >, std::vector<til::numeric_array<std::size_t, 3> > > >
    Meshes;
  Meshes meshes(nlevel);
  meshes[0].first = getVertices(mesh);
  meshes[0].second = getFaceIndices(mesh);

  Meshes meshes2(nlevel);
  meshes2[0].first = getVertices(mesh2);
  meshes2[0].second = getFaceIndices(mesh2);

  for (int i = 1; i < nlevel; ++i)
  {
    std::size_t n = getVertices(mesh).size();
    shared_ptr<std::vector<std::vector<std::size_t> > > neigh = getNeighborIndices(mesh);
        
    std::cout << "Smoothing mesh..." << std::flush;
    {
      typedef til::xsr::Integer_access<std::vector<til::numeric_array<float,3> > >              VertexIndexing;
      typedef til::xsr::Integer_access<std::vector<std::vector<std::size_t> > >                 NeighborIndexing;
      typedef til::xsr::Iterator_access<std::vector<til::numeric_array<float,3> >::iterator>    VertexIndexing2;
      typedef til::LaplacianSmoothing<VertexIndexing, VertexIndexing2, NeighborIndexing>        Smoother;
      //typedef til::ScaleIndependantLaplacianSmoothing<VertexIndexing, NeighborIndexing, VertexIndexing2> Smoother;
      
      std::vector<til::numeric_array<float, 3> > vertices(n);
      VertexIndexing tmp1(getVertices(mesh));
      NeighborIndexing tmp2(*neigh);
      VertexIndexing2 tmp3;
      Smoother lambdaSmoothing(tmp1, tmp3, tmp2, lambda);
      VertexIndexing tmp4(vertices);
      NeighborIndexing tmp5(*neigh);
      Smoother muSmoothing(tmp4, tmp3, tmp5, -mu);

      for (int i = 0; i < niter; ++i)
      {
        lambdaSmoothing(0, n, vertices.begin());
        muSmoothing(0, n, getVertices(mesh).begin());
      }
    }
    std::cout << "OK" << std::endl;

    std::cout << "Subsampling mesh..." << std::flush;
    {
      til::quantizer(getVertices(mesh), getFaceIndices(mesh), std::size_t(getVertices(mesh).size() / 4.0));
    }
    std::cout << "OK" << std::endl;

    meshes[i].first = getVertices(mesh);
    meshes[i].second = getFaceIndices(mesh);
  }
  
  
  for (int i = 1; i < nlevel; ++i)
  {
    std::size_t n = getVertices(mesh2).size();
    shared_ptr<std::vector<std::vector<std::size_t> > > neigh = getNeighborIndices(mesh2);
        
    std::cout << "Smoothing mesh..." << std::flush;
    {
      typedef til::xsr::Integer_access<std::vector<til::numeric_array<float,3> > >              VertexIndexing;
      typedef til::xsr::Integer_access<std::vector<std::vector<std::size_t> > >                 NeighborIndexing;
      typedef til::xsr::Iterator_access<std::vector<til::numeric_array<float,3> >::iterator>    VertexIndexing2;
      typedef til::LaplacianSmoothing<VertexIndexing, VertexIndexing2, NeighborIndexing>        Smoother;
      //typedef til::ScaleIndependantLaplacianSmoothing<VertexIndexing, NeighborIndexing, VertexIndexing2> Smoother;
      
      std::vector<til::numeric_array<float, 3> > vertices(n);
      VertexIndexing tmp1(getVertices(mesh2));
      NeighborIndexing tmp2(*neigh);
      VertexIndexing2 tmp3;
      Smoother lambdaSmoothing(tmp1, tmp3, tmp2, lambda);
      VertexIndexing tmp4(vertices);
      NeighborIndexing tmp5(*neigh);
      Smoother muSmoothing(tmp4, tmp3, tmp5, -mu);

      for (int i = 0; i < niter; ++i)
      {
        lambdaSmoothing(0, n, vertices.begin());
        muSmoothing(0, n, getVertices(mesh2).begin());
      }
    }
    std::cout << "OK" << std::endl;

    std::cout << "Subsampling mesh..." << std::flush;
    {
      til::quantizer(getVertices(mesh2), getFaceIndices(mesh2), std::size_t(getVertices(mesh2).size() / 4.0));
    }
    std::cout << "OK" << std::endl;

    meshes2[i].first = getVertices(mesh2);
    meshes2[i].second = getFaceIndices(mesh2);
  }
    
  for (int i = 0; i < nlevel; ++i)
  {
    std::cout << "Writing mesh..." << std::flush;
    {
      char name[512];
      sprintf(name, "level%i-0.mesh", i);
      AimsSurfaceTriangle s;
      til::Mesh1 tmp;
      getVertices(tmp) = meshes[i].first;
      getFaceIndices(tmp) = meshes[i].second;
      til::convert(s, tmp);
      Writer<AimsSurfaceTriangle> w(name);
      w.write(s);
    }
    std::cout << "OK" << std::endl;
    std::cout << "Writing mesh..." << std::flush;
    {
      char name[512];
      sprintf(name, "level%i-1.mesh", i);
      AimsSurfaceTriangle s;
      til::Mesh1 tmp;
      getVertices(tmp) = meshes2[i].first;
      getFaceIndices(tmp) = meshes2[i].second;
      til::convert(s, tmp);
      Writer<AimsSurfaceTriangle> w(name);
      w.write(s);
    }
    std::cout << "OK" << std::endl;
  }
  

  // multiresolution registration
  
  /*
  for (int i = 0; i < nlevel; ++i)
  {
    sigma *= 2.0;
  }
  */
  
  //double sigma = 30.0;
  typedef std::vector<til::numeric_array<float,3> >::const_iterator iterator;
  shared_ptr<std::vector<std::vector<std::size_t> > > neighc;
  shared_ptr<std::vector<iterator> > quant;
  std::cout << "sigma = " << sigma << std::endl;
  for (int level = nlevel-1; level >= 0; --level)
  {
    std::cout << "Multiresolution level " << level << " (sigma = " << sigma << ")" << std::endl;
    
    //std::size_t n = meshes[level].first.size();
    til::math::IsotropicGaussianKernel<til::numeric_array<float,3>, double> gaussian(sigma);
    
    std::cout << "Computing neighborhoods..." << std::flush;
    neighc = circular_neighborhoods(meshes[level].first, meshes[level].second);
    std::cout << "OK" << std::endl;

    std::cout << "Computing Gonzalez clustering..." << std::flush;
    quant = til::gonzalez_clustering(meshes2[level].first, til::square(sigma / beta));
    std::cout << "OK" << std::endl;
    
    std::cout << "Computing normalized RBF with their support..." << std::flush;
    std::vector<std::vector<std::pair<std::size_t, double> > > defoBasis(quant->size());
    computingNormalizedRBFWithSupport(meshes2[level].first, quant, gaussian, defoBasis);
    std::cout << "OK" << std::endl;
    /*
    {
      std::cout << "Computing RBF trace..." << std::flush;
      std::vector<std::vector<double> > rbf(quant->size(), std::vector<double>(n));
      for (std::size_t j = 0; j < rbf.size(); ++j)
      {
        //std::transform(meshes[level].first.begin(), meshes[level].first.end(), rbf[j].begin(), bind2nd(gaussian, *(*quant)[j]));
        for (std::size_t i = 0; i < n; ++i)
        {
          rbf[j][i] = gaussian((*quant)[j]->data(), meshes2[level].first[i].data());
        }
        
      }
      std::cout << "OK" << std::endl;
  
      // It's probably better to threshold the RBF before normalization; first because thresholding does not
      // strictly preserve normalization (even though for mild thresholds it should not impact much); second
      // because before normalization we know the shape of the kernel and in particular its maximum of 1, meaning that
      // the threshold can be taken more confidently; last, because normalization does not undo thresholding -- zeroed
      // coefficients remain null.
      std::cout << "Thresholding RBF..." << std::flush;
      {
        for (std::size_t j = 0; j < rbf.size(); ++j)
        {
          for (std::size_t i = 0; i < n; ++i)
          {
            if (rbf[j][i] < 10e-4) rbf[j][i] = 0.0;
          }
        }
      }
      std::cout << "OK" << std::endl;
  
      std::cout << "Normalizing RBFs..." << std::flush;
      {
        std::vector<double> sum(n, 0.0);
        for (std::size_t j = 0; j < rbf.size(); ++j)
        {
          std::transform(rbf[j].begin(), rbf[j].end(), sum.begin(), sum.begin(), std::plus<double>());
          / *
          for (std::size_t i = 0; i < n; ++i)
          {
            sum[i] += rbf[j][i];
          }
          * /
        }
        for (std::size_t j = 0; j < rbf.size(); ++j)
        {
          std::transform(rbf[j].begin(), rbf[j].end(), sum.begin(), rbf[j].begin(), std::divides<double>());
          / *
          for (std::size_t i = 0; i < n; ++i)
          {
            rbf[j][i] /= sum[i];
          }
          * /
        }
      }
      std::cout << "OK" << std::endl;
      
      std::cout << "Writing defo basis..." << std::flush;
      {
        char name[256];
        for (std::size_t i = 0; i < min(std::size_t(10), rbf.size()); ++i)
        {
          sprintf(name, "defobase%i", i);
          til::writeTexture(rbf[i], name);
        }
      }
      std::cout << "OK" << std::endl;
      
      std::cout << "Converting RBF into compact form..." << std::endl;
      {
        for (std::size_t j = 0; j < rbf.size(); ++j)
        {
          for (std::size_t i = 0; i < n; ++i)
          {
            if (rbf[j][i] > 0.0)
            {
              defoBasis[j].push_back(std::make_pair(i, rbf[j][i]));
            }
          }
        }
      }
      std::cout << "OK" << std::endl;
    }
    */
    
    std::cout << "Generating KDTree..." << std::flush;
    typedef til::KDTree<std::size_t, til::MeshTraits<til::Mesh_N>::VertexCollection> MyKDTree;
    MyKDTree kdtree(meshes[level].first);
    makeKDTree(meshes[level].first, kdtree);
    std::cout << "OK" << std::endl;
    
    std::cout << "Starting minimization (" << 3*quant->size() << " degrees of freedom)" << std::endl;
    typedef til::Find_closest< double, MyKDTree > Finder;
    Finder finder(kdtree);
    typedef RBFReg2 <
      std::vector<til::numeric_array<float,3> >,
      std::vector<std::vector<std::size_t> >,
      Finder
    > Energy;
    //til::math::IsotropicGaussianKernel<til::Point<float,3> > gaussianK(sigma);
    Energy energy(meshes2[level].first, meshes[level].first, *neighc, finder, defoBasis);
    til::Powell<Energy> minalgo(energy);
    // initial scaling estimates of the parameters
    minalgo.initStd() = std::vector<float>(3*quant->size(), 2.0f);
    std::vector<float> params(3*quant->size(), 0.0f);
    params = minalgo(params);
    std::cout << "Minimization finished" << std::endl;
  
  
    
    std::cout << "Applying transformation..." << std::flush;
    {
      std::size_t iParam = 0;
      for (std::vector<std::vector<std::pair<std::size_t, double> > >::iterator iKernel = defoBasis.begin(); iKernel != defoBasis.end(); ++iKernel, iParam += 3)
      {
        // for each kernel, loop on its definition domain
        for (std::vector<std::pair<std::size_t, double> >::iterator iiKernel = iKernel->begin(); iiKernel != iKernel->end(); ++iiKernel)
        {
          for (std::size_t j = 0; j < 3; ++j)
          {
            meshes2[level].first[iiKernel->first][j] += iiKernel->second * params[iParam + j];
          }
        }
      }
    }
    std::cout << "OK" << std::endl;
  
    std::cout << "Writing result..." << std::flush;
    {
      char name[256];
      sprintf(name, "bord%i", level);
      til::write_mesh(meshes2[level].first, meshes2[level].second, name);
    }
    std::cout << "OK" << std::endl;
    
    /*
    std::cout << "Applying transformation..." << std::flush;
    {
      std::size_t n = meshes2[level].first.size();
      std::vector<til::Point<float,3> > newVertices(n, til::Point<float,3>(0,0,0));
      for (std::size_t i = 0; i < n; ++i)
      {
        //newVertices[i] = meshes2[level].first[i];
        double sum = 0.0;
        std::size_t ip = 0;
        for (std::size_t j = 0; j < quant->size(); ++j)
        {
          double coeff = gaussian((*quant)[j]->data(), meshes2[level].first[i].data());
          newVertices[i].data()[0] += params[ip++] * coeff;
          newVertices[i].data()[1] += params[ip++] * coeff;
          newVertices[i].data()[2] += params[ip++] * coeff;
          sum += coeff;
        }
        newVertices[i].data() *= 1/sum;
      }
      //std::transform(newVertices.begin(), newVertices.end(), meshes2[level].first.begin());
      for (std::size_t i = 0; i < n; ++i)
      {
        meshes2[level].first[i].data() += newVertices[i].data();
      }
    }
    std::cout << "OK" << std::endl;

    std::cout << "Writing result..." << std::flush;
    {
      char name[256];
      sprintf(name, "bordd%i", level);
      til::write_mesh(meshes2[level].first, meshes2[level].second, name);
    }
    std::cout << "OK" << std::endl;
    */
    
    //if (level == 0) continue;
    //--level;

    for (int mylevel = level-1; mylevel >= 0; --mylevel)
    {
      std::cout << "Applying transformation..." << std::flush;
      {
        std::size_t n = meshes2[mylevel].first.size();
        std::vector<til::numeric_array<float,3> > newVertices(n, til::numeric_array<float,3>(0,0,0));
        for (std::size_t i = 0; i < n; ++i)
        {
          //newVertices[i] = meshes2[level].first[i];
          double sum = 0.0;
          std::size_t ip = 0;
          for (std::size_t j = 0; j < quant->size(); ++j)
          {
            double coeff = gaussian(*(*quant)[j], meshes2[mylevel].first[i]);
            newVertices[i][0] += params[ip++] * coeff;
            newVertices[i][1] += params[ip++] * coeff;
            newVertices[i][2] += params[ip++] * coeff;
            sum += coeff;
          }
          newVertices[i] *= 1/sum;
        }
        //std::transform(newVertices.begin(), newVertices.end(), meshes2[level].first.begin());
        for (std::size_t i = 0; i < n; ++i)
        {
          meshes2[mylevel].first[i] += newVertices[i];
        }
      }
      std::cout << "OK" << std::endl;
    }

    std::cout << "Writing result..." << std::flush;
    {
      char name[256];
      sprintf(name, "bordd%i", level);
      til::write_mesh(meshes2[level].first, meshes2[level].second, name);
    }
    std::cout << "OK" << std::endl;
  }
}

void conv(AimsData<float> & res)
{
  til::numeric_array<std::size_t, 3> dim;
  dim[0] = res.dimX();
  dim[1] = res.dimY();
  dim[2] = res.dimZ();

  til::MultiDWTND::Direct<til::DWTND<til::DWTCubicConjugate<float, til::bc::Cyclic>::Direct,3>,3>()(&res[0], dim);
  til::multi_dwtND_shuffle(&res[0], dim);
}

void inverseDWT(const std::vector<float> & coeffs, std::size_t s, AimsData<float> & res)
{
  assert(coeffs.size() >= s*s*s);

  til::numeric_array<std::size_t, 3> dim;
  dim[0] = res.dimX();
  dim[1] = res.dimY();
  dim[2] = res.dimZ();
  
  // clear res
  for (int i = 0; i < res.dimX() * res.dimY() * res.dimZ(); ++i)
  {
    res[i] = 0;
  }
  
  // put coeffs inside res
  int n = 0;
  for (std::size_t k = 0; k < s; ++k)
  for (std::size_t j = 0; j < s; ++j)
  for (std::size_t i = 0; i < s; ++i)
  {
    res(i,j,k) = coeffs[n];
    //std::cout << i << " " << j << " " << k << " : " << res(i,j,k) << std::endl;
    ++n;
  }

  /* 
  std::cout << "Writing im " << num << std::flush;
  {
    static int num = 0;
    std::ostringstream os;
    os << "im" << num;
    Writer<AimsData<float> > w(os.str());
    w.write(res);
    ++num;
  }
  std::cout << "OK" << std::endl;
  */
  
  til::multi_dwtND_unshuffle(&res[0], dim);
  til::multi_idwtND_cubic(&res[0], dim);
}

void testWavConjugate()
{
  //AimsData<float> vx(57, 62, 39);
  AimsData<float> vx(64,64,64);
  std::size_t n = 4;
  {
    std::vector<float> coeff(n*n*n, 0.0f);
    inverseDWT(coeff, n, vx);
  }
  til::aimswrite(vx, "toto0");
  {
    std::vector<float> coeff(n*n*n, 0.0f);
    coeff[0] = 1;
    inverseDWT(coeff, n, vx);
  }
  til::aimswrite(vx, "toto1");
  {
    std::vector<float> coeff(n*n*n, 0.0f);
    coeff[7] = 1;
    inverseDWT(coeff, n, vx);
  }
  til::aimswrite(vx, "toto2");
  {
    std::vector<float> coeff(n*n*n, 0.0f);
    coeff[n*n*n-1] = 1;
    inverseDWT(coeff, n, vx);
  }
  til::aimswrite(vx, "toto3");
  {
    double res = 0.0;
    for (int i = 0; i < vx.dimX() * vx.dimY() * vx.dimZ(); ++i)
    {
      res += vx[i] * vx[i];
    }
    std::cout << "hand convol: " << res << std::endl;
  }
  
  til::numeric_array<std::size_t, 3> dim;
  dim[0] = vx.dimX();
  dim[1] = vx.dimY();
  dim[2] = vx.dimZ();
  til::MultiDWTND::Direct<til::DWTND<til::DWTCubicConjugate<float, til::bc::Cyclic>::Direct,3>,3>()(&vx[0], dim);
  til::multi_dwtND_shuffle(&vx[0], dim);
  
  std::cout << "auto convol: " << vx(n-1, n-1, n-1) << std::endl;
  til::aimswrite(vx, "toto4");
}

void generateStupidImages()
{
  {
    AimsData<short> im(64,64,64);  
    for (int k = 0; k < im.dimZ(); ++k)
    for (int j = 0; j < im.dimY(); ++j)
    for (int i = 0; i < im.dimX(); ++i)
    {
      if (til::norm2(i-25, j-25, k-25) < 15*15)
      {
        im(i,j,k) = 1;
      }
    }
    til::aimswrite(im, "ellipse1");
  }
  {
    AimsData<short> im(64,64,64);  
    for (int k = 0; k < im.dimZ(); ++k)
    for (int j = 0; j < im.dimY(); ++j)
    for (int i = 0; i < im.dimX(); ++i)
    {
      if (til::norm2(0.9*(i-25), 1.1*(j-35), 1.1*(k-35)) < 15*15)
      {
        im(i,j,k) = 1;
      }
    }
    til::aimswrite(im, "ellipse2");
  }
}


// /home/cathier/code/c++/myproject-main-linux-release/bin/toto -m /home/Panabase/pascal/fibers/jeanpierre/jeanpierre_Lwhite.clean.transfo.mesh -i /home/Panabase/pascal/fibers/jeanpierre/jeanpierre_9430_t2_diffusion.ima -o tmp3.dim
void meshInterior(int argc, char* argv[])
{
  typedef til::Mesh_N MyMesh;
  MyMesh mesh;
  typedef AimsData<short> Image;
  Image im;
  std::string oname;
  Writer<Image> wi;
  {
    aims::Reader<AimsSurfaceTriangle> rm;
    aims::Reader<Image> ri;
    aims::AimsApplication app(argc, aims_const_hack(argv), "meshInterior");
    app.addOption(rm, "-m", "input mesh");
    app.addOption(ri, "-i", "prototype image");
    app.addOption(wi, "-o", "output image");
    app.initialize();

    til::read_mesh(rm, mesh);
    ri.read(im);
  }
  // clear image
  for (int i = 0; i < im.dimX() * im.dimY() * im.dimZ(); ++i)
  {
    im[i] = 0;
  }
  
  {
    typedef std::vector<til::numeric_array<float, 3> > Vertices;
    typedef std::vector<til::numeric_array<std::size_t, 3> > Faces;
    Vertices & vertices = getVertices(mesh);
    Faces & faces = getFaceIndices(mesh);
    til::numeric_array<float, 3> p;
    std::size_t faceCount = 0, nFaces = faces.size();
    for (Faces::iterator iF = faces.begin(); iF != faces.end(); ++iF)
    {
      indicator(faceCount++, nFaces);
      for (float x = 0; x < 1; x += 0.1)
      for (float y = 0; y < 1-x; y += 0.1)
      {
        p = x * vertices[(*iF)[0]] + y * vertices[(*iF)[1]] + (1-x-y) * vertices[(*iF)[2]];
        im(int(p[0] / im.sizeX() + 0.5), 
           int(p[1] / im.sizeY() + 0.5), 
           int(p[2] / im.sizeZ() + 0.5)) = 1;
      }
    }
  }
  
  /*
  typedef til::KDTree<std::size_t, til::MeshTraits<MyMesh>::VertexCollection> MyKDTree;
  MyKDTree kdtree(getVertices(mesh));
  makeKDTree(getVertices(mesh), kdtree);  
  typedef til::Find_closest<double, MyKDTree> Finder;
  Finder finder(kdtree);
  shared_ptr<std::vector<std::vector<std::size_t> > > neighc = circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  til::numeric_array<double, 3> normal;
  for (int k = 0; k < im.dimZ(); ++k)
  {
    indicator(k, im.dimZ());
    for (int j = 0; j < im.dimY(); ++j)
    for (int i = 0; i < im.dimX(); ++i)
    {
      til::numeric_array<float, 3> tmp(i * im.sizeX(), j * im.sizeY(), k * im.sizeZ());
      std::size_t c = finder(tmp);
      //std::pair<std::size_t, std::size_t> p2 = til::closest_neighborhood<double>(tmp, c, getVertices(mesh), (*neighc)[c]);
      //til::geo::triangle_normal(getVertices(mesh)[c], getVertices(mesh)[p2.first], getVertices(mesh)[p2.second], normal);
      normal = til::closest_normal<double>(tmp, c, getVertices(mesh), (*neighc)[c]);
      if (til::dot(normal, getVertices(mesh)[c] - tmp, til::prec<double>()) >= 0)
      {
        im(i,j,k) = 1;
      }
    }
  }
  */
  wi.write(im);
}

void testSameSign()
{
  std::cout << std::boolalpha << boost::logic::tribool(false) << " " << boost::logic::tribool(true) << " " << boost::logic::tribool(boost::logic::indeterminate) << std::endl;

  std::cout << std::boolalpha << til::fuzzy::is_positive(1.0) << " " << bool(til::fuzzy::is_positive(1.0)) << std::endl;
  std::cout << std::boolalpha << til::fuzzy::is_positive(0.0) << " " << bool(til::fuzzy::is_positive(0.0)) << std::endl;
  std::cout << std::boolalpha << til::fuzzy::is_positive(-1.0) << " " << bool(til::fuzzy::is_positive(-1.0)) << std::endl;
  
  std::cout << std::boolalpha << til::fuzzy::same_sign(1.0,1.0) << std::endl;
  std::cout << std::boolalpha << til::fuzzy::same_sign(-1.0,-1.0) << std::endl;
  std::cout << std::boolalpha << til::fuzzy::same_sign(1.0,-1.0) << std::endl;
  std::cout << std::boolalpha << til::fuzzy::same_sign(-1.0,1.0) << std::endl;
  std::cout << std::boolalpha << til::fuzzy::same_sign(0.0,-1.0) << std::endl;
  std::cout << std::boolalpha << til::fuzzy::same_sign(0.0,1.0) << std::endl;
  std::cout << std::boolalpha << til::fuzzy::same_sign(-1.0,0.0) << std::endl;
  std::cout << std::boolalpha << til::fuzzy::same_sign(1.0,0.0) << std::endl;
}

void defimage(int argc, char * argv[])
{
  typedef float prec_type;
  typedef AimsData<short> AimsImage;
  typedef AimsData<prec_type> AimsFloat;

  AimsFloat im1;
  std::string name_out;
  double step, ampl;
  {
    std::string name_in;
    AimsApplication app(argc, aims_const_hack(argv), "defimage");
    app.addOption(name_in, "-i", "input image" );
    app.addOption(name_out, "-o", "output image" );
    app.addOption(step, "-step", "");
    app.addOption(ampl, "-ampl", "");
    app.initialize();

    std::cout << "Reading image..." << std::flush;
    aims::AnyTypeReader<AimsData<float> > r(name_in);
    r.read(im1);
    std::cout << "OK" << std::endl;
  }
  Writer<AimsFloat> w(name_out);

  AimsFloat im2(im1.dimX(), im1.dimY(), im1.dimZ());
  carto::rc_ptr<aims::Interpolator> interp = getLinearInterpolator(im1);
  
  //aims::LinearInterpolator<float> interp2(&im1, false);
  //interp2._image = &im1;

  float sizeX = im1.sizeX();
  float sizeY = im1.sizeY();
  float sizeZ = im1.sizeZ();
  
  for (int k = 0; k < im1.dimZ(); ++k)
  {
  std::cout << k << "..." << std::flush;
  for (int j = 0; j < im1.dimY(); ++j)
  for (int i = 0; i < im1.dimX(); ++i)
  {
    //im2(i,j,k) = interp2.trilin(i*sizeX + ampl*sin(2*M_PI*j/step), j*sizeY + ampl*sin(2*M_PI*i/step), k*sizeZ);
    //im2(i,j,k) = interp2(i*im1.sizeX() + ampl*sin(2*M_PI*j/step), j*im1.sizeY() + ampl*sin(2*M_PI*i/step), k*im1.sizeZ());
    im2(i,j,k) = (*interp)(i*sizeX + ampl*sin(2*M_PI*j/step), j*sizeY + ampl*sin(2*M_PI*i/step), k*sizeZ);
  }
  }

  std::cout << "Writing transformed image..." << std::flush;
  w.write(im2);
  std::cout << "OK" << std::endl;
  
  std::cout << "Writing transformations..." << std::flush;
  for (int k = 0; k < im1.dimZ(); ++k)
  for (int j = 0; j < im1.dimY(); ++j)
  for (int i = 0; i < im1.dimX(); ++i)
  {
    im2(i,j,k) = ampl*sin(2*M_PI*j/step);
  }
  til::aimswrite(im2, name_out + ".x");
  for (int k = 0; k < im1.dimZ(); ++k)
  for (int j = 0; j < im1.dimY(); ++j)
  for (int i = 0; i < im1.dimX(); ++i)
  {
    im2(i,j,k) = ampl*sin(2*M_PI*i/step);
  }
  til::aimswrite(im2, name_out + ".y");
  for (int k = 0; k < im1.dimZ(); ++k)
  for (int j = 0; j < im1.dimY(); ++j)
  for (int i = 0; i < im1.dimX(); ++i)
  {
    im2(i,j,k) = 0;
  }
  til::aimswrite(im2, name_out + ".z");
  std::cout << "OK" << std::endl;
}


void testPower()
{
  //std::cout << til::functor::Square<double>(0) << " " << til::functor::Square<double>(1) << " " << til::functor::Square<double>(2) << " " << til::functor::Square<double>(3) << std::endl;
  
  try
  {
    til::functor::Pow<double, 0>()(0);
  }
  catch (til::functor::Pow<double, 0>::InvalidArgument &)
  {
    std::cout << "Caught OK" << std::endl;
  }
  
  std::cout <<                                             til::functor::Pow<double, 0>()(1) << " " << til::functor::Pow<double, 0>()(2) << " " << til::functor::Pow<double, 0>()(3) << std::endl;
  std::cout << til::functor::Pow<double, 1>()(0) << " " << til::functor::Pow<double, 1>()(1) << " " << til::functor::Pow<double, 1>()(2) << " " << til::functor::Pow<double, 1>()(3) << std::endl;
  std::cout << til::functor::Pow<double, 2>()(0) << " " << til::functor::Pow<double, 2>()(1) << " " << til::functor::Pow<double, 2>()(2) << " " << til::functor::Pow<double, 2>()(3) << std::endl;
  std::cout << til::functor::Pow<double, 3>()(0) << " " << til::functor::Pow<double, 3>()(1) << " " << til::functor::Pow<double, 3>()(2) << " " << til::functor::Pow<double, 3>()(3) << std::endl;
  std::cout << til::functor::Pow<double, 4>()(0) << " " << til::functor::Pow<double, 4>()(1) << " " << til::functor::Pow<double, 4>()(2) << " " << til::functor::Pow<double, 4>()(3) << std::endl;
  std::cout << til::functor::Pow<double, 5>()(0) << " " << til::functor::Pow<double, 5>()(1) << " " << til::functor::Pow<double, 5>()(2) << " " << til::functor::Pow<double, 5>()(3) << std::endl;
  std::cout << til::functor::Pow<double, 6>()(0) << " " << til::functor::Pow<double, 6>()(1) << " " << til::functor::Pow<double, 6>()(2) << " " << til::functor::Pow<double, 6>()(3) << std::endl;
  std::cout << til::functor::Pow<double, 13>()(0) << " " << til::functor::Pow<double, 13>()(1) << " " << til::functor::Pow<double, 13>()(2) << " " << til::functor::Pow<double, 13>()(3) << std::endl;
}

void testConvert2()
{
  til::numeric_array<float, 3> a;
  til::numeric_array<double, 3> b;
  a[0] = 1.29;
  a[1] = 9.75e7;
  a[2] = -0.0007256;
  std::cout << a << std::endl;
  convert2(a).into(b);
  std::cout << b << std::endl;
}

void testSphericalCartesian()
{
  std::cout << "testing spherical coords..." << std::flush;
  for (int i = 0; i < 1000; ++i)
  {
    double x = til::rand(-1.0, 1.0);
    double y = til::rand(-1.0, 1.0);
    double z = til::rand(-1.0, 1.0);

    double x2, y2, z2;
    {    
      double theta, phi, rho;
      til::geo::cartesian2spherical(x, y, z, theta, phi, rho);      
      til::geo::spherical2cartesian(theta, phi, rho, x2, y2, z2);
    }
    if (std::fabs(x-x2) > 1e-6 ||
        std::fabs(y-y2) > 1e-6 ||
        std::fabs(z-z2) > 1e-6)
    {
      std::cout << "error!" << std::endl;
    }
    else std::cout << "." << std::flush;
  }
  std::cout << "OK" << std::endl;
}

void testNNeigh()
{
  typedef til::Mesh_N MyMesh;
  MyMesh mesh;
  {
    Reader<AimsSurfaceTriangle> r("/home/cathier/volatile/data/twins/nih1412T_Lwhite.mesh");
    AimsSurfaceTriangle s;
    til::read_mesh(r, mesh);
  }
  typedef std::vector<std::vector<std::size_t> > Neighborhoods;
  shared_ptr<Neighborhoods> pneigh = circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  Neighborhoods neigh_n;
  til::get_n_neighborhood(*pneigh, neigh_n, 1);
  if (neigh_n.size() != pneigh->size())
    std::cout << "failed " << neigh_n.size() << " " << pneigh->size() << std::endl;
  else std::cout << "OK " << neigh_n.size() << " " << pneigh->size() << std::endl;
  for (std::size_t i = 0; i < neigh_n.size(); ++i)
  {
    if (neigh_n[i].size() != (*pneigh)[i].size())
      std::cout << "failed " << i << " " << neigh_n[i].size() << " " << (*pneigh)[i].size() << std::endl;
    else
      std::cout << "OK " << i << " " << neigh_n[i].size() << " " << (*pneigh)[i].size() << std::endl;
  }    
}

void testassert()
{
  std::cout << "true" << std::endl;
  assert(true);
  std::cout << "false" << std::endl;
  assert(false);
}

void testsort()
{
  for (int i = 0; i < 10000; ++i)
  {
    double a = til::rand(-1.0, 1.0);
    double b = til::rand(-1.0, 1.0);
    double c = til::rand(-1.0, 1.0);
    til::sort(a, b, c);
    if (a < b || b < c)
    {
      std::cout << "error" << std::flush;
    }
    else
    {
      std::cout << "." << std::flush;
    }
  }
}

void testnface()
{
  typedef til::Mesh_N MyMesh;
  MyMesh mesh;
  {
    Reader<AimsSurfaceTriangle> r("/home/cathier/volatile/data/twins/nih1412T_Lwhite.mesh");
    AimsSurfaceTriangle s;
    til::read_mesh(r, mesh);
  }
  typedef std::vector<std::vector<std::size_t> > Neighborhoods;
  shared_ptr<Neighborhoods> pneigh = circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  
  std::vector<std::vector<std::size_t> > nfaces;
  til::neighboring_faces(*pneigh, getFaceIndices(mesh), nfaces);
  
  for (std::size_t i = 0; i < pneigh->size(); ++i)
  {
    std::cout << i << "\n";
    for (std::size_t j = 0; j < (*pneigh)[i].size(); ++j)
    {
      std::cout << (*pneigh)[i][j] << " ";
    }
    std::cout << "\n";
    for (std::size_t j = 0; j < nfaces[i].size(); ++j)
    {
      std::cout << getFaceIndices(mesh)[nfaces[i][j]][0] << " " << getFaceIndices(mesh)[nfaces[i][j]][1] << " " << getFaceIndices(mesh)[nfaces[i][j]][2] << "\n";
    }
    std::cout << std::endl;
  }
}

int main( int argc, char* argv[] )
{
  try
  {
#ifdef __TIME__
    std::cout << "(Compiled at " << __TIME__ << " on " << __DATE__ << ")" << std::endl;
#endif
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
    
    //testBundles(argc, argv);
    //regBundles(argc, argv);
    
    /*
    testPowell();
    testCG();
    testCG2();
    testGD();
    */
    
    //testsqrt();
    //testApplyAffineTransform();
    //testMatrixDet();
    //geoinflate(argc, argv);
    //testSimpleRegistration();
    //testcurv();
    //testGraphConvertion();
    //testKDTree0();
    //totodebug(argc, argv);
    //testSimplification(argc, argv);
    //testSimplification2(argc, argv);
    //testQuantization(argc, argv);
    //testQuantization2(argc, argv);
    //powellWavelet();
    //writeEmptyFloatImage();
    //testUnsignedDiff();
    //testdwt();
    //testdwt2D();
    //dwtImage(argc, argv);
    //testQuant2(argc, argv);
    //testQuantization3(argc, argv);
    removingPointsWithHighCurvature(argc, argv);
    //clearBorders(argc, argv);
    //testGonzCluster(argc, argv);
    //testSelfIntersection(argc, argv);
    //testTriangleIntersection(argc, argv);
    //testCurvatureSmoothing(argc, argv);
    //testCurvatureSmoothing2(argc, argv);
    //testGaussPatchApprox(argc, argv);
    //testGaussPatchApprox2(argc, argv);
    //testLambdaMu(argc, argv);
    //testMultiscaleRegistration(argc, argv);
    //removingPointsWithHighCurvature2(argc, argv);
    //checkBundleTransfo(argc, argv);
    //testWavConjugate();
    //generateStupidImages();
    //meshInterior(argc, argv);
    //testSameSign();
    //defimage(argc, argv);
    //testPower();
    //testConvert2();
    //testSphericalCartesian();
    //testNNeigh();
    //testassert();
    //testsort();
//     testnface();
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
	catch( carto::user_interruption & )
	{
		// Exceptions thrown by command line parser (already handled, simply exit)
	}
	catch( exception & e )
	{
		cerr << e.what() << endl;
	}
}

