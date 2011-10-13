
// includes from STL
//#include <math.h>
#include <cmath>
//#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include <sys/times.h>
#include <sys/time.h>
#include <vector>

// includes from BOOST
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
//#include "aims/curve/bundles.h"
#include "aims/getopt/getopt2.h"
#include "aims/io/reader.h"
#include "aims/io/writer.h"
//#include "aims/mesh/inflate.h"
#include "aims/mesh/surface.h"
#include "aims/mesh/surfacegen.h"
#include "aims/mesh/texture.h"
//#include "aims/resampling/linearInterpolator.h"
//#include "aims/roi/maskIterator.h"
//#include "aims/utility/converter_volume.h"
#include "aims/vector/vector.h"

// includes from TIL
#include "til/til_common.h"
#include "til/AffineMap.h"
#include "til/TExpr.h"
#include "til/Vector3.h"

// includes from TIL
#include "aims_wrap.h"
#include "BinaryTree.h"
#if GCC_VERSION < 030404
//#include "connectivityMatrix.h"
#endif
#include "containers.h"
#include "domino.h"
#include "Forces.h"
#include "func_iterator.h"
#include "globalTraits.h"
#include "histogram.h"
#include "kdtree.h"
#include "MeshTraits.h"
#include "meshUtils.h"
#include "minTools.h"
#include "miscUtils.h"
#include "triangle_mesh_geodesic_map.h"
#include "SparseVector.h"

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


template < typename TMesh >
std::vector<til::sparse_vector<double> >
//std::vector<std::vector<double> >
getDistToNeighbors(const TMesh & mesh, double distance)
{
  std::cout << "Computing geomap (1)..." << std::flush;
  std::size_t n = size(getVertices(mesh));
  til::ghost::GMapStop_AboveThreshold<double> stopGhost(distance);
  //til::TriangleMeshGeodesicMap<TMesh, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage<til::sparse_vector, TMesh, double > > 
  til::GraphDistanceMap<TMesh, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage<til::sparse_vector, TMesh, double > > 
  //til::GraphDistanceMap<TMesh, double, til::ghost::GMapStop_AboveThreshold<double> > 
    geomap(mesh, stopGhost);
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
    geomap.compute();
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
  //til::TriangleMeshGeodesicMap<TMesh, double, Plugin, til::policy::GMap_DefaultStorage<til::sparse_vector, TMesh, double > > 
  til::GraphDistanceMap<TMesh, double, Plugin, til::policy::GMap_DefaultStorage<til::sparse_vector, TMesh, double > > 
    geomap(mesh);
  std::vector<std::size_t> startPoints(1);
  std::vector<double> dist(1, 0.0);
  //std::vector<double> res(n, 0.0);
  std::vector<shared_ptr<til::sparse_vector<double> > > res(n);
  for (std::size_t i=0; i < n; ++i)
  {
    startPoints[0] = i;
    geomap.init(startPoints, dist);
    geomap.stopGhost().init(nindex[i].begin(), nindex[i].end());
    geomap.compute();
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
  //til::TriangleMeshGeodesicMap<TMesh, double, Plugin, til::policy::GMap_DefaultStorage<til::sparse_vector, TMesh, double > > 
  til::GraphDistanceMap<TMesh, double, Plugin, til::policy::GMap_DefaultStorage<til::sparse_vector, TMesh, double > > 
    geomap(mesh, neighc);
  std::vector<std::size_t> startPoints(1);
  std::vector<double> dist(1, 0.0);
  //std::vector<double> res(n, 0.0);
  shared_ptr<til::sparse_vector<double> > res(n);
  //for (std::size_t i=0; i < n; ++i)
  {
    startPoints[0] = i;
    geomap.init(startPoints, dist);
    geomap.stopGhost().init(nindex[i].begin(), nindex[i].end());
    geomap.compute();
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
  //til::TriangleMeshGeodesicMap<TMesh, double, Plugin, til::policy::GMap_DefaultStorage<til::sparse_vector, TMesh, double > > 
  //til::GraphDistanceMap<TMesh, double, Plugin, til::policy::GMap_DefaultStorage<til::sparse_vector, TMesh, double > > 
  til::GraphDistanceMap<TMesh, double, Plugin  > 
    geomap(mesh, neighc);
  std::vector<std::size_t> startPoints(1);
  std::vector<double> dist(1, 0.0);
  //std::vector<double> res(n, 0.0);
  shared_ptr<std::vector<double> > res;
  //for (std::size_t i=0; i < n; ++i)
  {
    startPoints[0] = i;
    geomap.init(startPoints, dist);
    geomap.stopGhost().init(nindex[i].begin(), nindex[i].end());
    geomap.compute();
    res = geomap.distanceMap();
  }
  //std::cout << "OK" << std::endl;
  return res;
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
  shared_ptr<std::vector<std::vector<std::size_t> > > neighc = circular_neighborhoods(mesh);
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
  std::vector<til::Vector<float,3> > fcenter(n);
  std::vector<til::Vector<float,3> > fdist(n);
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
        Writer<AimsSurfaceTriangle> w(s);
        w.write(tmp); 
        
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
      std::vector<til::Vector<float,3> >::iterator iC = fcenter.begin();
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
    til::Vector<float,3> sumcf(0,0,0);
    std::cout << "Applying forces..." << std::flush;
    for (std::size_t i = 0; i < size(getVertices(mesh)); ++i)
    {
      sumcf += fcenter[i];
      tmp1 += til::norm<double>(fcenter[i]);
      tmp2 += til::norm<double>(fdist[i]);
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
    Writer<AimsSurfaceTriangle> w("/home/cathier/tmp/finalevel");
    w.write(s);
  } 
  {
    Writer<AimsSurfaceTriangle> w("/home/cathier/tmp/mevel");
    w.write(wmesh);
  }
}

int main( int argc, char* argv[] )
{
  try
  {
    std::cout << "main starts now!" << std::endl;
    geoinflate(argc, argv);
  	return 0;
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

