
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
//#include "containers.h"
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

#include "totodebug.h"

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


template < typename TMesh, typename TPrec, typename TConn, typename Finder >
class FiberDistance : public std::unary_function< til::Vector<float,12>, float >
{
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
    std::cout << "bbst2 " << getVertices(mesh2)[0] << " " << getVertices(mesh1)[0] << std::endl;
    std::cout << "bbst " << getVertices(m_mesh2)[0] << " " << getVertices(m_mesh1)[0] << std::endl;
   }

  float operator()(const til::Vector<float,12> & params)
  {
    til::AffineMap<float> a;
    a.transfo().setMatrix(
      til::Matrix3<float>(
        params[0],params[1],params[2],
        params[3],params[4],params[5],
        params[6],params[7],params[8]
        ));
    a.transfo().setTransl(til::Vector<float,3>(params[9],params[10],params[11]));
    double resdist = 0.0;
    double resconn = 0.0;

    //for (const_vertex_iterator iVertex2 = getVertices(mesh2).begin(); iVertex2 != getVertices(mesh2).end(); ++iVertex2)
    for (std::size_t i2 = 0; i2 < til::size(getVertices(m_mesh2)); ++i2)
    {
      std::size_t i1 = m_finder(a(getVertices(m_mesh2)[i2]));
      resdist += til::dist2<double>(getVertices(m_mesh2)[i2], getVertices(m_mesh1)[i1]);
      m_nn[i2] = i1;
      if (i2 == 0)
      {
        std::cout << "@ " << getVertices(m_mesh2)[i2] << " " << a(getVertices(m_mesh2)[i2]) << " " << getVertices(m_mesh1)[i1] << " " << getVertices(m_mesh1)[0] << std::endl;
      }
    }
    resdist /= std::sqrt(til::det(a.transfo().getMatrix())) + 128 * std::numeric_limits<double>::epsilon();

    for (std::size_t i2 = 0; i2 < til::size(getVertices(m_mesh2)); ++i2)
    {
      for (til::sparse_vector<double>::Map::const_iterator c2 = m_conn2[i2].getMap().begin();
        c2 != m_conn2[i2].getMap().end(); ++c2)
      {
        resconn += c2->second * m_conn1[m_nn[i2]][m_nn[c2->first]];
      }
    }

    double res = -m_alpha*resconn + (1-m_alpha)*resdist;

    std::cout << "Func called at " << params << " : " << resdist << " " << -resconn << " " << res << std::endl;
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


void testBundles(int argc, char * argv[])
{
  typedef til::Mesh_N MyMesh;
  MyMesh mesh;
  double alpha = 0.5;
  {
    Reader<AimsSurfaceTriangle> r;
    AimsApplication app( argc, aims_const_hack(argv), "testBundles" );
    app.addOption(r, "-i", "input mesh" );
    app.addOption(alpha, "-alpha", "alpha");
    app.initialize();

    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }

  std::cout << "First : " << getVertices(mesh)[0] << std::endl;
  
  const double DIST = 5.0;
  std::cout << "Computing geomap" << std::endl;
  til::ghost::GMapStop_AboveThreshold<double> stopGhost(DIST);
  til::TriangleMeshGeodesicMap<MyMesh, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage<til::sparse_vector, MyMesh, double > > geomap(mesh, stopGhost);
  std::vector<std::size_t> startPoints(1);
  std::vector<double> dist(1, 0.0);
  typedef std::vector< std::pair<std::size_t, double> > QuickMap;
  std::vector< QuickMap > res(til::size(getVertices(mesh)));
  //til::MeshTraits<til::Mesh2_N>::VertexCollection::const_iterator iVertex = getVertices(mesh).begin();
  //for(; iVertex != getVertices(mesh).end(); ++iVertex)
  std::cout << "wow" << std::endl;
  {
    til::MapHistogram<std::size_t> h;
    std::vector<std::size_t> nneigh(til::size(getVertices(mesh)));
    for (std::size_t i = 0; i < til::size(getVertices(mesh)); ++i)
    {
      //startPoints[0] = &(getVertices(mesh)[i]);
      startPoints[0] = i;
      geomap.init(startPoints, dist);
      geomap.compute();
      //std::cout << count << std::endl;
      shared_ptr<til::sparse_vector<double> > tmp = geomap.distanceMap();
      res[i].resize(til::size(tmp->getMap()));
      {
        using namespace til::expr;
        til::detail::loop_xx(castTo(*_1, *_2), res[i], tmp->getMap());
      }
      //til::loop(res[i], tmp.getMap(), Convert());
      //til::convert(res[i], tmp.getMap());
      //std::cout << count << std::endl;
      h.accumulate(til::size(res[i]));
      nneigh[i] = til::size(res[i]);
      //res[i] = geomap.distanceMap();
      //sres.setDefaultValue(0);
      //std::transform(sres.begin(), sres.end(), res.begin(), res.begin(), std::plus<double>());
    }
    std::cout << "Finishing" << std::endl;

    // writing number of neighbors in the DIST-cell
    {
      Texture1d t(1, til::size(nneigh));
      til::convert(t, nneigh);
      Writer<Texture1d> w("/volatile/guevara/tmp/nneigh");
      w.write(t);
    }
  }

  // writing distances in the DIST-cell of the first vertex
  {
    Texture1d t(1, til::size(getVertices(mesh)));
    for (std::size_t i = 0; i < til::size(res[0]); ++i)
    {
      t.item(res[0][i].first) = res[0][i].second;
    }
    Writer<Texture1d> w("/volatile/guevara/tmp/toubou");
    w.write(t);
  }

  // loading fibers
  boost::shared_ptr<std::vector<std::vector<til::Point<float,3> > > > pfibers;
  {
    aims::BundleReader r("/home/Panabase/pascal/fibers/res/old/res.bundles");
    til::BundleLoader loader;
    r.addBundleListener(loader);
    r.read();
    pfibers = loader.getFibers();
  }
  typedef std::vector<til::Point<float,3> > Fiber;
  typedef std::vector<Fiber> Fibers;
  Fibers & fibers = *pfibers;

  std::cout << "Number of fibers: " << til::size(fibers) << std::endl;
  /*
  for (std::size_t i = 0; i < til::size(fibers[0]); ++i)
  {
    std::cout << fibers[0][i] << std::endl;
  }
  */

  std::cout << "Generating kdtree" << std::endl;
  typedef til::KDTree<std::size_t, til::MeshTraits<MyMesh>::VertexCollection> MyKDTree;
  MyKDTree kdt(getVertices(mesh));
  makeKDTree(getVertices(mesh), kdt);

  std::cout << "Looking for closest points" << std::endl;
  typedef til::sparse_vector<double> Connectivity;
  typedef std::vector< Connectivity > Connectivities;
  Connectivities conn(til::size(getVertices(mesh)), til::sparse_vector<double>(til::size(getVertices(mesh))));
  for (Fibers::const_iterator iFiber = fibers.begin(); iFiber != fibers.end(); ++iFiber)
  {
    til::FindClosest< double, MyKDTree > fc(kdt);
    std::size_t A = fc(iFiber->front());
    std::size_t B = fc(iFiber->back());
    /*
    for (QuickMap::const_iterator i = res[A].begin(); i != res[A].end(); ++i)
    for (QuickMap::const_iterator j = res[B].begin(); j != res[B].end(); ++j)
    {
      double e1 = til::square(i->second);
      double e2 = til::square(j->second);
      double w = std::exp(- ( e1 + e2 ) / ( 2*DIST*DIST ));
      conn[i->first][j->first] += w;
      conn[j->first][i->first] += w;
    }
    */
    for (QuickMap::const_iterator j = res[B].begin(); j != res[B].end(); ++j)
    {
      double e2 = til::square(j->second);
      double w = std::exp( -e2  / ( 2*DIST*DIST ));
      conn[A][j->first] += w;
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
    Texture1d t(1, til::size(getVertices(mesh)));
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

  // Create another brain
  MyMesh mesh2(mesh);
  {
    til::AffineMap<float> a;
    a.transfo().setMatrix(til::Matrix3<float>(1.1, -0.1, 0.1, 0.1, 0.9, 0.1, -0.1, 0, 1));
    a.transfo().setTransl(til::Vector<float,3>(10, 10, 10));
    //a.transfo().setMatrix(til::Matrix3<float>(1,0,0,0,1,0,0,0,1));
    //a.transfo().setTransl(til::Vector<float,3>(0,0,0));
    using namespace til::expr;
    mesh2.vertices() = 
    boost::shared_ptr<til::MeshTraits<til::Mesh_N>::VertexCollection>(new til::MeshTraits<til::Mesh_N>::VertexCollection(size(getVertices(mesh))));
    //mesh2.vertices() = shared_ptr(new std::vector<til::Point<float,3> >(size(getVertices(mesh))));
    //til::detail::loop_xx(*_2 = bind(a,*_1), getVertices(mesh), getVertices(mesh2));
    til::detail::loop_xx(*_2 = func(a)(*_1), getVertices(mesh), getVertices(mesh2));
  }
  std::cout << "OK" << std::endl;
  std::cout << getVertices(mesh)[0] << " " << getVertices(mesh2)[0] << std::endl;

  {
    AimsSurfaceTriangle s;
    til::convert(s, mesh2);
    Writer<AimsSurfaceTriangle> w("/volatile/guevara/tmp/reginit");
    w.write(s);
  } 

  // Registration between brains
  {
    typedef til::FindClosest< double, MyKDTree > Finder;
    Finder finder(kdt);
    typedef FiberDistance<MyMesh, double, Connectivities, Finder > Energy;
    Energy dist(mesh, mesh2, conn, conn, alpha, finder);
    til::Powell<Energy> minalgo(dist);
    // initial scaling estimates of the parameters
    {
      std::vector<float> initStd(12, 1.0f);
      initStd[9] = initStd[10] = initStd[11] = 50;
      minalgo.initStd() = initStd;
    }
    
    // Initial starting point, corresponding to the identity transfom.
    til::Vector<float,12> params;
    params[0] = 1;
    params[4] = 1;
    params[8] = 1;    
    params = minalgo(params);
    std::cout << "Minimization finished: " << params << std::endl;
    std::cout << "First : " << getVertices(mesh2)[0] << std::endl;

    {
      using namespace til::expr;
      til::AffineMap<float> a;
      a.transfo().setMatrix(til::Matrix3<float>(params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8]));
      a.transfo().setTransl(til::Vector<float,3>(params[9],params[10],params[11]));
      for (std::size_t i = 0; i < size(getVertices(mesh2)); ++i)
      {
        getVertices(mesh2)[i] = a(getVertices(mesh2)[i]);
        //a(getVertices(mesh2)[i]);
      }
      //til::detail::loop_x(*_1 = func(a)(*_1), getVertices(mesh2));
    }
    std::cout << "First : " << getVertices(mesh2)[0] << std::endl;
    {
      AimsSurfaceTriangle s;
      til::convert(s, mesh2);
      Writer<AimsSurfaceTriangle> w("/volatile/guevara/tmp/regres");
      w.write(s);
    } 
  }  
}

int main( int argc, char* argv[] )
{
  try
  {
    std::cout << "main starts now!" << std::endl;
    testBundles(argc, argv);
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

