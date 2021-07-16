// includes from AIMS
#include "aims/getopt/getopt2.h"
#include "aims/io/reader.h"
#include "aims/io/writer.h"

// includes from TIL
#include <cathier/aims_wrap.h>
#include "cathier/mesh_decimation.h"

using namespace aims;
using namespace std;


void removingPointsWithHighCurvature(int argc, char * argv[])
{
  til::Mesh_N mesh;
  std::cout << "Reading mesh..." << std::flush;
  Writer<AimsSurfaceTriangle> w;
  cerr << "WARNING THIS COMMAND IS OBSOLETE. Please use AimsMeshCleaner instead." << endl;
  float maxCurv;
  {
    bool do_work = false;
    Reader<AimsSurfaceTriangle> r;
    AimsApplication app( argc, aims_const_hack(argv), "removingPointsWithHighCurvature");
    app.addOption(r, "-i", "input mesh" );
    app.addOption(w, "-o", "output mesh" );
    app.addOption(do_work, "--work", "really do something instead of just printing a deprecation warning." );
    app.addOption(maxCurv, "-maxCurv", "maximum unoriented Gaussian curvature allowed in the mesh");
    app.initialize();

    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);
  }
  std::size_t n =  getVertices(mesh).size();

  std::cout << "Get circular neighbor..." << std::flush;
  boost::shared_ptr<std::vector<std::vector<std::size_t> > > neighc = circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
  std::cout << "OK" << std::endl;

  std::cout << "Converting into graph..." << std::flush;
  typedef til::MeshVertexNodeX<std::size_t>  VertexNode;
  typedef til::MeshFaceNodeX<VertexNode>     FaceNode;
  std::list<VertexNode> graph_vertices;
  std::list<FaceNode> graph_faces;
  til::list2graph_mesh_conversion(getVertices(mesh), getFaceIndices(mesh), *neighc, graph_vertices, graph_faces);
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

  std::vector<float> unorientedGaussianCurvatures(n);
  typedef til::xsr::graph::Position<VertexNode::VertexIndex> A;
  typedef til::xsr::graph::Neighbors<VertexNode::VertexIndex> B;
  A a; B b; til::MeshCurvature2<A,B,float> mc(a,b);
  // NB: this hack doesn't work in gcc3.3 :(
  //til::MeshCurvature2<A,B,float> mc((A()),(B()));
  unsigned int count = 0;
  {
    til::Vertex_remover<VertexNode, FaceNode> vertexRemover(graph_vertices, graph_faces);
    bool flagFinished;
    unsigned int iter = 0;
    do
    {
      flagFinished = true;
      std::cout << "pass " << ++iter << "... " << std::flush;
      for (std::list<VertexNode>::iterator i = graph_vertices.begin(); i != graph_vertices.end(); )
      {
        mc.process(i);
        if (mc.unorientedGaussianCurvature() > maxCurv)
        {
          if (vertexRemover(i))
          {
            flagFinished = false;
            ++count;
            continue;
          }
        }
        ++i;
      }
      std::cout << " Removed " << count << " vertices" << std::endl;
    } while (!flagFinished);
  }
  std::cout << "Remove " << count << " points" << std::endl;

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


int main( int argc, char* argv[] )
{
  try
  {
#ifdef __TIME__
    std::cout << "(Compiled at " << __TIME__ << " on " << __DATE__ << ")" << std::endl;
#endif

    removingPointsWithHighCurvature( argc, argv );
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

