
#include <list>
#include <vector>

#include "aims/io/reader.h"
#include "aims/mesh/surface.h"

#include "til/til_common.h"

#include "aims_wrap.h"
#include "Mesh.h"
#include "meshUtils.h"

using namespace aims;
using namespace carto;

namespace
{
  struct Temp
  {
    double angle;
    std::size_t index;
    friend bool operator<(const Temp & x, const Temp & y) { return x.angle < y.angle; }
  };
  struct TempComp
  {
    bool operator()(const Temp & x, const Temp & y) { return x.angle < y.angle; }
  };


struct MeshVertexNode;

struct MeshFaceNode
{
  boost::array<std::list<MeshVertexNode>::iterator,3> face;
};

struct MeshVertexNode
{
  til::Point<float,3> pos;
  std::list<std::list<MeshVertexNode>::iterator> neighbors;
  std::list<std::list<MeshFaceNode>::iterator> faces;
};

struct GroComp
{
  typedef std::list<MeshVertexNode>::iterator T;

  GroComp(T a, T b, T c)
    : m_a(a)
    , m_b(b)
    , m_c(c)
  {}
  
  boost::array<T,3>
  helper(T a, T b, T c)
  {
    boost::array<T,3> tmp;
    tmp[0] = a;
    tmp[1] = b;
    tmp[2] = c;
    return tmp;
  }
  
  bool operator()(const std::list<MeshFaceNode>::iterator & f)
  {
    return
      f->face == helper(m_a, m_b, m_c) ||
      f->face == helper(m_c, m_a, m_b) ||
      f->face == helper(m_b, m_c, m_a) ||
      f->face == helper(m_b, m_a, m_c) ||
      f->face == helper(m_c, m_b, m_a) ||
      f->face == helper(m_a, m_c, m_b)
      ;
  }
  
  T m_a;
  T m_b;
  T m_c;
};


struct ItComp
{
  bool operator()
  (
   const std::list<MeshVertexNode>::const_iterator & i,
   const std::list<MeshVertexNode>::const_iterator & j
  )
  {
    return &*i < &*j;
  }
};

std::list<boost::array<std::size_t,3> >
simple_delaunay_triangulation(std::vector<til::Point<float, 2> > points)
{
  typedef boost::array<std::size_t,3> Face;
  typedef std::list<Face> FaceList;
  
  // Get size now before points are added
  std::size_t n = points.size();
  
  // Don't waste time calling this for less than 4 points.
  assert(til::size(points) >= 4);
  
  FaceList faces;

  // compute bounding triangle
  {
    float miny, k1, k2;
    k1 = k2 = -std::numeric_limits<float>::max();
    miny = std::numeric_limits<float>::max();
    for (std::size_t i = 0; i < til::size(points); ++i)
    {
      miny = til::min(miny, points[i][1]);
      k1 = til::max(k1, points[i][1] - std::sqrt(3.0f)*points[i][0]);
      k2 = til::max(k2, points[i][1] + std::sqrt(3.0f)*points[i][0]);
    }
    
    til::Point<float,2> M(
      (k2 - k1) / (2*std::sqrt(3.0f)),
      2*miny/3 + (k1+k2)/6
      );
      
    til::Point<float,2> A((miny - k1) * (1/std::sqrt(3.0f)), miny);
    til::Point<float,2> B((k2 - miny) * (1/std::sqrt(3.0f)), miny);
    til::Point<float,2> C( (k2 - k1) / (2*std::sqrt(3.0f)), (k1+k2)/2);

    A += 2.0f*(A-M);
    B += 2.0f*(B-M);
    C += 2.0f*(C-M);
    
    points.push_back(A);
    points.push_back(B);
    points.push_back(C);
    
    boost::array<std::size_t,3> face = { points.size()-1, points.size()-2, points.size()-3 };
    faces.push_back(face);
  }

  //Face face = {0,1,2};
  //faces.push_back(face);
  //for (std::size_t i = 3; i < til::size(points); ++i)
  for (std::size_t i = 0; i < n; ++i)
  {
    std::cout << "another point :" << points[i] << std::endl;

    // Remove triangles if we lie in their circumcircle, and collect their vertex.
    std::set<std::size_t> neighbors;
    //e.insert(0);
    //e.insert(i-1);
    {
      FaceList::iterator f = faces.begin();
      while (f != faces.end())
      {
        std::cout << "face " << points[(*f)[0]] << " " <<  points[(*f)[1]] << " " <<  points[(*f)[2]] << ":";
        if (til::is_in_circumcircle<double>(points[i].data(), points[(*f)[0]].data(), points[(*f)[1]].data(), points[(*f)[2]].data()))
        {
          std::cout << "inside" << std::endl;
          neighbors.insert((*f)[0]);
          neighbors.insert((*f)[1]);
          neighbors.insert((*f)[2]);
          f = faces.erase(f);
        }
        else
        {
          std::cout << "outside circle" << std::endl;
          ++f;
        }
      }
    }
    
    assert(neighbors.size() >= 3);
    
    // order neighbors according to the angle they make with the current point
    
    std::vector<Temp> tmp;
    tmp.reserve(neighbors.size());
    for (std::set<std::size_t>::iterator iNeighbor = neighbors.begin(); iNeighbor != neighbors.end(); ++iNeighbor)
    {
      Temp t;
      t.index = *iNeighbor;
      til::Vector<float,2> v = points[*iNeighbor]-points[i];
      t.angle = til::geo::angle(v[0], v[1]);
      tmp.push_back(t);
    }
    //std::sort(tmp.begin(), tmp.end(), TempComp());
    std::sort(tmp.begin(), tmp.end());
    
    // Form new faces
    for (std::size_t it = 0; it < tmp.size(); ++it)
    {
      std::cout << "Adding face " << i << " " << tmp[it].index << " " << tmp[ (it+1) % til::size(tmp) ].index << std::endl;
      Face face = {i, tmp[it].index, tmp[ (it+1) % til::size(tmp) ].index };
      faces.push_back(face);
    }
  }

    {
      std::cout << "faces" << std::endl;
      for (FaceList::iterator f = faces.begin(); f != faces.end(); ++f)
      {
        std::cout << (*f)[0] << " " << (*f)[1] << " " << (*f)[2] << std::endl;
      }
      std::cout << "endfaces" << std::endl;
    }


  // Remove extra faces
  {
    FaceList::iterator f = faces.begin();
    while (f != faces.end())
    {
      // Remove faces with supertriangle vertex
      if ((*f)[0] >= n ||
          (*f)[1] >= n ||
          (*f)[2] >= n)
      {
        f = faces.erase(f);
      }
      else
      {
        ++f;
      }
    }
  }


  // Check if we have to resort to normal testing to remove faces
  if (faces.size() != n-2)
  {
    std::cout << "NI!" << faces.size() << "-" << n << std::endl;
    {
      std::cout << "faces2" << std::endl;
      for (FaceList::iterator f = faces.begin(); f != faces.end(); ++f)
      {
        std::cout << (*f)[0] << " " << (*f)[1] << " " << (*f)[2] << std::endl;
      }
      std::cout << "endfaces" << std::endl;
    }
    /*
    FaceList::iterator f = faces.begin();
    while (f != faces.end())
    {
      // Remove faces with supertriangle vertex
      if (0)
      {
        f = faces.erase(f);
      }
      else
      {
        ++f;
      }
    }
    */
  }
  
  exit(0);
  /*
  std::cout << "nfaces: " << til::size(faces) << std::endl;
  for (FaceList::const_iterator i = faces.begin(); i != faces.end(); ++i)
  {
    std::cout << (*i)[0] << " " << (*i)[1] << " " << (*i)[2] << std::endl;
  }
  */  
  return faces;
}

std::vector<til::Point<float,2> >
simple_neighborhood_flattening
(
 const til::Point<float,3>                                & point,
// til::Point<float,3>                                point,
 const std::vector<std::list<MeshVertexNode>::iterator>   & neighbors
)
{
  std::cout << "nneighbors " << til::size(neighbors) << std::endl;
  std::cout << "point " << point << std::endl;
  std::cout << "neighbors " << std::endl;
  for (std::size_t i = 0; i < neighbors.size(); ++i) std::cout << neighbors[i]->pos << std::endl;
  

  std::vector<double> angles;
  std::vector<double> norms;
  angles.reserve(til::size(neighbors));
  norms.reserve(til::size(neighbors));
  
  std::vector<std::list<MeshVertexNode>::iterator>::const_iterator n = neighbors.begin();
  til::const_cyclic_iterator<std::vector<std::list<MeshVertexNode>::iterator> > n2(neighbors, ++neighbors.begin());
  //std::cout << "point: " << point << std::endl;
  for (; n != neighbors.end(); ++n, ++n2)
  {
    std::cout << "a " << point << std::endl;
    std::cout << (*n2)->pos << " " << (*n)->pos << std::endl;
    std::cout << "b " <<  point << std::endl;
    std::cout << (*n2)->pos-point<< " " << (*n)->pos-point << std::endl;
    std::cout << "c " <<  point << std::endl;
    norms.push_back(til::norm<double>((*n)->pos-point));
    std::cout << "d " <<  point << std::endl;
    double tmp4 = til::norm<double>((*n2)->pos-point);
    std::cout << "d1 " <<  point << std::endl;
    double tmp3 = til::norm<double>((*n)->pos-point);
    std::cout << "d2 " <<  point << std::endl;
    double tmp1 = tmp4 * tmp3;
    std::cout << "d3 " <<  point << std::endl;
    double tmp2 = til::dot((*n2)->pos - point, (*n)->pos - point);
    std::cout << "d4 " <<  point << std::endl;
    double tmp5 = std::acos( tmp2 / tmp1);
    std::cout << "d5 " <<  point << std::endl;
    //angles.push_back(std::acos(til::dot((*n2)->pos-point, (*n)->pos-point) / (til::norm<double>((*n2)->pos-point) * til::norm<double>((*n)->pos-point))));

    angles.push_back(tmp5);
    std::cout << angles.size() << " " << angles.capacity() << std::endl;
    std::cout << "e " <<  point << std::endl;
  }

  
  std::cout << "norms: ";
  std::copy(norms.begin(), norms.end(), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;

  std::cout << "angles: ";
  std::copy(angles.begin(), angles.end(), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;
  
  
  std::partial_sum(angles.begin(), angles.end(), angles.begin());

  
  std::cout << "angles summed: ";
  std::copy(angles.begin(), angles.end(), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;
  
  
  //double totalAngle = std::accumulate(angles.begin(), angles.end(), 0.0);
  std::transform(angles.begin(), angles.end(), angles.begin(), std::bind2nd(std::multiplies<double>(), 2*M_PI/angles.back()));

  
  std::cout << "angles normalized: ";
  std::copy(angles.begin(), angles.end(), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;
  
  
  std::vector<til::Point<float, 2> > res(til::size(neighbors));
  for (std::size_t i = 0; i < size(res); ++i)
  {
    res[i][0] = std::cos(angles[i]) * norms[i];
    res[i][1] = std::sin(angles[i]) * norms[i];
    //res[i][0] = std::cos(angles[i]);
    //res[i][1] = std::sin(angles[i]);
  }
  
  std::cout << "newpoints" << std::endl;
  for (std::size_t i = 0; i < res.size(); ++i)
  {
    std::cout << res[i] << std::endl;
  }
  std::cout << "end newpoints" << std::endl;
  
  
  return res;
}

std::list<MeshVertexNode>::iterator
remove_vertex(std::list<MeshVertexNode>::iterator i, std::list<MeshVertexNode> & graph_vertices, std::list<MeshFaceNode> & graph_faces)
{
  // remove faces point belongs to
  for (std::list<std::list<MeshFaceNode>::iterator>::iterator j = i->faces.begin(); j != i->faces.end(); ++j)
  {
    // Remove index to this face
    for (std::size_t k = 0; k < 3; ++k)
    {
      // don't remove face from facelist of point i itself -- because j is iterating on this very list,
      // plus, it will be deleted with i anyway.
      if ((*j)->face[k] == i) continue;
      //std::size_t tmp = til::size((*j)->face[k]->faces);
      (*j)->face[k]->faces.remove(*j);
      //assert(til::size((*j)->face[k]->faces) == tmp-1);
    }
    // remove face itself
    
    graph_faces.erase(*j);
    //assert(til::size(graph_faces) == tmp-1);
  }
  // remove point in neighbor's list
  for (std::list<std::list<MeshVertexNode>::iterator>::iterator j = i->neighbors.begin(); j != i->neighbors.end(); ++j)
  {
    (*j)->neighbors.remove(i);
    //assert(til::size((*j)->neighbors) == tmp-1);
  }
  // remove point
  return graph_vertices.erase(i);
}
}

void totodebug(int argc, char * argv[])
{
  til::Mesh_N mesh;
  
  std::cout << "Reading mesh..." << std::flush;
  {
    Reader<AimsSurfaceTriangle> r("/home/cathier/data/icbm100T_Lwhite.mesh");
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
  shared_ptr<std::vector<std::vector<std::size_t> > > neighc = getCircularNeighborIndices(mesh);
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
  
  std::cout << "Converting into graph..." << std::flush;
  std::list<MeshVertexNode> graph_vertices;
  std::list<MeshFaceNode> graph_faces;
  {
    std::vector<std::list<MeshVertexNode>::iterator> index_vertex(til::size(getVertices(mesh)));
    std::vector<std::list<MeshFaceNode>::iterator> index_face(til::size(getFaceIndices(mesh)));
    {
      // First pass: push all points and get an index2iterator translation
      for (std::size_t i = 0; i < n; ++i)
      {
        MeshVertexNode m;
        m.pos = getVertices(mesh)[i];
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
          index_vertex[i]->faces.push_back(index_face[*j]);
        }
        for (std::size_t j = 0; j < til::size((*neighc)[i]); ++j)
        {
          index_vertex[i]->neighbors.push_back(index_vertex[(*neighc)[i][j]]);
        }
      }
    }
  }
  std::cout << "OK" << std::endl;
  

  std::cout << "OK" << std::endl;
  {
    std::cout << "Removing 4-connected vertices..." << std::flush;
    int count = 0;
    bool done;
    do
    {
      done = true;
      std::list<MeshVertexNode>::iterator i = graph_vertices.begin();
      std::list<MeshVertexNode>::iterator itmp;
      int ind = 0;
      while (i != graph_vertices.end())
      {
        if (til::size(i->neighbors) != til::size(i->faces))
        {
          std::cout << "wF!" << til::size(i->neighbors) << "-" << til::size(i->faces) << "-";
        }
        
        if (til::size(i->neighbors) <= 4)
        {
          std::cout << "." << std::flush;
          done = false;
          ++count;

          // save list of neighbors
          std::vector<std::list<MeshVertexNode>::iterator> neighbors(til::size(i->neighbors));
          std::copy(i->neighbors.begin(), i->neighbors.end(), neighbors.begin());
          
          // remove point
          itmp = remove_vertex(i, graph_vertices, graph_faces);
          
          // Add new faces
          {
            // get delaunay faces
            std::list<boost::array<std::size_t,3> > tri;
            if (til::size(neighbors) >= 4)
            {
              tri = simple_delaunay_triangulation(simple_neighborhood_flattening(i->pos, neighbors));
            }
            else if (til::size(neighbors) == 3)
            {
              boost::array<std::size_t,3> tmp = {0,1,2};
              tri.push_back(tmp);
            }
            else
            {
              std::cout << "w2!";
            }
            
            if (til::size(tri) != til::size(neighbors) - 2)
            {
              std::cout << "wK!" << til::size(tri) << "-" << til::size(neighbors)-2;
            }
            
            // adding faces
            std::list<boost::array<std::size_t,3> >::iterator newf = tri.begin();
            for (; newf != tri.end(); ++newf)
            {
              // convert into graph faces
              MeshFaceNode f;
              for (std::size_t n = 0; n < 3; ++n)
              {
                f.face[n] = neighbors[(*newf)[n]];
              }
              // add face to graph
              std::list<MeshFaceNode>::iterator gf = graph_faces.insert(graph_faces.end(), f);
              // add face index to face points
              for (std::size_t n = 0; n < 3; ++n)
              {
                GroComp grocomp(f.face[0], f.face[1], f.face[2]);
                std::list<std::list<MeshFaceNode>::iterator>::iterator p = find_if(f.face[n]->faces.begin(), f.face[n]->faces.end(), grocomp);
                if (p != f.face[n]->faces.end())
                {
                  std::cout << "wu!" << &*f.face[0] << " " << &*f.face[1] << " " << &*f.face[2] << " ";
                  std::cout << "wu!" << &*(*p)->face[0] << " " << &*(*p)->face[1] << " " << &*(*p)->face[2] << std::endl;
                }
                f.face[n]->faces.push_back(gf);
              }
              // adding neighbors
              {
                for (std::size_t n = 0; n < 3; ++n)
                {
                  //int tmp = abs(int((*newf)[n]) - int((*newf)[(n+1)%3]));
                  //if (tmp != 1 && tmp != til::size(neighbors)-1)
                  {
                    if (find(f.face[n]->neighbors.begin(), f.face[n]->neighbors.end(), f.face[(n+1)%3]) == f.face[n]->neighbors.end())
                    {
                      if (find(f.face[(n+1)%3]->neighbors.begin(), f.face[(n+1)%3]->neighbors.end(), f.face[n]) != f.face[(n+1)%3]->neighbors.end())
                      {
                        std::cout << "wT!";
                      }
                      f.face[n]->neighbors.push_back(f.face[(n+1)%3]);
                      f.face[(n+1)%3]->neighbors.push_back(f.face[n]);
                    }
                  }
                }
              }
            }
          }
          i = itmp;
        }
        else
        {
          ++i;
        }
        ++ind;
      }
      std::cout << "pass: " << count << std::endl;
    //} while (!done);
    } while (0);
    std::cout << "OK" << std::endl;
    std::cout << "Removed " << count << " 3-connected points" << std::endl;
  }  
}


