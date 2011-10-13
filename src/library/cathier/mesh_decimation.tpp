#ifndef TIL_MESH_DECIMATION_TPP_
#define TIL_MESH_DECIMATION_TPP_

#include "til/miscTools.h"

namespace til
{
  
  //-------------------------------------------------------------------------------------------

  namespace
  {
    
    template < typename TFaces, typename TPoints, typename TNeighbors >
    void remove_delaunay_faces(TFaces & faces, TNeighbors & neighbors, const TPoints & points)
    {
      typename TFaces::iterator f = faces.begin();
      while (f != faces.end())
      {
        if (til::geo::is_in_circumcircle<double>(points[i], points[(*f)[0]], points[(*f)[1]], points[(*f)[2]]))
        {
          neighbors.insert((*f)[0]);
          neighbors.insert((*f)[1]);
          neighbors.insert((*f)[2]);
          f = faces.erase(f);
        }
        else
        {
          ++f;
        }
      }
    }
        
  } // namespace
  
  template < typename TFace >
  std::list<boost::array<std::size_t,3> >
  SimpleDelaunayTriangulation::operator()(std::vector<til::numeric_array<float, 2> > points)
  {
    typedef std::vector<til::numeric_array<float, 2> > Points;
    typedef TFace Face;
    typedef std::list<Face> FaceList;
        
    // Don't waste time calling this for less than 4 points.
    assert(points.size() >= 4);
    
    FaceList faces;
  
    // compute bounding triangle and add points to list
    bounding_triangle(points.begin(), points.end(), std::back_inserter(points));
    // add triangle to list
    boost::array<std::size_t,3> face = { points.size()-1, points.size()-2, points.size()-3 };
    faces.push_back(face);
  
    for (std::size_t i = 0; i < n; ++i)
    {
      // Remove triangles if we lie in their circumcircle, and collect their vertex.
      std::set<std::size_t> neighbors;
      remove_delaunay_faces(faces, neighbors, points);
      {
        FaceList::iterator f = faces.begin();
        while (f != faces.end())
        {
          if (til::geo::is_in_circumcircle<double>(points[i], points[(*f)[0]], points[(*f)[1]], points[(*f)[2]]))
          {
            neighbors.insert((*f)[0]);
            neighbors.insert((*f)[1]);
            neighbors.insert((*f)[2]);
            f = faces.erase(f);
          }
          else
          {
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
        til::numeric_array<float,2> v = points[*iNeighbor]-points[i];
        t.angle = til::geo::angle(v[0], v[1]);
        tmp.push_back(t);
      }
      //std::sort(tmp.begin(), tmp.end(), TempComp());
      std::sort(tmp.begin(), tmp.end());
      
      // Form new faces
      for (std::size_t it = 0; it < tmp.size(); ++it)
      {
        //std::cout << "Adding face " << i << " " << tmp[it].index << " " << tmp[ (it+1) % til::size(tmp) ].index << std::endl;
        Face face = {i, tmp[it].index, tmp[ (it+1) % til::size(tmp) ].index };
        faces.push_back(face);
      }
    }
    /*
      {
        std::cout << "faces" << std::endl;
        for (FaceList::iterator f = faces.begin(); f != faces.end(); ++f)
        {
          std::cout << (*f)[0] << " " << (*f)[1] << " " << (*f)[2] << std::endl;
        }
        std::cout << "endfaces" << std::endl;
      }
  */
  
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
  /*
      {
        FaceList::iterator f = faces.begin();
        std::cout << "cross2D: " << til::cross(points[(*f)[1]] - points[(*f)[0]], points[(*f)[2]] - points[(*f)[0]]) << std::endl;
      }
  */
  
    // Check that normals are ok
    {
      for (FaceList::iterator f = faces.begin(); f != faces.end(); ++f)
      {
        if (til::cross(points[(*f)[1]] - points[(*f)[0]], points[(*f)[2]] - points[(*f)[0]]) < 0)
        {
          std::cout << "wS!";
        }
      }
    }
  
  
    // Check if we have to resort to normal testing to remove faces
    if (faces.size() != n-2)
    {
      /*
      std::cout << "NI!" << faces.size() << "-" << n << std::endl;
      {
        std::cout << "faces2" << std::endl;
        for (FaceList::iterator f = faces.begin(); f != faces.end(); ++f)
        {
          std::cout << (*f)[0] << " " << (*f)[1] << " " << (*f)[2] << std::endl;
        }
        std::cout << "endfaces" << std::endl;
      }
      */
          
      // Remove faces with bad normal
      /*
      {
        FaceList::iterator f = faces.begin();
        while (f != faces.end())
        {
          if (til::cross(points[(*f)[1]] - points[(*f)[0]], points[(*f)[2]] - points[(*f)[0]]) < 0)
          {
            f = faces.erase(f);
          }
          else
          {
            ++f;
          }
        }
      }
      */
  
      {
        FaceList::iterator f = faces.begin();
        while (f != faces.end())
        {
          bool flagdel = false;
          int count = 0;
          for (int i = 0; i < 3; ++i)
          {
            if ((*f)[i] > (*f)[(i+1)%3])
            {
              if (++count > 1)
              {
                flagdel = true;
                //f = faces.erase(f);
                break;
              }
            }
          }
          if (flagdel) f = faces.erase(f);
          else ++f;
        }
      }
      
      // Check again that we have the expected number of faces
      if (faces.size() != n-2)
      {
        std::cout << "NI2!" << faces.size() << "-" << n << std::endl;
        {
          std::cout << "faces2" << std::endl;
          for (FaceList::iterator f = faces.begin(); f != faces.end(); ++f)
          {
            std::cout << (*f)[0] << " " << (*f)[1] << " " << (*f)[2] << std::endl;
          }
          std::cout << "endfaces" << std::endl;
          std::cout << "points" << std::endl;
          for (std::size_t i = 0; i < n; ++i) std::cout << points[i] << std::endl;
          std::cout << "end points" << std::endl;
        }
      }
    }
    
    // check that all outer edges is missing
    {
      std::set<boost::array<std::size_t,2> > outer_edges;
      for (FaceList::iterator f = faces.begin(); f != faces.end(); ++f)
      {
        for (std::size_t i = 0; i < 3; ++i)
        {
          std::size_t delta = std::abs(int((*f)[i]) - int((*f)[(i+1)%3]));
          if (delta == 1 || delta == (n-1))
          {
            boost::array<std::size_t,2> tmp = {min((*f)[i], (*f)[(i+1)%3]), max((*f)[i], (*f)[(i+1)%3])};
            outer_edges.insert(tmp);
          }
        }
      }
      if (outer_edges.size() != n)
      {
        std::cout << "wU!";
      }
    }
    
    //exit(0);
    
    /*
    std::cout << "nfaces: " << til::size(faces) << std::endl;
    for (FaceList::const_iterator i = faces.begin(); i != faces.end(); ++i)
    {
      std::cout << (*i)[0] << " " << (*i)[1] << " " << (*i)[2] << std::endl;
    }
    */  
    return faces;
  }
  
  
  
} // namespace til

#endif /*MESH_DECIMATION_TPP_*/
