#ifndef GEOMETRICS_H_
#define GEOMETRICS_H_

/// \file Collects code regarding basic geometrical operations.

// includes from TIL
#include "til/numeric_array.h"
#include "til/Matrix3.h"

// includes from TIL
#include "fuzzy_logic.h"
#include "miscUtils.h"


namespace til { namespace geo
{

    //-----------------------------------//
   //  Small handy geometric functions  //
  //-----------------------------------//

  //---------------------------------------------------------------------------

  /// Heron's formula.
  /// The Heron's formula yields the surface of a triangle given the length of
  /// its edges.
  template < typename T >
  inline T herons_formula(T a, T b, T c);

  //---------------------------------------------------------------------------

  /// Heron's formula.
  /// The Heron's formula yields the surface of a triangle given the length of 
  /// its edges.
  template < typename T >
  inline T herons_formula(numeric_array<T,3> const & lengths);

  //---------------------------------------------------------------------------

  /// Returns the cosine of angle A given the lengths of the triangle edges.
  template < typename T >
  inline T law_of_cosines(T a, T b, T c);

  //---------------------------------------------------------------------------

  /// Returns the tangent of angle A given the area and the lengths of the
  /// edges of the triangle.
  template < typename T >
  inline T law_of_tangents(T a, T b, T c, T area);

  //---------------------------------------------------------------------------

  /// Returns the tangent of angle A given the lengths of the edges of the
  /// triangle.
  template < typename T >
  inline T law_of_tangents(T a, T b, T c);

  //---------------------------------------------------------------------------
    
  /// Returns the angle between (x,y) and the x-axis, with the norm of (x,y) provided.
  /// NB: this version is clearly intended for fast computation and thus no checking on norm is done:
  /// in particular, you have to check yourself that the norm is not too small.
  template < typename T >
  inline T angle(T x, T y, T norm);

  //---------------------------------------------------------------------------

  /// Returns the angle between (x,y) and the x-axis.
  template < typename T >
  inline T angle(T x, T y);

  //---------------------------------------------------------------------------

  /// Converts cartesian coordinates into spherical coordinates.
  template < typename T >
  void cartesian2spherical(T x, T y, T z, T & theta, T & phi, T & rho);

  //---------------------------------------------------------------------------

  /// Converts spherical coordinates into cartesian coordinates.
  // NB: I put an inline here in the hope that when a rho=1 is passed as an argument, 
  // multiplications whith rho will be dropped.
  template < typename T >
  inline void spherical2cartesian(T theta, T phi, T rho, T & x, T & y, T & z);

  //---------------------------------------------------------------------------

    //--------------------//
   //  Triangle2Segment  //
  //--------------------//

  /// Approximate an ill-conditionned triangle by a segment.
  /// Basically, it figures out which triangle vertices to keep in order to approximate the triangle.
  /// NB: this doesn't include the condition test. Callers have to decide when to call this class.
  template < typename TArray >
  class Triangle2Segment
  {
  public: // typedef
    typedef typename TArray::value_type prec_type;
  
  public: // functions
  
    /// Return true if successfull
    bool operator()
    (
      const TArray & A,   ///< Input: first triangle vertex
      const TArray & B,   ///< Input: second triangle vertex
      const TArray & C,   ///< Input: third triangle vertex
      TArray & X,         ///< Output: first segment end
      TArray & Y          ///< Output: second segment end
    );

  private: // functions
  
    void doit
    (
      const TArray & A,   ///< Input: first triangle vertex
      const TArray & B,   ///< Input: second triangle vertex
      const TArray & C,   ///< Input: third triangle vertex
      TArray & X,         ///< Output: first segment end
      TArray & Y          ///< Output: second segment end
    );
    
  private: // data
    TArray m_N;
  };

  //---------------------------------------------------------------------------  
  
  template < typename TArray >
  bool triangle2segment(const TArray & A, const TArray & B, const TArray & C, TArray & X, TArray & Y)
  {
    return Triangle2Segment<TArray>()(A, B, C, X, Y);
  }

  //---------------------------------------------------------------------------
  
    //------------------//
   //  TriangleNormal  //
  //------------------//

  /// Compute the normal of a triangle.
  /// Basically, it makes sure that the computed normal is numerically stable.
  //template < typename TArray, typename TPrec, typename TResArray = typename change_precision<TArray, TPrec>::type >
  template < typename TArray, typename TResArray >
  class TriangleNormal
  {
  public: // exceptions
    /// This exeception is thrown whenever a reliable computation of the triangle normal cannot
    /// be achieved.
    class InvalidTriangle : public std::exception {};
  public: // typedefs
    //typedef typename TArray::value_type prec_type;
    //typedef TPrec prec_type;
    typedef typename TResArray::value_type prec_type;
    //typedef typename change_precision<TArray, TPrec>::type TPrecArray;
  public: // functions

    /// Compute a triangle's normal given its three vertices.
    /// Return true if successful.
    bool operator()
    (
      const TArray & A,   ///< Input: first triangle vertex
      const TArray & B,   ///< Input: second triangle vertex
      const TArray & C,   ///< Input: third triangle vertex
      TResArray & N       ///< Output: triangle normal
    );
  };

  
  //---------------------------------------------------------------------------

  /// Compute the normal of a triangle.
  /*
  template < typename TPrec, typename TArray >
  typename change_precision<TArray, TPrec>::type
  triangle_normal(const TArray & A, const TArray & B, const TArray & C)
  {
    return TriangleNormal<TArray,TPrec,typename change_precision<TArray, TPrec>::type>()(A,B,C);
  }
  */

  //---------------------------------------------------------------------------

  template < typename TArray, typename TResArray >
  //typename change_precision<TArray, TPrec>::type
  bool
  triangle_normal(const TArray & A, const TArray & B, const TArray & C, TResArray & result);

  //---------------------------------------------------------------------------

  namespace
  {
    template < typename TNormal, typename TPoint >
    inline void normal2D(const TPoint & A, const TPoint & B, TNormal & N)
    {
      N[0] = A[1] - B[1];
      N[1] = B[0] - A[0];
    }
  }


  //---------------------------------------------------------------------------

    //----------------//
   //  IsIncludedIn  //
  //----------------//

  /// A namespace structure collecting geometrical object inclusion algorithms
  struct IsIncludedIn
  {
    /// Checks if a 2D point lies within a 2D triangle.
    template < typename TArray >
    class PointInTriangle2D
    {
    public: // typedefs    
      typedef typename TArray::value_type prec_type;
    public: // functions
      bool operator()
      (
        const TArray & P,   ///< a point
        const TArray & A,   ///< first triangle vertex
        const TArray & B,   ///< second triangle vertex
        const TArray & C    ///< third triangle vertex
      )
      {
        normal2D(A, B, m_N);
        prec_type sAB = dot(P - A, m_N);
        prec_type s = dot(C - A, m_N);
        if (!fuzzy::same_sign(s, sAB)) return false;
        normal2D(B, C, m_N);
        prec_type sBC = dot(P - B, m_N);
        if (!fuzzy::same_sign(sBC, sAB)) return false;
        normal2D(C, A, m_N);
        prec_type sCA = dot(P - C, m_N);
        if (!fuzzy::same_sign(sCA, sAB)) return false;
        return true;
      }
    private: // data, internals
      numeric_array<prec_type,2> m_N;
    };

    /// Checks if a 2D point lies within the circumcircle of a triangle, whose vertices are given
    /// in a counterclockwise fashion.
    template < typename TPrec, typename TPoint >
    struct PointInCircumcircle2D_counterclockwise
    {
      bool                ///< True iff point lies within circumcircle.
      operator()
      (
       const TPoint & p,  ///< a point
       const TPoint & a,  ///< first triangle vertex
       const TPoint & b,  ///< second triangle vertex
       const TPoint & c   ///< third triangle vertex.
      )
      {
        TPrec d = 0.0;
        Matrix3<TPrec> m;
      
        m(0,0) = a[0];
        m(1,0) = a[1];
        m(2,0) = norm2(a, prec<TPrec>());
        m(0,1) = b[0];
        m(1,1) = b[1];
        m(2,1) = norm2(b, prec<TPrec>());
        m(0,2) = c[0];
        m(1,2) = c[1];
        m(2,2) = norm2(c, prec<TPrec>());
        d += det(m);
      
        m(0,2) = p[0];
        m(1,2) = p[1];
        m(2,2) = norm2(p, prec<TPrec>());
        d -= det(m);
      
        m(0,1) = c[0];
        m(1,1) = c[1];
        m(2,1) = norm2(c, prec<TPrec>());
        d += det(m);
      
        m(0,0) = b[0];
        m(1,0) = b[1];
        m(2,0) = norm2(b, prec<TPrec>());
        d -= det(m);
        
        return d > 0;
      }
    };

    /// Checks if a 2D point lies within the circumcircle of a 2D triangle.
    template < typename TPrec, typename TPoint >
    struct PointInCircumcircle2D
    {
      bool                ///< True iff point lies within circumcircle
      operator()
      (
       const TPoint & p,  ///< a point
       const TPoint & a,  ///< first triangle vertex
       const TPoint & b,  ///< second triangle vertex
       const TPoint & c   ///< third triangle vertex
      )
      {
        // ensure a, b, c lie counter-clockwise
        if ((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]) >= 0)
        {
          return IsIncludedIn::PointInCircumcircle2D_counterclockwise<TPrec, TPoint>()(p,a,b,c);
          //return is_in_circumcircle_counterclockwise<TPrec>(p,a,b,c);
        }
        else
        {
          return IsIncludedIn::PointInCircumcircle2D_counterclockwise<TPrec, TPoint>()(p,a,c,b);
          //return is_in_circumcircle_counterclockwise<TPrec>(p,a,c,b);
        }
      }
    };
  };

  //---------------------------------------------------------------------------

    //-----------//
   //  Project  //
  //-----------//

  /// A namespace structure collecting geometrical projection algorithms.
  struct Project
  {
    
    template < typename TArray, typename TPrecArray >
    class PointOnLine3D
    {
    public: // typedefs    
      typedef typename TPrecArray::value_type prec_type;
      //typedef typename TArray::value_type prec_type;
      //typedef typename change_precision<TArray, TPrec>::type TPrecArray;
    //private: // typedefs
      //typedef numeric_array<prec_type,2> Array2D;
    public: // functions
      void operator()
      (
        const TArray & P,
        const TArray & A, const TArray & B,
        TPrecArray & output
      )
      {
        using namespace til::expr;
        detail::loop_xxx(*_1 = cast<prec_type>(*_2 - *_3), m_N.begin(), m_N.end(), B.begin(), A.begin());
        //m_N = B - A;
        m_N *= 1/norm<prec_type>(m_N);
        prec_type lambda = dot(P-A, m_N, prec<typename TPrecArray::value_type>());
        detail::loop_xxx(*_1 = *_2 + *_3 * lambda, output.begin(), output.end(), A.begin(), m_N.begin());
        //output = A + m_N * dot(P-A, m_N, prec<typename TPrecArray::value_type>());
      }
    private:
      TPrecArray m_N;
    };
    
    template < typename TArray, typename TPrecArray >
    class PointOnSegment3D
    {
    public: // typedefs
      typedef typename TPrecArray::value_type prec_type;
    public: // functions
      void operator()
      (
        const TArray & P,
        const TArray & A, const TArray & B,
        TPrecArray & output
      )
      {
        using namespace til::expr;
        detail::loop_xxx(*_1 = cast<prec_type>(*_2 - *_3), m_N.begin(), m_N.end(), B.begin(), A.begin());
        prec_type lambda =  dot(P-A, m_N, prec<typename TPrecArray::value_type>());
        if (lambda < 0)
        {
          std::copy(A.begin(), A.end(), output.begin());
        }
        else if (lambda > norm2(m_N, prec<typename TPrecArray::value_type>()))
        {
          std::copy(B.begin(), B.end(), output.begin());
        }
        else
        {
          detail::loop_xxx(*_1 = *_2 + *_3 * lambda, output.begin(), output.end(), A.begin(), m_N.begin());
        }
      }
    private: // data
      TPrecArray m_N;      
    };
    
    /// Project a point onto a 3D triangle.
    /// TODO: Do we need these projections? For registration, I guess we could simply use the distance to the
    /// projected point, which may be more efficient.
    template < typename TArray, typename TOutputArray >
    class PointOnTriangle3D
    {
    public: // typedefs    
      typedef typename TOutputArray::value_type prec_type;
      //typedef typename TArray::value_type prec_type;
      //typedef typename change_precision<TArray, TPrec>::type TPrecArray;
    private: // typedefs
      typedef numeric_array<prec_type,2> Array2D;
    public: // functions
      void operator()
      (
        const TArray & P,
        const TArray & A, const TArray & B, const TArray & C,
        TOutputArray & output
      )
      {
        // Compute triangle normal
        if (!triangle_normal(A, B, C, m_N))
        {
          // If triangle is degenerate, compute projection on a segment.
          TArray X, Y;
          triangle2segment(A, B, C, X, Y);
          PointOnSegment3D<TArray, TOutputArray>()(P, X, Y, output);
          // TODO: projection on a point should follow...
          return;
        }
        output += dot(A-P, m_N, prec<typename TOutputArray::value_type>())*m_N;
        // project everybody on the closest canonical 2D plane
        std::size_t i = min_index(std::abs(m_N[0]),std::abs(m_N[1]),std::abs(m_N[2]));
        std::size_t j = (i+1)%3;
        std::size_t k = (i+2)%3;
        Array2D a; a[0] = A[j]; a[1] = A[k];
        Array2D b; b[0] = B[j]; b[1] = B[k];
        Array2D c; c[0] = C[j]; c[1] = C[k];
        Array2D p; p[0] = p[j]; p[1] = p[k];

        Array2D n;
        normal2D(a, b, n);
        prec_type sAB = dot(p - a, n);
        prec_type s = dot(c - a, n);
        normal2D(b, c, n);
        prec_type sBC = dot(p - b, n);
        normal2D(c, a, n);
        prec_type sCA = dot(p - c, n);
                
        if (!fuzzy::same_sign(s, sAB))
        {
          if (!fuzzy::same_sign(s, sBC))
          {
            std::copy(B.begin(), B.end(), output.begin());
            //output = B;
          }
          else if (!fuzzy::same_sign(s, sCA))
          {
            std::copy(A.begin(), A.end(), output.begin());
            //output = A;
          }
          else
          {
            m_projline(P, A, B, output);
          }
        }
        else if (!fuzzy::same_sign(s, sBC))
        {
          if (!fuzzy::same_sign(s, sCA))
          {
            std::copy(C.begin(), C.end(), output.begin());
            //output = C;
          }
          else
          {
            m_projline(P, B, C, output);
          }
        }
        else if (!fuzzy::same_sign(s, sCA))
        {
          m_projline(P, C, A, output);
        }
      }
    
    private: // data, internals
      TOutputArray m_N;
      Project::PointOnLine3D<TArray, TOutputArray> m_projline;
    };
  };

  //---------------------------------------------------------------------------

    //-------------------//
   //  AreIntersecting  //
  //-------------------//
  
  /// A namespace structure collecting geometrical object intersection algorithms
  struct AreIntersecting
  {
    /// Check whether two 2D segments are intersecting or not
    template < typename TArray >
    class Segments2D
    {
    public: // typedefs    
      typedef typename TArray::value_type prec_type;
    public: // functions
      bool operator()
      (
        const TArray & A1, const TArray & B1,
        const TArray & A2, const TArray & B2
      )
      {
        normal2D(A1, B1, m_N);
        if (fuzzy::same_sign(dot(A2 - A1, m_N), dot(B2 - A1, m_N)))
        {
          return false;
        }        
        normal2D(A2, B2, m_N);
        if (fuzzy::same_sign(dot(A1 - A2, m_N), dot(B1 - A2, m_N)))
        {
          return false;
        }
        return true;
      }
    private: // data, internals
      numeric_array<prec_type,2> m_N;
    };

    /// Check whether two 2D triangles are intersecting or not
    template < typename TArray >
    class Triangles2D
    {
    public: // functions
      bool operator()
      (
        const TArray & A1, const TArray & B1, const TArray & C1,
        const TArray & A2, const TArray & B2, const TArray & C2
      )
      {
        // Check if edges are intersecting
        if (m_interseg(A1,B1,A2,B2)) return true;
        if (m_interseg(A1,C1,A2,B2)) return true;
        if (m_interseg(B1,C1,A2,B2)) return true;

        if (m_interseg(A1,B1,A2,C2)) return true;
        if (m_interseg(A1,C1,A2,C2)) return true;
        if (m_interseg(B1,C1,A2,C2)) return true;

        if (m_interseg(A1,B1,B2,C2)) return true;
        if (m_interseg(A1,C1,B2,C2)) return true;
        if (m_interseg(B1,C1,B2,C2)) return true;
        
        // Check if one triangle is completely included into the other
        if (m_pointInTriangle(A1, A2, B2, C2)) return true;
        if (m_pointInTriangle(A2, A1, B1, C1)) return true;
        
        return false;
      }
      
    private: // data, internals
      // NB: these are put as private data, so that the user can avoid frequent reallocation of their own internals.
      AreIntersecting::Segments2D<TArray>      m_interseg;
      IsIncludedIn::PointInTriangle2D<TArray>  m_pointInTriangle;
    };


    /// Check whether two triangles are intersecting or not.
    /// TODO: I think I could use a 'strict' argument, that would decide on whether we want a strict or
    /// loose interesction, or more precisely, if we consider the triangles as topologically open or closed.
    /// Right now they are considered as open.
    template < typename TArray >
    class Triangles3D
    {  
    public: // typedefs
    
      typedef typename TArray::value_type prec_type;
      typedef numeric_array<prec_type,2> Array2D;
      //typedef TPrec prec_type;
      //typedef typename change_precision<TArray, TPrec>::type TPrecArray;

    public: // functions
          
      /// Check if two triangles in 3D are intersecting
      bool operator()
      (
        const TArray & A1, const TArray & B1, const TArray & C1,
        const TArray & A2, const TArray & B2, const TArray & C2
      )
      {
        const prec_type EPSILON = 128 * std::numeric_limits<prec_type>::epsilon();
  
        triangle_normal(A1, B1, C1, m_N1);
        m_N1 *= 1/norm<prec_type>(m_N1);
        // Checking if the vertices of the second triangle are on one side of the first
        prec_type dA2 = dot(m_N1, A2-A1);
        prec_type dB2 = dot(m_N1, B2-A1);
        prec_type dC2 = dot(m_N1, C2-A1);
        bool sAB2 = fuzzy::same_sign(dA2, dB2);
        bool sAC2 = fuzzy::same_sign(dA2, dC2);
        bool sBC2 = fuzzy::same_sign(dB2, dC2);
        // NB: one can be surprised to have three tests instead of two. However, it improves robustness, esp.
        // when one of the distance is exactly zero.
        if (sAB2 && sAC2 && sBC2) return false;
  
        triangle_normal(A2, B2, C2, m_N2);
        m_N2 *= 1/norm<prec_type>(m_N2);
        // Checking if the vertices of the first triangle are on one side of the second
        prec_type dA1 = dot(m_N2, A1-A2);
        prec_type dB1 = dot(m_N2, B1-A2);
        prec_type dC1 = dot(m_N2, C1-A2);
        bool sAB1 = fuzzy::same_sign(dA1, dB1);
        bool sAC1 = fuzzy::same_sign(dA1, dC1);
        bool sBC1 = fuzzy::same_sign(dB1, dC1);
        if (sAB1 && sAC1 && sBC1) return false;  
  
        m_D = cross<prec_type>(m_N1, m_N2);
        prec_type D0 = std::abs(m_D[0]);
        prec_type D1 = std::abs(m_D[1]);
        prec_type D2 = std::abs(m_D[2]);
        if (max(D0, D1, D2) < EPSILON)
        {  // coplanar, 2D case

          // Project points on a 2D plane
          // Instead of projecting on the plane with normal D, we project on a canonical plane
          // maximizing their areas -- this avoid many multiplications.
          std::size_t i = min_index(D0,D1,D2);
          std::size_t j = (i+1)%3;
          std::size_t k = (i+2)%3;
          
          Array2D a1; a1[0] = A1[j]; a1[1] = A1[k];
          Array2D b1; b1[0] = B1[j]; b1[1] = B1[k];
          Array2D c1; c1[0] = C1[j]; c1[1] = C1[k];
  
          Array2D a2; a2[0] = A2[j]; a2[1] = A2[k];
          Array2D b2; b2[0] = B2[j]; b2[1] = B2[k];
          Array2D c2; c2[0] = C2[j]; c2[1] = C2[k];

          return m_intersectTriangle2D(a1,b1,c1,a2,b2,c2);
        }
        else
        {
          dA1 = std::abs(dA1);
          dB1 = std::abs(dB1);
          dC1 = std::abs(dC1);
          dA2 = std::abs(dA2);
          dB2 = std::abs(dB2);
          dC2 = std::abs(dC2);
          
          /*
          // NB: numerically, this doesn't seem to work, even though moeller pretends it does
          std::size_t i = max_index(D0, D1, D2);
          prec_type a1 = A1[i];
          prec_type b1 = B1[i];
          prec_type c1 = C1[i];
          prec_type a2 = A2[i];
          prec_type b2 = B2[i];
          prec_type c2 = C2[i];
          */
  
          prec_type a1 = dot<prec_type>(A1, m_D);
          prec_type b1 = dot<prec_type>(B1, m_D);
          prec_type c1 = dot<prec_type>(C1, m_D);
          prec_type a2 = dot<prec_type>(A2, m_D);
          prec_type b2 = dot<prec_type>(B2, m_D);
          prec_type c2 = dot<prec_type>(C2, m_D);
          
                  
          prec_type emin, emax;
          if (!sAB1 && !sAC1)
          {
            emin = line_coord(a1, b1, dA1, dB1);
            emax = line_coord(a1, c1, dA1, dC1);
          }
          else if (!sAB1)
          {
            emin = line_coord(b1, a1, dB1, dA1);
            emax = line_coord(b1, c1, dB1, dC1);
          }
          else
          {
            emin = line_coord(c1, a1, dC1, dA1);
            emax = line_coord(c1, b1, dC1, dB1);
          }
          if (emin > emax) std::swap(emin, emax);
          
  
          prec_type fmin, fmax;
          if (!sAB1 && !sAC1)
          {
            fmin = line_coord(a2, b2, dA2, dB2);
            fmax = line_coord(a2, c2, dA2, dC2);
          }
          else if (!sAB1)
          {
            fmin = line_coord(b2, a2, dB2, dA2);
            fmax = line_coord(b2, c2, dB2, dC2);
          }
          else
          {
            fmin = line_coord(c2, a2, dC2, dA2);
            fmax = line_coord(c2, b2, dC2, dB2);
          }
          if (fmin > fmax) std::swap(fmin, fmax);
          
          if (fmin > emax || fmax < emin) return false;
  
          /*
          std::cout << A1 << " " << B1 << " " << C1 << std::endl;
          std::cout << A2 << " " << B2 << " " << C2 << std::endl;
    
          std::cout << m_N1 << std::endl;
          std::cout << A2-A1 << " " << B2-A1 << " " << C2-A1 << std::endl;
          std::cout << dA2 << " " << dB2 << " " << dC2 << std::endl;
          //std::cout << sAB2 << " " << sAC2 << std::endl;
          std::cout << m_N2 << std::endl;
          std::cout << A1-A2 << " " << B1-A2 << " " << C1-A2 << std::endl;
          std::cout << dA1 << " " << dB1 << " " << dC1 << std::endl;
  
          std::cout << m_D << std::endl;
          std::cout << D0 << " " << D1 << " " << D2 << std::endl;
          std::cout << i << std::endl;        
  
          std::cout << emin << " " << emax << " " << fmin << " " << fmax << std::endl;
          */
          
          return true;
        }
      }
      
    private: // functions
  
      prec_type line_coord(prec_type a, prec_type b, prec_type da, prec_type db)
      { return a + (b-a) * da / (da + db); }
  
    
    private: // data, internals
    
      AreIntersecting::Triangles2D<Array2D> m_intersectTriangle2D;
      /*  
      TPrecArray m_N1;
      TPrecArray m_N2;
      TPrecArray m_D;
      */
      ///*
      TArray m_N1;
      TArray m_N2;
      TArray m_D;
      //*/
    };
  };  


  //---------------------------------------------------------------------------

    //--------------------//
   //  Helper functions  //
  //--------------------//

  /// Check whether a point lies within the circumcircle of a given triangle.
  template < typename TPrec, typename TPoint >
  inline bool is_in_circumcircle
  (
       const TPoint & p,  ///< a point
       const TPoint & a,  ///< first triangle vertex
       const TPoint & b,  ///< second triangle vertex
       const TPoint & c   ///< third triangle vertex
  )
  {
    return IsIncludedIn::PointInCircumcircle2D<TPrec, TPoint>()(p,a,b,c);
  }


}} // namespace til::geo

// package include
#include "geometrics.tpp"


#endif /*GEOMETRICS_H_*/
