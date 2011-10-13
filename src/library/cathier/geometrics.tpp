

namespace til { namespace geo
{

  //---------------------------------------------------------------------------------------------

  template < typename T >
  inline T herons_formula(T a, T b, T c)
  {
    T s = (a + b + c) / T(2.0);
    // NB: it actually happened that this product is negative even though a, b and c are positive.
    // This happens when one of the numbers is very small compared to the others.
    // So the max really is necessary.
    return std::sqrt(std::max(s*(s-a)*(s-b)*(s-c), T(0)));
  }

  //---------------------------------------------------------------------------------------------

  template < typename T >
  inline T herons_formula(numeric_array<T,3> const & lengths)
  {
    return herons_formula(lengths[0], lengths[1], lengths[2]);
  }

  //---------------------------------------------------------------------------------------------

  template < typename T >
  inline T law_of_cosines(T a, T b, T c)
  {
    return (b*b + c*c - a*a) / (2*b*c);
  }

  //---------------------------------------------------------------------------------------------

  template < typename T >
  inline T law_of_tangents(T a, T b, T c, T area)
  {
    return 4 * area / (b*b + c*c - a*a);
  }

  //---------------------------------------------------------------------------------------------

  template < typename T >
  inline T law_of_tangents(T a, T b, T c)
  {
    return law_of_tangents(a, b, c, herons_formula(a, b, c));
  }

  //---------------------------------------------------------------------------------------------
    
  // TODO: Add a policy to return the angle between 0 and 2PI, or between -PI and PI.
  template < typename T >
  inline T angle(T x, T y, T norm)
  {
    if (x < 0)
    {
      if (y < 0)  return -M_PI - 2.0 * std::atan( y / (norm - x) );
      else        return  M_PI - 2.0 * std::atan( y / (norm - x) );
    }
    else          return 2.0 * std::atan( y / (norm + x) );
  }

  //---------------------------------------------------------------------------------------------

  // TODO: a policy could be used for the case when the norm is too small.
  template < typename T >
  inline T angle(T x, T y)
  {
    T norm = x*x + y*y;
    if (norm < 128*std::numeric_limits<T>::epsilon())
    {
      return 0;
    }
    norm = std::sqrt(norm);
    return angle(x,y,norm);
  }

  //---------------------------------------------------------------------------------------------

  /// Converts cartesian coordinates into spherical coordinates
  template < typename T >
  void
  cartesian2spherical(T x, T y, T z, T & theta, T & phi, T & rho)
  {
    rho = norm(x, y, z);
    theta = angle(x, y);
    phi = std::acos(z / rho);
  }

  //---------------------------------------------------------------------------------------------

  template < typename T >
  inline
  void
  spherical2cartesian(T theta, T phi, T rho, T & x, T & y, T & z)
  {
    T rsinphi = rho * std::sin(phi);
    x = std::cos(theta) * rsinphi;
    y = std::sin(theta) * rsinphi;
    z = rho * std::cos(phi);
  }

  //---------------------------------------------------------------------------
  
  template < typename TArray >
  bool Triangle2Segment<TArray>::operator()
  (
    const TArray & A,   ///< Input: first triangle vertex
    const TArray & B,   ///< Input: second triangle vertex
    const TArray & C,   ///< Input: third triangle vertex
    TArray & X,         ///< Output: first segment end
    TArray & Y          ///< Output: second segment end
  )
  {
    const prec_type EPSILON = 128 * std::numeric_limits<prec_type>::epsilon();
    m_N = B - A;
    if (norm_inf(m_N) < EPSILON)
    {
      m_N = C - A;
      if (norm_inf(m_N) < EPSILON)
      {
        return false;
      }
      this->doit(A, C, B, X, Y);
    }
    else
    {
      this->doit(A, B, C, X, Y);
    }
    return true;
  }

  //---------------------------------------------------------------------------

  template < typename TArray >
  void Triangle2Segment<TArray>::doit
  (
    const TArray & A,   ///< Input: first triangle vertex
    const TArray & B,   ///< Input: second triangle vertex
    const TArray & C,   ///< Input: third triangle vertex
    TArray & X,         ///< Output: first segment end
    TArray & Y          ///< Output: second segment end
  )
  {
    prec_type lambda = dot(m_N, C-A);
    if (lambda < 0)
    {
      X = B;
      Y = C;
    }
    else if (lambda > norm2(m_N))
    {
      X = A;
      Y = C;
    }
    else
    {
      X = A;
      Y = B;
    }
  }

  //---------------------------------------------------------------------------

  template < typename TArray, typename TResArray >
  bool TriangleNormal<TArray, TResArray >::operator()
  (
    const TArray & A,   ///< Input: first triangle vertex
    const TArray & B,   ///< Input: second triangle vertex
    const TArray & C,   ///< Input: third triangle vertex
    TResArray & N       ///< Output: triangle normal
  )
  {
    const prec_type EPSILON = 128 * std::numeric_limits<prec_type>::epsilon();
    N = cross(B-A, C-A, prec<prec_type>());
    if (norm_inf(N) > EPSILON) return true;
    N = cross(C-B, A-B, prec<prec_type>());
    if (norm_inf(N) > EPSILON) return true;
    N = cross(A-C, B-C, prec<prec_type>());
    if (norm_inf(N) > EPSILON) return true;
    return false;
    //throw InvalidTriangle();
  }

  //---------------------------------------------------------------------------

  template < typename TArray, typename TResArray >
  bool
  triangle_normal(const TArray & A, const TArray & B, const TArray & C, TResArray & result)
  {
    return TriangleNormal<TArray, TResArray>()(A,B,C, result);
  }

}} // namespace til::geo

