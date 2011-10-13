#ifndef TIL_MATRIX3_TOOLS_H
#define TIL_MATRIX3_TOOLS_H

/// \file Belongs to matrix3 package -- do not include directly, include til/Matrix3.h instead


namespace til
{

  //---------------------------------------------------------------------------

  /// compute mat * in
  template < typename T1, typename T2, typename T3 >
  void
  apply
  (
    Matrix3<T1> const & mat
  , numeric_array<T2,3> const & vec
  , numeric_array<T3,3> & vecout
  )
  {
    for (std::size_t j = 0; j < 3; ++j)
    {
      vecout[j] = mat(0,j) * vec[0];
    }
    for (std::size_t i = 1; i < 3; ++i)
    {
      T2 veci = vec[i];
      for (std::size_t j = 0; j < 3; ++j)
      {
        vecout[j] += mat(i,j) * veci ;
      }
    }
  }

  //---------------------------------------------------------------------------

  /// compute transpose(mat) * in
  template < typename T1, typename T2, typename T3 >
  void
  apply_t
  (
    Matrix3<T1> const & mat
  , numeric_array<T2,3> const & vec
  , numeric_array<T3,3> & vecout
  )
  {
    for (std::size_t j = 0; j < 3; ++j)
    {
      vecout[j] = mat(j,0) * vec[0];
    }
    for (std::size_t i = 1; i < 3; ++i)
    {
      T2 veci = vec[i];
      for (std::size_t j = 0; j < 3; ++j)
      {
        vecout[j] += mat(j,i) * veci;
      }
    }
  }

  //---------------------------------------------------------------------------

  /*
  template < typename TMatrix, typename TVector >
  TMatrix diag(const TVector & v)
  {
    TMatrix m;
    redim(m, i, i);
    fill(m, typename value_type_of<TMatrix>::type(0));
    for (std::size_t i = 0; i < size(v); ++i)
    {
      m(i,i) = v(i);
    }
  }
  */

  template < typename T >
  T det(const Matrix3<T> & m)
  {
    return m(0,0)*m(1,1)*m(2,2)+m(1,0)*m(2,1)*m(0,2)+m(2,0)*m(0,1)*m(1,2) -
            m(0,2)*m(1,1)*m(2,0)-m(1,2)*m(2,1)*m(0,0)-m(2,2)*m(0,1)*m(1,0);
  }  

  /*
  /// Determinent of a square matrix.
  template < typename TPrec, typename TMatrix >
  TPrec det(const TMatrix &m)
  {
    assert(m.dims()[0] == m.dims()[1]);
    std::size_t n = m.dims()[0];
    // Actually we need to extract the first row, in the same format as the underlying data.
    typename TMatrix::row_type x;
    resize(x, n);
    std::copy(m.begin(), m.begin()+n, x.begin());
    // Actually we need a cyclic iterator here. Or a cyclic range for that matter
    for (std::size_t i = 1; i < n; ++i)
    {
      x[i] *= m( (i+1) % n, i);
    }
    TPrec res = x[n-1];
    for (int i = n-2; i >= 0; --i)
    {
      res = x[i] - res;
    }
    return res;
  }
  */
	//} // namespace linalg
} // namespace til

#endif

