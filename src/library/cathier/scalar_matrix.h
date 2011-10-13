#ifndef SCALARMATRIX_H_
#define SCALARMATRIX_H_

// includes from TIL
#include "til/numeric_array.h"

namespace til
{

  //---------------------------------------------------------------------------
  
  /// A scalar matrix.
  /// A scalar matrix is a super simple matrix. It is a diagonal matrix, whose diagonal elements are
  /// all equal.
  template < typename T >
  class Scalar_matrix
  {
  public: // typedefs
    typedef T value_type;
  public: // constructors
    explicit ScalarMatrix(prec_type value) : m_value(value) {}
  public: // operators
    value_type operator()(const numeric_array<std::size_t,2> & pos)
    {
      return (pos[0] == pos[1] ? m_value : 0);
    }
  public: // set & get
    value_type   value() const { return m_value; }
    value_type & value()       { return m_value; }
  private: // data
    value_type m_value;
  };
  
  //---------------------------------------------------------------------------

  /// Compute x^T.M.x
  template < typename T, typename TArray >
  T matdot(ScalarMatrix<T> m, const TArray & x)
  {
    return norm2<T>(n) * m.value();
  }

  //---------------------------------------------------------------------------

}


#endif /*SCALARMATRIX_H_*/
