#ifndef TIL_MATRIX3_H
#define TIL_MATRIX3_H

// includes from BOOST library
#include "boost/array.hpp"

// includes from TIL library
#include "til/til_common.h"
#include "til/multi_array.h"
#include "til/numeric_array.h"


// Ignore specific warnings
#ifdef _MSC_VER
#pragma warning (push)
#pragma warning (disable : 4231) // nonstandard extension used : 'extern' before template explicit instantiation
#endif


// TODO: This is crap, this has to change entirely...

namespace til
{
  //---------------------------------------------------------------------------

  template < typename TStorage >
  class matrix
  {
  public: // typedefs
    typedef matrix<TStorage>                      Self;
    typedef typename TStorage::value_type         value_type;
    typedef typename TStorage::reference          reference;
    typedef typename TStorage::const_reference    const_reference;
    typedef typename TStorage::range              range;
    typedef typename TStorage::const_range        const_range;
  public: // constructors
    /// Default contructor.
    matrix() : m_data() {}
    /// Initialization with size.
    matrix(numeric_array<std::size_t,2> dims) : m_data(dims) {}
    /// Initialization with size and fill value.
    matrix(numeric_array<std::size_t,2> dims, value_type v) : m_data(dims, v) {}
    /// Initialization from container.
    matrix(const TStorage & data) : m_data(data) {}
  public: // set & get
    const TStorage & data() const { return m_data; }
    TStorage & data() { return m_data; }
  public: // functions
    /// Returns matrix size.
    numeric_array<std::size_t,2> dims() const { return m_data.dims(); }
    range whole_range() { return m_data.whole_range(); }
    const_range whole_range() const { return m_data.whole_range(); }
  public: // operators
    reference operator()(const numeric_array<std::size_t,2> & pos)
    {
      return m_data(pos);
    }
    const_reference operator()(const numeric_array<std::size_t,2> & pos) const
    {
      return m_data(pos);
    }

  private: // data
    TStorage m_data;
  };

  //---------------------------------------------------------------------------

  template < typename TStorage >
  inline std::ostream & 
  operator<< (std::ostream & os, const matrix<TStorage> & m)
  {
    return m.data();
  }

  //---------------------------------------------------------------------------

  /// A full, square matrix whose size is known at compile time.
  template < std::size_t D, typename TPrec >
  struct simple_matrix
  {
    typedef matrix<multi_array<numeric_array<TPrec, D*D>, D> > type;
  };

  //---------------------------------------------------------------------------

  /// A mathematical matrix.
  template < typename T>
  class Matrix3
  {
  public: // typedefs
    typedef Matrix3<T> Self;
    typedef T value_type;
    typedef boost::array<T,3> row_type;

    typedef boost::array<T,3*3> Data;
    typedef typename Data::iterator iterator;
    typedef typename Data::const_iterator const_iterator;

  public: // constructors & destructor

    /// Empty matrix of zeros.
    Matrix3() { this->reset(); }

    /// No initialization.
    Matrix3(NoInit) {}

    /// Initialize with elements.
	  Matrix3(T xx, T xy, T xz, T yx, T yy, T yz, T zx, T zy, T zz)
	  {
		  (*this)(0,0) = xx;
		  (*this)(1,0) = xy;
		  (*this)(2,0) = xz;
		  (*this)(0,1) = yx;
		  (*this)(1,1) = yy;
		  (*this)(2,1) = yz;
		  (*this)(0,2) = zx;
		  (*this)(1,2) = zy;
		  (*this)(2,2) = zz;
	  }

    // TODO: actually, maybe all we need is a range.
    //row_type row(std::

    /*virtual*/ ~Matrix3() {} // Denis: removed virtual to avoid link errors on Mac

  public: // set & get

    boost::array<std::size_t,2> dims() const
    {
      boost::array<std::size_t,2> d = { {3, 3} };
      return d;
    }

    const T & operator()(std::size_t i, std::size_t j) const
    {
      assert(i<3);
      assert(j<3);
      return m_values[i+3*j];
    }

    T & operator()(std::size_t i, std::size_t j)
    {
      assert(i<3);
      assert(j<3);
      return m_values[i+3*j];
    }

	  // Reset the matrix to default values
    void reset() { for (std::size_t i=0; i<3*3; ++i) m_values[i] = T(); }

    iterator begin() { return m_values.begin(); }
    iterator end() { return m_values.end(); }
    const_iterator begin() const { return m_values.begin(); }
    const_iterator end() const { return m_values.end(); }
  	
  public: // operators
  
    void operator+=(Self const & mat)
    {
      for (int j = 0; j < 3; ++j)
      for (int i = 0; i < 3; ++i)
      {
        (*this)(i,j) += mat(i,j);
      }
    }

    void operator-=(Self const & mat)
    {
      for (int j = 0; j < 3; ++j)
      for (int i = 0; i < 3; ++i)
      {
        (*this)(i,j) -= mat(i,j);
      }
    }
    
  private: // data

	  boost::array<T,3*3> m_values;
  };

  //---------------------------------------------------------------------------

  template < typename T >
  inline Matrix3<T>
  operator-(Matrix3<T> const & mat)
  {
    Matrix3<T> res;
    for (int j = 0; j < 3; ++j)
    for (int i = 0; i < 3; ++i)
    {
      res(i,j) = -mat(i,j);
    }
    return res;
  }

  //---------------------------------------------------------------------------

  // EXPIMP_TEMPLATE template class TIL_API Matrix3<double>;

  //---------------------------------------------------------------------------

  template < typename T1, typename T2 >
  inline numeric_array<typename combine<T1,T2>::type,3>
  operator*(const Matrix3<T1> & m, const numeric_array<T2,3> & v)
  {
    numeric_array<typename combine<T1,T2>::type, 3> res;
    for (std::size_t j = 0; j < 3; ++j)
    {
      res[j] = m(0,j) * v[0];
      for (std::size_t i = 1; i < 3; ++i)
      {
        res[j] += m(i,j) * v[i];
      }
    }
    return res;
  }

  //---------------------------------------------------------------------------

  template < typename T1, typename T2 >
  inline T2 operator*(const Matrix3<T1> & m, const T2 & v)
  {
    T2 res;
    for (std::size_t j = 0; j < 3; ++j)
    {
      res[j] = m(0,j) * v[0];
      for (std::size_t i = 1; i < 3; ++i)
      {
        res[j] += m(i,j) * v[i];
      }
    }
    return res;
  }

  //---------------------------------------------------------------------------

  template < typename T >
  inline std::ostream & 
  operator<< (std::ostream & os, const Matrix3<T> & m)
  {
    for (int i = 0; i < 3; ++i)
    {
      os << (i == 0? "[ " : "  ");
      for (int j = 0; j < 3; ++j)
      {
        os << m(j,i) << " ";
      }
      os << (i == 2? "]\n":"\n");
    }
    return os;
  }

  //---------------------------------------------------------------------------

} // namespace til


#ifdef _MSC_VER
#pragma warning(pop)
#endif

// package includes
#include "til/matrix3Tools.h"

#endif
