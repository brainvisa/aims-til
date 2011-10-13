#ifndef TIL_SYMMATRIX3_H
#define TIL_SYMMATRIX3_H

// Local includes
#include "til/til_common.h"
#include "til/til_declarations.h"
#include "til/Matrix3.h"
#include "til/value_proxy.h"


// Ignore specific warnings
#ifdef _MSC_VER
#pragma warning (push)
#pragma warning (disable : 4231) // nonstandard extension used : 'extern' before template explicit instantiation
#endif

// Namespace 
namespace til {
	//namespace linalg {

  /// A class to store a 3*3 symetric matrix
  template < typename T >
  class SymMatrix3 : public Matrix3<T>
  {
  public: // typedefs
    
    typedef SymMatrix3<T>     Self;
    typedef Matrix3<T>        Base;
    // TODO: this should actually be const T &. Go and check Matrix3
    typedef T                 const_reference;
    typedef const T *         const_pointer;
  private: // class

    class ValueProxy;

  public: // Constructors and destructor

    /// Empty matrix full of zeros
    SymMatrix3() : Base() { this->reset(); }
    
    /// No initialization
    SymMatrix3(NoInit) {}

    /// Initialize with elements
	  SymMatrix3(T xx, T xy, T xz, T yy, T yz, T zz)
	  {
		  (*this)(0,0) = xx;
		  (*this)(1,0) = xy;
		  (*this)(2,0) = xz;
		  (*this)(1,1) = yy;
		  (*this)(2,1) = yz;
		  (*this)(2,2) = zz;
	  }

  public: // operators

    /// get value at position (i,j).
    T operator()(std::size_t i, std::size_t j) const
    { return this->Base::operator()(i,j); }

    /// get non-const access to value at position (i,j).
    ValueProxy operator()(std::size_t i, std::size_t j)
    { return ValueProxy(i,j,*this); }

	  /// Set a value to the matrix.
    /// This function is used by the value proxy
    void set(const std::pair<std::size_t, std::size_t> & i, T value)
    {
      this->Base::operator()(i.first, i.second) = value;
      this->Base::operator()(i.second, i.first) = value;
    }

    T get(const std::pair<std::size_t, std::size_t> & i)
    { return this->Base::operator()(i.first, i.second); }
  };

  template < typename T>
  class SymMatrix3<T>::ValueProxy : public value_proxy<SymMatrix3<T>, std::pair<std::size_t,std::size_t> >
  {
  public: // typedefs
    typedef value_proxy<SymMatrix3<T>, std::pair<std::size_t,std::size_t> > Base;
    typedef typename Base::value_type value_type;
  public: // constructors
    ValueProxy(std::size_t i, std::size_t j, SymMatrix3<T> & m) :
          Base(std::make_pair(i,j), m)
    {      
    }
    void operator=(typename boost::call_traits<value_type>::param_type value)
    {
      this->Base::operator=(value);
    }
  };


  EXPIMP_TEMPLATE template class TIL_API SymMatrix3<double>;


//} // namespace linalg
} // namespace til

#ifdef _MSC_VER
#pragma warning(pop)
#endif

// Package includes
#include "til/symmatrix3tools.h"

#endif

