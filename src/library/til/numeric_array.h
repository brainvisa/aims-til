#ifndef TIL_NUMERIC_ARRAY_H
#define TIL_NUMERIC_ARRAY_H

// includes from BOOST
#include "boost/array.hpp"
#include "boost/call_traits.hpp"

// includes from TIL
#include "til/basic_iterator.h"
#include "til/basic_range.h"
#include "til/cat2type.h"
#include "til/labels.h"
#include "til/std_wrap.h"
#include "til/traits.h"

// Self mental note: the idea here is to have an extra layer over
// standard containers. In this layer, we define all the operators
// we want. The idea here is that choices of a best iterating method
// for a given storage type and a given mathematical operation are done
// here, so that at the upper levels, no iterators are used.

// this may seem odd, but actually I've come to think that I come to
// a point where genericity has a limit, and this limit is in the flexibility.
// Clearly, sparse containers must have a standard face, that can fill in
// the gaps and produce a value for every index. Still, these iterators
// must certainly not be used when operating on them. The efficient way 
// would be to use a sparse iterator -- and then, we would not take
// long before realizing that all the for loops, ++i1, ++i2, that are
// at the core of stl algorithms, won't do. This would be true for
// RLE as well.

// In short: for any storage strategy, we are left redesigning not only
// the iterators, but also the loops, because more often than not, we
// cannot simply ++i1, ++i2 when operating on two containers. This means
// that unfortunately, there is no generic algorithms that would actually
// do for us.

// This becomes even worse when storage types are mixed.
// Consider a sparse matrix and a non-sparse matrix: what iterators to
// use to compute an operation between them? Well, that depends: if
// it is a multiplication, probably


namespace til
{

  //--------------------------------------------------------------------------

    //----------------------//
   //  numeric_array_impl  //
  //----------------------//

  namespace detail
  {
    /// An implementation of a numeric array.
    template < typename T, std::size_t D >
    class numeric_array_impl
      : public boost::array<T,D>
      , public numeric_container_label
    {
    public: // typedefs
    
      typedef numeric_array_impl<T,D>                 Self;
      typedef boost::array<T,D>                       Base;
      typedef typename Base::reverse_iterator         reverse_iterator;
      typedef typename Base::const_reverse_iterator   const_reverse_iterator;

    public: // classes

      class iterator;
      class const_iterator;
      class range;
      class const_range;

    public: // constructors & destructor
       
      /// Default constructor
      numeric_array_impl();
      
      /// Standard constructors.
      /// NB: obviously, the size is irrelevant here; however, we need these for genericity
      //explicit numeric_array_impl(std::size_t d);
      
      //numeric_array_impl(std::size_t d, typename boost::call_traits<T>::param_type t);

      /*
      /// Conversion from another type.
      template < typename X >
      explicit numeric_array_impl(const X & v) { functor::CastTo<Self,X>()(*this, v); }
      */

    public: // initialization

      /*
      template < typename X >
      inline void operator=(const X & v)
      {
        functor::CastTo<X,Self>()(v,*this);
      }
      */

    public: // standard interface

      /// Return container size.
      std::size_t size() const { return D; }

    public: // iterators

      iterator begin() { return iterator(this->Base::begin()); }
      iterator end()   { return iterator(this->Base::end()); }
      const_iterator begin() const { return const_iterator(this->Base::begin()); }
      const_iterator end()   const { return const_iterator(this->Base::end()); }

      reverse_iterator rbegin() { return this->Base::rbegin(); }
      reverse_iterator rend()   { return this->Base::rend(); }
      const_reverse_iterator rbegin() const { return this->Base::rbegin(); }
      const_reverse_iterator rend()   const { return this->Base::rend(); }

    public: // range

      range whole_range() { return range(this->begin(), this->end()); }
      const_range whole_range() const { return const_range(this->begin(), this->end()); }

    public: // set & get
      
      /// Read access to n-th value.
      inline const T & operator[](std::size_t n) const
      {
        assert(n < D);
        return this->elems[n];
      }
      /// Read-write access to n-th value.
      inline T & operator[](std::size_t n)
      {
        assert(n < D);
        return this->elems[n];
      }

    public: // artithmetic operators

      // arithmetic with other vectors
  #define TIL_DEFINE(op)                                                \
      template < typename T2 >                                          \
      void operator op (const numeric_array_impl<T2,D> & v)             \
      { for (std::size_t i = 0; i < D; ++i) (*this)[i] op v[i]; }       \
      
      TIL_DEFINE(+=)
      TIL_DEFINE(-=)
      TIL_DEFINE(*=)
      TIL_DEFINE(/=)

  #undef TIL_DEFINE
      
      // arithmetic with constants
      // deliberately not templated!
      void operator+=(const T & x)
      { for (std::size_t i = 0; i < D; ++i) (*this)[i] += x; }
      void operator-=(const T & x)
      { for (std::size_t i = 0; i < D; ++i) (*this)[i] -= x; }
      void operator*=(const T & x)
      { for (std::size_t i = 0; i < D; ++i) (*this)[i] *= x; }

      void operator/=(const T & x)
      { this->divide(x, cat2type<is_floating_point, T>()); }
      //{ this->operator*=(1/x); }

    private: // functions
    
      /// optimization for floating points: we multiply with the inverse rather than dividing.
      template < typename X >
      void divide(const X & x, label::Passed<is_floating_point>)
      { this->operator*=(1/x); }

      /// for non-floating point values, we apply the regular division.
      template < typename X >
      void divide(const X & x, label::Failed<is_floating_point>)
      { for (std::size_t i = 0; i < D; ++i) (*this)[i] /= x; }


      /*      
      void operator/=(const T & v)
      { this->divide(v, cat2type<is_integer,T>()); }
      
    private: // functions

    
      void divide(const T & x, label::Passed<is_integer>)
      {
      }

      void divide(const T & x, label::Failed<is_integer>)
      {
      }
      */
    };


    //--------------------------------------------------------------------------

    /*
    template < typename T, std::size_D >
    inline typename numeric_array<T,D>::iterator
    begin(numeric_array<T,D> & x)
    { return typename numeric_array<T,D>::iterator(x.begin()); }

    template < typename T, std::size_D >
    inline typename numeric_array<T,D>::const_iterator
    begin(const numeric_array<T,D> & x)
    { return typename numeric_array<T,D>::iterator(x.begin()); }
    */

    //--------------------------------------------------------------------------


    template < typename T, std::size_t D >
    class numeric_array_impl<T,D>::const_iterator
      : public basic_iterator<const T>
    {
    public: // typedefs
      typedef basic_iterator<const T> Base;
    public: // constructors & destructor
      const_iterator() : Base() {}
      const_iterator(const T * p) : Base(p) {}
      const_iterator(iterator i) : Base(i) {}
    };


    //--------------------------------------------------------------------------


    template < typename T, std::size_t D >
    class numeric_array_impl<T,D>::iterator
      : public basic_iterator<T>
    {
    public: // typedefs
      typedef basic_iterator<T> Base;
    public: // constructors & destructor
      iterator() : Base() {}
      iterator(T* p) : Base(p) {}
    };

    //-------------------------------------------------------------------------

    template < typename T, std::size_t D >
    class numeric_array_impl<T,D>::range
      : public basic_range<typename numeric_array_impl<T,D>::iterator>
    {
    public: // typedefs
      typedef typename numeric_array_impl<T,D>::iterator iterator;
      typedef basic_range<iterator> Base;
    public: // constructors
      range(iterator begin, iterator end) : Base(begin, end) {}
    };

    template < typename T, std::size_t D >
    class numeric_array_impl<T,D>::const_range : public basic_range<const_iterator>
    {
    public: // typedefs
      typedef basic_range<const_iterator> Base;
    public: // constructors
      const_range(const_iterator begin, const_iterator end) : Base(begin, end) {}
    };

    //-------------------------------------------------------------------------
  
  } // namespace detail


  //---------------------------------------------------------------------------

  template < typename T, std::size_t D >
  class numeric_array
    : public detail::numeric_array_impl<T,D>
  {
  public: // typedefs
    typedef detail::numeric_array_impl<T,D> Base;
  public: // constructors
    numeric_array() : Base() {}
    //explicit numeric_array(std::size_t d) : Base(d) {}
    //numeric_array(std::size_t d, T t) : Base(d,t) {}
    numeric_array(const numeric_array<T,D> & x) : Base(x) {}
    //template < typename X >
    //numeric_array(const X & x) : Base(x) {}
  /*
  public: // operators
    template < typename X >
    void operator=(const X & x) { this->Base::operator=(x); }
    */
  };

  //---------------------------------------------------------------------------

  template < typename T >
  class numeric_array<T,3>
    : public detail::numeric_array_impl<T,3>
  {
  public: // typedefs
    typedef detail::numeric_array_impl<T,3> Base;
  public: // constructors

    // TODO: find a way that would enable us to initialize like a({2,3}).
    // That would avoid overloading problems, esp. for 2D vectors, with
    // the "standard" form COnstructor(size, value).

    numeric_array() : Base() {}
    //explicit numeric_array(std::size_t d) : Base(d) {}
    //numeric_array(std::size_t d, T t) : Base(d,t) {}
    //template < typename X >
    //numeric_array(const X & x) : Base(x) {}
    numeric_array(T x, T y, T z) : Base() { this->init(x,y,z); }
    numeric_array(const numeric_array<T,3> & x) : Base(x) {}
  public: // initialization
    void init(T x, T y, T z) { (*this)[0] = x; (*this)[1] = y; (*this)[2] = z; }
  /*
  public: // operators
    template < typename X >
    void operator=(const X & x) { this->Base::operator=(x); }
    */
  };

  //---------------------------------------------------------------------------

} // namespace til

// Package includes
#include "til/numeric_array.tpp"
#include "til/numeric_array_operators.h"
#include "til/numeric_array_tools.h"


#endif

