#ifndef GHOST_ARRAY_H_
#define GHOST_ARRAY_H_

// include from BOOST
#include "boost/array.hpp"

namespace til
{
  
  template < typename T, std::size_t D >
  class ghost_array
  {
  public: // typedefs
  
    typedef T           value_type;
    typedef T &         reference;
    typedef const T &   const_reference;
    typedef T *         iterator;
    typedef const T *   const_iterator;
    
  public: //
    ghost_array(T * t) : m_(t) {}
    
    iterator begin() { return m_; }
    iterator end() { return m_ + D; }
    const_iterator begin() const { return m_; }
    const_iterator end() const { return m_ + D; }
    
    const_reference operator[](std::size_t i) const { return *(m_+i); }
    reference operator[](std::size_t i) { return *(m_+i); }
    
    operator boost::array<T,D> () const
    {
      boost::array<T,D> res;
      std::copy(this->begin(), this->end(), res.begin());
      return res;
    }
    
  private: // data
    T * m_;
  };
  
}

#endif /*GHOST_ARRAY_H_*/
