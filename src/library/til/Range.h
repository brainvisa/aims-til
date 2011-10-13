#ifndef TIL_RANGE_H
#define TIL_RANGE_H

// includes from TIL library
#include "til/til_declarations.h"
#include "til/Box.h"

namespace til
{

  template < typename T, std::size_t D >
  class rectangular_range
  {
  private: // typedefs
    typedef numeric_array<T,D> Vector;
  public: // constructors
    rectangular_range() : m_min(), m_max() {}
    rectangular_range(const Vector & min, const Vector & max)
      : m_min(min), m_max(max) {}
  public: // set & get
    Vector & min() { return m_min; }
    Vector & max() { return m_max; }
    const Vector & min() const { return m_min; }
    const Vector & max() const { return m_max; }
  private: // data
    Vector m_min;
    Vector m_max;
  };

  template < typename T, std::size_t D, typename TNDIterator >
  class rectangular_range_indicator
    : public rectangular_range<T,D>
  {
  public: // typedef
    typedef rectangular_range<T,D> Base;
    typedef typename Base::Vector Vector;
  public: // constructors
    rectangular_range_indicator(const Base & rec) : Base(rec), m_pos() {}
    rectangular_range_indicator(const Base & rec, const Vector & pos) : Base(rec), m_pos(pos) {}
    rectangular_range_indicator(const Vector & min, const Vector & max) : Base(min, max), m_pos() {}
    rectangular_range_indicator(const Vector & min, const Vector & max, const Vector & pos) : Base(min, max), m_pos(pos) {}
  public: // set & get
    Vector & pos() { return m_pos; }
  public: // operator
    void operator++()
    {
      if (++m_pos[0] < this->max()[0])
      {
        ++m_i;
        return;
      }
      m_pos[0] = this->min()[0];
      for (std::size_t i = 1; i < D; ++i)
      {
        if (++m_pos[i] < this->max()[i]) break;
        m_pos[i] = this->min()[i];
        set_coord(m_i, i, m_pos[i]);
      }
      m_i.from_pos(m_pos);
    }
  private: // data
    Vector m_pos;
    TNDIterator m_i;
  };


  /// Describe integer cube ranges (e.g. describing subimage).
  /// NB: Template parameter T has to be an integer.
  // TODO: rename as BoxRange or RectangularRange
  template < typename T, std::size_t D >
  class Range
    : public Box<T,D>
  {
  public: // typedefs
    typedef Range<T,D> Self;
    typedef Box<T, D> Base;

  public: // constructors & destructor

	  /// Default constructor
	  Range() : Base() {};

 	  /// Set range bounds in all dimensions
    Range(const numeric_array<T,D> & minBounds, const numeric_array<T,D> & maxBounds)
		  : Base(minBounds, maxBounds) {};

    /// Set max bounds, assuming min bounds are zeros.
    Range(const numeric_array<T, D> & maxBounds) : Base(numeric_array<T,D>(0,0,0), maxBounds) {}

  public: // set & get

    /*
	  /// Set range using its center and half sizes.
	  void setCenterAndHalfSizes(const numeric_array<T,D> & center, const numeric_array<T,D> & halfSizes);

	  /// Set range using its center and its size.
	  void setCenterAndSizes(const numeric_array<T,D> & center, const numeric_array<T,D> & halfSizes);
    */

    /// Get range size
	  numeric_array<T,D> dims() const
    { return this->max_bounds() - this->min_bounds() + 1; }
  };

  /*
  template <typename T, std::size_t D>
  void Range<T,D>::setCenterAndHalfSizes(const numeric_array<T,D> & center, const numeric_array<T,D> & halfSize)
  {
    this->set_bounds(center - halfSize, center + halfSize);
  }

  template <typename T, std::size_t D>
  void Range<T,D>::setCenterAndSizes(const numeric_array<T,D> & center, const numeric_array<T,D> & size)
  {
    for (int i = 0; i < D; ++i) { assert(size[i]%2); }
    this->set_bounds(center - size/2, center + size/2);
  }
  */

} // namespace til


// package include
#include "til/range_tools.h"

#endif

