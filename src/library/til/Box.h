#ifndef TIL_BOX_H
#define TIL_BOX_H

// Standard library includes
#include <stdexcept>

// Local includes
#include "til/til_common.h"
#include "til/til_declarations.h"
#include "til/numeric_array.h"

// Ignore specific warnings
#ifdef _MSC_VER
#pragma warning (push)
#pragma warning (disable : 4231) // nonstandard extension used : 'extern' before template explicit instantiation
#endif // _MSC_VER

// Namespace 

namespace til {


  /// A 3D box parallel to canonical axes
  template < typename T, std::size_t D >
  class Box
  {
  public: // typedefs
    typedef Box<T,D> Self;

  public: // constructors & destructor

	  /// Default constructor. All coordinates are set to zero.
	  Box();

	  /// Set min and max bounds for all dimensions.
    Box(const numeric_array<T,D> & posMin, const numeric_array<T,D> & posMax);

  public: // functions
  	
    /// Get min bounds
    const numeric_array<T,D> & min_bounds() const { return m_minBounds; }
    /// Get max bounds
	  const numeric_array<T,D> & max_bounds() const { return m_maxBounds; }

	  // Set bounds
	  // NB: box bounds have to be consistent at all time, e.g. one cannot
	  // set a min x that is greater than the current max x.
    // This is why there is not direct writable access to bounds, we have to go through
    // those set functions
    
    /// Set min bounds.
    void set_min_bounds(const numeric_array<T,D> & min);
    /// Set max bounds.
    void set_max_bounds(const numeric_array<T,D> & min);
    /// Set min bound on axis i.
    void set_min_bound(std::size_t i, T value);
    /// Set max bound on axis i.
    void set_max_bound(std::size_t i, T value);
	  /// Set min and max bounds.
	  void set_bounds(const numeric_array<T,D> & minBounds, const numeric_array<T,D> & maxBounds);
    /// Set min and max bounds for axis i.
    void set_bounds(std::size_t i, T min, T max);

  private: // data

	  numeric_array<T,D> m_minBounds;
	  numeric_array<T,D> m_maxBounds;
  };

#ifdef TIL_EXPORT_SOME_CLASSES
  EXPIMP_TEMPLATE template class TIL_API Box<int,3>;
#endif

  template < typename T, std::size_t D >
  // This is naturally assuming that default constructor initialize stuff
  Box<T,D>::Box() : m_minBounds(), m_maxBounds()
  {
  }

  template < typename T, std::size_t D >
  Box<T,D>::Box(const numeric_array<T,D> & minBounds, const numeric_array<T,D> & maxBounds)
  {
      this->set_bounds(minBounds, maxBounds);
  }

  template < typename T, std::size_t D >
  void Box<T,D>::set_bounds(const numeric_array<T,D> & minBounds, const numeric_array<T,D> & maxBounds)
  {
    for (std::size_t i = 0; i < D; ++i) assert(minBounds[i] <= maxBounds[i]);
    m_minBounds = minBounds;
    m_maxBounds = maxBounds;
  }

  template < typename T, std::size_t D >
  void Box<T,D>::set_bounds(std::size_t i, T minBound, T maxBound)
  {
    assert(minBound <= maxBound);
    m_minBounds[i] = minBound;
    m_maxBounds[i] = maxBound;
  }

  template < typename T, std::size_t D >
  void Box<T,D>::set_min_bounds(const numeric_array<T,D> & minBounds)
  {
    for (int i = 0; i < D; ++i) assert(minBounds[i] <= m_maxBounds[i]);
    m_minBounds = minBounds;
  }

  template < typename T, std::size_t D >
  void Box<T,D>::set_max_bounds(const numeric_array<T,D> & maxBounds)
  {
    for (int i = 0; i < D; ++i) assert(maxBounds[i] >= m_minBounds[i]);
    m_maxBounds = maxBounds;
  }

  template < typename T, std::size_t D >
  void Box<T,D>::set_min_bound(std::size_t i, T value)
  {
    assert(value <= m_maxBounds[i]);
    m_minBounds[i] = value;
  }

  template < typename T, std::size_t D >
  void Box<T,D>::set_max_bound(std::size_t i, T value)
  {
    assert(value >= m_minBounds[i]);
    m_maxBounds[i] = value;
  }

} // namespace til

#ifdef _MSC_VER
#pragma warning (pop)
#endif

// Package includes
#include "til/boxTools.h"

#endif

