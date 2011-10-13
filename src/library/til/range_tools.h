#ifndef TIL_RANGE_TOOLS_H
#define TIL_RANGE_TOOLS_H

/// \file Belongs to Range package. Do not include directly, include "til/Range.h" instead

// includes from STL
#include <cassert>

namespace til
{
  /// Set range using its center and half sizes.
  template <typename T, std::size_t D>
  void setCenterAndHalfSizes
  (
    Range<T,D> & range,
    const numeric_array<T,D> & center,
    const numeric_array<T,D> & halfSize
  )
  {
    range.set_bounds(center - halfSize, center + halfSize);
  }

 /// Set range using its center and its size.
  template <typename T, std::size_t D>
  void setCenterAndSizes
  (
    Range<T,D> & range,
    const numeric_array<T,D> & center,
    const numeric_array<T,D> & size
  )
  {
    for (int i = 0; i < D; ++i) { assert(size[i]%2); }
    range.set_bounds(center - size/2, center + size/2);
  }


} // namespace til

#endif

