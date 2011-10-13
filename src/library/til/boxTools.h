#ifndef TIL_BOXTOOLS_H
#define TIL_BOXTOOLS_H

/// \File Belongs to package Box
/// Do not include directly, include til/Box.h instead.


namespace til
{

  /// Check whether a point lies within box
  template < typename T, typename TStorage, std::size_t D >
  inline 
  typename boost::enable_if<is_numeric_container<TStorage>, bool>::type
  contains(const Box<T,D> & box, const TStorage & v)
  {
    return all_greater_equal(v, box.min_bounds()) && all_less_equal(v, box.max_bounds());
  }

  /// Check whether box1 contains box2
  template < typename T1, typename T2, std::size_t D >
  inline bool contains(const Box<T1,D> & boxOutside, const Box<T2,D> & boxInside)
  {
	  return (contains(boxOutside, boxInside.min_bounds()) &&
			  contains(boxOutside, boxInside.max_bounds()));
  }


  /// Computes the intersection of two boxes. The result is stored in
  /// the first input box.
  // NB: the range checking is done by the box object itself; if bounds are
  // not consistent, box throws an exception
  template < typename T, std::size_t D >
  Box<T,D> intersection(const Box<T,D> & b1, const Box<T,D> & b2)
  {
	  // NB: we first collect the numbers, and then assign them to the box,
	  // because using Box.setX/Y/ZMin/Max individually will most
	  // of the time break the box consistency rule (having min's < max's)
	  // and hence throw an exception
    return Box<T,D>(
      max(b1.min_bounds(), b2.min_bounds()),    
      min(b1.max_bounds(), b2.max_bounds())
      );
  }

  /// Standard operator<< to print the coordinates of a Box object.
  template < typename T, std::size_t D >
  std::ostream& operator<<(std::ostream& os, const Box<T,D> &box)
  {
	  return os << "(" __ box.min_bounds()  __ "," __ box.max_bounds() << " )" << std::endl;
  }


  /// Return the size of a box
  template < typename T, std::size_t D>
  numeric_array<T,D> size(const Box<T,D> & box)
  {
	  return box.max_bounds() - box.min_bounds();
  }

  /// Return the volume of a box
  // TODO: return type might not be appropriate (e.g. char)
  template < typename T, std::size_t D >
  T volume(const Box<T,D> &box)
  {
    T res = 1;
    for (std::size_t i = 0; i < D; ++i)
    {
      res *= box.max_bounds()[i] - box.min_bounds()[i];
    }
    return res;
  }

} // namespace til


#endif

