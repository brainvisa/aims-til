#ifndef TIL_CONSTIMAGECVOLUMETRICITERATOR_H
#define TIL_CONSTIMAGECVOLUMETRICITERATOR_H

// includes from TIL library
#include "til/labels.h"
#include "til/numeric_array.h"
#include "til/Range.h"
#include "til/miscTools.h"


// Namespace 
namespace til {


  // Predeclaration
  template < typename T > class ImageC;

  // A class to iterate through all values of an image while knowing the
  // image coordinates of the current pixel
  // Allows to read but not to write to image
  template < typename T >
  //class ConstImageCVolumetricIterator : public VolumetricImageIterator_label
  class ConstVolumetricIterator<ImageC<T> > : public VolumetricImageIterator_label
  {
  public: // typdefs

	  typedef T			      value_type;
	  typedef const T &	  reference;
	  typedef ImageC<T>	  TImage;

  public: // constuctors & destructors

	  // default constructor
    ConstVolumetricIterator<ImageC<T> >() { m_index = 0; m_end = 0; }

	  // Iterates on the whole image
    ConstVolumetricIterator<ImageC<T> >(const ImageC<T> & im) : m_im(const_cast<ImageC<T>&>(im)) { this->init(); }

	  // Iterates only on the region of interest specified by the box
	  ConstVolumetricIterator<ImageC<T> >(const ImageC<T> & im, const Range<int,3> & roi) : m_im(const_cast<ImageC<T>&>(im)) { this->init(roi); }

	  virtual ~ConstVolumetricIterator<ImageC<T> > () {}


  public: // initialization

	  virtual void init();
	  virtual void init(const Range<int,3> &box);


  public: // set & get

	  // set current position
	  //void setPos(int x, int y, int z);
	  void set_pos(const numeric_array<int,3> &pos);// { this->setPos(EXPAND_VECTOR(pos)); }

	  //void setUnsafePos(int x, int y, int z);
	  void setUnsafePos(const numeric_array<int,3> &pos);// { this->setUnsafePos(EXPAND_VECTOR(pos)); }

    /*
	  // get current position
	  int getX() const { return m_pos.getX(); }
	  int getY() const { return m_pos.getY(); }
	  int getZ() const { return m_pos.getZ(); }
    */
	  //int getPos(int i) const { return m_pos.get(i); }
	  const numeric_array<int,3> & pos() const { return m_pos; }

	  // Get the region of interest
	  const Range<int,3> & roi() const { return m_roi;}

	  // Get current image
	  const ImageC<T> & image() const { return m_im; }


  public: // functions

	  /// Get value of a neighbor WITHOUT RANGE CHECKING.
	  /// NB: the offset is passed, not the actual position
	  INLINE T getUnsafeValue(int offsetx, int offsety, int offsetz) const
	  {
		  return *(m_index + offsetx + m_im.dim()[0]*(offsety + m_im.dim()[1]*offsetz));
		  //return m_im.getUnsafeValue(
		  //	this->getX() + offsetx,
		  //	this->getY() + offsety,
		  //	this->getZ() + offsetz);
	  }

	  /// Get value of a neighbor WITHOUT RANGE CHECKING.
	  /// NB: the offset is passed as a template parameter
	  template < int offsetx, int offsety, int offsetz >
	  INLINE T getUnsafeValue() const
	  {
		  return *(m_index + offsetx + m_im.dim()[0]*(offsety + m_im.dim()[1]*offsetz));
	  }

	  /// Get value of a neighbor WITHOUT RANGE CHECKING.
	  /// NB: the offset is passed, not the actual position
	  INLINE T getUnsafeValue(const numeric_array<int,3> &offset) const
	  {
		  //return this->getUnsafeValue(EXPAND_VECTOR(offset));
      // TODO: a vector-vector functor could actually take care of this
      return *(m_index + offset[0] + m_im.dim()[0]*(offset[1] + m_im.dim()[1]*offset[2]));
	  }

    
	  template < class Extrapolator, int offsetx, int offsety, int offsetz >
	  INLINE T getValue() const
	  {
		  if (containsNeighbor<offsetx, offsety, offsetz>(*this))
		  {
			  return this->getUnsafeValue<offsetx, offsety, offsetz>();
		  }
		  else
		  {
        return Extrapolator::getExtrapolatedValue(m_im, this->pos() + til::numeric_array<int,3>(offsetx, offsety, offsetz));
/*				  this->pos()[0] + offsetx,
				  this->pos()[1] + offsety,
				  this->pos()[2] + offsetz);*/
		  }
	  }

	  /// Get value of a neighbor.
	  /// Extrapolation is done if the position lies beyond image range.
	  /// NB: the offset is passed, not the actual position
    /*
    template < class Extrapolator >
	  T getValue(int offsetx, int offsety, int offsetz) const
	  {
		  // NB: here we keep this test and do not call directly
		  // Extrapolator::getValue because we can access the neighbors
		  // faster here than using im.getValue
		  if (contains(m_im,
			  this->getX() + offsetx,
			  this->getY() + offsety,
			  this->getZ() + offsetz))
		  {
			  return this->getUnsafeValue(offsetx, offsety, offsetz);
		  }
		  else
		  {
			  return Extrapolator::getExtrapolatedValue(m_im, 
				  this->getX() + offsetx,
				  this->getY() + offsety,
				  this->getZ() + offsetz);
		  }
	  }
    */

	  /// Get value of a neighbor.
	  /// Extrapolation is done if the position lies beyond image range.
	  /// NB: the offset is passed, not the actual position
	  template < class Extrapolator >
	  T getValue(const numeric_array<int,3> & offset) const
	  {
		  // NB: here we keep this test and do not call directly
		  // Extrapolator::getValue because we can access the neighbors
		  // faster here than using im.getValue
		  if (contains(m_im, m_pos + offset))
		  {
			  return this->getUnsafeValue(offset);
		  }
		  else
		  {
			  return Extrapolator::getExtrapolatedValue(m_im, m_pos + offset);
		  }
		  //return this->getValue<Extrapolator>(EXPAND_VECTOR(offset));
	  }

	  // Get value of a neighbor
	  // NB: the offset is passed, not the actual position
	  // TODO: hey, why do we need this operator for, with all the
	  // previous methods?
    /*
	  T operator()(int offsetx, int offsety, int offsetz) const
	  {
		  return m_im.getValue(
			  this->getX() + offsetx,
			  this->getY() + offsety,
			  this->getZ() + offsetz);
	  }
    */

	  INLINE T operator()(const numeric_array<int,3> &offset) const
	  {
		  //return this->operator()(EXPAND_VECTOR(offset));
      return m_im.getValue(m_pos + offset);
	  }

	  /// Go to next image element
	  INLINE void operator++();

	  /// Go to next image element, return whether operation succeeded
	  INLINE bool next();

	  // Print some info about the object
	  // friend void printInfo FRIEND_TEMPLATE_NO_ARG (ConstImageCVolumetricIterator<T> &);

	  /// Go to next element in the following direction
	  INLINE void next(ImageAxis axis);

	  /// Test whether the iterator has reached the end of the image or not
	  INLINE bool isAtEnd() const { return (m_index == 0); }

	  /// Return the value of the current element
	  reference operator*() const { return *m_index; }


  public: // set & get
  //protected: // set & get

	  T * getIndex() const { return m_index; }


  private: // functions

	  void _next();

  private: // data

	  // Current position in the volume
    numeric_array<int,3> m_pos;

	  // Precomputed offsets to move in each direction
    numeric_array<int,3> m_offset;
  	
	  // ImageC
	  ImageC<T> & m_im;

	  // Region of interest in the image
	  Range<int,3> m_roi;

	  // Pointer to current element	
	  T* m_index;

	  // Pointer on the last element
	  T* m_end;
  };

  template < typename T >
  void ConstVolumetricIterator<ImageC<T> > ::init(const Range<int,3> &box)
  {
	  /*
	  if (!isAllocated(im))
	  {
		  m_index = 0;
		  m_end = 0;
		  return;
	  }
	  */
  	
	  if (!contains(getRange(m_im), box))
	  {
		  throw std::domain_error("ROI lies outside image range");
	  }

	  m_roi = box;
	  m_pos = box.min_bounds();

	  m_index = m_im.getUnsafePointerOf(m_pos);
	  m_end = m_im.getUnsafePointerOf(numeric_array<int,3>(box.max_bounds()));

	  m_offset[0] = 1;
	  m_offset[1] = m_im.dim()[0];
	  m_offset[2] = m_im.dim()[0]*m_im.dim()[1];
  }


  template < typename T >
  void ConstVolumetricIterator<ImageC<T> > ::init()
  {
	  /*
	  if (!isAllocated(im))
	  {
		  m_index = 0;
		  return;
	  }
	  */
	  this->init(getRange(m_im));
  }


  template < typename T >
  void ConstVolumetricIterator<ImageC<T> >::_next()
  {
	  m_pos[0] = m_roi.min_bounds()[0];
  	
	  if (++(m_pos[1]) > m_roi.max_bounds()[1])
	  {
		  m_pos[1] = m_roi.min_bounds()[1];
  		
		  if (++(m_pos[2]) > m_roi.max_bounds()[2])
		  {
			  m_index = 0;
			  //return *this;
			  return;
		  }
	  }
  	
	  // TODO: to be faster, use offsets to jump in y and z directions
	  m_index = m_im.getUnsafePointerOf(m_pos);
  }

  template < typename T >
  //INLINE ConstImageCVolumetricIterator<T> & 
  INLINE void
  ConstVolumetricIterator<ImageC<T> > ::operator++()
  {
	  // We hit a border

	  if (++(m_pos[0]) > m_roi.max_bounds()[0])
	  //if (++(m_pos.m_values[0]) > m_roi.m_posMax.m_values[0])
	  //if (++(m_pos.m_values[0]) > m_roi.getXMax())
	  {
		  this->_next();
	  }
	  // We're still in the volume

	  else
	  {
		  ++m_index;
	  }

	  //return *this;
  }

  template < typename T >
  void ConstVolumetricIterator<ImageC<T> >::set_pos(const numeric_array<int,3> & pos)
  {
    if (!contains(m_roi, pos))
      throw std::out_of_range("Point does not lie within iterator range");
    this->setUnsafePos(pos);
  }

  /*
  template < typename T >
  void ConstVolumetricIterator<ImageC<T> > ::setPos(int x, int y, int z)
  {
	  if (!contains(m_roi, x, y, z))
	  {
		  throw std::out_of_range("Point does not lie within iterator range");
	  }

	  this->setUnsafePos(x,y,z);
  }
  */

  template < typename T >
  inline void ConstVolumetricIterator<ImageC<T> >::setUnsafePos(const numeric_array<int,3> & pos)
  {
    m_pos = pos;
    m_index = m_im.getUnsafePointerOf(m_pos);
  }

  /*
  template < typename T >
  inline void ConstVolumetricIterator<ImageC<T> > ::setUnsafePos(int x, int y, int z)
  {
	  m_pos.setX(x);
	  m_pos.setY(y);
	  m_pos.setZ(z);

	  m_index = m_im.getUnsafePointerOf(m_pos);
  }
  */

  template < typename T >
  INLINE bool
  ConstVolumetricIterator<ImageC<T> >::next()
  {
	  if (++(m_pos[0]) > m_roi.max_bounds()[0])
	  {
		  m_pos[0] = m_roi.min_bounds()[0];
		  if (++(m_pos[1]) > m_roi.max_bounds()[1])
		  {
			  m_pos[1] = m_roi.min_bounds()[1];
			  if (++(m_pos[2]) > m_roi.max_bounds()[2])
			  {
				  m_index = 0;
				  return false;
			  }
		  }
		  // TODO: to be faster, use offsets to jump in y and z directions
		  m_index = m_im.getUnsafePointerOf(m_pos);
	  }
	  else
	  {
		  ++m_index;
	  }
	  return true;
  }

  template < typename T >
  INLINE void
  ConstVolumetricIterator<ImageC<T> >::next(ImageAxis axis)
  {
	  if (++m_pos[axis] > m_roi.max_bounds()[axis])
	  {
		  // We hit a border
		  m_index = 0;
	  }
	  else
	  {
		  // We're still in the volume
		  m_index += m_offset[axis];
	  }
  }


  template < typename T >
  void printInfo(ConstVolumetricIterator<ImageC<T> >  &it)
  {
	  std::cout << "Pointer: " << it.m_index << std::endl;
	  std::cout << "Position: " << it.m_pos << std::endl;
	  std::cout << "Offset: " << it.m_offset << std::endl;
	  std::cout << "Range: " << it.m_roi.min_bounds() __ it.m_roi.max_bounds() << std::endl;	
  }

} // namespace til

#endif
