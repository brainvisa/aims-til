#ifndef TIL_IMAGE_BASE_H
#define TIL_IMAGE_BASE_H

// include from TIL library
#include "til/til_common.h"
#include "til/labels.h"
#include "til/numeric_array.h"

namespace til 
{

	/// type of voxel size
	typedef float t_voxsize;

	
	/// Collects common code accross all image classes.
	/// Namely, everything related to dimension, voxel size, space origin...
	/// Useless by itself
	class ImageBase : public Image_label
	{

	public: // set & get

		// Get image dimensions
    /*
		INLINE int dim()[0] const { return m_dim.getX();}		///< get X dimension
		INLINE int dim()[1] const { return m_dim.getY();}		///< get Y dimension
		INLINE int dim()[2] const { return m_dim.getZ();}		///< get Z dimension
    
		int getDim(int i) const { return m_dim.get(i); }		///< get i-th dimension
		const Vector<int,3> & getDim() const { return m_dim; }	///< get image dimension
    */

    const numeric_array<int,3> & dim() const { return m_dim; }	 ///< get image dimension


		/// Get total number of elements in image.
		int size() const { return this->dim()[0]*this->dim()[1]*this->dim()[2]; }


		// Get voxel size
    /*
		t_voxsize vdim()[0] const { return m_vDim.getX(); }				///< get X voxel size
		t_voxsize vdim()[1] const { return m_vDim.getY(); }				///< get Y voxel size
		t_voxsize vdim()[2] const { return m_vDim.getZ(); }				///< get Z voxel size
    
		t_voxsize vdim(int i) const { return m_vDim.get(i); }		///< get i-th voxel size
    */
		const numeric_array<t_voxsize,3> & vdim() const { return m_vDim; }	///< get voxel size


		// space origin    
		//const Vector<double,3> & origin() const { return m_ori; }		///< Get space origin.
		//Vector<double,3> & origin() { return m_ori; }		///< Get space origin.

    // TODO: shouldn't this be private?
		/// Set the voxel coordinates
		void set_vdim(const numeric_array<t_voxsize,3> &vdim)
		{
      assert(all_greater_equal(vdim, t_voxsize(0)));
      m_vDim = vdim;
		}

	protected: // methods

		void set_dim(const numeric_array<int,3> & dim)
		{
      assert(all_greater_equal(dim,0));
			m_dim = dim;
		}

		// Test whether these coordinates fall into image range
		/*
    bool contains(int i, int j, int k) const
		{
			return (
        i >= 0 && 
        j >= 0 && 
        k >= 0 && 
				i < m_dim[0] && 
				j < m_dim[1] &&
				k < m_dim[2]
        );
		}
    */
    bool contains(const numeric_array<int,3> & p) const
    {
      return all_greater_equal(p, 0) && all_less(p, m_dim);
    }

	private: // data

		// Image dimensions
		numeric_array<int,3> m_dim;

    // Voxel size
		numeric_array<t_voxsize,3> m_vDim;

    // TODO: this is crap actually. I mean, either you have a full matrix, that combines
    // voxel size and orientation and everything, or you don't. This half-assed origin
    // is useless because it doesn't account for axes flipping.
    // Origin
		//Vector<double,3> m_ori;
	};


	/// A trait class to assign iterators to image types
	template < class TImage >
	struct Iterator {};
	
} // namespace


#endif

