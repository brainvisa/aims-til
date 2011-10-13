#ifndef TIL_AFFINE_H
#define TIL_AFFINE_H

// include from TIL library
#include "til/til_common.h"
#include "til/labels.h"
#include "til/Matrix3.h"
#include "til/numeric_array.h"
#include "til/traits.h"

// Namespace 
namespace til
{
	//namespace math {

	/// A affine mathematical object.
	/// Note that this object represents a mathematical object, *not* a
	/// coordinate transform between images. Indeed, a mathematical affine
	/// transform is completely independant of images and knows nothing about
	/// them. On the other hand, an affine coordinate transform between images
	/// has to be able to deal with such things as world/volume coordinates, 
	/// and have standardized API with other types of mapping to be able to 
	/// template over such types.
	/// For affine coordinate transforms, see AffineMap
  
  // TODO: I don't see any reason to get stuck in 3D here
  
	template < typename T >
	class Affine
	{

	public: // constructors & destructor

		Affine() {};
		Affine(const Matrix3<T> & m, const numeric_array<T,3> & transl)
			: m_mat(m), m_transl(transl) {};


	public: // set & get

		const Matrix3<T> & getMatrix() const { return m_mat; }
		const numeric_array<T,3> & getTransl() const { return m_transl; }
		Matrix3<T> & getMatrix() { return m_mat; }
		numeric_array<T,3> & getTransl() { return m_transl; }

		void setMatrix(const Matrix3<T> &m) { m_mat = m; }
		void setTransl(const numeric_array<T,3> &transl) { m_transl = transl; }


	public: // functions

		void reset()
		{
			m_mat.reset();
			fill(m_transl, 0);
		}


	private: // data 

		Matrix3<T> m_mat;
		numeric_array<T,3> m_transl;
	};


  namespace functor
  {
    template < typename T1, typename T2 >
    struct Mul<Affine<T1>, numeric_array<T2,3> >
     : public std::binary_function<const Affine<T1> &, const numeric_array<T2,3>, numeric_array<typename combine<T1, T2>::type, 3> >
    {
      typedef std::binary_function<const Affine<T1> &, const numeric_array<T2,3>, numeric_array<typename combine<T1, T2>::type, 3> > Base;
      typedef typename Base::first_argument_type         first_argument_type;
      typedef typename Base::second_argument_type        second_argument_type;
      typedef typename Base::result_type                 result_type;
  
      Mul() : Base() {}
  
      result_type operator()(first_argument_type a, second_argument_type v) const
      {
        return a.getMatrix() * v + a.getTransl();
      }
    };
  }


  /// Multiplication between a 3D affine transform and a 3D vector
	template <typename T1, typename T2>
  INLINE
	numeric_array<typename combine<T1, T2>::type, 3>
	operator*(const Affine<T1> & a, const numeric_array<T2,3> & v)
	{
    return functor::Mul<Affine<T1>, numeric_array<T2,3> >()(a,v);
    //return a.getMatrix() * v + a.getTransl();
	}

  template <typename T1, typename T2>
  inline T2 operator*(const Affine<T1> & a, const T2 & v)
  {
    return a.getMatrix() * v + a.getTransl();
    //return a.getMatrix() * v + a.getTransl();
  }

  /*
  /// Mutliplication between a 3D affine transform and a 4D vector
  template <typename T1, typename T2>
  INLINE
  numeric_array<typename combine<T1, T2>::type, 4>
  operator*(const Affine<T1> &a, const numeric_array<T2, 4> &v)
  {
    return a.getMatrix() * v + a.getTransl();
  }
  */

	//} // namespace math
} // namespace til

// Package includes
#include "til/affineTransformTools.h"

#endif

