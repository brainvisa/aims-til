#ifndef TIL_AFFINE_MAP_H
#define TIL_AFFINE_MAP_H

// includes from TIL library
#include "til/til_common.h"
#include "til/AffineTransform.h"
#include "til/labels.h"
#include "til/Matrix3.h"
#include "til/numeric_array.h"
#include "til/traits.h"

namespace til
{

	/// Defines an affine coordinate transform.
	// NB: It is not clear yet how to disambiguate between world coord and
	// voxel coords.
	template < typename T >
	class AffineMap : public Mapping_label, public detemplated_functor_label
	{
  public: // typedefs

    template < typename X >
    struct TypeStruct
    {
      typedef X Type;
    };

    typedef typename std::unary_function<numeric_array<T,3>, numeric_array<T,3> >::argument_type argument_type;
    typedef typename std::unary_function<numeric_array<T,3>, numeric_array<T,3> >::result_type result_type;
  
	public: // constructors & destructor
		AffineMap() {};
		AffineMap(const Matrix3<T> &m, const numeric_array<T,3> &transl)
			: m_data(m, transl) {};

  public: // set & get
  
    const Affine<T> & transfo() const { return m_data; }
    Affine<T> & transfo() { return m_data; }

	public: // functions

    // removing this, because I want to have a real functor here... otherwise it
    // starts being a mess in template expression... unless I decide to have the equivalent
    // of std::unary/binary_function for detemplated functors...
    /*
    template < typename V >
    inline
    Vector<typename combine<T,V>::type, 3>
    operator()(const Vector<V,3> &v) const
    {
      return m_data*v;
			//return (*static_cast<Affine<T>*>(this))*v;
			//return (*this)*v;
    }
    */
    
    template < typename X >
    inline X
    operator()(const X & v) const { return m_data * v; }
    
    //inline Vector<T,3> operator()(const Vector<T,3> & v) const { return m_data*v; }
    //inline Point<T,3> operator()(const Point<T,3> & v) const { return m_data*v; }
    
		/*
		INLINE
		Vector<double,3>
		operator()(const Vector<int,3> &v) const
		{
			return static_cast<Affine<T> >(*this)*v;
		}
		*/
    
    private: // data
    
      Affine<T> m_data;
	};

}


#endif

