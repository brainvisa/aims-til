#ifndef TIL_SCALING_H
#define TIL_SCALING_H

// includes from TIL library
#include "til/til_common.h"
#include "til/til_declarations.h"
#include "til/traits.h"

namespace til {

	/// Namespace for all linear algebra related stuff
	//namespace linalg {

		/// Represent a mathematical 3D scaling, i.e. a multiplication along
	/// all canonical axis, plus a translation
	template < typename T >
	class Scaling
	{
	public: // constructors & destructor

		Scaling() {}
		Scaling(const numeric_array<T,3> &scale, const numeric_array<T,3> &transl)
			: m_scale(scale), m_transl(transl) {}


	public: // set & get

		const numeric_array<T,3> & getScale() const { return m_scale; }
		const numeric_array<T,3> & getTransl() const { return m_transl; }
		numeric_array<T,3> & getScale() { return m_scale; }
		numeric_array<T,3> & getTransl() { return m_transl; }

		void setScale(const numeric_array<T,3> &scale) { m_scale = scale; }
		void setTransl(const numeric_array<T,3> &transl) { m_transl = transl; }


	private: // data

		numeric_array<T,3> m_scale;
		numeric_array<T,3> m_transl;
	};


	template < typename T1, typename T2 >
	INLINE
	numeric_array<T1,3>
	operator*(const Scaling<T1> &s, const numeric_array<T2,3> &v)
	{
		return mul<numeric_array<T1,3> >(s.getScale(), v) + s.getTransl();
	}

	//} // namespace linalg
} // namespace til

#endif

