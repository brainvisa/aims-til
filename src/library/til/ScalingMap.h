#ifndef TIL_SCALING_MAP_H
#define TIL_SCALING_MAP_H

// includes from TIL library
#include "til/numeric_array.h"
#include "til/scaling.h"

namespace til
{

	template < typename T >
	class ScalingMap : public Scaling<T>, public Mapping_label
	{
	public: // constructors & destructor

		ScalingMap() {}
		ScalingMap(const numeric_array<T,3> & scale, const numeric_array<T,3> & transl)
			: Scaling<T>(scale, transl) {}


	public: // functions

		template < typename V >
		INLINE
		numeric_array<typename combine<T,V>::type, 3>
		operator()(const numeric_array<V,3> &v) const
		{
      // TODO: EXTREMELY INEFFICIENT, I WANT TO THROW ALL THIS VECTOR/POINT SHIT
      // OUT OF THE WINDOW!!!!!!!!!!
      return numeric_array<typename combine<T,V>::type, 3>(static_cast<Scaling<T> >(*this) * v);
			//return Point<typename combine<T,V>::type, 3>(static_cast<Scaling<T> >(*this)*v.data());
			//return (*this)*v;
		}
	};

} // namespace til

// Package includes
#include "scalingMapTools.h"

#endif


