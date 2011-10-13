#ifndef TIL_RECURSIVE_FILTER_COMPOSITION_H
#define TIL_RECURSIVE_FILTER_COMPOSITION_H

// Local includes

#include "til/til_common.h"

namespace til {


	// Applies filter1, then filter2

	template < TRecursiveFilter1, TRecursiveFilter2 >
	class RecursiveFilterComposition
	{
	public: // constuctors & destructor

		RecursiveFilterComposition(const TRecursiveFilter1 &filter1, const TRecursiveFilter2 &filter2)
		{
			m_filter1 = filter1;
			m_filter2 = filter2;
		}

	public: // methods

		void apply(const T *in, T *out, int length) const
		{
			filter1.apply(in, out);
			filter2.apply(out, 
		}


	private: // data

		TRecursiveFilter1 m_filter1;
		TRecursiveFilter2 m_filter2;
	};


} // namespace

#endif

