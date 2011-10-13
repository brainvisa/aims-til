#ifndef TIL_CONDITIONAL_ITERATOR_H
#define TIL_CONDITIONAL_ITERATOR_H

// includes from TIL library
#include "til/til_common.h"

// namespace
namespace til
{
	/// Conditional iterator for images.
	/// As opposed to a standard iterator, a conditional iterator takes a
	/// boolean functor as a parameter, and will iteratate for those elements
	/// that returns true.
	/// The ConditionalIterator class actually contains only this mechanism;
	/// all the iterator mechanism itself should be provided by an iterator
	/// from which ConditionalIterator derives, and which is a parameter of
	/// templation. So ConditionalIterator just provides this functionality
	/// to already existing iterators.
	/// Example: itlin(im, _1 >0) provides a linear iterator for image im
	/// for voxels that are strictly positive.

	template < class TIterator, class BoolFunctor >
	class ConditionalIterator : public TIterator
	{
	public: // constructors & destructor
		
		/// Constructor for const iterators
		ConditionalIterator(const typename TIterator::TImage &im, const BoolFunctor &boolFunctor = BoolFunctor()) : TIterator(im), m_boolFunctor(boolFunctor) {};
		
		/// Constructor for non-const iterators
		ConditionalIterator(typename TIterator::TImage &im, const BoolFunctor &boolFunctor = BoolFunctor()) : TIterator(im), m_boolFunctor(boolFunctor) {};

	public: // functions

		/// Go to the next element if possible.
		/// Note that with conditional iterators the concept of a for loop, with
		/// increase and test at a different place, is awkward.
		bool next()
		{
			do
			{
				// iterate if you can, but stop if you gotta stop!
				if (!this->TIterator::next()) return false;
			}
			// Not a positive point? Then iterate again
			while (!m_boolFunctor((*this)));
			
			// Okay, we found a positive point, and we are not at the end
			return true;
		}

	private:
		BoolFunctor m_boolFunctor;
	};
}

#endif

