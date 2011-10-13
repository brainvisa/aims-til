#ifndef TIL_EQUIVALENCECHAIN_H
#define TIL_EQUIVALENCECHAIN_H


// Standard library includes

#include <iostream>
#include <stdexcept>
#include <time.h>
#include <vector>


// Local includes
#include "til/til_common.h"



// Namespace 

namespace til {



// This class stores equivalences between labels, coded with type T,
// typically used in segmented images (of type T).
// Equivalence is e.g. saying "object with label X and object with
// label Y are equivalent" (the meaning of equivalency is of course
// left to the user, e.g. equivalence = "exact same object")

// Equivalence is stored a compact way, via a chain: a label points to
// an equivalent label.
// Memory usage is therefore very low (just N, the total # of labels)
// Update is reasonably quick (to add an equivalence between two labels,
// one has to link the end of the chain) if in practice
// equivalence are added between labels that are near the end of their
// respective chains

// However some operations are in general slow, e.g. "are label X and Y
// equivalent?" need to go to the end of the chain for both labels (much
// slower than matrix look-up of course).

// Label 0 is special and always represent the background. It cannot be
// part of any chain.


class TIL_API EquivalenceChain
{

public:

	// Allocate equivalence chain
	// The maximum number of labels cannot be higher than the
	// maximum of type T minus 1
	EquivalenceChain(int maxNumberOfLabels);


	// Destructor
	~EquivalenceChain();


	// clear the chain (everything points to background zero)
	void reset();


	// clear the chain, with 'nLabels' valid and unique labels
	// pointing to themselves
	void reset(int nLabels);


	// Set equivalence between label1 and label2
	void setEquivalence(int label1, int label2);


	// Create a new label value
	// Returns 0 (background) if impossible (chain is already full)
	int getNewLabel();


	// All equivalent labels now point to the same new label
	// Those new labels are compact, i.e. they form a sequence from
	// 1 to N (N being the number returned by the function)
	int mergeLabels();


	// Get label pointed by label n
	int operator[] (int n) { return m_chain[n]; }


	// Get the end of the chain which label belongs to
	int getLastLabel(int label)
	{
		if (label != m_chain[label])
		{
			return m_chain[label] = getLastLabel(m_chain[label]);
		}
		return label;
	}

/*
	int getLastLabel(int label)
	{
		//int count = 0;

		if (label < 0 || label > m_maxLabel)
		{
			throw std::out_of_range("Label out of range");
		}

		while (label != m_chain[label])
		{
			//++count;
			label = m_chain[label];
		}
		return label;
	}
*/

	int getGreatestLabel() { return m_maxLabel; }

public: // friends

	friend void print(const EquivalenceChain &);


private: // functions

	void allPointToLastLabel();
	int _fillWithLastLabel(int i);


private: // data

	// The current maximum label used
	int m_maxLabel;

	// The size of the equivalence chain
	int m_size;

	// The equivalence chain
	int * m_chain;
};




INLINE int EquivalenceChain::getNewLabel()
{
	if (m_maxLabel == m_size) return 0;
	++m_maxLabel;
	m_chain[m_maxLabel] = m_maxLabel;
	return m_maxLabel;
}


INLINE void EquivalenceChain::setEquivalence(int label1, int label2)
{

	/*
	if (label1 != label2)
	{
		std::vector<int> v1;
		this->getChain(label1, v1);
		std::vector<int> v2;
		this->getChain(label2, v2);

		int last = *(v1.end());

		std::vector<int>::iterator i;
		for (i=v1.begin(); i != v1.end(); ++i)
		{
			*i = last;
		}
		for (i=v2.begin(); i != v2.end(); ++i)
		{
			*i = last;
		}
	}
	*/
	if (label1 != label2)
	{
		int last1 = this->getLastLabel(label1);
		int last2 = this->getLastLabel(label2);

		if (last1 != last2)
		{
			m_chain[last1] = last2;
		}
	}
}


} // namespace

#endif

