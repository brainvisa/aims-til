#include "til/EquivalenceChain.h"




// Namespace 

namespace til {


EquivalenceChain::EquivalenceChain(int maxNLabels)
{
	if (maxNLabels < 1)
	{
		throw std::invalid_argument("Chain too small");
	}
	m_size = maxNLabels;
	m_chain = new int[m_size+1];
	this->reset();
}



EquivalenceChain::~EquivalenceChain() { delete [] m_chain; }



void EquivalenceChain::reset()
{
	for (int i = 0; i <= m_size; ++i) m_chain[i] = 0;
	m_maxLabel = 0;
}



void EquivalenceChain::reset(int nLabels)
{
    int i;
	for (i = 0; i <= nLabels; ++i) m_chain[i] = i;
	for (     ; i <= m_size ; ++i) m_chain[i] = 0;
	m_maxLabel = nLabels;
}



// Fill the chain starting from i with the last element of this chain
// so that all elements are now pointing to the last element


int EquivalenceChain::_fillWithLastLabel(int i)
{
	int pointedElem = m_chain[i];

	if (pointedElem == m_chain[pointedElem])
	{
		// The pointed element is the last element of the chain:
		// do nothing and return the value
		return pointedElem;
	}
	else 
	{
		// The element is not the last element
		// go to the next element of the chain
		// and assigne the current element to the last element

		return (m_chain[i] = _fillWithLastLabel(m_chain[pointedElem]));
	}
}


// Instead of pointing to the next element, each label points to the last element
// of its equivalence chain
// Equivalent of doing a i = findLastLabel(list, i), but (supposedly) faster


void EquivalenceChain::allPointToLastLabel()
{
	for (int i=m_maxLabel; i>=1; --i)
	//for (int i=1; i<=m_maxLabel; ++i)
	{
		_fillWithLastLabel(i);
	}
}



int EquivalenceChain::mergeLabels()
{
	
	this->allPointToLastLabel();
	
	int i, j, alt, neu;
	
	// Allocate new chain

	int* newChain = new int[m_size+1];

    for (i = 0; i <= m_size; ++i) newChain[i] = 0;

    neu = 1;
    for (i = 1; i <= m_maxLabel; ++i)
    {
        alt = m_chain[i];

        if (alt != 0)
        {
            for (j = i; j <= m_maxLabel; ++j)
			{
                if (m_chain[j] == alt)
                {
                    newChain[j] = neu;
					m_chain[j] = 0;
                }
			}
		++neu;
        }
    }
	
	delete [] m_chain;
	m_chain = newChain;

    return neu-1;
}



void print(const EquivalenceChain & ec)
{
	for (int i = 0; i <= ec.m_size; ++i)
	{
		std::cout << i << ":" << ec.m_chain[i] << " " << std::flush;
	}
}

} // namespace

