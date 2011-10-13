#ifndef TIL_KDTREE_TPP_
#define TIL_KDTREE_TPP_

namespace til
{
  
  //----------------------------------------------------------------------------------
  
    //-----------------//
   //  KDTreeFactory  //
  //-----------------//
  
  template < typename TKDTree >
  void KDTreeFactory<TKDTree>::build(TKDTree & res)
  {
    // Create a vector containing the indices of the collection elements
    makeIndexVector(m_v, m_vIndex);

    // Set kdtree as a global variable inside the class
    m_pKdt = &res;

    // Start kdtree building iterations
    typename std::vector<index_type>::iterator left = m_vIndex.begin();
    typename std::vector<index_type>::iterator right = m_vIndex.end()-1;
    genKDIter(left, right, m_pKdt->addChild(0, 0u), 0);
  }
  
  
  template < typename TKDTree >
  void KDTreeFactory<TKDTree>::genKDIter
  (
    typename std::vector<index_type>::iterator left,
    typename std::vector<index_type>::iterator right,
    typename TKDTree::iterator iRes,
    int dim
  )
  {
    // We want the median to have a mean index.
    typename std::vector<index_type>::iterator median = left+(right-left)/2;
        
    // sort and iterate
    // NB: The right+1 is because STL algorithms expect to have an element out
    // of range as the "end" element.
    m_dim = dim;
    std::nth_element(left, median, right+1, m_compare);
    iRes->value() = *median;
    iRes->dim() = dim;

    // Iterate right and left, if range permits it
    if (left < median)
    {
      genKDIter(left, median-1,  m_pKdt->addChild(iRes, TKDTree::LEFT), (dim+1)%3);
    }
    if (median < right)
    {
      genKDIter(median+1, right, m_pKdt->addChild(iRes, TKDTree::RIGHT), (dim+1)%3);
    }
  }

  template < typename TKDTree >
  void KDTreeFactory<TKDTree>::makeIndexVector(const Collection & v, std::vector<index_type> & vIndex)
  {
    vIndex.resize(size(v));
    typename std::vector<indexed_type>::const_iterator iV = v.begin();
    typename std::vector<index_type>::iterator iVIndex = vIndex.begin();
    for (; iV != v.end(); ++iV, ++iVIndex)
    {
      getIndex(iV, *iVIndex);
      //*iVIndex = til::getIndex<index_type>(v, iV);
    }
  }
  
  
  //----------------------------------------------------------------------------------
  
    //----------------//
   //  Find_closest  //
  //----------------//
  
    
  template < typename TPrecision, typename TKDTree >
  typename Find_closest<TPrecision,TKDTree>::index_type
  Find_closest<TPrecision,TKDTree>::operator()(const indexed_type & vec)
  {
    this->init();
    m_vec = vec;
    lookMax(m_pKdt->root());
    return m_minIndex;
  }
  
  template < typename TPrecision, typename TKDTree >
  void Find_closest<TPrecision,TKDTree>::init()
  {
    m_min = std::numeric_limits<TPrecision>::max();
    m_minIndex = 0;
    m_niter = 0;
  }
  
  template < typename TPrecision, typename TKDTree >
  void Find_closest<TPrecision,TKDTree>::lookMax(typename TKDTree::const_iterator iNode)
  {
    ++m_niter;
    // Stop here if we reached end of tree
    // TODO: this test might be put below, to avoid an extra function call
    // everytime we reach an end point.
    if (iNode == 0) return;
    index_type node = *iNode;
    //NB: this was used when node is a pointer. Dunno if its removal changes something. hope not.
    //if (node == 0) return;
    
    // Compute the distance to the current node point    
    TPrecision dist2point = dist2(m_pKdt->get(node), m_vec, prec<TPrecision>());
    //std::cout << dist2point << "*" << std::flush;
    if (dist2point < m_min)
    {
      m_min = dist2point;
      m_minIndex = node;
    }

    // Continue looking downwards, starting from the lower half or the higher
    // half, depending on which the vector belongs to.
    TPrecision diffWall = m_vec[iNode->dim()] - m_pKdt->get(node)[iNode->dim()];

    if (diffWall < 0)
    {
      // Point belong to lower half: process this half first
      lookMax(iNode.child(TKDTree::LEFT));
      // We process the other half only if the min distance found so far
      // is larger than the distance to the halfplane border
      if (m_min > square(diffWall))
      {
        lookMax(iNode.child(TKDTree::RIGHT));
      }
    }
    else
    {
      // Point belong to higher half: process this half first
      lookMax(iNode.child(TKDTree::RIGHT));
      // We process the other half only if the min distance found so far
      // is larger than the distance to the halfplane border
      if (m_min > square(diffWall))
      {
        lookMax(iNode.child(TKDTree::LEFT));
      }
    }    
  }
  

} // namespace til

#endif /*KDTREE_TPP_*/
