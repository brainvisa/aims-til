

namespace til
{
  //---------------------------------------------------------------------------
  
  template < typename TImage, typename TNeighborhood, typename TExtrapolator >
  int connectedComponents(TImage &im, const TNeighborhood &nh)
  {
    // typdefs
    typedef typename TImage::value_type value_type;
  
    // the "-1" is to be able to write tests for loops like
    // i<=maxLabel that can actually be false
    const int MAX_LABEL = (int) min(100000.0, double(std::numeric_limits<value_type>::max()) - 1);
  
    /*
    // Set first slice to background
    {
      Range range = getRange(im); range.setZMax(0);
      typename Iterator<TImage>::Volumetric iIm(im, range);
      for(; !iIm.isAtEnd(); ++iIm) { *iIm = 0; }
    }
    // Set top border to background
    {
      Range range = getRange(im); range.setYMax(0);
      typename Iterator<TImage>::Volumetric iIm(im, range);
      for(; !iIm.isAtEnd(); ++iIm) { *iIm = 0; }
    }
    // Set left border to background
    {
      Range range = getRange(im); range.setXMax(0);
      typename Iterator<TImage>::Volumetric iIm(im, range);
      for(; !iIm.isAtEnd(); ++iIm) { *iIm = 0; }
    }
    */
  
    value_type neighborLabels[3*3*3]; // contains the labels of neighbors
    for (int n=0; n<3*3*3; ++n) neighborLabels[n] = value_type(0);
    int nLabeledNeighbors;      // counts the number of labelled neighbors
    value_type currentLabel = 0;
  
    // Allocate equivalence chain
    EquivalenceChain lut(MAX_LABEL);
  
    // Allocate variables used in the loop
    value_type minLabel;
  #ifdef VERBOSE
    int k, oldk = -1;
  #endif
    
  
    // For all pixels
  
    typename Iterator<TImage>::Volumetric iIm(im);
    for (; !iIm.isAtEnd(); ++iIm)
    {
  
  #ifdef VERBOSE
      k = iIm.getZ();
      if (oldk != k)
      {
        oldk = k;
        std::cout << "Completed " << (int)((100.0 * k) / (im.getZ()-1) + 0.5) << "%     \r" << std::flush;
      }
  #endif
  
      // If pixel is not background, it has to be labeled
      if (*iIm != 0)
      {
        // Get and count labels of neighbors already processed      
        // Process only backward neighbors (13 out of 26)
        nLabeledNeighbors = 0;
        
        SPAN_NEIGHBORS(-1,-1,-1)
        SPAN_NEIGHBORS( 0,-1,-1)
        SPAN_NEIGHBORS(+1,-1,-1)
  
        SPAN_NEIGHBORS(-1, 0,-1)
        SPAN_NEIGHBORS( 0, 0,-1)
        SPAN_NEIGHBORS(+1, 0,-1)
  
        SPAN_NEIGHBORS(-1,+1,-1)
        SPAN_NEIGHBORS( 0,+1,-1)
        SPAN_NEIGHBORS(+1,+1,-1)
  
        SPAN_NEIGHBORS(-1,-1, 0)
        SPAN_NEIGHBORS( 0,-1, 0)
        SPAN_NEIGHBORS(+1,-1, 0)
  
        SPAN_NEIGHBORS(-1, 0, 0)
  
  
        // If no neighbor is labeled : create new label
        if (nLabeledNeighbors == 0)
        {
          currentLabel = lut.getNewLabel();
  
          // Is label table full?
          if (currentLabel == 0)
          {
            //std::cout << "Labels compressed" << std::endl;
            value_type numberOfLabels = lut.mergeLabels();
  
            typename Iterator<TImage>::Linear iTmp(im);
            iTmp.setEnd(iIm.getIndex()-1);
            for (; !iTmp.isAtEnd(); ++iTmp)
            {
              if (*iTmp)
              {
                *iTmp = lut[*iTmp];
              }
            }
  
            lut.reset(numberOfLabels);
            currentLabel = lut.getNewLabel();
  
            // Table is still full after label merging: error
            if (currentLabel == 0)
            {
              throw std::overflow_error("Too many regions");
            }
          }
  
          *iIm = currentLabel;
        }
        else
        {
          // Take the lowest label
          minLabel = std::numeric_limits<value_type>::max();
          for (int n=0; n<13; ++n)
          {
            if (neighborLabels[n] > 0 && neighborLabels[n] < minLabel)
            {
              minLabel = neighborLabels[n];
            }
          }
  
          *iIm = minLabel;
  
          // If there were more than one label, modify the look-up table
          // to ensure that all labels correspond to the same region
          if (nLabeledNeighbors >= 2)
          {
            int lastLabelCurrent = lut.getLastLabel(*iIm);
  
            for (int n=0; n<13; ++n)
            {
              if (neighborLabels[n] != 0)
              {
                lut.setEquivalence(neighborLabels[n], lastLabelCurrent);
              }
            }
          }
        }
      }
    }
  
  
    // Merge all equivalent labels
      lut.mergeLabels();
  
    // Update label image accordingly
    for (typename Iterator<TImage>::Linear iTmp(im); !iTmp.isAtEnd(); ++iTmp)
    {
      if (*iTmp) *iTmp = lut[*iTmp];
    }
  
    return lut.getGreatestLabel();
  }

  //---------------------------------------------------------------------------

  template < typename TImage, typename TNeighborhood >
  int connectedComponents(TImage &im, const TNeighborhood &nh)
  {
    return connectedComponents<TImage, TNeighborhood, ZeroExtrapolator>(im,nh);
  }

  //---------------------------------------------------------------------------

} // namespace til

