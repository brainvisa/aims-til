#ifndef TIL_CONNECTEDCOMPONENTS3D_H
#define TIL_CONNECTEDCOMPONENTS3D_H

/// \file
/// Connected components labeling

// Standard library includes
#include <iostream>

// Local includes
#include "til/til_common.h"
#include "til/EquivalenceChain.h"
#include "til/ImageExtrapolator.h"
#include "til/Neighborhood.h"
#include "til/neighborhoodTools.h"


// macro used inside connected components
// It puts the neighbor label inside neighborLabels and increase
// nLabeledNeighbors accordingly
// To be undefined at the end of the file
#define SPAN_NEIGHBORS(i,j,k)                                                                   \
  if (nh.template isNeighbor<(i),(j),(k)>() &&                                                  \
  ((neighborLabels[COFFSET(i,j,k)] = iIm.template getValue<TExtrapolator,(i),(j),(k)>()) != 0)) \
    ++nLabeledNeighbors;

namespace til
{

  //---------------------------------------------------------------------------

  /// Label the connected components of image im.
  /// Returns the number of connected components.
  template < typename TImage, typename TNeighborhood, typename TExtrapolator >
  int connectedComponents(TImage &im, const TNeighborhood &nh);

  //---------------------------------------------------------------------------

  /// Label the connected components of image im.
  /// Returns the number of connected components.
  /// Default extrapolation of zero
  // TODO: this has to change: extrapolation of zero by default??? Remove #include "extrapol"
  template < typename TImage, typename TNeighborhood >
  int connectedComponents(TImage &im, const TNeighborhood &nh)
  {
  	return connectedComponents<TImage, TNeighborhood, ZeroExtrapolator>(im,nh);
  }

  //---------------------------------------------------------------------------

} // namespace til

// Undefine local macros
#undef SPAN_NEIGHBORS

// package include
#include "connected_components.tpp"

#endif

