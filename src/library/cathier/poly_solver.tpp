
namespace til { namespace math
{
  
  //---------------------------------------------------------------------------
  
  template < typename TPrec, typename TInfinitySolutionsPolicy >
  void 
  PolySolver_real< TPrec, TInfinitySolutionsPolicy >::
  solve(TPrec a, TPrec b)
  {
    const TPrec EPSILON = 128*std::numeric_limits<TPrec>::epsilon();
    
    if (std::abs(a) < EPSILON)
    {
      if (std::abs(b) < EPSILON)
      {
        m_infSolPolicy(*this);
      }
      else
      {
        m_nSols = 0;
      }
    }
    else
    {
      m_nSols = 1;
      m_sols[0] = -b / a;
    }
  }
  
  //---------------------------------------------------------------------------

  template < typename TPrec, typename TInfinitySolutionsPolicy >
  void 
  PolySolver_real< TPrec, TInfinitySolutionsPolicy >::
  solve(TPrec a, TPrec b, TPrec c)
  {
    const TPrec EPSILON = 128*std::numeric_limits<TPrec>::epsilon();
    
    // Check that a is not too small -- otherwise we'll have numerical problems, that
    // are solved here by considering that a small a makes the polynomial first degree
    if (std::abs(a) < EPSILON)
    {
      this->solve(b,c);
      return;
    }
    
    TPrec delta = b*b-4*a*c;
    
    // Start with this case, as it is probably the most frequent case
    if (delta > EPSILON)
    {
      m_nSols = 2;
      delta = std::sqrt(delta);
      a *= 2;
      m_sols[0] = (-b-delta)/a;
      m_sols[1] = (-b+delta)/a;
    }
    // TODO: I am not sure whether this should not be delta < 0 instead...
    else if (delta < -EPSILON)
    {
      m_nSols = 0;
    }
    else
    {
      m_nSols = 1;
      m_sols[0] = -b / (2*a);
    }
  }
  
  //---------------------------------------------------------------------------

}} // namespace til::math

