


namespace til { namespace math {


  //---------------------------------------------------------------------------

  namespace detail
  {
    //-------------------------------------------------------------------------

    template < typename T >
    T ahuber(T ax, T K)
    {
      // NB: testing K >= 0 is actually not needed.
      if (ax <  K)
      {
        return ax*ax;
      }
      else
      {
        return 2 * K * (ax - K);
      }
    }

    //-------------------------------------------------------------------------

  } // namespace detail

  //---------------------------------------------------------------------------

  template < typename T >
  T huber(T x, T K)
  {
    return detail::ahuber(std::abs(x), K);
  }
  
  //---------------------------------------------------------------------------

  template < typename T >
  T d_huber(T x, T K)
  {
    if (std::abs(x) < K)
    {
      return 2 * x;
    }
    else if (x > 0)
    {
      return 2 * K;
    }
    else
    {
      return (-2) * K;
    }
  }

  //---------------------------------------------------------------------------

  template < typename T >
  T shrink(T x, T K)
  {
    assert(K >= 0);
    if (std::abs(x) < K)
    {
      return 0;
    }
    else if (x > 0)
    {
      return x - K;
    }
    else
    {
      return x + K;
    }
  }

  //---------------------------------------------------------------------------

  /// Constructor with Gaussian standard deviation.
  template < typename TPrec >
  Gaussian<TPrec>::Gaussian(prec_type sigma)
  {
    this->setSigma(sigma);
  }

  //---------------------------------------------------------------------------

  template < typename TPrec >
  void
  Gaussian<TPrec>::setSigma(prec_type sigma)
  {
    assert(sigma > 0);
    m_2var = 2 * sigma * sigma;
    if (m_2var < 128 * std::numeric_limits<prec_type>::epsilon())
    {
      throw StandardDeviationTooSmall();
    }
  }

  //---------------------------------------------------------------------------

  template < typename TPrec >
  typename Gaussian<TPrec>::prec_type
  Gaussian<TPrec>::operator()(prec_type x) const
  {
    return squared_input(x*x);
  }
    
  //---------------------------------------------------------------------------

  template < typename TPrec >
  typename Gaussian<TPrec>::prec_type
  Gaussian<TPrec>::squared_input(prec_type x2) const
  {
    return std::exp(-x2 / m_2var);
  }

  //---------------------------------------------------------------------------
  
}} // namespace til::math

