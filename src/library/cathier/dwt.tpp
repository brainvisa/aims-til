
namespace til
{

  //---------------------------------------------------------------------------

  /// namespace for boundary conditions
  // TODO: could actually be in a separate file... included from here.
  namespace bc
  {
    struct ThinMirror
    {
      std::size_t operator()(int i, std::size_t n)
      {
        assert(i >= -1);
        if (i>=0) { assert(std::size_t(i) <= n); }
        if (i == -1) return 3;
        else if (i == int(n)) return n-4;
        //else if (i == int(n)) return ( n%2 ? n-2 : n-4 );
        else return i;
      }
    };

    struct Cyclic
    {
      std::size_t operator()(int i, std::size_t n)
      {
        assert(i >= -1);
        if (i>=0) { assert(std::size_t(i) <= n); }
        if (i == -1) return (n-1);
        else if (i == int(n)) return 0;
        else return i;
      }
    };
  }
  


  //---------------------------------------------------------------------------
   
    //----------------------------------------------//
   //  Various anonymous helper functions for DWT  //
  //----------------------------------------------//
  
  namespace
  {
    /// Computes the number of elements reached when hopping in an array of size length with a step of size s.
    struct SteppedSize : public std::binary_function<std::size_t, std::size_t, std::size_t>
    {
      std::size_t operator()(std::size_t length, std::size_t step)
      { return 1 + (length-1) / step; }
    };

    /// Initial step size -- it should be the largest power of two x so that 1+(4-1)*x is smaller than n.
    struct StepInit : public std::unary_function<std::size_t, std::size_t>
    {
      std::size_t operator()(std::size_t n)
      {
        // This loop may actually be more efficient than a closed form formula for reasonable n.
        std::size_t res;
        for (res = 1; 1+(4-1)*res <= n; res *= 2);
        return res/2;
      }
    };

    inline std::size_t index_offset(
    std::size_t x, std::size_t y, std::size_t z,
    std::size_t dimx, std::size_t dimy)
    {
      return x + dimx * ( y + dimy * z);
    }
    
    inline std::size_t index_offset(
    std::size_t x, std::size_t y,
    std::size_t dimx)
    {
      return x + dimx * y;
    }
    

    /*
    inline std::size_t wrap(int i, std::size_t n)
    {
      assert(i >= -1);
      if (i>=0) { assert(std::size_t(i) <= n); }
      if (i == -1) return 3;
      else if (i == int(n)) return n-4;
      //else if (i == int(n)) return ( n%2 ? n-2 : n-4 );
      else return i;
    }
    */
  }

  //---------------------------------------------------------------------------

  template < typename T, typename BC >
  void DWTCubic<T, BC>::Direct::
  operator()(T * s, std::size_t n, std::size_t step)
  {
    for (std::size_t i = 0; i < n; i+=2) s[i*step] -= ( s[BC()(i+1,n)*step] + s[BC()(i-1,n)*step] ) * T(1.0/4.0);
    for (std::size_t i = 1; i < n; i+=2) s[i*step] -= ( s[BC()(i+1,n)*step] + s[BC()(i-1,n)*step] );
    for (std::size_t i = 0; i < n; i+=2) s[i*step] += ( s[BC()(i+1,n)*step] + s[BC()(i-1,n)*step] ) * T(3.0/16.0);

    // This is the "constant L1 norm" normalization
    //for (std::size_t i = 0; i < n; i+=2) s[i*step] *= T(4);

    // This is the usual normalization
    for (std::size_t i = 0; i < n; i+=2) s[i*step] *= T(2);
    for (std::size_t i = 1; i < n; i+=2) s[i*step] *= T(1.0/2.0);
  }
  //---------------------------------------------------------------------------

  template < typename T, typename BC >
  void DWTCubic<T,BC>::Inverse::
  operator()(T * s, std::size_t n, std::size_t step)
  {
    // This is the "constant L1 norm" normalization
    //for (std::size_t i = 0; i < n; i+=2) s[i*step] *= T(1.0/4.0);
    
    std::size_t step2 = 2 * step;
    std::size_t ns = n * step;
    // This is the usual normalization
    for (std::size_t i = 0;    i < ns; i += step2) s[i] *= T(1.0/2.0);
    for (std::size_t i = step; i < ns; i += step2) s[i] *= T(2);
    
    for (std::size_t i = 0, j = 0;    i < n; i+=2, j += step2) s[j] -= ( s[BC()(i+1,n)*step] + s[BC()(i-1,n)*step] ) * T(3.0/16.0);
    for (std::size_t i = 1, j = step; i < n; i+=2, j += step2) s[j] += ( s[BC()(i+1,n)*step] + s[BC()(i-1,n)*step] );
    for (std::size_t i = 0, j = 0;    i < n; i+=2, j += step2) s[j] += ( s[BC()(i+1,n)*step] + s[BC()(i-1,n)*step] ) * T(1.0/4.0);
  }

  //---------------------------------------------------------------------------

  template < typename T, typename BC >
  void DWTCubicConjugate< T, BC >::Direct::
  operator()(T * s, std::size_t n, std::size_t step)
  {
    for (std::size_t i = 1; i < n; i+=2) s[i*step] += ( s[BC()(i+1,n)*step] + s[BC()(i-1,n)*step] ) * T(1.0/4.0);
    for (std::size_t i = 0; i < n; i+=2) s[i*step] += ( s[BC()(i+1,n)*step] + s[BC()(i-1,n)*step] );
    for (std::size_t i = 1; i < n; i+=2) s[i*step] -= ( s[BC()(i+1,n)*step] + s[BC()(i-1,n)*step] ) * T(3.0/16.0);

    //for (std::size_t i = 0; i < n; i+=2) s[i*step] *= T(1.0/(4.0));

    //for (std::size_t i = 0; i < n; i+=2) s[i*step] *= T(std::sqrt(2.0)/4);
    //for (std::size_t i = 1; i < n; i+=2) s[i*step] *= T(std::sqrt(2.0));

    for (std::size_t i = 0; i < n; i+=2) s[i*step] *= T(1.0/2.0);
    for (std::size_t i = 1; i < n; i+=2) s[i*step] *= T(2);
  }
    
  //---------------------------------------------------------------------------

  template < typename T, typename BC >
  void DWTCubicConjugate< T, BC >::Inverse::
  operator()(T * s, std::size_t n, std::size_t step)
  {
    for (std::size_t i = 0; i < n; i+=2) s[i*step] *= T(2);
    for (std::size_t i = 1; i < n; i+=2) s[i*step] *= T(1.0/2.0);
    for (std::size_t i = 1; i < n; i+=2) s[i*step] += ( s[BC()(i+1,n)*step] + s[BC()(i-1,n)*step] ) * T(3.0/16.0);
    for (std::size_t i = 0; i < n; i+=2) s[i*step] -= ( s[BC()(i+1,n)*step] + s[BC()(i-1,n)*step] );
    for (std::size_t i = 1; i < n; i+=2) s[i*step] -= ( s[BC()(i+1,n)*step] + s[BC()(i-1,n)*step] ) * T(1.0/4.0);
  }

  //---------------------------------------------------------------------------

  template < typename DWT >
  void DWTND<DWT, 2>::operator()
  (
    typename DWT::prec_type * im
  , numeric_array<std::size_t,2> dim
  , int dir
  , numeric_array<std::size_t,2> step
  , numeric_array<std::size_t,2> real_dim
  )
  {
    assert(dir >= 0);
    assert(dir < 2);
    if (dir == 0)
    {
      for (std::size_t j = 0; j < dim[1]; ++j)
      {
        m_dwt(im+j*step[1]*real_dim[0], dim[0], step[0]);
      }
    }
    else if (dir == 1)
    {
      for (std::size_t i = 0; i < dim[0]; ++i)
      {
        m_dwt(im+i*step[0], dim[1], real_dim[0]*step[1]);
      }
    }
  }
  
  //---------------------------------------------------------------------------
  
  template < typename DWT >
  void DWTND<DWT,3>::operator()(prec_type * im, numeric_array<std::size_t,3> dim, int dir, numeric_array<std::size_t,3> step, numeric_array<std::size_t,3> real_dim)
  {
    assert(dir >= 0);
    assert(dir < 3);
    if (dir == 0)
    {
      for (std::size_t k = 0; k < dim[2]; ++k)
      for (std::size_t j = 0; j < dim[1]; ++j)
      {
        m_dwt(im+index_offset(0, j*step[1], k*step[2], real_dim[0], real_dim[1]), dim[0], step[0]);
      }
    }
    else if (dir == 1)
    {
      for (std::size_t k = 0; k < dim[2]; ++k)
      for (std::size_t i = 0; i < dim[0]; ++i)
      {
        m_dwt(im+index_offset(i*step[0], 0, k*step[2], real_dim[0], real_dim[1]), dim[1], real_dim[0]*step[1]);
      }
    }
    else if (dir == 2)
    {
      for (std::size_t j = 0; j < dim[1]; ++j)
      for (std::size_t i = 0; i < dim[0]; ++i)
      {
        m_dwt(im+index_offset(i*step[0], j*step[1], 0, real_dim[0], real_dim[1]), dim[2], real_dim[0]*real_dim[1]*step[2]);
      }
    }
  }

  //---------------------------------------------------------------------------

  template < typename DWT >
  void MultiDWT::Direct<DWT>::operator()
  (
    prec_type * s,    ///< Input array. This also contains the output, the computation being in-place.
    std::size_t n     ///< Length of the input.
  )
  {
    std::size_t sm = 1;
    while (n >= 4)
    {
      m_dwt(s, n, sm);
      n = SteppedSize()(n, 2);
      sm *= 2;
    }
  }

  //---------------------------------------------------------------------------

  template < typename IDWT >
  void MultiDWT::Inverse<IDWT>::operator()
  (
    prec_type * s,  ///< The input array. Also contains the output, since the computation is in-place.
    std::size_t n   ///< The input length.
  )
  {
    for (std::size_t sm = StepInit()(n); sm > 0; sm /= 2)
    {
      m_idwt(s, SteppedSize()(n, sm), sm);
    }
  }

  //---------------------------------------------------------------------------

  template < typename DWT, std::size_t D >
  void MultiDWTND::Direct<DWT,D>::operator()
  (
    prec_type * im
  , const numeric_array<std::size_t,D> dim
  )
  {
    numeric_array<std::size_t,D> sm;
    std::fill(sm.begin(), sm.end(), 1);
    numeric_array<std::size_t,D> sm_bis = sm;
    numeric_array<std::size_t,D> ssize = dim;
    numeric_array<std::size_t,D> ssize_bis = dim;

    numeric_array<bool,D> f;
    for (std::size_t i = 0; i < D; ++i) f[i] = (dim[i] >= 4);
    
    while (or_all(f.begin(), f.end()))
    {
      // loop in all dimensions, starting from the lowest.
      for (std::size_t i = 0; i < D; ++i)
      {
        if (f[i])
        {
          //std::cout << "dir " << i << " step " << sm[0] << " " << sm[1] << " " << sm[2] << " dim " <<  ssize[0] << " " << ssize[1] << " " << ssize[2] << std::endl;
          m_dwt(im, ssize, i, sm, dim);
          if (ssize_bis[i] >= (4-1)*2+1)
          {
            //ssize_bis[i] = 1 + (ssize_bis[i]-1) / 2;
            ssize_bis[i] = SteppedSize()(ssize_bis[i], 2);
            sm_bis[i] *= 2;
          }
          else f[i] = false;
        }
      }
      sm = sm_bis;
      ssize = ssize_bis;
    }
  }

  //---------------------------------------------------------------------------

  template < typename IDWT, std::size_t D >
  void
  MultiDWTND::Inverse< IDWT, D >::
  operator()(prec_type * im, const numeric_array<std::size_t,D> dim)
  {
    numeric_array<std::size_t,D> sm;
    numeric_array<std::size_t,D> ssize;
    std::transform(dim.begin(), dim.end(), sm.begin(), StepInit());
    // NB: sm_bis has to be initialized here, because it might not be initialized in the
    // following loop.
    numeric_array<std::size_t,D> sm_bis = sm;
    
    do
    {
      // loop in all dimensions, starting with the highest
      for (int i = D-1; i >= 0; --i)
      {
        if (sm[i] >= 1 && greater_than_others(sm, i))
        {
          std::transform(dim.begin(), dim.end(), sm.begin(), ssize.begin(), SteppedSize());
          //for (std::size_t j = 0; j < D; ++j) ssize[j] = SteppedSize()(dim[j], sm[j]);
          //std::cout << "dir " << i << " step " << sm[0] << " " << sm[1] << " dim " <<  ssize[0] << " " << ssize[1] << " rdim " <<  dim[0] << " " << dim[1] << std::endl;
          //std::cout << "dir " << i << " step " << sm[0] << " " << sm[1] << " " << sm[2] << " dim " <<  ssize[0] << " " << ssize[1] << " " << ssize[2] << std::endl;
          m_idwt(im, ssize, i, sm, dim);
          sm_bis[i] = sm[i] / 2;
        }
      }
      sm = sm_bis;
    } while (!all_less(sm, 1));
  }

  //---------------------------------------------------------------------------

  template < typename T >
  void MultiDWTPower< T >::
  operator()(T * s, std::size_t length, std::size_t hop)
  {
    for (std::size_t step = 2; SteppedSize()(length, step) >= 4; step *= 2)
    {
      for (std::size_t i = 0; i < length; i += step)
      {
        s[i*hop] *= 2;
      }
    }
  }

  //---------------------------------------------------------------------------

  template < typename T >
  void MultiDWTNDPower<T,3>::
  operator()(T * im, numeric_array<std::size_t,3> dim)
  {
    MultiDWTPower<T> power;
    for (std::size_t k = 0; k < dim[2]; ++k)
    for (std::size_t j = 0; j < dim[1]; ++j)
    {
      power(im+(k*dim[1]+j)*dim[0], dim[0], 1);
    }
    for (std::size_t k = 0; k < dim[2]; ++k)
    for (std::size_t i = 0; i < dim[0]; ++i)
    {
      power(im+i+ k*dim[0]*dim[1], dim[1], dim[0]);
    }
    for (std::size_t j = 0; j < dim[1]; ++j)
    for (std::size_t i = 0; i < dim[0]; ++i)
    {
      power(im+i+j*dim[0], dim[2], dim[0]*dim[1]);
    }     
  }

  //---------------------------------------------------------------------------
  
  namespace
  {
    /// Reorder sequence by putting the even index first, then the odd.
    template < typename T >
    void shuffle_step(T * s, std::size_t length, T * buffer, std::size_t jump)
    {
      for (std::size_t i = 0; i < length; ++i)
      {
        buffer[i] = s[i*jump];
      }
      for (T * p = buffer; p < buffer+length; p+=2, s+=jump)
      {
        *s = *p;
      }
      for (T * p = buffer+1; p < buffer+length; p+=2, s+=jump)
      {
        *s = *p;
      }
    }
    
    /// Reorder sequence by dispatching first half on even indexes, second half on odd.
    template < typename T >
    void unshuffle_step(T * s, std::size_t length, T * buffer, std::size_t jump)
    {
      for (std::size_t i = 0, j = 0; i < length; ++i, j += jump)
      {
        buffer[i] = s[j];
      }
      std::size_t jump2 = 2*jump;
      std::size_t jl = jump * length;
      for (T * p = s; p < s + jl; p += jump2)
      {
        *p = *(buffer++);
      }
      for (T * p = s+jump; p < s + jl; p += jump2)
      {
        *p = *(buffer++);
      }
    }
  } // namespace

  //---------------------------------------------------------------------------

  template < typename T >
  void MultiDWTShuffle<T>::shuffle(T * s, std::size_t length, std::size_t jump)
  {
    if (m_buffer.size() < length) m_buffer.resize(length);
    for (std::size_t step = 1; SteppedSize()(length, step) >= 4; step *= 2)
    {
      //std::cout << "length " << SteppedSize()(step, length) << " step " << step << std::endl;
      shuffle_step(s, SteppedSize()(length, step), &m_buffer.front(), jump);
    }
  }
    
  //---------------------------------------------------------------------------

  template < typename T >
  void MultiDWTShuffle<T>::unshuffle(T * s, std::size_t length, std::size_t jump)
  {
    if (m_buffer.size() < length) m_buffer.resize(length);
    std::size_t step;
    for (step = 1; SteppedSize()(length, step) >= 4; step *= 2);
    do
    {
      step/=2;
      unshuffle_step(s, SteppedSize()(length, step), &m_buffer.front(), jump);
    } while (step > 1);
  }
  
  //---------------------------------------------------------------------------

  template < typename T >
  void MultiDWTNDShuffle<T,2>::shuffle(T * im, numeric_array<std::size_t,2> dim)
  {
    MultiDWTShuffle<T> shuffle;
    for (std::size_t j = 0; j < dim[1]; ++j)
    {
      shuffle.shuffle(im+j*dim[0], dim[0], 1);
    }
    for (std::size_t i = 0; i < dim[0]; ++i)
    {
      shuffle.shuffle(im+i, dim[1], dim[0]);
    }
  }
    
  //---------------------------------------------------------------------------

  template < typename T >
  void MultiDWTNDShuffle<T,2>::unshuffle(T * im, numeric_array<std::size_t,2> dim)
  {
    MultiDWTShuffle<T> shuffle;
    for (std::size_t j = 0; j < dim[1]; ++j)
    {
      shuffle.unshuffle(im+j*dim[0], dim[0], 1);
    }
    for (std::size_t i = 0; i < dim[0]; ++i)
    {
      shuffle.unshuffle(im+i, dim[1], dim[0]);
    }
  }
  
  //---------------------------------------------------------------------------

  template < typename T >
  void MultiDWTNDShuffle<T,3>::shuffle(T * im, numeric_array<std::size_t,3> dim)
  {
    MultiDWTShuffle<T> shuffle;
    for (std::size_t k = 0; k < dim[2]; ++k)
    for (std::size_t j = 0; j < dim[1]; ++j)
    {
      shuffle.shuffle(im+(k*dim[1]+j)*dim[0], dim[0], 1);
    }
    for (std::size_t k = 0; k < dim[2]; ++k)
    for (std::size_t i = 0; i < dim[0]; ++i)
    {
      shuffle.shuffle(im+i+ k*dim[0]*dim[1], dim[1], dim[0]);
    }
    for (std::size_t j = 0; j < dim[1]; ++j)
    for (std::size_t i = 0; i < dim[0]; ++i)
    {
      shuffle.shuffle(im+i+j*dim[0], dim[2], dim[0]*dim[1]);
    }
  }
  
  //---------------------------------------------------------------------------

  template < typename T >
  void MultiDWTNDShuffle<T,3>::unshuffle(T * im, numeric_array<std::size_t,3> dim)
  {
    MultiDWTShuffle<T> shuffle;
    for (std::size_t k = 0; k < dim[2]; ++k)
    for (std::size_t j = 0; j < dim[1]; ++j)
    {
      shuffle.unshuffle(im+(k*dim[1]+j)*dim[0], dim[0], 1);
    }
    for (std::size_t k = 0; k < dim[2]; ++k)
    for (std::size_t i = 0; i < dim[0]; ++i)
    {
      shuffle.unshuffle(im+i+ k*dim[0]*dim[1], dim[1], dim[0]);
    }
    for (std::size_t j = 0; j < dim[1]; ++j)
    for (std::size_t i = 0; i < dim[0]; ++i)
    {
      shuffle.unshuffle(im+i+j*dim[0], dim[2], dim[0]*dim[1]);
    }
  }

  //---------------------------------------------------------------------------

  /// In-place (4-2) DWT via lifting.
  template < typename T, typename BC >
  inline void dwt_cubic(T * s, std::size_t n, std::size_t step, BC)
  {
    typename DWTCubic<T, BC>::Direct()(s, n, step);
  }
  /// Same as above with cyclic bc by default
  template < typename T >
  inline void dwt_cubic(T * s, std::size_t n, std::size_t step)
  {
    dwt_cubic(s, n, step, bc::Cyclic());
  }

  /// In-place inverse (4-2) DWT via lifting.
  template < typename T, typename BC >
  inline void idwt_cubic(T * s, std::size_t n, std::size_t step, BC)
  {
    typename DWTCubic<T, BC>::Inverse()(s, n, step);
  }
  /// Same as above with cyclic bc by default
  template < typename T >
  inline void idwt_cubic(T * s, std::size_t n, std::size_t step)
  {
    idwt_cubic(s, n, step, bc::Cyclic());
  }
  
  template < typename T, std::size_t D >
  inline void dwtND_cubic(T * im, numeric_array<std::size_t, D> dim)
  {
    numeric_array<std::size_t, D> step;
    for (std::size_t i = 0; i < D; ++i)
    {
      step[i] = 1;
    }
    
    for (std::size_t i = 0; i < D; ++i)
    {
      DWTND<typename DWTCubic<T, bc::Cyclic>::Direct,D>()(im, dim, i, step, dim);
    }
  }

  template < typename T, std::size_t D >
  inline void dwtND_cubicConjugate(T * im, numeric_array<std::size_t, D> dim)
  {
    numeric_array<std::size_t, D> step;
    for (std::size_t i = 0; i < D; ++i)
    {
      step[i] = 1;
    }
    
    for (std::size_t i = 0; i < D; ++i)
    {
      DWTND<typename DWTCubicConjugate<T, bc::Cyclic>::Direct,D>()(im, dim, i, step, dim);
    }
  }

  
  /*
  /// In-place multidimensional (4-2) DWT along a given direction.
  template < typename T, std::size_t D >
  inline void dwtND_cubic(T * im, numeric_array<std::size_t,D> dim, int dir, numeric_array<std::size_t,D> step, numeric_array<std::size_t,D> real_dim)
  {
    DWTND<DWTCubic<T>,D>()(im, dim, dir, step, real_dim);
  }

  /// In-place multidimensional inverse (4-2) DWT along a given direction.
  template < typename T, std::size_t D >
  inline void idwtND_cubic(T * im, numeric_array<std::size_t,D> dim, int dir, numeric_array<std::size_t,D> step, numeric_array<std::size_t,D> real_dim)
  {
    DWTND<IDWTCubic<T>,D>()(im, dim, dir, step, real_dim);
  }
  */
    
  template < typename T, typename BC >
  void multi_dwt_cubic(T * s, std::size_t dim, BC)
  {
    MultiDWT::Direct<typename DWTCubic<T, BC>::Direct>()(s, dim);
  }
  template < typename T >
  void multi_dwt_cubic(T * s, std::size_t dim)
  {
    multi_dwt_cubic(s, dim, bc::Cyclic());
  }

  template < typename T, typename BC >
  void multi_idwt_cubic(T * s, std::size_t dim, BC)
  {
    MultiDWT::Inverse<typename DWTCubic<T, BC>::Inverse>()(s, dim);
  }
  template < typename T >
  void multi_idwt_cubic(T * s, std::size_t dim)
  {
    multi_idwt_cubic(s, dim, bc::Cyclic());
  }
  
  template < typename T, std::size_t D, typename BC >
  void multi_dwtND_cubic(T * im, numeric_array<std::size_t,D> dim, BC)
  {
    MultiDWTND::Direct<DWTND<typename DWTCubic<T, BC>::Direct,D>,D>()(im, dim);
  }
  template < typename T, std::size_t D >
  void multi_dwtND_cubic(T * im, numeric_array<std::size_t,D> dim)
  {
    multi_dwtND_cubic(im, dim, bc::Cyclic());
  }

  template < typename T, std::size_t D, typename BC >
  void multi_idwtND_cubic(T * im, numeric_array<std::size_t,D> dim, BC)
  {
    MultiDWTND::Inverse<DWTND<typename DWTCubic<T, BC>::Inverse,D>,D>()(im, dim);
  }
  template < typename T, std::size_t D >
  void multi_idwtND_cubic(T * im, numeric_array<std::size_t,D> dim)
  {
    multi_idwtND_cubic(im, dim, bc::Cyclic());
  }
  
  template < typename T >
  void multi_dwt_shuffle(T * s, std::size_t n)
  {
    MultiDWTShuffle<T>().shuffle(s,n);
  }
  
  template < typename T >
  void multi_dwt_unshuffle(T * s, std::size_t n)
  {
    MultiDWTShuffle<T>().unshuffle(s, n);
  }
  
  template < typename T, std::size_t D >
  void multi_dwtND_shuffle(T * im, numeric_array<std::size_t, D> dim)
  {
    MultiDWTNDShuffle<T,D>().shuffle(im, dim);
  }
  
  template < typename T, std::size_t D >
  void multi_dwtND_unshuffle(T * im, numeric_array<std::size_t, D> dim)
  {
    MultiDWTNDShuffle<T,D>().unshuffle(im, dim);
  }
  
  template < typename T, std::size_t D >
  void multi_dwtND_power(T * im, numeric_array<std::size_t, D> dim)
  {
    MultiDWTNDPower<T,D>()(im, dim);
  }
  
} // namespace til

