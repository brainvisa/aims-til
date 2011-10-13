#ifndef TIL_DWT_H_
#define TIL_DWT_H_

/// \file Stuff related to discrete wavelet transforms.

// includes from TIL
#include "til/numeric_array.h"

namespace til
{
  
  //---------------------------------------------------------------------------
  
    //------------//
   //  DWTCubic  //
  //------------//

  /// In-place (4-2) DWT via lifting.
  template < typename T, typename BC >
  struct DWTCubic
  {
    /// Direct DWT.
    struct Direct
    {
      typedef T prec_type;
      void operator()(T * s, std::size_t n, std::size_t step);
    };
    
    /// Inverse DWT.
    struct Inverse
    {
      typedef T prec_type;
      void operator()(T * s, std::size_t n, std::size_t step);
    };
  };
  
  //------------------------------------------------------------------------------------------------

    //---------------------//
   //  DWTCubicConjugate  //
  //---------------------//
  
  /// Transposed filter of (4-2).
  template < typename T, typename BC >
  struct DWTCubicConjugate
  {
    /// Direct DWT.
    struct Direct
    {
      typedef T prec_type;
      void operator()(T * s, std::size_t n, std::size_t step);
    };
    
    /// Inverse DWT.
    struct Inverse
    {
      typedef T prec_type;
      void operator()(T * s, std::size_t n, std::size_t step);
    };
  };
  
  //---------------------------------------------------------------------------
   
    //---------//
   //  DWTND  //
  //---------//

  /// Extension of DWT to multidimensional arrays.
  template < typename DWT, std::size_t D >
  class DWTND;

  //---------------------------------------------------------------------------

  template < typename DWT >
  class DWTND< DWT, 2 >
  {
  public: // typedefs
    typedef typename DWT::prec_type prec_type;
    void operator()
    (
      prec_type * im
    , numeric_array<std::size_t,2> dim
    , int dir
    , numeric_array<std::size_t,2> step
    , numeric_array<std::size_t,2> real_dim
    );
  private: // data
    DWT m_dwt;
  };
  
  //---------------------------------------------------------------------------

  template < typename DWT >
  class DWTND< DWT, 3 >
  {
  public: // typedefs
    typedef typename DWT::prec_type prec_type;
    void operator()
    (
      prec_type * im
    , numeric_array<std::size_t,3> dim
    , int dir
    , numeric_array<std::size_t,3> step
    , numeric_array<std::size_t,3> real_dim
    );
  private: // data
    DWT m_dwt;
  };

  /*
  //---------------------------------------------------------------------------

    //--------------//
   //  DWTND (2D)  //
  //--------------//

  /// Specialization in two dimensions.
  template < typename DWT >
  class DWTND<DWT, 2>
  {
  public: // typedefs
    typedef typename DWT::prec_type prec_type;

  public: // operators
    void operator()(prec_type * im, numeric_array<std::size_t,2> dim, int dir, numeric_array<std::size_t,2> step, numeric_array<std::size_t,2> real_dim)
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
    
  private: // data
    DWT m_dwt;
  };
  
    //--------------//
   //  DWTND (3D)  //
  //--------------//
  
  /// Specialization in three dimensions.
  template < typename DWT >
  class DWTND<DWT,3>
  {
  public: // typedefs
    typedef typename DWT::prec_type prec_type;

  public: // operators
    void operator()(prec_type * im, numeric_array<std::size_t,3> dim, int dir, numeric_array<std::size_t,3> step, numeric_array<std::size_t,3> real_dim);
    
  private: // data
    DWT m_dwt;
  };
  //------------------------------------------------------------------------------------------------
  */
   
    //------------//
   //  MultiDWT  //
  //------------//

  /// Multi-scale in-place DWT.
  /// Applies successive DWT on the base coefficients.
  struct MultiDWT
  {
    /// Direct multiscale DWT.
    template < typename DWT >
    class Direct
    {
    public: // typedefs
      typedef typename DWT::prec_type prec_type;

    public: // operators
      void operator()
      (
        prec_type * s,    ///< Input array. This also contains the output, the computation being in-place.
        std::size_t n     ///< Length of the input.
      );

    private: // data
      DWT m_dwt;
    };
    
    /// Inverse multiscale DWT.
    template < typename IDWT >
    class Inverse
    {
    public: // typedefs
      typedef typename IDWT::prec_type prec_type;

    public: // operators
      void operator()
      (
        prec_type * s,  ///< The input array. Also contains the output, since the computation is in-place.
        std::size_t n   ///< The input length.
      );

    private: // data
      IDWT m_idwt;
    };
  };

  //------------------------------------------------------------------------------------------------

    //--------------//
   //  MultiDWTND  //
  //--------------//

  /// Multi-scale in-place DWT in N-D.
  /// Applies successive DWT on the base coefficients along the different axes.
  struct MultiDWTND
  {
    template < typename DWT, std::size_t D >
    class Direct
    {
    public: // typedefs
      typedef typename DWT::prec_type prec_type;
    public: // operator
      void operator()(prec_type * im, const numeric_array<std::size_t,D> dim);
    
    private: // function
    
      template < typename TIterator >
      static bool or_all(TIterator begin, TIterator end)
      {
        for (; begin != end; ++begin)
        {
          if (*begin) return true;
        }
        return false;
      }
  
    private: // data
      DWT m_dwt;
    };
    

    template < typename IDWT, std::size_t D >
    class Inverse
    {
    public: // typedefs
      typedef typename IDWT::prec_type prec_type;

    public: // operator    
      void operator()(prec_type * im, const numeric_array<std::size_t,D> dim);
      
    private: // static functions
    
      static inline bool greater_than_others(numeric_array<std::size_t,D> sm, int i0)
      {
        for (std::size_t i = 0; i < D; ++i)
        {
          // NB: I removed this test, because I think, in general, it slows down the code. However it works
          // only if the test is always true for two equal values.
          //if (i == i0) continue;
          if (sm[i] > sm[i0]) return false;
        }
        return true;
      }
  
    private: // data
      IDWT m_idwt;
    };
  };

  //------------------------------------------------------------------------------------------------
    
  template < typename T >
  struct MultiDWTPower
  {
    void operator()(T * s, std::size_t length, std::size_t hop);
  };

  //------------------------------------------------------------------------------------------------

  template < typename T, std::size_t D >
  struct MultiDWTNDPower;

  template < typename T >
  struct MultiDWTNDPower< T, 3 >
  {
    void operator()(T * im, numeric_array<std::size_t,3> dim);
  };
  
  //------------------------------------------------------------------------------------------------

    //--------------//
   //  DWTShuffle  //
  //--------------//

  template < typename T >
  class MultiDWTShuffle
  {
  public: // constructors
    MultiDWTShuffle() : m_buffer() {}
    MultiDWTShuffle(std::size_t maxlength) : m_buffer(maxlength) {}
  public: // functions
    void shuffle(T * s, std::size_t length, std::size_t jump = 1);    
    void unshuffle(T * s, std::size_t length, std::size_t jump = 1);
  private: // data  
    std::vector<T> m_buffer;
  };

  //------------------------------------------------------------------------------------------------

  template < typename T, std::size_t D >
  struct MultiDWTNDShuffle;
  
  //------------------------------------------------------------------------------------------------

  template < typename T >
  struct MultiDWTNDShuffle<T,2>
  {
    void shuffle(T * im, numeric_array<std::size_t,2> dim);
    void unshuffle(T * im, numeric_array<std::size_t,2> dim);
  };

  //------------------------------------------------------------------------------------------------

  template < typename T >
  struct MultiDWTNDShuffle<T,3>
  {
    void shuffle(T * im, numeric_array<std::size_t,3> dim);
    void unshuffle(T * im, numeric_array<std::size_t,3> dim);
  };

  //------------------------------------------------------------------------------------------------

    //---------------------------------//
   //  Helper functions for all this  //
  //---------------------------------//

  /// In-place (4-2) DWT via lifting.
  template < typename T, typename BC >
  inline void dwt_cubic(T * s, std::size_t n, std::size_t step, BC);

  /// Same as above with cyclic bc by default
  template < typename T >
  inline void dwt_cubic(T * s, std::size_t n, std::size_t step);

  /// In-place inverse (4-2) DWT via lifting.
  template < typename T, typename BC >
  inline void idwt_cubic(T * s, std::size_t n, std::size_t step, BC);

  /// Same as above with cyclic bc by default
  template < typename T >
  inline void idwt_cubic(T * s, std::size_t n, std::size_t step);
  
  template < typename T, std::size_t D >
  inline void dwtND_cubic(T * im, numeric_array<std::size_t, D> dim);

  template < typename T, std::size_t D >
  inline void dwtND_cubicConjugate(T * im, numeric_array<std::size_t, D> dim);
  
  template < typename T, typename BC >
  void multi_dwt_cubic(T * s, std::size_t dim, BC);

  template < typename T >
  void multi_dwt_cubic(T * s, std::size_t dim);

  template < typename T, typename BC >
  void multi_idwt_cubic(T * s, std::size_t dim, BC);

  template < typename T >
  void multi_idwt_cubic(T * s, std::size_t dim);
  
  template < typename T, std::size_t D, typename BC >
  void multi_dwtND_cubic(T * im, numeric_array<std::size_t,D> dim, BC);

  template < typename T, std::size_t D >
  void multi_dwtND_cubic(T * im, numeric_array<std::size_t,D> dim);

  template < typename T, std::size_t D, typename BC >
  void multi_idwtND_cubic(T * im, numeric_array<std::size_t,D> dim, BC);

  template < typename T, std::size_t D >
  void multi_idwtND_cubic(T * im, numeric_array<std::size_t,D> dim);
  
  template < typename T >
  void multi_dwt_shuffle(T * s, std::size_t n);
  
  template < typename T >
  void multi_dwt_unshuffle(T * s, std::size_t n);
  
  template < typename T, std::size_t D >
  void multi_dwtND_shuffle(T * im, numeric_array<std::size_t, D> dim);
  
  template < typename T, std::size_t D >
  void multi_dwtND_unshuffle(T * im, numeric_array<std::size_t, D> dim);
  
  template < typename T, std::size_t D >
  void multi_dwtND_power(T * im, numeric_array<std::size_t, D> dim);
  
} // namespace til

// package include
#include "dwt.tpp"

#endif /*DWT_H_*/
