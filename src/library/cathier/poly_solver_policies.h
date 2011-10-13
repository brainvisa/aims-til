#ifndef TIL_POLY_SOLVER_POLICIES_H
#define TIL_POLY_SOLVER_POLICIES_H


namespace til { namespace math { namespace policy
{    

  //---------------------------------------------------------------------------

  /// If the polynomial has an infinity of solutions, set a flag to true that
  /// can be queried
  class InfinitySolutions_Advertise
  {
  public: // typedefs
  
    typedef InfinitySolutions_Advertise Self;
  
  public: // constructors & destructor
    InfinitySolutions_Advertise() : m_hasInfSols(false) {}
    
  public: // functions

    /// Called by solver
    template < typename TPrec >
    void operator()(PolySolver_real<TPrec, Self> & solver)
    {
      m_hasInfSols = true;
    }
      
    /// Query if equation has an infinity of solutions
    bool hasInfSols() { return m_hasInfSols; }

  private: // data
  
    bool m_hasInfSols;
  };

  //---------------------------------------------------------------------------

  /// If the polynomial has an infinity of solutions, ignore all of them
  struct InfinitySolutions_None
  {
    typedef InfinitySolutions_None Self;
    /// Called by solver
    template < typename TPrec >
    void operator()(PolySolver_real<TPrec, Self> & solver)
    {
      solver.set_nsols(0);
    }
  };

  //---------------------------------------------------------------------------

  /// If the polynomial has an infinity of solutions, throw an exception
  struct InfinitySolutions_Throw
  {
    class InfinityOfSolutions : public std::exception {};
    
    typedef InfinitySolutions_Throw Self;
    /// Called by solver
    template < typename TPrec >
    void operator()(PolySolver_real<TPrec, Self> &)
    {
      throw InfinityOfSolutions();
    }
  };

  //---------------------------------------------------------------------------

}}} // namespace til::math::policy


#endif

