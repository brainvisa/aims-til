#ifndef POLYSOLVER_H_
#define POLYSOLVER_H_


namespace til { namespace math {
  
  
  namespace policy { class InfinitySolutions_None; }
    
    
  //---------------------------------------------------------------------------
  
    //-------------------//
   //  PolySolver_real  //
  //-------------------//
    
  /// Find real roots of a polynomial.
  /// The Infinity Solutions Policy tells what to do in case there are an 
  /// infinity of solutions. There is no default policy, because this is an 
  /// important choice that can be too easily ignored.
  template < typename TPrec, typename TInfinitySolutionsPolicy >
  class PolySolver_real
  {  
  public: // constructors -----------------------------------------------------

    PolySolver_real() : m_infSolPolicy() {};

  public: // friends ----------------------------------------------------------

    // unfortunately, one cannot have a template parameter as a friend...
    // TODO: This stuff has to be resolved, because obvisoulsy I don't want
    // to add a friend to this class everytime a new policy is designed.
    friend class policy::InfinitySolutions_None;
    
  public: // set & get --------------------------------------------------------

    /// Return the number of solutions.
    int nsols() const { return m_nSols; }
    /// Return an array containing the solutions.
    const TPrec * sols() const { return m_sols; }
    /// Return the "infinity of solutions" policy.
    TInfinitySolutionsPolicy getInfSolPolicy() const { return m_infSolPolicy; }
    
  public: // functions --------------------------------------------------------
  
    ///  solve first order polynomial equation [ a X  + b = 0]
    void solve(TPrec a, TPrec b);
    
    /// solve 2nd order polynomial equation [ a X^2  + b X  + C = 0]
    void solve(TPrec a, TPrec b, TPrec c);
  
  private: // set & get -------------------------------------------------------

    void set_nsols(int n) { m_nSols = n; }
    TPrec * get_sols() { return m_sols; }

  private: // data ------------------------------------------------------------
    // number of solutions
    // NB: we don't use std::vector because the solver has to be as fast as possible, AND the maximum size
    // of the array is known and very small.
    int m_nSols;
    // solutions
    // Well, we can't solve polynomial equations with a degree higher than 5, so a max of 5 solution will do.
    TPrec m_sols[5];
    // policy for infinity of solutions 
    TInfinitySolutionsPolicy m_infSolPolicy;
  };

  //---------------------------------------------------------------------------



}} // namespace til::math

// package include
#include "poly_solver.tpp"
#include "poly_solver_policies.h"

#endif /*POLYSOLVER_H_*/
