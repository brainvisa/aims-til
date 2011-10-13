#ifndef MINTOOLS_H_
#define MINTOOLS_H_

/// \file Provide multidimensional minimization tools, massively inspired from the Numerical Recipies.

// includes from STL
#include <algorithm>
#include <cmath>
#include <functional>

// includes from TIL
//#include "til/til_common.h"
#include "til/labelClasses.h"
#include "til/TExpr.h"
#include "til/traits.h"

// includes from TIL
#include "globalTraits.h"
#include "miscUtils.h"


namespace til
{
  
  //-----------------------------------------------------------------------------------------------------
  
    //----------//
   //  Mnbrak  //
  //----------//
  
  /// A class to find a bracketing triplet around a function minimum.
  template < typename TFunctor >
  class Mnbrak
  {
  public: // typedefs
  
    typedef typename TFunctor::argument_type            input_type;
    typedef typename value_type_of<input_type>::type    input_prec;
    typedef typename TFunctor::result_type              output_prec;

  public: // constructors & destructor
  
    Mnbrak(TFunctor functor)
      : m_functor(functor)
      , m_verbose(false)
    {}

  public: // set & get
  
    bool & verbose() { return m_verbose; }
  
  public: // functions
  
  void operator()(input_prec & ax, input_prec & bx, input_prec & cx, input_type point, input_type dir);

  private: // constants
  
    static const input_prec GOLD;
    static const input_prec EPSILON;
    // maximum magnification step
    static const input_prec GLIMIT;
  
  private: // data

    TFunctor m_functor;
    output_prec m_fa;
    output_prec m_fb;
    output_prec m_fc;    

  private: // data, internals

    input_type m_buf;
    bool m_verbose;
  };
    

  //-----------------------------------------------------------------------------------------------------

    //---------//
   //  Brent  //
  //---------//
  
  /// A class to minimize a functional within a bracketing triplet.
  template < typename TFunctor >
  class Brent
  {
  public: // typedefs
  
    typedef typename TFunctor::argument_type                          input_type;
    typedef typename value_type_of<input_type>::type                  input_prec;
    typedef typename TFunctor::result_type                            output_prec;
    typedef typename combine<input_prec, output_prec>::type           max_prec;

  public: // constructors & destructor

    Brent();
    Brent(TFunctor functor);

  public: // init
    
    /// Somewhat meaningful default values
    void init()
    {
      m_itmax = 100;
      m_tol = std::max<max_prec>(std::sqrt(std::numeric_limits<input_prec>::epsilon()), 
                                 std::sqrt(std::numeric_limits<output_prec>::epsilon()));
      m_verbose = false;
    }
    
  public: // set & get
  
    void setMaxIterations(int itmax) { m_itmax = itmax; }

    max_prec tol() const { return m_tol; }
    void setTol(input_prec tol) { m_tol = tol; }
    
    // Return the argmin
    input_prec xmin() const { return m_xmin; }
    // Return the function value at argmin
    output_prec fmin() const { return m_fmin; }

    bool & verbose() { return m_verbose; }
    
  public: // operators
    
    /// Minimize functor starting from p and in the direction dir
    /// [a,m,b] is a bracketing triplet of the minimum
    void operator()(input_prec a, input_prec m, input_prec b, input_type point, input_type dir);
        
  private: // checks
  
    typedef check< ! std::numeric_limits<input_prec >::is_integer > Check_input_is_floating_point;
    typedef check< ! std::numeric_limits<output_prec>::is_integer > Check_output_is_floating_point;
   
  private: // static constants
  
    static const input_prec CGOLD;
   
  private: // data
    
    // The functor to minimize
    TFunctor m_functor;
    // The min position
    input_prec m_xmin;
    // The min value
    output_prec m_fmin;
    // imprecision tolerence on our result
    output_prec m_tol;
    // maximum number of iterations
    int m_itmax;
    bool m_verbose;
    
  private: // data, internals
    input_type m_buf;
  };
  
  //-----------------------------------------------------------------------------------------------------

  template < typename TFunctor >
  class LMLike
    : public line_minimization_algorithm_label
  {
  public: // typedefs
    typedef typename TFunctor::argument_type            input_type;
    typedef typename value_type_of<input_type>::type    input_prec;
    typedef typename TFunctor::result_type              output_prec;
  
  public: // constructors & destructor
  
    LMLike(TFunctor functor)
      : m_functor(functor)
    {
      this->init();
    }

    LMLike(TFunctor functor, input_prec startfactor, input_prec coeff)
      : m_functor(functor)
      , m_factor(startfactor)
      , m_coeff(coeff)
    {
    }


  public: // initialization

    void init()
    {
      m_factor = 1;
      m_coeff = 2;
    }

  public: // set & get
  
    TFunctor & functor() { return m_functor; }
    const TFunctor & functor() const { return m_functor; }

    input_prec & factor() { return m_factor; }
    const input_prec & factor() const { return m_factor; }
    
    const input_prec & coeff() const { return m_coeff; }
    void setCoeff(input_prec & coeff)
    {
      assert(coeff > 1.0);
      m_coeff = coeff;
    }

  public: // operators
  
    output_prec operator()(input_type & p, const input_type & dir)
    {
      using namespace expr;   
      assert(p.size() == dir.size());
      resize(m_buf, p.size());
      
      // Compute the value at current point
      output_prec f0 = m_functor(m_buf);
      // Compute the value at estimated point
      detail::loop_xxx(*_1 = *_2 + m_factor * *_3, m_buf.begin(), m_buf.end(), p.begin(), dir.begin());
      output_prec fv = m_functor(m_buf);
      // TODO: Could use amjito rule...
      if (fv < f0)
      {
        detail::loop_xx(*_1 += m_factor **_2 , p.begin(), p.end(), dir.begin());
        m_factor *= m_coeff;
      }
      else
      {
        do
        {
          m_factor /= m_coeff;
          if (m_factor < 128 * std::numeric_limits<input_prec>::epsilon())
          {
            std::cerr << "Bad descent direction" << std::endl;
            break;
          }
          detail::loop_xxx(*_1 = *_2 + m_factor * *_3, m_buf.begin(), m_buf.end(), p.begin(), dir.begin());
          fv = m_functor(m_buf);  
        } while (fv >= f0);
        detail::loop_xx(*_1 += m_factor **_2 , p.begin(), p.end(), dir.begin());
      }
      //std::cout << "Factor: " << m_factor << std::endl;
      return fv;
    }

  private: // data, input
  
    TFunctor m_functor;
    input_prec m_factor;
    input_prec m_coeff;
    
  private: // data, internal
  
    input_type m_buf;
  };


  //-----------------------------------------------------------------------------------------------------

    //-------------//
   //  FixedStep  //
  //-------------//

  /// Use a fixed step during the descent.
  /// This does not check if we are actually optimizing the criterion. Yet this can be useful in some
  /// cases, e.g. when doing alternate minimization, for which we know that it's not really worth to super-fine
  /// one part of the arguments while the other is fixed. So step should be small. OR, the size of the
  /// direction argument has actually already been chosen wisely. In this case a fixed step (say, of 1)
  /// could be used.
  template < typename TFunctor >
  class FixedStep
    : public line_minimization_algorithm_label
  {
  public: // typedefs  
    typedef typename TFunctor::argument_type            input_type;
    typedef typename value_type_of<input_type>::type    input_prec;
    typedef typename TFunctor::result_type              output_prec;

  public: // constructors & destructor
  
    explicit FixedStep(input_prec stepsize)
      : m_stepsize(stepsize)
    {}

  public: // operators
  
    output_prec operator()(input_type & p, const input_type & dir)
    {
      using namespace expr;
      detail::loop_xx(*_1 += *_2 * m_stepsize, p.begin(), p.end(), dir.begin());    
      //TODO: not sure if I should actually return sth. I guess I should change the design so that
      // operator() does not return anything, and we can get the value using getValue or something, 
      // that would not be implemented for FixedStep.
      return 0;
    }
    
  private: // data, input
  
    input_prec m_stepsize;
  };

  //-----------------------------------------------------------------------------------------------------

    //------------------//
   //  SmartFixedStep  //
  //------------------//

  /// Try a fix step, but if it doesn't minimize the criteria, use only a fraction of the step and 
  /// start all over again.
  template < typename TFunctor >
  class SmartFixedStep
    : public line_minimization_algorithm_label
  {
  public: // typedefs  
    typedef typename TFunctor::argument_type            input_type;
    typedef typename value_type_of<input_type>::type    input_prec;
    typedef typename TFunctor::result_type              output_prec;

  public: // constructors & destructor
  
    SmartFixedStep(TFunctor functor, input_prec stepsize, input_prec coeff = input_prec(0.5))
      : m_functor(functor)
      , m_stepsize(stepsize)
      , m_coeff(coeff)
    {}

  public: // operators
  
    output_prec operator()(input_type & p, const input_type & dir)
    {
      using namespace expr;

      resize(m_buf, p.size());
      // Compute the value at current point
      output_prec f0 = m_functor(m_buf);
      // Compute the value at estimated point
      input_prec mystep = m_stepsize;
      detail::loop_xxx(*_1 = *_2 + mystep * *_3, m_buf.begin(), m_buf.end(), p.begin(), dir.begin());
      output_prec fv = m_functor(m_buf);
      while (fv >= f0)
      {
        mystep *= m_coeff;
        std::cout << "mystep = " << mystep << std::endl;
        if (mystep < 128 * std::numeric_limits<input_prec>::epsilon())
        {
          std::cerr << "Bad descent direction" << std::endl;
          break;
        }
        detail::loop_xxx(*_1 = *_2 + mystep * *_3, m_buf.begin(), m_buf.end(), p.begin(), dir.begin());
        fv = m_functor(m_buf);
      }
      detail::loop_xx(*_1 += *_2 * mystep, p.begin(), p.end(), dir.begin());
      return fv;
    }
    
  private: // data, input  
    TFunctor m_functor;
    input_prec m_stepsize;
    input_prec m_coeff;
    
  private: // data, internal
    input_type m_buf;
  };

  //-----------------------------------------------------------------------------------------------------

    //-----------//
   //  LineMin  //
  //-----------//

  /// A class to minimize a multidimensional functional along a line.
  template < typename TFunctor >
  class LineMin
    : public line_minimization_algorithm_label
  {
  public: // typedefs  
    typedef typename TFunctor::argument_type            input_type;
    typedef typename value_type_of<input_type>::type    input_prec;
    typedef typename TFunctor::result_type              output_prec;
  
    typedef Mnbrak<TFunctor> Bracketer;
    typedef Brent<TFunctor> Minimizer;


  public: // constructors & destructor
  
    LineMin(TFunctor functor) :
        m_mnbrak(functor)
      , m_brent(functor)
      , m_verbose(false)
    {}

  public: // set & get
  
    /// Return minimizer
    Minimizer & minimizer() { return m_brent; }
    
    /// Return bracketer
    Bracketer & bracketer() { return m_mnbrak; }

    bool & verbose() { return m_verbose; }

  public: // operators
  
    output_prec operator()(input_type & p, const input_type & dir);
    
  private: // functors
    Mnbrak<TFunctor>  m_mnbrak;
    Brent<TFunctor>   m_brent;
    bool m_verbose;
  };


  //-----------------------------------------------------------------------------------------------------

  namespace detail
  {
    /// Simple code factoring.
    template < typename TFunctor >
    class IterativeMininizationAlgorithm_basis
    {
    public: // typedefs
    
      typedef typename TFunctor::argument_type            input_type;
      typedef typename value_type_of<input_type>::type    input_prec;
      typedef typename TFunctor::result_type              output_prec;
      
    public: // constructors

      explicit IterativeMininizationAlgorithm_basis(TFunctor functor)
        : m_functor(functor)
        , m_nitermax()
        , m_niter()
        , m_minstep()
      {
        m_ftol = std::max<output_prec>(std::sqrt(std::numeric_limits<input_prec>::epsilon()),
                                       std::sqrt(std::numeric_limits<output_prec>::epsilon()));
        m_minstep = std::sqrt(std::numeric_limits<input_prec>::epsilon());
      }
      
    public: // set & get
    
      // error tolerence on the function minimum
      output_prec       & ftol()       { return m_ftol; }
      output_prec const & ftol() const { return m_ftol; }

      // maximum number of iterations allowed
      unsigned int &       maxIter()        { return m_nitermax; }
      unsigned int const & maxIter() const  { return m_nitermax; }
  
      // number of iterations done
      unsigned int &       nIter()          { return m_niter; }
      unsigned int const & nIter() const    { return m_niter; }

      TFunctor &       functor()       { return m_functor; }
      TFunctor const & functor() const { return m_functor; }

      input_prec       & min_step()       { return m_minstep; }
      input_prec const & min_step() const { return m_minstep; }

    private: // data
    
      TFunctor m_functor;
      unsigned int m_nitermax;
      unsigned int m_niter;
      // error tolerence on the function minimum
      output_prec m_ftol;
      input_prec m_minstep;
    };


    template < typename TFunctor, typename TDFunctor >
    class IterativeGradMininizationAlgorithm_basis
      : public IterativeMininizationAlgorithm_basis<TFunctor>
    {
    public: // typedefs
      typedef IterativeMininizationAlgorithm_basis<TFunctor>  Basis;
      typedef typename Basis::input_type                      input_type;
      typedef typename Basis::input_prec                      input_prec;
      typedef typename Basis::output_prec                     output_prec;
    public: // constructors
      IterativeGradMininizationAlgorithm_basis(TFunctor functor, TDFunctor dfunctor)
        : Basis(functor)
        , m_dfunctor(dfunctor)
      {}
    public: // set & get
      TDFunctor       & dfunctor()       { return m_dfunctor; }
      TDFunctor const & dfunctor() const { return m_dfunctor; }
    private: // data
      TDFunctor m_dfunctor;
    };
  }

  //-----------------------------------------------------------------------------------------------------

    //----------//
   //  Powell  //
  //----------//
  
  /// Powell multidimensional minimization.
  /// NB: TFunctor should be a stl-compliant unary functor. Furthermore, it is assumed here
  /// that its argument is a random-access array of floating point numbers.
  /// E.g. arrays of 2D or 3D vectors, or arrays of complex numbers, are invalid argument types here.
  template < typename TFunctor, typename TLineMin = LineMin<TFunctor> >
  class Powell : public detail::IterativeMininizationAlgorithm_basis<TFunctor>
  {
  public: // typedefs
  
    typedef detail::IterativeMininizationAlgorithm_basis<TFunctor>    Basis;
    typedef typename Basis::input_type                                input_type;
    typedef typename Basis::input_prec                                input_prec;
    typedef typename Basis::output_prec                               output_prec;
  
  public: // constructors & destructor
  
    explicit Powell(TFunctor functor)
      : Basis(functor)
      , m_linemin(functor)
      , m_initStd()
      , m_dirs()
      , m_pmin()
    {
      this->init();
    }

  public: // init
  
    /// Default initialization.
    /// Set max number of iterations and and error tolerance to arbitrary but reasonable values.
    void init()
    {
      this->maxIter() = 100;
    }
    
  public: // set & get
  
    // intial estimate of range search along all directions
    std::vector<input_prec> & initStd() { return m_initStd; }

  public: // operators
  
    /// Minimize starting from p.
    // NB: the starting point is passed by value because we must have a copy anyway.
    input_type operator()(input_type p);
    
  private: // functions
    
    /// Set directions to the canonical n-basis
    void initDirections(std::size_t n);

    /// Minimize over all directions, and return the direction yielding the
    /// biggest decrease, together with the decrease
    void minOverAllDirections(typename std::vector<input_type>::iterator & iDirBig, output_prec & delta);
  
  private: // data, input
  
    // Line minimization algorithm
    TLineMin m_linemin;
    // Possible initialization of direction vector norm to take into account different scalings in 
    // input parameters
    std::vector<input_prec> m_initStd;
    
  private: // data, internals
    // current directions for minimization
    std::vector<input_type> m_dirs;
    // current estimate of function minimum
    output_prec m_fmin;
    // current estimate of function minimum parameters
    input_type m_pmin;
  };


  //-----------------------------------------------------------------------------------------------------

  /// First derivative estimator of a monodimensional functional.
  /// NB: Works only for real functions. TODO: could it be elegantly generalized other (e.g. complex) functions?
  template < typename TFunctor >
  class DerivativeEstimator
  {
  public: // typedefs
    typedef typename TFunctor::argument_type input_prec;
    typedef typename TFunctor::result_type output_prec;
  public: // constructors
    DerivativeEstimator(TFunctor functor) : m_functor(functor) {}
  public: // set & get

    // Set delta step in estimation formula
    void setDelta(input_prec delta)
    {
      assert(delta > this->mindelta());
      m_delta = delta;
      m_h = output_prec(1) / (2 * m_delta);
    }
    // Delta step
    const input_prec & delta() const { return m_delta; }
        
  public: // operator
    /// Estimate gradient at position x by computing ( f(x+delta) - f(x-delta) )/ (2*delta)
    output_prec operator()(input_prec x)
    { return m_h * (m_functor(x + m_delta) - m_functor(x - m_delta)); }
    
  private: // functions
    
    input_prec mindelta()
    {
      return 128 * std::max<input_prec>( std::numeric_limits< input_prec  >::epsilon(),
                                         std::numeric_limits< output_prec >::epsilon() );
    }
    
  private: // data, input
    TFunctor m_functor;
    input_prec m_delta;
    output_prec m_h;
  };

  /// Gradient estimator of a multidimensional functional.
  template < typename TFunctor >
  class GradientEstimator
  {
  public: // typedefs
  
    typedef typename TFunctor::argument_type            input_type;
    typedef typename value_type_of<input_type>::type    input_prec;
    typedef typename TFunctor::result_type              output_prec;

  public: // constructors
  
    
    GradientEstimator(TFunctor functor, input_prec delta)
      : m_functor(functor)
    {
      this->setDelta(delta);
    }
    
  public: // set & get
  
    // Set delta step in estimation formula
    void setDelta(input_prec delta)
    {
      assert(delta > this->mindelta());
      m_delta = delta;
      m_h = output_prec(1) / (2 * m_delta);
    }
    // Delta step
    const input_prec & delta() const { return m_delta; }
    
  public: // operator
  
    void operator()(input_type p, input_type & grad)
    {
      //std::cout << "[begin] Gradient estimation" << std::endl;
      std::size_t n = p.size();
      for (std::size_t i = 0; i < n; ++i)
      {
        input_type p1(p), p2(p);
        p1[i] += m_delta;
        p2[i] -= m_delta;
        grad[i] = m_h * ( m_functor(p1) - m_functor(p2) );
      }
      //std::cout << "[end] Gradient estimation" << std::endl;
    }
    
  private: // functions
    
    input_prec mindelta()
    {
      return 128 * std::max<input_prec>( std::numeric_limits< input_prec  >::epsilon(),
                                         std::numeric_limits< output_prec >::epsilon() );
    }

  private: // data, input
  
    TFunctor m_functor;
    input_prec m_delta;
    output_prec m_h;
  };


  //-----------------------------------------------------------------------------------------------------

    //-------------------//
   //  GradientDescent  //
  //-------------------//

  /// Gradient descent minimization.
  template < typename TFunctor, typename TGradFunctor, typename TLineMin = LineMin<TFunctor> >
  class GradientDescent
    : public detail::IterativeGradMininizationAlgorithm_basis<TFunctor, TGradFunctor>
  {
  public: // typedefs
  
    typedef detail::IterativeGradMininizationAlgorithm_basis<TFunctor, TGradFunctor>   Basis;
    typedef typename Basis::input_type                      input_type;
    typedef typename Basis::input_prec                      input_prec;
    typedef typename Basis::output_prec                     output_prec;

  public: // constructors
  
    GradientDescent(TFunctor functor, TGradFunctor gradfunctor, TLineMin linemin);
    
  public: // operators
  
    input_type operator()(input_type p);
    
  private: // functions
  
    void init()
    {
      this->maxIter() = 100;
    }
    
  private: // data
  
    // Line minimization algorithm
    TLineMin m_linemin;
    // current estimate of function minimum
    output_prec m_fmin;    
    input_type m_grad;
  };

  //-----------------------------------------------------------------------------------------------------

    //-----------------------//
   //  PRConjugateGradient  //
  //-----------------------//

  /// Polak-Ribiere conjugate gradient minimization.
  /// NB: right now, the default linemin functor is not optimal for general CG, as it doesn't use the gradient.
  /// Although it's a suitable choice for some (my:-) settings.
  // TODO: for all gradient-based algorithm, I think we may want to add a check on the (max) norm of the gradient
  // to stop the algorithm if it is too small.
  // TODO: un truc qui pourrait etre utile, c'est pour le bracket, partir avec des valeurs qui prennent en compte les resultats
  // precedents. Par exemple, si on a trouve lambda, ben on bracket a partir de 20*max des lambdas des 3 dernieres iterations.
  // Et en plus on pourrait en profiter pour changer les criteres d'arrets pour etre sur qu'on s'arrete pas uniquement parce
  // que le bracket initial est trop grand par rapport au lambda optimal, comme ca arrive parfois.
  template < typename TFunctor, typename TGradFunctor, typename TLineMin = LineMin<TFunctor> >
  class PRConjugateGradient
    : public detail::IterativeGradMininizationAlgorithm_basis<TFunctor, TGradFunctor>
  {
  public: // typedefs
  
    typedef detail::IterativeGradMininizationAlgorithm_basis<TFunctor, TGradFunctor>   Basis;
    typedef typename Basis::input_type                      input_type;
    typedef typename Basis::input_prec                      input_prec;
    typedef typename Basis::output_prec                     output_prec;
  
  public: // constructors
  
    PRConjugateGradient(TFunctor functor, TGradFunctor gradfunctor, TLineMin linemin)
      : Basis(functor, gradfunctor)
      , m_linemin(linemin)
    {
      this->init();
    }

  public: // operators
  
    input_type operator()(input_type p);
    
  private: // functions
  
    void init() { this->maxIter() = 100; }
    
  private: // data
  
    // Line minimization algorithm
    TLineMin m_linemin;
    // current estimate of function minimum
    output_prec m_fmin;
    input_type m_grad;
    input_type m_conjgrad;
  };
  
  
  //-----------------------------------------------------------------------------------------------------

  namespace detail
  {
    template < class T, template < typename > class C >
    class linear_map_base
    {
    public: // typedefs
      typedef C < T > Container;
    public: // constructors
      linear_map_base(Container & data) : m_data(data) {}
    public: // set & get
      Container & get() { return m_data; }
      Container const & get() const { return m_data; }
    private: // data
      Container & m_data;
    };
  }

  /// Linear mapping of a container.
  /// Useful for optimization
  /// Note that this should not be the only way to pass arguments. E.g one could imagine to minimize
  /// on a vector and a matrix and thus needing help for mapping even though their is no
  /// outer container.
  template < class T, template < typename > class C >
  class linear_map
    : public detail::linear_map_base<T,C>
  {
  public: // typedefs
    typedef detail::linear_map_base<T,C> Base;
    typedef typename Base::Container Container;
  public: // operators
    T operator[](std::size_t i) { return this->get()[i]; }
  public: // constructors
    linear_map(Container & data) : Base(data) {}
  };
  
  template < class T, std::size_t D, template < typename > class C >
  class linear_map< til::numeric_array<T, D>, C >
    : public detail::linear_map_base<til::numeric_array<T, D>, C>
  {
  public: // typedefs
    typedef detail::linear_map_base<til::numeric_array<T, D>, C> Base;
    typedef typename Base::Container Container;
  public: // constructors
    linear_map(Container const & data) : Base(data) {}
  public: // functions
    T operator[](std::size_t i) { return this->get()[i/D][i%D]; }
  };
}


// Package includes
#include "minTools.tpp"

#endif /*MINTOOLS_H_*/
