#ifndef TIL_MINTOOLS_TPP
#define TIL_MINTOOLS_TPP

namespace til
{
  
  
  //---------------------------------------------------------------------------
  
  /// A class to find a bracketing triplet around a function minimum.
  template < typename TFunctor >
  void Mnbrak<TFunctor>::operator()(input_prec & ax, input_prec & bx, input_prec & cx, input_type point, input_type dir)
  {
    if (m_verbose) std::cout << "[ BEGIN Bracketing ]" << std::endl;
    
    using namespace expr;
    // check that input dimension are equal
    assert(size(point) == size(dir));
    // resize buffer accordingly
    resize(m_buf, size(point));
    
    input_prec q, r, u, ulim;
    output_prec fu;
    
    detail::loop_xxx(*_1 = *_2 + ax * *_3, m_buf.begin(), m_buf.end(), point.begin(), dir.begin());
    if (m_verbose) std::cout << "lambda = " << ax << std::endl;
    m_fa = m_functor(m_buf);

    //m_fa = m_functor(point + ax * dir);
    detail::loop_xxx(*_1 = *_2 + bx * *_3, m_buf.begin(), m_buf.end(), point.begin(), dir.begin());
    if (m_verbose) std::cout << "lambda = " << bx << std::endl;
    m_fb = m_functor(m_buf);
    //m_fb = m_functor(point + bx * dir);

    // swap A and B so that we go downhill from A to B
    if (m_fb > m_fa)
    {
      std::swap(ax, bx);
      std::swap(m_fa, m_fb);
    }
    // First guess for c
    cx = bx + GOLD * (bx - ax);

    detail::loop_xxx(*_1 = *_2 + cx * *_3, m_buf.begin(), m_buf.end(), point.begin(), dir.begin());
    if (m_verbose) std::cout << "lambda = " << cx << std::endl;
    m_fc = m_functor(m_buf);
    //m_fc = m_functor(point + cx * dir);
    // Keep estimating until we have a bracketing triplet
    while (m_fb > m_fc)
    {
      r = ( bx - ax ) * ( m_fb - m_fc );
      q = ( bx - cx ) * ( m_fb - m_fa );
      // parabolic fit estimate
      u = bx - ((bx-cx)*q - (bx-ax)*r) / (2.0 * sign(max(std::abs(q-r), EPSILON), q-r));
      ulim = bx + GLIMIT * ( cx - bx );
      // Check whether u is between b and c
      if (same_sign(bx-u, u-cx))
      {
        detail::loop_xxx(*_1 = *_2 + u * *_3, m_buf.begin(), m_buf.end(), point.begin(), dir.begin());
        if (m_verbose) std::cout << "lambda = " << u << std::endl;
        fu = m_functor(m_buf);
        //fu = m_functor(point + u * dir);
        // If the value is smaller than fc, then [b,u,c] is a bracket
        if (fu < m_fc)
        {
          // rename [b u] as [a b], so that we return bracket [a b c]
          shift_right(u, bx, ax);
          shift_right(fu, m_fb, m_fa);
          if (m_verbose) std::cout << ax << " " << bx << " " << cx << std::endl;
          if (m_verbose) std::cout << "[ END Bracketing ]" << std::endl;
          return;
        }
        // If the value is bigger than fb, then [a,b,u] is a bracket
        // NB: This probably happens very seldom, since u is the minimum parabolic fit.
        else if (fu > m_fb)
        {
          // rename u as b so that we return bracket [a b c ]
          cx = u;
          m_fc = fu;
          if (m_verbose) std::cout << ax << " " << bx << " " << cx << std::endl;
          if (m_verbose) std::cout << "[ END Bracketing ]" << std::endl;
          return;
        }
        
        // NB: in all other cases, it will be assumed that [b c u] will be our next triplet,
        // be it a bracket or not, therefore we don't have special treatments and returns.
        
        // If parabolic fit did not bracket, do a golden ratio [b - c -- *]
        u = cx + GOLD * (cx - bx);
        detail::loop_xxx(*_1 = *_2 + u * *_3, m_buf.begin(), m_buf.end(), point.begin(), dir.begin());
        if (m_verbose) std::cout << "lambda = " << u << std::endl;
        fu = m_functor(m_buf);
        //fu = m_functor(point + u * dir);
      }
      // Check that u, the parabolic trial, is outside c but still before ulim
      else if (same_sign(cx-u, u-ulim))
      {
        detail::loop_xxx(*_1 = *_2 + u * *_3, m_buf.begin(), m_buf.end(), point.begin(), dir.begin());
        if (m_verbose) std::cout << "lambda = " << u << std::endl;
        fu = m_functor(m_buf);
        //fu = m_functor(point + u * dir);
        // if fu < fc, do a golden ration [c - u -- *], and  get rid of b.
        if (fu < m_fc)
        {
          shift_right(u + GOLD * (u-cx), u, cx, bx);
          detail::loop_xxx(*_1 = *_2 + u * *_3, m_buf.begin(), m_buf.end(), point.begin(), dir.begin());
          if (m_verbose) std::cout << "lambda = " << u << std::endl;
          shift_right(m_functor(m_buf), fu, m_fc, m_fb);
          //shift_right(m_functor(point + u * dir), fu, m_fc, m_fb);
        }
        // (else, we have a bracket)
      }
      // Check that u is beyond ulim and cx
      else if (same_sign(u - ulim, ulim - cx))
      {
        u = ulim;
        detail::loop_xxx(*_1 = *_2 + u * *_3, m_buf.begin(), m_buf.end(), point.begin(), dir.begin());
        if (m_verbose) std::cout << "lambda = " << u << std::endl;
        fu = m_functor(m_buf);
        //fu = m_functor(point + u * dir);
      }
      // If we are here, it means that u is before b, which does not make sense:
      // do a golden search.
      else
      {
        u = cx + GOLD*(cx-bx);
        detail::loop_xxx(*_1 = *_2 + u * *_3, m_buf.begin(), m_buf.end(), point.begin(), dir.begin());
        if (m_verbose) std::cout << "lambda = " << u << std::endl;
        fu = m_functor(m_buf);
        //fu = m_functor(point + u * dir);
      }
      
      // rename [b c u] into [a b c]
      
      shift_right(u, cx, bx, ax);
      shift_right(fu, m_fc, m_fb, m_fa);
    }
    if (m_verbose) std::cout << ax << " " << bx << " " << cx << std::endl;
    if (m_verbose) std::cout << "[ END Bracketing ]" << std::endl;
  }

  //---------------------------------------------------------------------------
  
  template < typename TFunctor >
  const typename Mnbrak<TFunctor>::input_prec Mnbrak<TFunctor>::GOLD = (std::sqrt(5.0) + 1.0) / 2.0;
  template < typename TFunctor >
  const typename Mnbrak<TFunctor>::input_prec Mnbrak<TFunctor>::EPSILON = 128 * std::numeric_limits<input_prec>::epsilon();
  template < typename TFunctor >
  const typename Mnbrak<TFunctor>::input_prec Mnbrak<TFunctor>::GLIMIT = 100;
  
  //---------------------------------------------------------------------------

  template < typename TFunctor >
  Brent<TFunctor>::Brent()
    : m_functor()
    , m_xmin()
    , m_fmin()
    , m_tol()
    , m_itmax()
  {
    this->init();
  }

  //---------------------------------------------------------------------------
  
  template < typename TFunctor >
  Brent<TFunctor>::Brent(TFunctor functor)
    : m_functor(functor)
    , m_xmin()
    , m_fmin()
    , m_tol()
    , m_itmax()
  {
    this->init();
  }
  
  //---------------------------------------------------------------------------
  
  /// A class to minimize a functional within a bracketing triplet.
  template < typename TFunctor >
  void Brent<TFunctor>::operator()
  (
    input_prec a
  , input_prec m
  , input_prec b
  , input_type point
  , input_type dir
  )
  {
    if (m_verbose) std::cout << "[ BEGIN Parabolic linemin ]" << std::endl;
    using namespace expr;
    typedef typename combine<input_prec, output_prec>::type max_prec;
    resize(m_buf, size(point));
    // The last and before last estimates of the minimum
    input_prec w, v;
    // Their corresponding function values
    output_prec fw, fv, fu;
    // The center of [a, b]
    input_prec x_med;
    // Distance moved on the last step last
    input_prec stepLength = 0.0;
    // Distance moved on the step before last
    input_prec stepLength2 = 0.0;
    input_prec u;
    max_prec r, p, q, tmp;
    max_prec tol1, tol2;
    // Swap bracket if necessary
    if (a > b) std::swap(a,b);
    // Initialize all estimated minima at the bracket point m
    m_xmin = w = v = m;
    detail::loop_xxx(*_1 = *_2 + m_xmin * *_3, m_buf.begin(), m_buf.end(), point.begin(), dir.begin());
    if (m_verbose) std::cout << "lambda = " << m_xmin << std::endl;
    fw = fv = m_fmin = m_functor(m_buf);
    //fw = fv = m_fmin = m_functor(point + m_xmin * dir);

    // main minimization loop
    for (int iter = 0; iter < m_itmax; ++iter)
    {
      x_med = 0.5 * ( a + b );
      //tol1 = m_tol * std::abs(m_xmin) + ZEPS;
      tol1 = m_tol;
      tol2 = 2 * tol1;
      // Convergence test
      //if (std::abs(m_xmin-x_med) <= tol2)
      if (std::abs(m_xmin-x_med) <= (tol2-input_prec(0.5)*(b-a)))
      {
        if (m_verbose) std::cout << "[ END Parabolic linemin ]" << std::endl;
        return;
      }

      // Parabolic fit
      if (std::abs(stepLength2) > tol1)
      {
        r = (m_xmin-w) * (m_fmin-fv);
        q = (m_xmin-v) * (m_fmin-fw);
        p = (m_xmin-v) * q - (m_xmin-w)*r;
        q = max_prec(2.0) * (q-r);
        if (q > 0) p = -p;
        q = std::abs(q);
        tmp = stepLength2;
        stepLength2 = stepLength;
        // Is parabolic fit acceptable?
        if (std::abs(p) >= std::abs(max_prec(0.5)*q*tmp) || p <= q * (a-m_xmin) || p >= q * (b-m_xmin))
        {
          // If not, take golden section step instead
          //stepLength2 = ( m_xmin >= x_med ? a - m_xmin : b - m_xmin );
          //stepLength = CGOLD * stepLength2;
          stepLength = CGOLD * ( m_xmin >= x_med ? a - m_xmin : b - m_xmin );
        }
        else
        {
          // Accept parabolic fit
          stepLength = p / q;
          /*
          u = m_xmin + stepLength;
          if (u-a < tol2 || b-u < tol2)
            stepLength = SIGN(tol1, x_med-m_xmin);
          */
        }
      }
      else
      {
          //stepLength2 = ( m_xmin >= x_med ? a - m_xmin : b - m_xmin );
          //stepLength = CGOLD * stepLength2;
          stepLength = CGOLD * ( m_xmin >= x_med ? a - m_xmin : b - m_xmin );
      }
      
      u = (std::abs(stepLength) >= tol1 ? m_xmin+stepLength : m_xmin+sign(tol1,stepLength));
      // This is the one function evaluation per iteration
      detail::loop_xxx(*_1 = *_2 + u * *_3, m_buf.begin(), m_buf.end(), point.begin(), dir.begin());
      if (m_verbose) std::cout << "lambda = " << u << std::endl;
      fu = m_functor(m_buf);
      //fu = m_functor(point + u * dir);
      // Now decide what to do with out function evaluation
      if (fu <= m_fmin)
      {
        if (u >= m_xmin) a = m_xmin; else b = m_xmin;
        shift_right(u, m_xmin, w, v);
        shift_right(fu, m_fmin, fw, fv);
      }
      else
      {
        if (u < m_xmin) a = u; else b = u;
        if (fu <= fw || w == m_xmin)
        {
          shift_right(u,w,v);
          shift_right(fu,fw,fv);
        }
        else if (fu <= fv || v == m_xmin || v == w)
        {
          v = u;
          fv = fu;
        }
      }
    }
    std::cerr << "Warning: Brent reached maximum iteration limit of " << m_itmax << std::endl;
    if (m_verbose) std::cout << "[ END Parabolic linemin ]" << std::endl;
    //return m_xmin;
  }

  //---------------------------------------------------------------------------
  
  template < typename TFunctor > 
  const typename Brent<TFunctor>::input_prec Brent<TFunctor>::CGOLD = (3.0 - std::sqrt(5.0)) / 2.0;
  
  //---------------------------------------------------------------------------

  template < typename TFunctor >
  typename LineMin<TFunctor>::output_prec
  LineMin<TFunctor>::operator()
  (
    input_type & p,
    const input_type & dir
  )
  {
    using namespace expr;
    input_prec a, m, b;
    a = input_prec(0);
    m = input_prec(1);
    m_mnbrak(a, m, b, p, dir);
    m_brent(a, m, b, p, dir);
    if (m_verbose) std::cout << "linemin: " << m_brent.xmin() << " " << m_brent.fmin() << std::endl;
    detail::loop_xx(*_1 += *_2 * m_brent.xmin(), p.begin(), p.end(), dir.begin());
    //p += dir * m_brent.xmin();
    return m_brent.fmin();
  }


  //---------------------------------------------------------------------------

  template < typename TFunctor, typename TLineMin >
  typename Powell<TFunctor, TLineMin>::input_type Powell<TFunctor, TLineMin>::operator()(input_type p)
  {
    using namespace expr;
    
    const output_prec EPSILON = 128 * std::numeric_limits<output_prec>::epsilon();
    
    // function value at point p, before line min
    output_prec fp;
    output_prec delta, fpext, tmp;

    input_type pt, pext, xt, dir;
    resize(pt, size(p));
    resize(pext, size(p));
    resize(xt, size(p));
    resize(dir, size(p));
    
    resize(m_pmin, size(p));
    
    typename std::vector<input_type>::iterator iDirBig;
    
    //int ibig;
    std::copy(p.begin(), p.end(), m_pmin.begin());
    //m_pmin = p;
    m_fmin = this->functor()(m_pmin);

    int n = size(m_pmin);
    
    this->initDirections(n);
    //resize(pext, n);
    
    //pt = m_pmin;

    for (this->nIter() = 0; this->nIter() < this->maxIter(); ++this->nIter())
    {
      //std::cout << "Iteration " << m_niter << " start from " << m_fmin << " : " << m_pmin << std::endl;
      std::cout << "Iteration " << this->nIter() << " start from " << m_fmin << std::endl;
      
      // save information from previous iteration
      fp = m_fmin;
      pt = m_pmin;
      
      // minimize over all directions  
      this->minOverAllDirections(iDirBig, delta);

      // Stopping criterion:
      // If we didn't minimize that much, return
      // if (fp - m_fmin < ftol) return m_pmin;
        if (2.0*(fp-m_fmin) <= this->ftol() * (std::abs(fp) + std::abs(m_fmin)) + EPSILON)
        {
          return m_pmin;
        }
 
        // Get information at an extrapolated point, twice as far as m_pmin from previous position.
      // This will help us to decide whether we want to update the set of directions or not
      detail::loop_xxx(*_1 = input_prec(2) * *_2 - *_3, pext.begin(), pext.end(), m_pmin.begin(), pt.begin());
      //pext = input_prec(2) * m_pmin - pt;
      fpext = this->functor()(pext);

      // Update directions only if minimization along the best direction is not "played out".
      if (fpext < fp)
      {
        // another magic criterion, that detect either if the biggest direction indeed contributed
        // much more than the others, or that we are near the bottom of the parabole in this
        // direction
        tmp = output_prec(2) * ( fp - output_prec(2) * m_fmin + fpext ) * square( fp - m_fmin - delta ) - delta * square( fp - fpext );
        if (tmp < 0)
        {
          detail::loop_xxx( *_1 = *_2 - *_3, dir.begin(), dir.end(), m_pmin.begin(), pt.begin());
          //dir = m_pmin - pt;
          m_fmin = m_linemin(m_pmin, dir);
          /*
          m_dirs[ibig] = m_dirs[n];
          m_dirs[n] = dir;
          */
          *iDirBig = dir;
          //m_dirs[ibig] = dir;
        }
      }
    } while(1);

    std::cerr << "Warning: Powell reached maximum iteration limit of " << this->maxIter() << std::endl;
    return m_pmin;
  }

  //---------------------------------------------------------------------------
  
  template < typename TFunctor, typename TLineMin >
  void Powell<TFunctor, TLineMin>::initDirections(std::size_t n)
  {
    m_dirs.resize(n);
    if (m_initStd.size())
    {
      for (std::size_t i = 0; i != n; ++i)
      {
        resize(m_dirs[i], n);
        std::fill(m_dirs[i].begin(), m_dirs[i].end(), input_prec(0));
        m_dirs[i][i] = m_initStd[i];
      }
    }
    else
    {
      for (std::size_t i = 0; i != n; ++i)
      {
        resize(m_dirs[i], n);
        std::fill(m_dirs[i].begin(), m_dirs[i].end(), input_prec(0));
        m_dirs[i][i] = input_prec(1);
      }
    }
  }

  //---------------------------------------------------------------------------

  template < typename TFunctor, typename TLineMin >
  void Powell<TFunctor, TLineMin>::minOverAllDirections(typename std::vector<input_type>::iterator & iDirBig, output_prec & delta)
  {
    delta = output_prec(-1);
    output_prec fpext;
    int count = 0;
    for (typename std::vector<input_type>::iterator iDir = m_dirs.begin(); iDir != m_dirs.end() ; ++iDir)
    {
      fpext = m_fmin;
      m_fmin = m_linemin(m_pmin, *iDir);
      //std::cout << "dirmin  " << count++ << " " << m_fmin << " : " << m_pmin << std::endl;
      std::cout << "dirmin  " << count++ << " " << m_fmin << std::endl;
      // Keep this direction if it was the biggest decrease
      if (fpext - m_fmin > delta)
      {
        delta = fpext - m_fmin;
        //ibig = i
        iDirBig = iDir;
      }
    }
  }

  //---------------------------------------------------------------------------
  
  // TODO: rewrite all these iterative minimization as steps only, let the user loop himself, that
  // allows for incredible flexibility, for printing, debugging, whatever.
  template < typename TFunctor, typename TGradFunctor, typename TLineMin >
  typename PRConjugateGradient<TFunctor, TGradFunctor, TLineMin>::input_type
  PRConjugateGradient<TFunctor, TGradFunctor, TLineMin>::operator()(input_type p)
  {
    const output_prec EPSILON = 128 * std::numeric_limits<output_prec>::epsilon();
    
    std::size_t n = p.size();
    
    // We use the function version of resize, and not the member function, because input_type might
    // not have this member function.
    resize(m_grad, n);
    resize(m_conjgrad, n);
    
    // temporary variable to hold gradient
    // TODO: right now, gradtmp is not allocated, this task is left to the gradient
    // algo. I think it sucks.
    input_type gradtmp;
    output_prec gg, dgg, gam;
    output_prec fp;
    output_prec dgg2;

    // initialization
    fp = this->functor()(p);
    this->dfunctor()(p, m_grad);
    std::for_each(m_grad.begin(), m_grad.end(), std::negate<input_prec>());
    std::copy(m_grad.begin(), m_grad.end(), m_conjgrad.begin());
    
    for (this->nIter() = 0; this->nIter() < this->maxIter(); ++this->nIter())
    {
      std::cout << "CG iter " << this->nIter(); 
            //<< " : ";
      //std::copy(p.begin(), p.end(), std::ostream_iterator<double>(std::cout, " "));
      std::cout << std::endl;
      
      input_type p_tmp = p;
      // line minimization along conjugate direction
      m_fmin = m_linemin(p, m_conjgrad);
      if (2*std::abs(fp - m_fmin) <= this->ftol() * (std::abs(m_fmin) + std::abs(fp)) + EPSILON)
      {
        //return m_pmin;
        return p;
      }
      // Stopping criterion based on maximum motion
      if (this->min_step() > 0.0)
      {
        input_prec mx = 0;
        for (std::size_t i = 0; i < n; ++i)
        {
          if (mx < std::abs(p[i] - p_tmp[i]))
          {
            mx = std::abs(m_conjgrad[i]);
          }
        }
        std::cout << "Max motion: " << mx << std::endl;
        if (mx < this->min_step())
        {
          return p;
        }
      }
      fp = m_fmin;
      this->dfunctor()(p, gradtmp);
      gg = dgg = 0;
      dgg2 = 0;
      for (std::size_t i = 0; i < n; ++i)
      {
        gg += square(m_grad[i]);
        dgg += (gradtmp[i]+m_grad[i]) * gradtmp[i];
        dgg2 += square(gradtmp[i]);
      }
      if (gg < EPSILON)
      {
        //return m_pmin;
        return p;
      }
      //std::cout << "dgg: " << dgg << " , dgg2: " << dgg2 << std::endl;
      gam = dgg / gg;
      for (std::size_t i = 0; i < n; ++i)
      {
        m_grad[i] = -gradtmp[i];
        m_conjgrad[i] = m_grad[i] + gam * m_conjgrad[i];
      }
    }
    std::cout << "max iter reached" << std::endl;
    return p;
  }
  

  //---------------------------------------------------------------------------

  template < typename TFunctor, typename TGradFunctor, typename TLineMin >
  GradientDescent<TFunctor, TGradFunctor, TLineMin>::GradientDescent
  (
    TFunctor      functor
  , TGradFunctor  gradfunctor
  , TLineMin      linemin
  )
    : Basis(functor, gradfunctor)
    , m_linemin(linemin)
  {
    this->init();
  }

  //---------------------------------------------------------------------------
  
  template < typename TFunctor, typename TGradFunctor, typename TLineMin >
  typename GradientDescent<TFunctor, TGradFunctor, TLineMin>::input_type
  GradientDescent<TFunctor, TGradFunctor, TLineMin>::operator()(input_type p)
  {
    const output_prec EPSILON = 128 * std::numeric_limits<output_prec>::epsilon();
    
    // temporary variable to hold gradient
    input_type gradtmp;
    output_prec fp;
  
    // initialization
    fp = this->functor()(p);
    
    for (this->nIter() = 0; this->nIter() < this->maxIter(); ++this->nIter())
    {
      std::cout << "GD iter " << this->nIter() << std::endl; // << " : " << p << std::endl;
      // compute gradient
      this->dfunctor()(p, m_grad);
      // negate gradient to point toward smaller values
      std::for_each(m_grad.begin(), m_grad.end(), std::negate<input_prec>());
      // line minimization along gradient direction
      m_fmin = m_linemin(p, m_grad);
      
      //if(0)
      if (2*std::abs(fp - m_fmin) <= this->ftol() * (std::abs(m_fmin) + std::abs(fp)) + EPSILON)
      {
        //return m_pmin;
        return p;
      }
      fp = m_fmin;
    }
    return p;
  }

  //---------------------------------------------------------------------------
    
} // namespace til

#endif /*MINTOOLS_H_*/
