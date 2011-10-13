#ifndef TIL_PROGRESS_INDICATOR_H_
#define TIL_PROGRESS_INDICATOR_H_

// includes from STL
#include <ctime>
#include <iostream>

// includes from Linux
//#include <sys/time.h> // gettimeofday

namespace til
{
  

  namespace callback
  {

    //----------------------------------------------------------------------------------------------


    /// ProgressIndicator callback that simply prints percentage.
    struct PI_Print
    {
      void operator()(std::size_t current, std::size_t total)
      {
        std::cout << 100 * current / total << "%..." << std::flush;
      }
    };


    //----------------------------------------------------------------------------------------------

      //------------//
     //  PI_Timer  //
    //------------//

    /// ProgressIndicator callback that prints percentage if computation is taking time
    template < typename TPICallback >
    class PI_Timer
    {
    public: // constructors

      PI_Timer()
        : m_start()
        , m_finish()
      {
        this->setDelay(2.0);
        this->init();
      }

      explicit PI_Timer(double delay)
        : m_delay(delay)
        , m_start()
        , m_finish()
      { this->init(); }

    public: // initialization
    
      void init()
      {
        m_active = false;
      }

    public: // set & get

      /// Set delay, in seconds
      void setDelay(double s) { m_delay = s * CLOCKS_PER_SEC; }
      /// Access callback object
      TPICallback & callback() { return m_callback; }

    public: // operators

      void operator()(std::size_t current, std::size_t total)
      {
        if (current == 0)
        {
          //gettimeofday(&start, &tz);
          m_start = std::clock();
          return;
        }
        else if (current == 1)
        {
          //gettimeofday(&finish, &tz);
          //double t = (finish.tv_sec-start.tv_sec) * 1000000L + (finish.tv_usec-start.tv_usec);
          m_finish = std::clock();
          double t = double(m_finish - m_start);
          if (t > m_delay)
          {
            m_active = true;
          }
        }

        if (m_active)
        {
          m_callback(current, total);
        }
      }
    
    private: // data, input

      double m_delay;
      TPICallback m_callback;

    private: // data, internal;

      //struct timeval m_start, m_finish;
      //struct timezone m_tz;
      std::clock_t m_start, m_finish;
      bool m_active;
    };
  } // namespace callback


  //----------------------------------------------------------------------------------------------
  

    //---------------------//
   //  ProgressIndicator  //
  //---------------------//

  /// A simple base class used to give progression feedback.
  /// You give the number of iterations, and the number of parts you want to split these iterations in. Then, 
  /// your callback is called when a sufficient number of iterations have been done.
  template < typename TPICallback >
  class ProgressIndicator
  {
  public: // constructors

    /// Give the total number of iterations, and the number of parts you want to divide this total by.
    ProgressIndicator(std::size_t parts, std::size_t total);
  
  public: // initialization
  
    void init();
    
  public: //set & get
  
    std::size_t & total() { return m_total; }
    std::size_t & parts() { return m_parts; }
    TPICallback & callback() { return m_callback; }
    
  public: // operators
    
    /// Call the callback if the iteration number i has reached the next part.
    void operator()(std::size_t i);

  private: // functions  
    void nextPart();

  private: // data, input
    std::size_t m_total;
    std::size_t m_parts;
    TPICallback m_callback;

  private: // data, internal
    double m_mark;
    std::size_t m_currentpart;
  };
  
  
  typedef ProgressIndicator<callback::PI_Timer<callback::PI_Print> > DelayedProgressIndicator;


  //----------------------------------------------------------------------------------------------

    //---------------//
   //  NoIndicator  //
  //---------------//
  
  /// Empty indicator class.
  struct NoIndicator
  {
    void operator()(std::size_t) {}
  };
  
} // namespace til

// package include
#include "progress_indicator.tpp"

#endif /*PROGRESS_INDICATOR_H_*/
