#ifndef TIL_ACCUMULATOR_H_
#define TIL_ACCUMULATOR_H_

// includes from TIL
#include "til/cat2type.h"
#include "til/functors.h"
//#include "globalTraits.h"


namespace til
{

  //-----------------------------------------------------------------------------------------------------

    //---------------------------------//
   //  Standard accumulator policies  //
  //---------------------------------//

  namespace policy
  {
    /// Record nothing while accumulating values.
    struct AccumulatorRecord_None
    {
      template < typename T >
      void value(T) {}
   
      template < typename TIterator >
      void range(TIterator, TIterator) {}
   
      template < typename TContainer >
      void container(TContainer &) {}
    };

    /// This Record policy counts the number of elements that have been added.
    /// Usefull for computing statistics.
    class AccumulatorRecord_Sum
    {
    public: // constructors & destructor
    
      AccumulatorRecord_Sum() { m_count = 0; }
    
    public: // set & get
    
      unsigned int get() const { return m_count; }
    
    public: // functions
      
      template < typename T >
      void value(T)
      { ++m_count; }
   
      template < typename TIterator >
      void range(const TIterator & start, const TIterator & finish)
      { m_count += std::distance(start, finish); }
   
      template < typename TContainer >
      void container(const TContainer & c)
      { m_count += size(c); }
      
    private:
      unsigned int m_count;
    };
  }


  //-----------------------------------------------------------------------------------------------------

    //---------------//
   //  Accumulator  //
  //---------------//
  
  /// A class to accumulate value.
  /// What is exactly meant by accumulate is actually left to the user, via the AccumulationPolicy.
  /// Note that the interface of this policy is actually that of a inplace functor (i.e.
  /// a functor that do an inplace operation and returns void, such as AddTo -- this is
  /// necessary for objects, because things like A = A+B cannot be simplified for say
  /// vectors, and it is thus highly inefficient.
  /// So the use of AddTo as the accumulation trait should be enough to make it
  /// a sum accumulator.
  template < typename T, typename TAccumulation, typename AccumulationPolicy, typename RecordPolicy = policy::AccumulatorRecord_None>
  class Accumulator
  {
  public: // constructors
  
    Accumulator() : 
        m_accumulation()
      , m_accumulationPolicy()
      , m_recordPolicy()
    {}

  public: // set & get
    
    /// Get accumulation result.
    const TAccumulation & get() const { return m_accumulation; }
    /// Get record policy.
    const RecordPolicy & recordPolicy() const { return m_recordPolicy; }

    /// clear accumulated values
    void clear() { m_accumulation.clear(); }

  public: // functions
    
    /// Accumulate a value.
    void accumulate(typename boost::call_traits<T>::param_type value)
    { m_accumulationPolicy(m_accumulation, value); }
    
    /// Accumulate values spanned by given range.
    template < typename TIterator >
    // ensure that the iterator points on the same type as the histo type
    typename boost::enable_if<boost::is_same<typename value_type_of<TIterator>::type, T> >::type
    accumulate(TIterator begin, TIterator end)
    {
      m_recordPolicy.range(begin, end);
      for (; begin != end; ++begin)
      {
        // TODO: We could actually propose to do something when the containee is also a container, also
        // it is not necessarily a good idea.
        //this->_accumulateValue(*begin);
        this->accumulate(*begin);
      }
    }
    
  private: // data
  
    TAccumulation         m_accumulation;
    AccumulationPolicy    m_accumulationPolicy;
    RecordPolicy          m_recordPolicy;
  };
  

  //-----------------------------------------------------------------------------------------------------

    //-------------------//
   //  MeanAccumulator  //
  //-------------------//
  
  /// A class to accumulate values and return their mean.
  template < typename T, typename TAccumulation >
  class MeanAccumulator : 
    public Accumulator<T, TAccumulation, functor::AddTo<TAccumulation, T>, policy::AccumulatorRecord_Sum >
  {
  public: // typedefs 

    typedef MeanAccumulator<T, TAccumulation> Self;
    typedef Accumulator<T, TAccumulation, functor::AddTo<TAccumulation, T>, policy::AccumulatorRecord_Sum > Base;

  public: // set & get
  
    /// Get mean of accumulated values.
    TAccumulation get()
    {
      return this->Base::get() / this->recordPolicy().get();
    }
  };
  
} // namespace til

#endif /*ACCUMULATOR_H_*/
