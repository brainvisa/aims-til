
namespace til
{
  
  //---------------------------------------------------------------------------

  template < typename T, typename TCount >
  const T & 
  rle_array<T,TCount>::
  operator[](std::size_t i) const
  {
    assert(i <= this->size());
    std::size_t current_i = 0;
    typename Data::const_iterator iData = m_data.begin();
    while (i >= iData->length())
    {
      i -= iData->length();
      ++iData;
    }
    return iData->value();
  }

  //---------------------------------------------------------------------------

  template < typename T, typename TCount >
  typename rle_array<T,TCount>::ValueProxy
  rle_array<T,TCount>::
  operator[](std::size_t i)
  {
    assert(i <= this->size());
    typename Data::iterator iData = m_data.begin();
    while (i >= iData->length())
    {
      i -= iData->length();
      ++iData;
    }
    return ValueProxy(std::make_pair(iData, i), *this);
  }

  //---------------------------------------------------------------------------


  template < typename T, typename TCount >
  void 
  rle_array<T,TCount>::
  set
  (
    std::pair<typename Data::iterator, TCount> & index
  , const T & value
  )
  {
    // Do nothing if nothing to do
    if (value == index.first->value()) return;

    // Is the current segment made of only one element?
    if (index.first->length() == 1)
    {
      assert(index.second == 0);

      typename Data::iterator iPreviousSegment = index.first;
      bool atBeginning = (index.first == m_data.begin());
      if (!atBeginning) --iPreviousSegment;
      typename Data::iterator iNextSegment = index.first;
      bool atEnd = (++iNextSegment == m_data.end());

      if (!atBeginning && iPreviousSegment->value() == value)
      {
        index.second = (iPreviousSegment->length())++;
        if (!atEnd && iNextSegment->value() == value)
        {
          iPreviousSegment->length() += iNextSegment->length();
          m_data.erase(iNextSegment);
        }
        m_data.erase(index.first);
        index.first = iPreviousSegment;
      }
      else if (!atEnd && iNextSegment->value() == value)
      {
        ++(iNextSegment->length());
        index.first = m_data.erase(index.first);        
      }
      // Then we remain a single point : simply update the value.
      else
      {
        // Simply overwrite the value
        index.first->value() = value;
        //index.first->setValue(value);
      }
    }
    // Are we at the very first of current repeated value?
    else if (index.second == 0)
    {
      typename Data::iterator iPreviousSegment = index.first;
      bool atBeginning = (index.first == m_data.begin());
      if (!atBeginning) --iPreviousSegment;
      
      // If so, decrease current length number by one
      --(index.first->length());

      // Is the previous value the same?
      if (!atBeginning && iPreviousSegment->value() == value)
      {
        // Then just increase its length number
        index.second = (iPreviousSegment->length())++;
        index.first = iPreviousSegment;
      }
      else
      {
        // else insert a new length value
        index.first = m_data.insert(index.first, Segment(1, value));
        index.second = 0;
      }
    }
    // Are we exactly at the end of current segment?
    else if (index.second == index.first->length()-1)
    { 
      typename Data::iterator iNextSegment = index.first;
      bool atEnd = (++iNextSegment == m_data.end());
            
      // Then, decrease current length number by one
      --(index.first->length());
      index.second = 0;
      
      if (!atEnd && iNextSegment->value() == value)
      {
        // If the next sequence has the same value, just increase
        // its length number by one
        ++(iNextSegment->length());
        index.first = iNextSegment;
      }
      else
      {
        // Insert a new length value
        index.first = m_data.insert(iNextSegment, Segment(1, value));
      }
    }
    // we are in the middle of a repeated value
    else
    {
      assert(index.second > 0);
      // First half
      m_data.insert(index.first, Segment(index.second, index.first->value()));
      // New element
      m_data.insert(index.first, Segment(1, value));
      // last half
      index.first->length() -= index.second+1;
      --index.first;
      index.second = 0;
    }
  }

/*
  template < typename T, typename TCount >
  class rle_array<T,TCount>::ValueProxy : public value_proxy<rle_array<T,TCount>, std::pair<typename Data::iterator, TCount> >
  {
  public: // typedefs
    typedef value_proxy<rle_array<T,TCount>, std::pair<typename Data::iterator, TCount> > Base;
    typedef typename Base::index_type index_type;
    typedef typename Base::container  container;
  public: // constructors & destructor    
      ValueProxy(index_type i, container * p) : Base(i, p) {}
  public: // operators
      // unfortunately we have to redefine this function by hand :(
      void operator=(typename boost::call_traits<T>::param_type value) { this->Base::operator=(value); }
  };
*/

  //---------------------------------------------------------------------------


  
  
} // namespace til


