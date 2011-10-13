
namespace til
{

  //---------------------------------------------------------------------------
    
  template < typename TPICallback >
  ProgressIndicator<TPICallback>::ProgressIndicator
  (
    std::size_t parts
  , std::size_t total
  )
    : m_total(total)
    , m_parts(parts)
  {
    this->init();
  }
  
  //---------------------------------------------------------------------------

  template < typename TPICallback >
  void ProgressIndicator<TPICallback>::init() { m_currentpart = 0; m_mark = 0.0; }
    
  //---------------------------------------------------------------------------

  template < typename TPICallback >
  void ProgressIndicator<TPICallback>::operator()(std::size_t i)
  {
    if ( i >= m_mark && i < m_mark + 1 )
    {
      //std::cout << i << " " << m_parts << " " << m_total << " " << m_currentpart << " " << m_mark << std::endl;
      m_callback(m_currentpart, m_parts);
      this->nextPart();
    }
  }

  //---------------------------------------------------------------------------

  template < typename TPICallback >
  void ProgressIndicator<TPICallback>::nextPart()
  {
    ++m_currentpart;
    m_mark = m_currentpart * m_total / double(m_parts);
  }

  //---------------------------------------------------------------------------  
  
} // namespace til


