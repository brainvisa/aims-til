#include <stdlib.h>

namespace til { namespace io
{

  //---------------------------------------------------------------------------

  template < typename TIterator, typename TProgressIndicator >
  void 
  SimpleLoader< TIterator, TProgressIndicator >::
  operator()(const char * filename, TIterator begin)
  {
    std::ifstream f(filename);
    if (!f) throw CannotOpenFile();
    std::string tmp;
    std::size_t i = 0;
    while (!std::getline(f, m_buffer).eof())
    {
      m_indicator(i);
      std::istringstream linestream(m_buffer);
      while (getline(linestream, tmp, ' '))
      {
        *begin = strtof(tmp.c_str(), 0);
        ++begin;
        ++i;
      }
    }    
  }

  //---------------------------------------------------------------------------

  template < typename TIterator, typename TProgressIndicator >
  void 
  SimpleLoader< TIterator, TProgressIndicator >::
  operator()(const std::string & filename, TIterator begin)
  {
    (*this)(filename.c_str(), begin);
  }

  //---------------------------------------------------------------------------

}} // namespace til::io


