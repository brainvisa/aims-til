
// includes from TIL
#include "io.h"

namespace til { namespace io
{

  //---------------------------------------------------------------------------

  GetArraySize::GetArraySize(std::size_t maxlinelength)
    : m_buffer(maxlinelength, '\0')
  {
  }
  
  //---------------------------------------------------------------------------

  void GetArraySize::operator()(const char * filename)
  {
    // open file
    std::ifstream f(filename);
    if (!f) throw CannotOpenFile();

    m_linecount = 0;
    m_columncount = 0;

    bool done = false;
    std::string tmp;
    while (!std::getline(f, m_buffer).eof())
    {
      if (!done)
      {
        std::istringstream linestream(m_buffer);
        while (std::getline(linestream, tmp, ' '))
        {
          ++m_columncount;
        }  
        done = true;
      }
      ++m_linecount;
    }
  }

  //---------------------------------------------------------------------------

  /// Returns the number of row and columns in an array file.
  std::pair<std::size_t, std::size_t>
  get_array_size(const char * filename, std::size_t maxLineLength)
  {
    GetArraySize f(maxLineLength);
    f(filename);
    return std::make_pair(f.ncolumns(), f.nlines());
  }

  //---------------------------------------------------------------------------

}} // namespace til::io


