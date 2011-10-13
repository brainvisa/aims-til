#ifndef TIL_IO_H_
#define TIL_IO_H_

// includes from STL
#include <exception>
#include <fstream>
#include <sstream>  // std::istringstream
#include <utility>  // std::pair


// includes from TIL
#include "progress_indicator.h"

namespace til { namespace io
{
  
  //---------------------------------------------------------------------------
  
    //---------------------//
   //  GetMatrixFileSize  //
  //---------------------//
  
  /// Return the number of row and columns of an array in a text file.
  class GetArraySize
  {
  public: // exceptions
  
    class CannotOpenFile : public std::exception {};
  
  public: // constructors
  
    explicit GetArraySize(std::size_t maxlinelength);

  public: // set & get
  
    void setMaxLineLength(std::size_t mll) { m_buffer.resize(mll); }
    /// Returns the number of lines.
    std::size_t nlines() const { return m_linecount; }
    /// Returns the number of columns.
    std::size_t ncolumns() const { return m_columncount; }    
    
  public: // operators
  
    /// Scan file to extract size information.
    /// Throws CannotOpenFile when, well, it cannot open the file.
    void operator()(const char * filename);
  
    /// Scan file to extract size information.
    void operator()(const std::string & filename) { (*this)(filename.c_str()); }
  
  private: // data, output
    std::size_t m_linecount, m_columncount;
    
  private: // data, internal
    std::string m_buffer;
  };
  
  //---------------------------------------------------------------------------
  
  /// Returns the number of row and columns in an array file.
  std::pair<std::size_t, std::size_t>
  get_array_size(const char * filename, std::size_t maxLineLength = 65536);
  
  //---------------------------------------------------------------------------

    //----------------//
   //  SimpleLoader  //
  //----------------//

  template < typename TIterator, typename TProgressIndicator = NoIndicator >
  class SimpleLoader
  {
  public: // exceptions

    class CannotOpenFile : public std::exception {};

  public: // constructors

    /// Construct object, giving maximum length, in characters, of a line on file.
    // TODO: somehow we should check when this max number is reached...
    explicit SimpleLoader(std::size_t maxlinelength)
    {
      m_buffer.resize(maxlinelength);
    }

  public: // set & get
  
    TProgressIndicator & progress_indicator() { return m_indicator; }
    
  public: // operators
  
    /// Load array from file, and store it from begin onwards.
    /// Throws CannotOpenFile if file cannot be opened.
    void operator()(const char * filename, TIterator begin);

    /// Load array from file, and store it from begin onwards.
    /// Throws CannotOpenFile if file cannot be opened.
    void operator()(const std::string & filename, TIterator begin);

  private: // data, input

    TProgressIndicator m_indicator;

  private: // data, internal

    std::string m_buffer;
  };  


  template < typename TIterator >
  void simple_load(const char * filename, TIterator begin, std::size_t maxlinelength = 65536)
  {
    SimpleLoader<TIterator> loader(maxlinelength);
    loader(filename, begin);
  }


  //----------------------------------------------------------------------------------------


    //---------------//
   //  ArrayWriter  //
  //---------------//

  /*
  class ArrayWriter
  {
  public: // operators
    
    void operator()(const char * filename)
  };

  {
    // convert from list to vector
    std::vector<std::size_t> res(linecount);
    std::copy(sets.front().begin(), sets.front().end(), res.begin());
    
    std::ofstream f(matoutname.c_str());
    if (!f) std::cerr << "Cannot open file " << matoutname << std::endl;
    else
    {
      til::DelayedProgressIndicator indicator(10, linecount);
      for (std::size_t i = 0; i < linecount; ++i)
      {
        indicator(i);
        for (std::size_t j = 0; j < linecount; ++j)
        {
          f << mat[res[i] * linecount + res[j]] << " ";
        }
        f << std::endl;
      }
    }
  }
  */

  
}} // namespace til::io

// package include
#include "io.tpp"

#endif /*IO_H_*/
