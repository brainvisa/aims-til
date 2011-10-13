// includes from STL
#include <map>
#include <iostream>
#include <fstream>

// includes from BOOST
 #include <boost/call_traits.hpp>
 #include <boost/type_traits.hpp>
 #include <boost/utility/enable_if.hpp>

// includes from TIL
#include "globalTraits.h"

// includes from CATHIER
#include "cathier/histogram.h"

namespace til
{
  // create printf_def to replace print(const MapHistogram_def & h) method
  void print_def(const MapHistogram_def & h)
  {
    MapHistogram_def::Map::const_iterator iH = h.get().begin();
    for (; iH != h.get().end(); ++iH)
    {
      std::cout << iH->first << " : " << iH->second << std::endl;
    }
  }

  void save_def(const MapHistogram_def & h, std::string filename, float valueFactor)
  {
    {
      std::fstream f;
      f.open(filename.c_str(), std::fstream::out);
      MapHistogram_def::Map::const_iterator iH = h.get().begin();
      for (; iH != h.get().end(); ++iH)
      {
        float value = iH->first*valueFactor;
        f << value << " : " << iH->second << std::endl;
      }
    }
  }


  MapHistogram_float::MapHistogram_float(float range)  : _range(range)
  {

  }


  void MapHistogram_float::accumulate(float value)
  {
    int hist_value = int(value / _range);
    h.accumulate(hist_value);
  }

  void MapHistogram_float::print()
  {
    print_def(h);
  }

  void MapHistogram_float::save(std::string filename,float valueFactor)
    {
      save_def(h, filename, valueFactor);
    }

} // namespace til


