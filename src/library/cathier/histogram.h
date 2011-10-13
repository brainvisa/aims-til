#ifndef _TIL_HISTOGRAM_H_
#define _TIL_HISTOGRAM_H_

// includes from STL
#include <map>

// includes from BOOST
#include <boost/call_traits.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

// includes from TIL
#include "til/Accumulator.h"

// includes from TIL
//#include "cat2type.h"
#include "globalTraits.h"

namespace til
{
  namespace policy
  {

    /// Histogram accumulation policy for Accumulator.
    /// This policy assumes that the histogram is computed by means of a map. Accumulation of a value
    /// then consists simply in increase the correct map entry by one
/*
    //template < typename T, typename TCount, template <typename,typename> class TMap = std::map >
    struct Accumulation_MapHistogram
    {
      typedef TMap<TCount,T> Map;
      
      void operator()(Map & m, const T & value)
      {
        // look if current value is already present in the histogram
        typename Map::iterator i = m.find(value);
        if (i == m.end())
        {
          // If not, add it in the histogram, and set count to one.
          m[value] = TCount(1);
        }
        else
        {
          // Otherwise, simply increase its count.
          ++(i->second);
        }
      }
    };
  }
*/
    // create Accumulation_MapHistogram_def to avoid Accumulation_MapHistogram<T,TCount,TMap> calls
    struct Accumulation_MapHistogram_def
    {
      typedef std::map<std::size_t, std::size_t > Map;
      
      void operator()(Map & m, const std::size_t & value)
      {
        // look if current value is already present in the histogram
        Map::iterator i = m.find(value);
        if (i == m.end())
        {
          // If not, add it in the histogram, and set count to one.
          m[value] = 1;
        }
        else
        {
          // Otherwise, simply increase its count.
          ++(i->second);
        }
      }
    };
    
 }

  /// Compute a histogram by means of a map.
  /// There's no binning (yet). I guess this is temporary until I have a need
  /// for a better, binning histogram.
  /// Note that binning policy should be independant from storage policy. I could have binning and still
  /// using maps, or not.
/*
  template < typename T, typename TCount = unsigned int, template <typename,typename> class TMap = std::map >
  class MapHistogram : public Accumulator<T, TMap<TCount, T>, policy::Accumulation_MapHistogram<T,TCount,TMap> >
  {
  public: // typedefs

    typedef TMap<T, TCount> Map;
  };
*/
  // create MapHistogram_def to avoid MapHistogram<std::size_t> calls
  // not allowed in new gcc version
  class MapHistogram_def : public Accumulator<std::size_t, std::map<std::size_t, std::size_t>, policy::Accumulation_MapHistogram_def >
  {
  public: // typedefs

    typedef std::map<std::size_t, std::size_t> Map;
  };
 

  
/*
  template < typename T, typename TCount, template <typename,typename> class TMap >
  void print(const MapHistogram_def & h)
  {
    typename MapHistogram<T,TCount,TMap>::Map::const_iterator iH = h.get().begin();
    for (; iH != h.get().end(); ++iH)
    {
      std::cout << iH->first << " : " << iH->second << std::endl;
    }
  }
*/
  // create printf_def to replace print(const MapHistogram_def & h) method
  void print_def(const MapHistogram_def & h);

  void save_def(const MapHistogram_def & h, std::string filename, float valueFactor);

  class MapHistogram_float
  {
      public:

      MapHistogram_float(float range);
      ~MapHistogram_float() {}
      void accumulate(float value);
      void print();
      void save(std::string filename, float valueFactor);
      protected:

      MapHistogram_def h;
      float _range;
  };

} // namespace

#endif /*_TIL_HISTOGRAM_H_*/
