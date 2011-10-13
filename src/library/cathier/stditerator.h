#ifndef _STDITERATOR_H_
#define _STDITERATOR_H_

#include "til/templateTools.h"


template < typename T, typename T2 = T >
struct stditerator 
  : public til::type_if<boost::is_const<T>::value, typename T2::const_iterator, typename T2::iterator>
{};
  

#endif //_STDITERATOR_H_
