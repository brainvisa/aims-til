#ifndef TIL_META_H
#define TIL_META_H

///\file Basic meta-programming tools

namespace til { namespace meta /// < Namespace for meta-programming tools
{

  /// Type representing an integer
  template < int i > struct int_type
  {
    enum { value = i };
  };

  template < typename T1, typename T2 > struct add_type;
  
  /// Addition of two integer types
  template < int i, int j>
  struct add_type<int_type<i>, int_type<j> > 
  {
    typedef int_type<i+j> type;
  };

  

}}// namespace til::meta

#endif

