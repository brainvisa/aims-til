#ifndef STD_WRAP_H_
#define STD_WRAP_H_

///\file Wrap of STL functions to take objects as input.
// I think it is important to pass objects rather than iterators to these functions.
// I mean, from what I understand, the reasons why STL is passing iterators around is for technical
// reasons (and I mean technical ease, here) -- one of them being compatibility with C-style
// arrays, others being I think const and reference type attribute headache.

// Okay, actually it's more suble: the iterator can contain part of the algorithm itself,
// e.g. reverse or insert or whatever. Such iterators actually extends the power of the
// algorithm to degrees the algorithm is not even aware of. That could be taken care of
// here by using extra template arguments reflecting the iterator types to use in the algo.

// The pb is that iterators generally know almost nothing about the objects they are pointing to.
// STL iterators are quite plain and thus manipulating iterators is okay. But when we try to do
// more complicated stuff, like compressed arrays, this hurts.
// I'm perfectly aware that precisely, compressed arrays and stuff should rather not be STL compliant,
// because if they were, it would be too easy to misuse them.
// Still: if we want to have a function that is generic and would accept both, we need these functions.

namespace til
{
  /// Wrapper of std::fill
  template < typename TContainer >
  inline void fill(TContainer & c, typename boost::call_traits<typename TContainer::value_type>::param_type value)
  {
    std::fill(c.begin(), c.end(), value);
  }
}

#endif /*STD_WRAP_H_*/
