
namespace til { namespace fuzzy
{

  //---------------------------------------------------------------------------

  template < typename T >
  inline
  boost::logic::tribool
  is_positive(T x, T delta)
  {
    assert(delta >= 0);
    if (x >= delta) return true;
    else if (x <= -delta) return false;
    //else return boost::logic::tribool(boost::logic::indeterminate);
    else return boost::logic::indeterminate;
  }

  //---------------------------------------------------------------------------

  template < typename T >
  inline 
  boost::logic::tribool
  same_sign(T x, T y, T delta)
  {
    // NB: I am not sure that this solution is actually much worse than what could be
    // achieved by explicit checking, regarding the number of tests to make.
    return fuzzy::is_positive(x, delta) == fuzzy::is_positive(y, delta);
  }

  //---------------------------------------------------------------------------

}} // namespace til::fuzzy

