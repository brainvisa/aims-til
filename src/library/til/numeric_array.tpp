

namespace til
{

  namespace detail
  {
    //-------------------------------------------------------------------------

    template < typename T, std::size_t D >
    numeric_array_impl<T,D>::numeric_array_impl()
      : Base()
    {
      for (std::size_t i = 0; i < D; ++i) (*this)[i] = 0;
    }
    
    //-------------------------------------------------------------------------

  /*
    template < typename T, std::size_t D >
    numeric_array_impl<T,D>::numeric_array_impl
#ifdef NDEBUG
    (std::size_t)
#else
    (std::size_t d)
#endif
      : Base()
    {
      assert(d == D);
      for (std::size_t i = 0; i < D; ++i) (*this)[i] = 0;
    }

    //-------------------------------------------------------------------------
    
    template < typename T, std::size_t D >
    numeric_array_impl<T,D>::numeric_array_impl(std::size_t d, typename boost::call_traits<T>::param_type t)
      : Base()
    {
      assert(d == D);
      til::fill(*this, t);
    }
*/


  } // namespace detail

  
} // namespace til


