#ifndef FUNC_ITERATOR_H_
#define FUNC_ITERATOR_H_

namespace til
{
  namespace detail
  {
    template < typename TIterator, typename TFunctor >
    class func_iterator_impl : public TIterator
    {
    public: // typedefs
      typedef TIterator Base;
      typedef typename TFunctor::result_type value_type;
      
    public: // constructors & destructor
      func_iterator_impl()                                   : TIterator( ), m_f( ) {}
      explicit func_iterator_impl(const TIterator & i)       : TIterator(i), m_f( ) {}
      explicit func_iterator_impl(TFunctor f)                : TIterator( ), m_f(f) {}
      func_iterator_impl(const TIterator & i, TFunctor f)    : TIterator(i), m_f(f) {}
      
    public: // functions
      value_type operator*()
      {
        return m_f(this->Base::operator*());
      }
  
    private: // data
      TFunctor m_f;
    };
  }
  
  /// A functor iterator is an iterator that overload the operator* to a call
  /// to a functor on the pointed value.
  template < typename TContainer, typename TFunctor >
  class func_iterator 
    : public detail::func_iterator_impl<typename TContainer::iterator, TFunctor>
  {
  public: // typedefs
    typedef detail::func_iterator_impl<typename TContainer::iterator, TFunctor> Base;
  
  public: // constructors & destructor
    func_iterator() : Base() {}
    explicit func_iterator(const typename TContainer::iterator & i) : Base(i) {}
    explicit func_iterator(TFunctor f) : Base(f) {}
    func_iterator(const typename TContainer::iterator & i, TFunctor f) : Base(i,f) {}
  };
  
  
  /// A functor iterator is an iterator that overload the operator* to a call
  /// to a functor on the pointed value.
  template < typename TContainer, typename TFunctor >
  class const_func_iterator
    : public detail::func_iterator_impl<typename TContainer::const_iterator, TFunctor>
  {
  public: // typedefs
    typedef detail::func_iterator_impl<typename TContainer::const_iterator, TFunctor> Base;
  
  public: // constructors & destructor
    const_func_iterator() : Base() {}
    explicit const_func_iterator(const typename TContainer::const_iterator & i) : Base(i) {}
    explicit const_func_iterator(TFunctor f) : Base(f) {}
    const_func_iterator(const typename TContainer::const_iterator & i, TFunctor f) : Base(i, f) {}
  };
  
  
  template < typename TFunctor, typename TIterator >
  detail::func_iterator_impl<TIterator, TFunctor>
  func_it(TIterator it, TFunctor f = TFunctor())
  {
    return detail::func_iterator_impl<TIterator, TFunctor>(it, f);
  }
  
  
  /// A wrapper on stl containers that returns func iterators.
  /// It does not bring any functionality, but it is needed, because it is somewhat
  /// dangerous to have iterators and containers that are not consistent (e.g.
  /// having different value_type).
  template < typename TContainer, typename TFunctor >
  class func_container : public TContainer
  {
  public: // typedefs
  
    typedef TContainer Base;
  
    typedef typename TFunctor::result_type              value_type;
    typedef func_iterator<TContainer, TFunctor>         iterator;
    typedef const_func_iterator<TContainer, TFunctor>   const_iterator;
  
  public: // constructors & destructor
    func_container(TContainer & c) : m_c(c), m_f() {}
  
  public: // functions
  
    iterator begin() { return iterator(this->Base::begin(), m_f); }
    iterator end()   { return iterator(this->Base::end(), m_f); }
    const_iterator begin() const { return const_iterator(this->Base::begin(), m_f); }
    const_iterator end()   const { return const_iterator(this->Base::end(), m_f); }
  
    value_type operator[](int i) { return m_f(this->Base::operator[](i)); }
  
  private: // data
  
    TContainer & m_c;
    TFunctor m_f;
  };
  
} // namespace til

#endif /*FUNC_ITERATOR_H_*/
