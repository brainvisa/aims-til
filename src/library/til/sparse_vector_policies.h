#ifndef TIL_SPARSE_VECTOR_POLICIES_H
#define TIL_SPARSE_VECTOR_POLICIES_H

/// \file Baseline policies for sparse vectors.
/// These policies answers the following question: when no value is stored at the
/// current location, what value should be returned?
/// Although not really developped yet, baseline policies on sparse vectors can
/// do wild stuff. For exemple, one could decide to return the value of another
/// container: in this case, the sparse vector could help storing only differences
/// between containers.

#include "til/miscTools.h"

namespace til { namespace policy
{

  /// Take the default value of type T as the baseline value of a sparse vector.
  template < typename T >
  struct SVBaseline_DefaultValue
  {
    struct Exception : public std::runtime_error
    { Exception(const std::string & msg) : std::runtime_error(msg) {} };

    /// Get default value at current position.
    T operator()(std::size_t) const { return T(); }

    void setValue(const T & value)
    {
      if (value != T())
      {
        std::cerr << "Can't change baseline value to " << value << std::endl;
        throw Exception("Can't change baseline value");
      }
    }
  };


  template < typename T >
  struct SVBaseline_Value
  {
    struct Exception : public std::runtime_error
    { Exception(const std::string & msg) : std::runtime_error(msg) {} };

    
    SVBaseline_Value() : m_value() {}
    explicit SVBaseline_Value(const T & value) : m_value(value) {}

    /// Get default value at current position.
    T operator()(std::size_t) const { return m_value; }

    void setValue(const T & value)
    {
      m_value = value;
    }
    
    T m_value;
  };


}} // namespace til::policy

#endif

