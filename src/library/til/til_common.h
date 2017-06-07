#ifndef TIL_TIL_COMMON_H
#define TIL_TIL_COMMON_H

/// \file
/// General macros, definitions and functions.

#include <cmath> // sqrt
#include <iostream>

#ifdef _MSC_VER
#pragma inline_depth(30)
#endif

// A definition used to declare friend template functions inside classes
// and still being to compile both with MSVC and GCC
#ifdef _MSC_VER
#define FRIEND_TEMPLATE_NO_ARG
#else
#define FRIEND_TEMPLATE_NO_ARG <>
#endif

// Define INLINE
#ifdef _MSC_VER
#define INLINE __forceinline
#else
#define INLINE inline
#endif


#if defined( _WIN32 ) && !defined( _GNUC )
// A definition to export/import objects and declarations, depending on whether
// these are used inside or outide the library
#ifdef aimstil_EXPORTS
#define TIL_API  __declspec(dllexport)
#define EXPIMP_TEMPLATE
#else
#define TIL_API  __declspec(dllimport)
#define EXPIMP_TEMPLATE extern
#endif
#else // WIN32
// On linux, we don't need this crap apparently
#define TIL_API
#define EXPIMP_TEMPLATE
#endif

// Define M_PI if it is not defined already
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


// A macro useful to write wrapper functions
// TODO: rename into TIL_EXPAND_VECTOR3D
#define EXPAND_VECTOR(v) v[0], v[1], v[2]



/// Convenient macro used in cout lines to separate two numbers
#define __ <<" "<<

// Signal windows.h not to define 'min' and 'max' macros. This should helfully tackle the 
// belowmentionned problem.
#ifndef NOMINMAX
#define NOMINMAX 1
#endif
/*
// Undefine 'min' and 'max' macros, 
// defined in one of the windows include files
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif
*/

namespace til
{
  // Convenient typedefs on unsigned types
  typedef unsigned char uchar;
  typedef unsigned short ushort;

  /// Library warning.
  inline void warning(const char * msg) { std::cout << "TIL warning: " << msg << std::endl; }

  /// A type representing the absence of a type.
  struct None {};

  /// A type used in some construtors to disable initialization.
  struct NoInit {};
  const NoInit no_init = NoInit();

} // namespace til


// Package include
// TODO: not sure whether this should stay that way...
#include "til/misc_scalar_functions.h"

#endif

