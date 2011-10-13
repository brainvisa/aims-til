#ifndef TIL_TEXPR_MACROS_H
#define TIL_TEXPR_MACROS_H

/// \file
/// Some macros to ease the otherwise tedious and unreadable declaration
/// of template expression classes


//#define EXPRFUNC_1ARG(name, body) EXPRFUNC_1ARG_RET(name, body, typename TypeStruct<Iterator1>::Type)
//#define EXPRFUNC_2ARG(name, body) EXPRFUNC_2ARG_RET(name, body, typename TypeStruct<Iterator1 TIL_COMMA Iterator2>::Type)
//#define EXPRFUNC_3ARG(name, body) EXPRFUNC_3ARG_RET(name, body, typename TypeStruct<Iterator1 TIL_COMMA Iterator2 TIL_COMMA Iterator3>::Type)
#define EXPRFUNC_1ARG(name, body) EXPRFUNC_1ARG_ARG(name, body, i1)
#define EXPRFUNC_2ARG(name, body) EXPRFUNC_2ARG_ARG(name, body, i1, i2)
#define EXPRFUNC_3ARG(name, body) EXPRFUNC_3ARG_ARG(name, body, i1, i2, i3)

#define TIL_COMMA ,

// TODO: prefix TIL_ to those macro names.

#define EXPRFUNC_1ARG_RET(name, body, ret)  EXPRFUNC_1ARG_FULL(name, body, ret, i1)
#define EXPRFUNC_2ARG_RET(name, body, ret)  EXPRFUNC_2ARG_FULL(name, body, ret, i1, i2)
#define EXPRFUNC_3ARG_RET(name, body, ret)  EXPRFUNC_3ARG_FULL(name, body, ret, i1, i2, i3)

#define EXPRFUNC_1ARG_ARG(name, body, arg1)  EXPRFUNC_1ARG_FULL(name, body, typename TypeStruct<Iterator1>::Type, arg1)
#define EXPRFUNC_2ARG_ARG(name, body, arg1, arg2)  EXPRFUNC_2ARG_FULL(name, body, typename TypeStruct<Iterator1 TIL_COMMA Iterator2>::Type, arg1, arg2)
#define EXPRFUNC_3ARG_ARG(name, body, arg1, arg2, arg3)  EXPRFUNC_3ARG_FULL(name, body, typename TypeStruct<Iterator1 TIL_COMMA Iterator2 TIL_COMMA Iterator3>::Type, arg1, arg2, arg3)


#define EXPRFUNC_1ARG_FULL(name, body, ret, arg1) \
template < class Iterator1 >                      \
ret                                               \
name (Iterator1 & arg1)                           \
{                                                 \
body                                              \
}                                                 \

#define EXPRFUNC_2ARG_FULL(name, body, ret, arg1, arg2) \
template < class Iterator1, class Iterator2>            \
ret                                                     \
name (Iterator1 & arg1, Iterator2 & arg2)               \
{                                                       \
body                                                    \
}                                                       \

#define EXPRFUNC_3ARG_FULL(name, body, ret, arg1, arg2, arg3)     \
template < class Iterator1, class Iterator2, class Iterator3 >    \
ret                                                               \
name (Iterator1 & arg1, Iterator2 & arg2, Iterator3 & arg3)       \
{                                                                 \
body                                                              \
}                                                                 \



/// Used to define the return type of template expresions.
/// The return type of a template expression is accessible via
/// X::TypeStruct<Iterator1,Iterator2,Iterator3>::Type
#define EXPR_RESULT_TYPE(rettype)														\
template < class Iterator1, class Iterator2 = Iterator1, class Iterator3 = Iterator1>	\
struct TypeStruct																		\
{																						\
typedef rettype Type;																	\
}																						\

#define EXPR_ITERATORTYPE(rettype)														\
template < class Iterator1, class Iterator2 = Iterator1, class Iterator3 = Iterator1>	\
struct IteratorStruct																	\
{																						\
typedef rettype Type;																	\
}																						\


#endif

