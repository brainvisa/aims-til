#ifndef LAPACK_H_
#define LAPACK_H_


/// \file This is a hand-made header for the C interface of LAPACK functions, since I could not find it.
/// NB: This is C; use extern "C" to include in a C++ program.

void dposv_(char *, int *, int *, double *, int *, double *, int *, int *);
void dgesv_(int *, int *, double *, int *, int *, double *, int *, int *);
void dsyev_(char *, char *, int *, double *, int *, double *, double *, int *, int *);


#endif /*LAPACK_H_*/
