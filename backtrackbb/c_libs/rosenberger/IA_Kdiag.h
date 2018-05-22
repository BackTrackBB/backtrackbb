/*
 * Kdiag.h
 * Routines to diagonalize a 3x3 upper triangular matrix
 *
 * Revision 1.4  2018/05/22 Claudio Satriano <satriano@ipgp.fr>
 *               - compiling under Win10
 * Revision 1.3  2016/03/15 Claudio Satriano <satriano@ipgp.fr>
 *               - subtrack2 now only returns the polarization filter
 * Revision 1.2  2014/05/11 Claudio Satriano <satriano@ipgp.fr>
 *               - Added rl_filter option to subtrack2
 * Revision 1.1  2011/01/31 22:56:05  andreas
 * Initial revision
 */
#ifndef KDIAG
#define KDIAG

#define IA_DEBUG

#include <stdio.h>
#include <stdlib.h>
#ifdef _MSC_VER
  #include <io.h>
#else
  #include <unistd.h>
#endif
#include <string.h>
#include <math.h>

#define KOG_EPS 1.0e-10  /* Termination-cond ist KOG_EPS*KOG_EPS */
#define KOG_MAXIT 10

/** data vector: {ZZ,NS,EW} */
typedef struct {
    double z;   /* vz-component */
    double n;   /* hn-component */
    double e;   /* he-component */
} IA_Dvect;

typedef struct {
    /* all 3 x 2 matrices */
    double **U;
    double **S;
} USV32_struct;

USV32_struct *New_USV32_struct(void);

void Init_USV32_struct(double **X, USV32_struct *P);

void Destroy_USV32_struct(USV32_struct *P);

void printmat(double **M, int m, int n);

void printvect(double *d, int m);

void subtrack2(const IA_Dvect *x, double *pol,
               const float lambda, const float delta,
               char proj, char rl_filter,
               USV32_struct *P,char r2u_init);

void mat_mult(double **A, double **B, int a_rows, int a_columns, int b_columns, double **C);

double **mat_transp(double **A,  int n_rows, int n_columns);

/* compute norm squared */
double norm2(double *a, int n_rows);

void KogExact(double a00, double a01, double a11, double **Q, double **A, double **E);
/*  -----------------------------------------------------------------------
 * Exact variant of Kogbetlianz's algorithm after
 *
 * Charlier, J. P. and Vanbegin, M. and Dooren, P. Van, "On Efficient
 * Implementations of Kogbetliantz's Algorithm for Computing the
 * Singular Value Decomposition", Numerische Mathematik,
 * vol. 52, pp. 279-300, 1988
 *
 * Computes left and right rotation matrices Q and E as well as
 * diagonal matrix A from upper triangular 2x2 matrix with elements
 * a_ij
 * -----------------------------------------------------------------------*/

void KogApx(double a00, double a01, double a11, double **Q, double **A, double **E);
/*  -----------------------------------------------------------------------
 * Approximation variant of Kogbetlianz's algorithm after
 *
 * Charlier, J. P. and Vanbegin, M. and Dooren, P. Van, "On Efficient
 * Implementations of Kogbetliantz's Algorithm for Computing the
 * Singular Value Decomposition", Numerische Mathematik,
 * vol. 52, pp. 279-300, 1988
 *
 * Computes left and right rotation matrices Q and E as well as
 * diagonal matrix A from upper triangular 2x2 matrix with elements
 * a_ij
 * -----------------------------------------------------------------------*/

void KDiag2by2(double **A, int i, int j, double **R, double **Y, double **T, char apx);
/* -----------------------------------------------------------------------
   Diagonalize 2x2 upper triangular sub-matrix A_ii,A_ij,A_jj of 3x3 matrix
   A; return A's SVD in 3x3 matrices R,Y,T.
   Use approximate Kb algorithm if apx == 1
   R,Y and T must be pre-allocated
   -----------------------------------------------------------------------*/

void BA3_SVD(double **A, double **U, double **S, double **V, int *n);
/*  -----------------------------------------------------------------------
 * Singular value Decomposition of a 3 x 3 upper triangular matrix A.
 * On output n is number of iterations
 * Matrix A is destroyed !!
 *  ----------------------------------------------------------------------- */

void *IA_Ealloc(size_t size, const char *fail_msg);
void **alloc2 (size_t n1, size_t n2, size_t size);
float **alloc2float(size_t n1, size_t n2);
void free2float(float **p);
double **alloc2double(size_t n1, size_t n2);
void free2double(double **p);

void free2 (void **p);
float **ealloc2float(size_t n1, size_t n2);
double **ealloc2double(size_t n1, size_t n2);
#endif
