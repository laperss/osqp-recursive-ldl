#ifndef CS_ADDON_H
# define CS_ADDON_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

# include "types.h"   // CSC matrix type

/*****************************************************************************
* Create and free CSC Matrices                                              *
*****************************************************************************/

/**
 * Create Compressed-Column-Sparse matrix from existing arrays
    (no MALLOC to create inner arrays x, i, p)
 * @param  m     First dimension
 * @param  n     Second dimension
 * @param  nzmax Maximum number of nonzero elements
 * @param  x     Vector of data
 * @param  i     Vector of row indices
 * @param  p     Vector of column pointers
 * @return       New matrix pointer
 */
//csc* csc_matrix(c_int    m,
//                c_int    n,
//                c_int    nzmax,
//                c_float *x,
//                c_int   *i,
//                c_int   *p);


/*****************************************************************************
* Copy Matrices                                                             *
*****************************************************************************/

/**
 *  Copy sparse CSC matrix A to B (B is preallocated, NO MALOC)
 */
    
void copy_csc_matrix(const csc *A, csc *B);

void copy_csc_transpose(const csc *A, c_float * data);

void copy_csc_plus_sigma(const csc *A, csc *B, c_float offset);

/*****************************************************************************
* Printing (used only in debugging)                                          *
*****************************************************************************/

void print_csc_matrix(csc *M, const char *name);


/*****************************************************************************
* Matrices Conversion                                                       *
*****************************************************************************/

/**
 * C = compressed-column CSC from matrix T in triplet form
 *
 * csr_to_csc: convert a matrix in csc format to csr
 *
 */

void csr_to_csc(const c_int  n_row, const c_int  n_col, 
                const c_int * Ap, const c_int * Aj,  const c_float * Ax,
                c_int *  Bp,  c_int * Bi, c_float * Bx);

void csc_to_csr(const c_int n_row, const c_int n_col, 
                const c_int * Ap, const c_int * Ai, const c_float * Ax,
                c_int * Bp, c_int * Bj, c_float * Bx);

/*****************************************************************************
* Extra operations                                                          *
*****************************************************************************/

/**
 * A * B + rho*I
 * where only the upper diagonal is included.
 */

int A_times_B_plus_rho(c_int ncol, c_int values,
	                   c_float *Bx, c_int *Bi, c_int*Bp,
	                   c_float *Ax, c_int *Ai, c_int*Ap,
	                   c_float *Cx, c_int *Ci, c_int*Cp,
	                   c_int nzmax,
	                   c_float * rhoinv);

/**
 * C = A * B
 */
int A_times_B(c_int nrow, c_int ncol, c_int values,
	          c_float *Bx, c_int *Bi, c_int*Bp,
	          c_float *Ax, c_int *Ai, c_int*Ap,
	          c_float *Cx, c_int *Ci, c_int*Cp,
	          c_int nzmax);

/**
 * C = A * (B + I)
 */
int A_times_B_plus_I (c_int nrow, c_int ncol, c_int start_col, c_int values,
	                 c_float *Bx, c_int *Bi, c_int*Bp,
	                 c_float *Ax, c_int *Ai, c_int*Ap,
	                 c_float *Cx, c_int *Ci, c_int*Cp,
		             c_int nzmax, c_float scale);


void A_minus_B(c_int start_col, c_int stop_col, c_int n, 
               c_int * Ap, c_int * Ai, c_float * Ax,
               c_int * Bp, c_int * Bi, c_float * Bx );


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef CS_ADDON_H
