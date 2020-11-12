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

void copy_mat_offset(csc *A, csc *B, c_int num_col, c_int column_offset, c_int * count);
void copy_csc_mat_in_place(int An, const c_float *Ax, const c_int *Ai, const c_int*Ap,
                           c_float *Bx, c_int *Bi, c_int*Bp);
csc* copy_csc_mat_cols(const csc *A, csc *B,
                       c_int row_start, c_int col_start,
                       c_int row_stop, c_int col_stop);

void copy_csc_mat_transpose(const csc *A, c_float * data,
                                 c_int row_start, c_int col_start,
                                 c_int row_stop, c_int col_stop);

void copy_csc_mat_cols_ident(const csc *A, csc *B,
                             c_int row_start, c_int col_start,
                             c_int row_stop,   c_int col_stop, c_float offset);

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

int amub_col_plus_rho_upper_diag (int nrow, int ncol, int values,
				                  c_float *Bx, c_int *Bi, c_int*Bp,
				                  c_float *Ax, c_int *Ai, c_int*Ap,
				                  c_float *Cx, c_int *Ci, c_int*Cp,
				                  int nzmax,
				                  c_float * rhoinv);

/**
 * C = A * B
 */
int amub_col (int nrow, int ncol, int values,
	          c_float *Bx, c_int *Bi, c_int*Bp,
	          c_float *Ax, c_int *Ai, c_int*Ap,
	          c_float *Cx, c_int *Ci, c_int*Cp,
	          int nzmax);

/**
 * C = A * (B + I)
 */
int amub_col_plus_I (c_int nrow, c_int ncol, c_int start_col, c_int values,
	                 c_float *Bx, c_int *Bi, c_int*Bp,
	                 c_float *Ax, c_int *Ai, c_int*Ap,
	                 c_float *Cx, c_int *Ci, c_int*Cp,
		             c_int nzmax, double scale);
/**
 * C = A * B (row major)
 */
int amub_row (int nrow, int ncol, int values,
	          c_float *Ax, c_int *Ai, c_int*Ap,
	          c_float *Bx, c_int *Bi, c_int*Bp,
	          c_float *Cx, c_int *Ci, c_int*Cp,
	          int nzmax);

void subtract_csc_matrix_cols(int start_col, int stop_col, int n, 
                              c_int * Ap, c_int * Ai, c_float * Ax,
                              c_int * Bp, c_int * Bi, c_float * Bx );


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef CS_ADDON_H
