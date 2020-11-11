#ifndef CS_ADDON_H
# define CS_ADDON_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

# include "types.h"   // CSC matrix type
# include "lin_alg.h" // Vector copy operations

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
//void prea_copy_csc_mat(const csc *A,
//                       csc       *B);
void array_to_csc(csc* C, c_float* matrix);

void copy_mat_offset(csc *A, csc *B, c_int num_col, c_int column_offset, c_int * count);
void copy_csc_mat_in_place(int An, const c_float *Ax, const c_int *Ai, const c_int*Ap,
                           c_float *Bx, c_int *Bi, c_int*Bp);
csc* copy_csc_mat_cols(const csc *A, csc *B,
                       c_int row_start, c_int col_start,
                       c_int row_stop, c_int col_stop);

c_float *  copy_csc_submat(const csc *A,
                       c_int row_start, c_int col_start,
                       c_int row_stop, c_int col_stop);

void copy_csc_mat_transpose(const csc *A, c_float * data,
                                 c_int row_start, c_int col_start,
                                 c_int row_stop, c_int col_stop);

void copy_csc_mat_cols_ident(const csc *A, csc *B,
                             c_int row_start, c_int col_start,
                             c_int row_stop,   c_int col_stop, c_float offset);

/*****************************************************************************
* Printing                                                                   *
*****************************************************************************/

void print_csc_matrix(csc *M, const char *name);
void print_csc_matrix_row_major(csc *M, const char *name);

/**
 *  Print only selected columns of CSC matrix/rows of CSR matrix
 */
void print_csc_matrix_cols(csc *M, const char *name, int start_col, int stop_col);
/**
 *  Print CSC matrix as a dense matrix
 */
void print_csc_matrix_dense(csc *M, const char *name, int cols, int rows);

/*****************************************************************************
* Matrices Conversion                                                       *
*****************************************************************************/


/**
 * C = compressed-column CSC from matrix T in triplet form
 *
 * TtoC stores the vector of indices from T to C
 *  -> C[TtoC[i]] = T[i]
 *
 * @param  T    matrix in triplet format
 * @param  TtoC vector of indices from triplet to CSC format
 * @return      matrix in CSC format
 */
//csc* triplet_to_csc(const csc *T,
//                    c_int     *TtoC);

void csr_tocsc(const c_int  n_row, const c_int  n_col, 
           const c_int * Ap, const c_int * Aj,  const c_float * Ax,
           c_int *  Bp,  c_int * Bi, c_float * Bx);

void csc_tocsr(const c_int n_row, const c_int n_col, 
               const c_int * Ap, const c_int * Ai, const c_float * Ax,
                     c_int * Bp, c_int * Bj, c_float * Bx);

/*****************************************************************************
* Extra operations                                                          *
*****************************************************************************/

/**
 * A * B + rho*I
 */

int amub_col_plus_rho_upper_diag (int nrow, int ncol, int values,
				                  c_float *Bx, c_int *Bi, c_int*Bp,
				                  c_float *Ax, c_int *Ai, c_int*Ap,
				                  c_float *Cx, c_int *Ci, c_int*Cp,
				                  int nzmax,
				                  c_float * rhoinv);

/**
 * C = A * B + rho*I
 */
int amub_col_plus_rho (int nrow, int ncol, int values,
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
 * C = A * B + I
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

void csc_to_csr_with_diag(const c_int n_col,
                          const c_int n_row, 
                          const c_int * Ap, 
                          const c_int * Aj, 
                          const c_float * Ax,
                          c_int * Bp,
                          c_int * Bi,
                          c_float * Bx);

void subtract_csc_matrix_cols(int start_col, int stop_col, int n, 
                              c_int * Ap, c_int * Ai, c_float * Ax,
                              c_int * Bp, c_int * Bi, c_float * Bx );

/**
 * A + diag
 */
void apldia (int nrow, int job,
             c_float * Ax, c_int * Ai, c_int * Ap, c_int * diag,
             c_float * Bx, c_int * Bi, c_int *Bp);

void diapos(c_int n, c_int *ja, c_int * ia, c_int *idiag);

int aplb (int nrow,int ncol,int values,
          c_float *Ax, c_int *Ai, c_int*Ap,
          c_float *Bx, c_int *Bi, c_int*Bp,
          c_float *Cx, c_int *Ci, c_int*Cp,
          int nzmax,int * iw);


csc* transpose(const csc *A);

csc* cholesky_solve_t(const csc *L, csc *A, double* x);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef CS_ADDON_H
