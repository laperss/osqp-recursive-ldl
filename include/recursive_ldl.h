#ifndef RECURSIVE_LDL_H
# define RECURSIVE_LDL_H


# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


# include "types.h"
# include "constants.h"


/**
 * Recursive LDL data structure
 */
typedef struct {
  OSQPData data;
  c_int  Nmax;   ///< Maximum horizon
  c_int     N;   ///< Current base horizon
  c_int    N0;   ///< Current base horizon
  c_int    N_recursive_temp;   ///< Current base horizon
  c_int    Nx;   ///< Current base horizon (of matrix X)
  c_int Nplus;   ///< Current additional horizon
  c_int    nx;   ///< Number of states
  c_int    ny;   ///< Number of constraints
  c_int    nu;   ///< Number of inputs
  c_int    nt;   ///< Number of terminal constraints
  c_int nnz[2];  ///< Constant and linear factor for computing number of nonzeros in L
  csc     *A0;   ///< A0 in R(ny+nx * nu)
  csc     *Ai;   ///< Ai in R(ny+nx * nx+nu)
  csc    *Aij;   ///< A0 in R(ny+nx * nx+nu)
  csc     *AN;   ///< A0 in R(nt * nx)
  csc     *Q0;   ///< Q0 in R(nu * nu)
  csc     *Qi;   ///< Qi in R(nx+nu * nx+nu)
  csc     *QN;   ///< Qi in R(nx * nx)
  csc **Amats;   ///< Contains pointers to all cost matrices
  csc **Qmats;   ///< Contains pointers to all cost matrices
  csc **X_even;   ///< Contains pointers to earlier pivot points
  double *Ytemp;   ///< The bottom rows U'*P*Linv*Dinv 
  csc   *Yhat;   ///< The bottom rows U'*P*Linv*Dinv 
  csc   *Uhat;   ///< The bottom rows U'*P*Linv*Dinv 
  csc   *Vhat;   ///< The bottom rows V*Q*Minv'*Einv 
  csc      *Lz;   ///< The lower triangular matriz Q'ZQ=MEM'
  c_float  *Dzinv;///< The (inverted) diagonal matrix E
  c_int    *Pz;   ///< The permutation matrix Q
  csc      *Ly;   ///< The lower triangular matriz Yhat=NFN'
  c_float  *Dyinv;   ///< The (inverted) diagonal matrix F
  c_int    *P0;   ///< The permutation matrix P0
} OSQPDataRLDL;

void cleanup_rldl(OSQPDataRLDL * data);

/**
 * Solve quadratic program
 *
 * The final solver information is stored in the \a work->info  structure
 *
 * The solution is stored in the  \a work->solution  structure
 *
 * If the problem is primal infeasible, the certificate is stored
 * in \a work->delta_y
 *
 * If the problem is dual infeasible, the certificate is stored in \a
 * work->delta_x
 *
 * @param  work Workspace allocated
 * @return      Exitflag for errors
 */

enum new_linsys_solver_type {  QDLDL_HORIZON_SOLVER=10 };
c_int osqp_update_recursive(OSQPWorkspace* work, OSQPDataRLDL *data_rldl, c_int N);
c_int osqp_update_X_horizon(OSQPWorkspace* work, OSQPDataRLDL *data_rldl, c_int deltaN);
c_int osqp_update_Z_horizon(OSQPWorkspace* work, OSQPDataRLDL *data_rldl, c_int deltaN);
c_int osqp_update_combine_recursive(OSQPWorkspace* work, OSQPDataRLDL *data_rldl, c_int N);
c_int osqp_solve_recursive(OSQPWorkspace *work);
c_int osqp_solve_combine_recursive(OSQPWorkspace *work);
c_int osqp_setup_recursive(OSQPWorkspace** workp, OSQPDataRLDL* data, const OSQPSettings* settings, c_int Nmax, c_int N, c_int nx, c_int nu, c_int ny, c_int nt);
c_int osqp_setup_combine_recursive(OSQPWorkspace** workp, OSQPDataRLDL* data, const OSQPSettings* settings, c_int Nmax, c_int N, c_int Nplus, c_int nx, c_int nu, c_int ny, c_int nt);
c_int osqp_allocate_combine_recursive(OSQPWorkspace** workp, OSQPDataRLDL* data,
				      const OSQPSettings* settings, c_int Nmax,
				      c_int nx, c_int nu, c_int ny, c_int nt);


c_int osqp_partial_update_bounds(OSQPWorkspace *work, c_int start, c_int stop,
			                             const c_float *l_new,
			                             const c_float *u_new);

//c_int init_linsys_solver_qdldl_recursive(qdldl_solver ** sp, c_int n, c_int m, c_int P_nnz_max, c_int A_nnz_max, const csc * L, c_float sigma, const c_float * rho_vec, c_int polish);

void get_L_dimensions(OSQPWorkspace* work, c_int *n, c_int *m, c_int *nnz);
void get_permutation_matrix(OSQPWorkspace* work);
# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef RECURSIVE_LDL_H 
