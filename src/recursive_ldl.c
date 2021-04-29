#include "recursive_ldl.h"
#include "osqp.h"
#include "auxil.h"
#include "util.h"
#include "scaling.h"
#include "error.h"
#include "qdldl_interface.h"
#include "qdldl.h"

#include "cs_addon.h"
#include "kkt.h"



void cleanup_rldl(OSQPDataRLDL * data){
    if (data){
	for (int i=0;i<data->Nmax;i++){
	    if (data->X_even[i] != OSQP_NULL) {
		csc_spfree(data->X_even[i]); 
	    }
	}
	free(data->X_even);
	free(data->Amats);
	free(data->Qmats);

	if (data->Ytemp) {
	    c_free(data->Ytemp);
	}
	if (data->Uhat) {
	    csc_spfree(data->Uhat);
	}
	if (data->Vhat) {
	    csc_spfree(data->Vhat);
	}
	if (data->Yhat) {
	    csc_spfree(data->Yhat);
	}
	if (data->Lz) {
	    csc_spfree(data->Lz);
	}
	if (data->Pz) {
	    free(data->Pz);
	}
	if (data->Dzinv) {
	    free(data->Dzinv);
	}
	if (data->Ly) {
	    csc_spfree(data->Ly);
	}
	if (data->Dyinv) {
	    free(data->Dyinv);
	}
	if (data->P0) {
	    free(data->P0);
	}
    }
    c_free(data);    

}

// Solves LX = A
void QDLDL_Lsolve_mat(const QDLDL_int    n,
                      const QDLDL_int*   Lp,
                      const QDLDL_int*   Li,
                      const QDLDL_float* Lx,
                      const QDLDL_int    m,
                      QDLDL_float*       A){

    QDLDL_int i,j, k;
    for(i = 0; i < n; i++){
        for(j = Lp[i]; j < Lp[i+1]; j++){
            for (k=0; k < m; k++){

                A[Li[j] + k*n] -= Lx[j]*A[i + k*n];
            }
        }
    }
}


void QDLDL_Ltsolve_mat(const QDLDL_int        n,
                       const QDLDL_int*      Lp,
                       const QDLDL_int*      Li,
                       const QDLDL_float*    Lx,
                       const QDLDL_int        m,
                       QDLDL_float *          A){

    QDLDL_int i,j,k;
    for(i = n-1; i>=0; i--){
        for(j = Lp[i]; j < Lp[i+1]; j++){
            for (k=0; k < m; k++){
                A[i + k*n] -= Lx[j]*A[Li[j] + k*n];
            }
        }
    }
}

void QDLDL_solve_mat(const QDLDL_int        n,
                     const QDLDL_int*      Lp,
                     const QDLDL_int*      Li,
                     const QDLDL_float*    Lx,
                     const QDLDL_float*  Dinv,
                     const QDLDL_int        m,
                     QDLDL_float *          A){

    QDLDL_int i, j, k;
    QDLDL_Lsolve_mat(n, Lp, Li, Lx, m, A);
    
    for(i = 0; i < n; i++) {
        for (k=0; k < m; k++){
            A[i + k*n] *= Dinv[i];
        }
    }
    QDLDL_Ltsolve_mat(n,Lp,Li,Lx,m,A);

}


c_int osqp_partial_update_bounds(OSQPWorkspace *work, c_int start, c_int stop,
				 const c_float *l_new,
				 const c_float *u_new) {
    c_int i, exitflag = 0;

    // Check if workspace has been initialized
    if (!work) return osqp_error(OSQP_WORKSPACE_NOT_INIT_ERROR);
    if (start>=stop) return osqp_error(OSQP_DATA_VALIDATION_ERROR);

#ifdef PROFILING
    if (work->clear_update_time == 1) {
	work->clear_update_time = 0;
	work->info->update_time = 0.0;
    }
    osqp_tic(work->timer); // Start timer
#endif /* ifdef PROFILING */

    // Check if lower bound is smaller than upper bound
    for (i = 0; i < stop-start; i++) {
	if (l_new[i] > u_new[i]) {
#ifdef PRINTING
	    printf("lower bound must be lower than or equal to upper bound\n");
#endif /* ifdef PRINTING */
	    return 1;
	}
    }

    // Replace l and u by the new vectors
    prea_vec_copy(l_new, &work->data->l[start], stop-start);
    prea_vec_copy(u_new, &work->data->u[start], stop-start);
    
    // Scaling
    if (work->settings->scaling) {
	vec_ew_prod(work->scaling->E, work->data->l, work->data->l, work->data->m);
	vec_ew_prod(work->scaling->E, work->data->u, work->data->u, work->data->m);
    }
    // Reset solver information
    reset_info(work->info);

#if EMBEDDED != 1
    // Update rho_vec and refactor if constraints type changes
    int constr_type_changed = 0;

    for (i = start; i < stop; i++) {
	if ((work->data->l[i] < -OSQP_INFTY * MIN_SCALING) &&
	    (work->data->u[i] > OSQP_INFTY * MIN_SCALING)) {
	    // Loose bounds
	    if (work->constr_type[i] != -1) {
		work->constr_type[i] = -1;
		work->rho_vec[i]     = RHO_MIN;
		work->rho_inv_vec[i] = 1. / RHO_MIN;
		constr_type_changed  = 1;
	    }
	} else if (work->data->u[i] - work->data->l[i] < RHO_TOL) {
	    // Equality constraints
	    if (work->constr_type[i] != 1) {
		work->constr_type[i] = 1;
		work->rho_vec[i]     = RHO_EQ_OVER_RHO_INEQ * work->settings->rho;
		work->rho_inv_vec[i] = 1. / work->rho_vec[i];
		constr_type_changed  = 1;
	    }
	} else {
	    // Inequality constraints
	    if (work->constr_type[i] != 0) {
		work->constr_type[i] = 0;
		work->rho_vec[i]     = work->settings->rho;
		work->rho_inv_vec[i] = 1. / work->settings->rho;
		constr_type_changed  = 1;
	    }
	}
    }

#endif // EMBEDDED != 1

#ifdef PROFILING

    work->info->update_time += osqp_toc(work->timer);
#endif /* ifdef PROFILING */
    


    return exitflag;
}


c_int update_problem_size(qdldl_solver * s, c_int n, c_int m){
    s->L->n = n+m;
    s->L->m = n+m;
    s->n = n;
    s->m = m;
    return 0;
}   


c_int compute_Px(qdldl_solver * s, OSQPDataRLDL * data){
    // TODO: Clean up this function...
    c_int ncols = s->L->n;

    c_int nx_ny = data->nx+data->ny;
    c_int nx_nu = data->nx+data->nu;

    c_int * Px0 = (c_int *)c_malloc(sizeof(c_int)*ncols);
    c_int * Px0temp = (c_int *)c_malloc(sizeof(c_int)*ncols);
    c_int Xshift =  (data->Nx-1)*(nx_nu) + data->nu;

    c_int P_counter = 0;
    for (int j=0; j<data->nu;j++) Px0temp[j] = P_counter++;
    for (int i=0; i<data->Nx-1;i++){
	for (int j=0; j<nx_ny; j++)  Px0temp[Xshift   + i*nx_ny + j]  = P_counter++;
	for (int j=0; j<nx_nu; j++)  Px0temp[data->nu + i*nx_nu + j]  = P_counter++;	
    }

    
    for (int i=0;i<ncols;i++){
	Px0[i] = Px0temp[s->P[i]];
    }

    
    for (int i=0;i<(data->Nx-1)*(nx_nu+nx_ny)+data->nu;i++){
	s->P[i] = Px0[i];
    }

    free(Px0);
    free(Px0temp);

}

/*compute_Vhat: computes V*Q*M^(-T)*E^(-1)
 * L    = Lower triangular factorization matrix s.t. P'XP = LDL'
 * Dinv = Inverse of diagonal factorization matrix s.t. P'XP = LDL'
 * P    = Permutation matrix
 * Aij  = Part of V = [Aij, 0]
 * 
 */
static c_int compute_Vhat(csc* L, c_float * Dinv, c_int *P, csc * Aij, csc * Vhat,
			  csc * Yhat, double * Ytemp){
    c_int ncols = L->n;
    c_int nrows = Aij->m;
    c_float* Vtemp= (c_float *)c_calloc(nrows*ncols, sizeof(c_float));
    c_float* Ytemp_copy= (c_float *)c_malloc(nrows*nrows*sizeof(c_float));
    for (int i=0;i<nrows*nrows;i++){ 
	Ytemp_copy[i] = Ytemp[i];
    }

    // Compute V*Q: will have nx elements (assuming Aij=[0,0;-1,0]]
    for (int i=0;i<Aij->n;i++){
        for (int p=Aij->p[i]; p<Aij->p[i+1]; p++)
            Vtemp[P[i] + ncols*Aij->i[p]] = Aij->x[p];
    }

    // Invert by L
    for (int i=0; i<ncols; i++){
        for(int p = L->p[i]; p < L->p[i+1]; p++){
	    for (int k=0; k < nrows; k++){
		if (Vtemp[i + k*ncols] != 0.0)
		    Vtemp[L->i[p] + k*ncols] -= L->x[p]*Vtemp[i + k*ncols];
	    }
        }
    }
    
    // Invert by D and add to Vhat
    int nnz = 0;
    Vhat->p[0] = nnz;
    for (int i=0;i<ncols; i++){
        for (int j=0;j<nrows;j++){
            if (Vtemp[i + j*ncols] != 0){
                Vhat->i[nnz] = j;
                Vhat->x[nnz++] = Vtemp[i + j*ncols]*Dinv[i];
            }
        }
        Vhat->p[i+1] = nnz; 
    }
    Vhat->n = ncols;
    Vhat->m = nrows;

    // Continue with computing Yhat: Invert by L'
    for (int i=0;i<ncols;i++){	
        for (int p=Vhat->p[i];p<Vhat->p[i+1];p++){ // For each row
	    int row = Vhat->i[p];
	    for (int mcolumn=0;mcolumn<nrows;mcolumn++){ // For each column in new matrix
		Ytemp_copy[mcolumn + row*nrows] -=  Vhat->x[p]*Vtemp[i + mcolumn*ncols];
	    }
	}
    }

    nnz = 0;
    Yhat->p[0] = nnz;	
    for (int i=0;i<nrows;i++){ // For every column i in Y
	for(int j=0;j<=i;j++){ // For every row j in current column
	    Yhat->x[nnz] = Ytemp_copy[i + j*nrows];
	    Yhat->i[nnz++] = j;
	}
        Yhat->p[i+1] = nnz;	
    }

    free(Vtemp);
    free(Ytemp_copy);
}

static c_int compute_Uhat(csc* L, c_float * Dinv, c_int *P, csc * Ai, csc * Uhat, c_float * Yhat){
    
    c_int ncols = L->n;
    c_int nrows = Ai->m;
    c_int nz_max = Ai->p[Ai->n];
    c_int Pinv[ncols];

    c_float* Utemp= (c_float *)c_calloc(nrows*ncols, sizeof(c_float));
    
    for (int i=0;i<ncols;i++){
	Pinv[P[i]] = i;
    }
    
    for (int i=0;i<Ai->n;i++){
        for (int p=Ai->p[i]; p<Ai->p[i+1]; p++){
            Utemp[Pinv[ncols - Ai->n + i] + ncols*Ai->i[p]] = Ai->x[p];
        }
    }

    

    // Invert by L
    for (int i=0;i<ncols; i++){
	
        for(int p = L->p[i]; p < L->p[i+1]; p++){ 
            for (int k=0; k < nrows; k++){
		if (Utemp[i + k*ncols] != 0.0){
		    Utemp[L->i[p] + k*ncols] -= L->x[p]*Utemp[i + k*ncols];
		    nz_max++;
		}
            }
        }
    }


    Uhat->n = ncols;
    Uhat->m = nrows;
    Uhat->p = (c_int *)c_malloc((ncols+1)* sizeof(c_int));
    Uhat->i = (c_int *)c_malloc(sizeof(c_int)*nz_max);
    Uhat->x = (c_float *)c_malloc(sizeof(c_float)*nz_max);


    // Invert by D
    int nnz = 0;
    Uhat->p[0] = nnz;

    for (int i=0;i<ncols; i++){
        for (int j=0;j<nrows;j++){
            if (fabs(Utemp[i + j*ncols])>1e-25){
                Uhat->i[nnz] = j;
                Uhat->x[nnz++] = Utemp[i + j*ncols]*Dinv[i];
                //Utemp[i + j*ncols] *= Dinv[i];
            }
        }
        Uhat->p[i+1] = nnz; 
    }

    for (int i=0;i<ncols;i++){	
        for (int p=Uhat->p[i];p<Uhat->p[i+1];p++){
	    int row = Uhat->i[p];
	    for (int mcolumn=0;mcolumn<nrows;mcolumn++){
		if (fabs(Utemp[i + mcolumn*ncols])>1e-25){
		    Yhat[mcolumn + row*nrows] -= Uhat->x[p]*Utemp[i + mcolumn*ncols];
		}
	    }
	}
    }
    

    free(Utemp);

}

static c_int compute_Uhat_solver(qdldl_solver * s, csc * Ai, csc * Uhat, csc * Yhat){
    return compute_Uhat(s->L, s->Dinv, s->P, Ai, Uhat, Yhat); 

}

void update_L_matrix(qdldl_solver * s, OSQPDataRLDL * data, c_int deltaN){
    c_int nx_ny = data->nx+data->ny;
    c_int nx_nu = data->nx+data->nu;
    c_int Nnew = data->N + deltaN;
    // number of columns in the X matrix
    c_int ncols_x = (data->Nx - 1)*(nx_nu+nx_ny)+data->nu;
    
    // Compute final permutation matrix
    // First shift the Px matrix, then add Pz and Py
    c_int P_count = 0;
    for (int j=0;j<ncols_x;j++){
	if (s->P[P_count] >= (data->N)*nx_nu)
	    s->P[P_count++] +=  deltaN*(nx_nu);
	else
	    P_count++;
    }
    for (int j=0;j<data->Lz->n;j++) s->P[P_count++] = data->P0[ncols_x + nx_ny + j];
    for (int j=0;j<nx_ny;j++)       s->P[P_count++] = data->P0[ncols_x + j];
        
    // Lmat: Start by shifting Uhat by nx_nu
    int nnz =0;
    s->L->p[0] = 0;
    for (int i=0; i<ncols_x; i++){
	for (int p=s->L->p[i]; p<s->L->p[i+1];p++){
	    if (s->L->i[p] > ncols_x)
		s->L->i[p] += deltaN*(nx_nu + nx_ny);
	}
    }

    
    // Lmat: next, add the Z matrices: This could be left and just extended: TODO
    nnz = s->L->p[ncols_x];
    for (int i=0; i<data->Lz->n; i++){
	// Add Lz
	s->Dinv[i + ncols_x] = data->Dzinv[i];	
	for (int p=data->Lz->p[i]; p<data->Lz->p[i+1];p++){
	    s->L->x[nnz] = data->Lz->x[p];
	    s->L->i[nnz++] = data->Lz->i[p] + ncols_x;
	}
	// Add Vhat
	for (int p=data->Vhat->p[i]; p<data->Vhat->p[i+1]; p++){
	    s->L->x[nnz] = data->Vhat->x[p];
	    s->L->i[nnz++] = data->Vhat->i[p] + ncols_x + data->Lz->m;
	}
	s->L->p[ncols_x + i+1] = nnz;
    }
    
    // Lmat: Add the Yhat matrix to factorization
    for (int i=0; i<data->Ly->n; i++){
	s->Dinv[i + ncols_x + data->Lz->n] = data->Dyinv[i];
	for (int p=data->Ly->p[i]; p<data->Ly->p[i+1];p++){
	    s->L->x[nnz] = data->Ly->x[p];
	    s->L->i[nnz++] = data->Ly->i[p] + ncols_x + data->Lz->m;
	}
	s->L->p[ncols_x + data->Lz->n + i+1] = nnz;
    }
    
}



void compute_L_matrix(qdldl_solver * s, OSQPDataRLDL * data){
    // Fill in the rest of the L and D matrices:
    // L = [L, 0, 0; 0, M, 0; Uhat, Vhat, Ly]
    // D = diag(D, E, Dx)

    c_int nx_ny = data->nx+data->ny;
    c_int nx_nu = data->nx+data->nu;
    c_int n = data->Nmax*(data->nx+data->ny+data->nx+data->nu) + data->nt;
    c_int deltaN_max = data->Nmax - data->Nx;
    c_int nzmax = s->L->p[s->L->n] + data->Uhat->p[data->Uhat->n]
	+ data->Vhat->p[data->Vhat->n]
	+ data->Ly->p[data->Ly->n]
	+ deltaN_max*data->nnz[1] + data->nnz[0] + 6000;

    // Copy the Lx matrix from the solver.
    // This is needed so that Uhat can be added below L. 
    csc * Lx = s->L;
    c_float * Dx = s->Dinv;
    c_int * Px = s->P;

    // New L matrix will contain the complete factorization
    s->L = csc_spalloc(n,  n, nzmax, 1, 0); // Zero out...

    // Compute P
    // First nu + (N-1)*(nx+ny+nx+nu) variables belong to X
    c_int P_count = 0;
    for (int j=0;j<Lx->n;j++){
	s->P[P_count++] = data->P0[Px[j]];
    }
    for (int j=0;j<data->Lz->n;j++){
	s->P[P_count++] = data->P0[j + Lx->n + nx_ny];
    }
    for (int j=0;j<nx_ny;j++){
	s->P[P_count++] = data->P0[Lx->n+j];
    }
    
    int nnz =0;
    s->L->p[0] = 0;
    for (int i=0; i<Lx->n; i++){
	// Add Dx and Px
	//s->Dinv[i] = Dx[i]; // Already done...
	//s->P[i] = Px[i];
	
	// Add Lx
	for (int p=Lx->p[i]; p<Lx->p[i+1];p++){
	    s->L->x[nnz] = Lx->x[p];
	    s->L->i[nnz++] = Lx->i[p];
	}
	// Add Uhat
	for (int p=data->Uhat->p[i]; p<data->Uhat->p[i+1]; p++){
	    s->L->x[nnz] = data->Uhat->x[p];
	    s->L->i[nnz++] = data->Uhat->i[p] + Lx->m + data->Lz->m;
	}
	s->L->p[i+1] = nnz;
    }

    for (int i=0; i<data->Lz->n; i++){
	// Add Dz and Pz
	s->Dinv[i + Lx->n] = data->Dzinv[i]; // Is this the inverse?
	//s->P[i + Lx->n] = data->Pz[i] + Lx->m;
	
	// Add Lz
	for (int p=data->Lz->p[i]; p<data->Lz->p[i+1];p++){
	    s->L->x[nnz] = data->Lz->x[p];
	    s->L->i[nnz++] = data->Lz->i[p] + Lx->m;
	}
	// Add Vhat
	for (int p=data->Vhat->p[i]; p<data->Vhat->p[i+1]; p++){
	    s->L->x[nnz] = data->Vhat->x[p];
	    s->L->i[nnz++] = data->Vhat->i[p] + Lx->m + data->Lz->m;
	}
	s->L->p[Lx->n + i+1] = nnz;
    }

    for (int i=0; i<data->Ly->n; i++){
	// Add Dy and Py
	s->Dinv[i + Lx->n + data->Lz->n] = data->Dyinv[i];
	//s->P[i + Lx->n + data->Lz->n] = data->Py[i] + Lx->m + data->Lz->m;
	
	for (int p=data->Ly->p[i]; p<data->Ly->p[i+1];p++){
	    s->L->x[nnz] = data->Ly->x[p];
	    s->L->i[nnz++] = data->Ly->i[p] + Lx->m + data->Lz->m;
	}
	s->L->p[Lx->n + data->Lz->n + i+1] = nnz;
    }

    //Lx->n + data->Lz->n + data->Ly->n,
//	   Lx->n + data->Lz->n + data->Ly->n);
    s->L->n = Lx->n + data->Lz->n + data->Ly->n;
    s->L->m = Lx->n + data->Lz->n + data->Ly->n;
    csc_spfree(Lx);
    
}


// This is for pivoting when the RHS is [0,-I; 0, 0]
// Then we can use that the final matrix has the same zero structure
static c_int pivot_odd( const csc *X, //csc *L,
			csc *Ybar_T,
			csc *YbartL0,
			csc *L_csc,
			c_float *Dinv_arr, 
			c_int size1, // number of rows in Aij == L->n == L->m == X->n == X->m
			c_int size2, // number of columns in Aij == Ybar_T->m == Ybar_TL0->m
			c_int *count_L,
			c_int start_idx,
			c_int nx, c_int nu, c_int ny){
    // Allocate memory for the L variable
    c_int i, j, k;
    c_int sum_Lnz;
    c_int factor_status;
    c_int Lp[X->n+1];
    c_int Li[size1*size1];
    c_float Lx[size1*size1];

    c_float Dinv[size1];
    c_float D[size1];

    c_int etree[size1];
    c_int Lnz[size1];    
    c_int  iwork[size1*3];
    c_int  bwork[size1];
    c_float  fwork[size1];
    
    
    c_float * Aij = (c_float*) c_calloc(size1*size2, sizeof(c_float));
    memset(Aij, 0, (size1*size2)*sizeof(c_float));
    for (i=0; i<nx;i++){
        Aij[ny + i + i*(size1)] = -1.0;
    }

    // LDL factorization
    sum_Lnz = QDLDL_etree(X->n, X->p, X->i, iwork, Lnz, etree);
    factor_status = QDLDL_factor(X->n, X->p, X->i, X->x,
				 Lp, Li, Lx,
				 D, Dinv, Lnz,
				 etree, bwork, iwork, fwork);
    
    if (!factor_status){
	printf("ERROR IN FACTORIZATION: %i\n", factor_status);
	return factor_status;
    }
    
    // SOLVE MATRIX EQUATION for Ybar_T --------------
    QDLDL_Lsolve_mat(size1, Lp, Li, Lx, size2, Aij);
    
    for(i = 0; i < size1; i++) {
        for (k=0; k < size2; k++)
            Aij[i + k*size1] *= Dinv[i];
    }
    // Before final solve: Save L'*Ybar for later use, but transposed
    k = 0;
    c_int nz_count = 0;
    YbartL0->p[0] = 0;
    // Current Aij is L'*Ybar (dimension = 22 rows, 12 columns)
    for (i=0; i<size1; i++){ // Loop over each row
	for (j=0; j<size2; j++){ // Loop over each column
	    if (fabs(Aij[i + j*size1]) > 1e-20){
		YbartL0->i[k] = j;
		YbartL0->x[k++] = -Aij[i + size1*j];
		nz_count++;
	    }
	}
	YbartL0->p[i+1] = nz_count;
    }
    
    // Continue solving
    QDLDL_Ltsolve_mat(size1,Lp,Li,Lx,size2,Aij);
    // -----------------------------------------------

    // Save Ybar_T
    k = 0;
    nz_count = 0;
    for (i=0; i<size1; i++){ // Loop over rows
	Ybar_T->p[i] = nz_count;
	for (j=0; j<nx; j++){ // With is size2, but zeros after nx, ignore these
	    if (fabs(Aij[i + j*size1]) > 1e-20){
		Ybar_T->i[k] = j;
		Ybar_T->x[k++] = Aij[i + j*size1];
		nz_count++;
	    }
	}
    }
    Ybar_T->p[size1] = nz_count;
    Ybar_T->n = size1;
    Ybar_T->m = size2;    

	
    // Ybar is given in row major format/ this is Ybar transpose!
    YbartL0->n = Ybar_T->n; // Same number of columns as L
    YbartL0->m = Ybar_T->m; // Same number of rows as Ybar
    // TODO: Remember that Ybar is supposed to be negative here...
    // nrows=nx (only consider first nx rows), ncols, values (bool)
    // Because of the triangular shape - only last nx rows are included
    A_times_B_plus_I (nx, size1, ny, 1,
		      Ybar_T->x, Ybar_T->i, Ybar_T->p,
		      Lx, Li, Lp,
		      YbartL0->x, YbartL0->i, YbartL0->p,
		      nx*size1, -1.0);

    // Add Ltemp and Y0*Ltemp to large L matrix
    for (i=0; i<size1; i++){
	for (j=Lp[i]; j<Lp[i+1]; j++){
	    if (fabs(Lx[j])>1e-10){
		L_csc->i[*count_L] = start_idx + Li[j];
		L_csc->x[(*count_L)++] = Lx[j];
	    }
	}
	for (j=YbartL0->p[i]; j<YbartL0->p[i+1]; j++){
	    if (fabs(YbartL0->x[j])>1e-10){
		L_csc->i[*count_L] = start_idx + size1 + YbartL0->i[j];
		L_csc->x[(*count_L)++] = YbartL0->x[j];
	    }
	}
	L_csc->p[start_idx+i+1] = *count_L;
    }
    for ( i = 0; i < size1; i++){    
	Dinv_arr[start_idx + i] = -Dinv[i]; // Negative
    }
    
    free(Aij);

    return 0;
}





static c_int pivot_even(const csc *X, //csc *L,
			csc *Ybar_T,
			csc *YbartL0,
			c_float *A_f,
			csc *L_csc,
			c_float *Dinv_arr, 
			c_int size1,
			c_int size2,
			c_int *count_L,
			int start_idx,
			c_int * P,
			c_int * perm_count,
			c_int * P_count){
    c_int factor_status;
    c_int i, j, k;
    c_int sum_Lnz;

    c_int Lp[X->n+1];
    c_int Li[size1*size1];
    c_float Lx[size1*size1];

    c_float Dinv[size1];
    c_float D[size1];

    c_int etree[size1];
    c_int Lnz[size1];

    c_int  iwork[size1*3];
    c_int  bwork[size1];
    c_float  fwork[size1];

    if(X->n != size1)
	printf("WARNING: X->n != size1\n");

    for(i=0;i<size1;i++) P[perm_count[0]++] = P_count[0]++;
    
    sum_Lnz = QDLDL_etree(X->n, X->p, X->i, iwork, Lnz, etree);
    factor_status = QDLDL_factor(X->n, X->p, X->i, X->x,
				 Lp, Li, Lx,
				 D, Dinv, Lnz,
				 etree, bwork, iwork, fwork);

    if (!factor_status){
	printf("ERROR IN FACTORIZATION: %i\n", factor_status);
	return factor_status;
    }

    
 
    // SOLVE MATRIX EQUATION --------------
    QDLDL_Lsolve_mat(size1, Lp, Li, Lx, size2, A_f);
    for(c_int i = 0; i < size1; i++) {
        for (k=0; k < size2; k++){
            A_f[i + k*size1] *= Dinv[i];
        }
    }
    // Before final solve: Save L'*Ybar for later use, but transposed
    k = 0;
    c_int nz_count = 0;
    YbartL0->p[0] = 0;
    for (i=0; i<size1; i++){ // Loop over each row
	for (j=0; j<size2; j++){ // Loop over each column
	    if (fabs(A_f[i + j*size1]) > 1e-10){
		YbartL0->i[k] = j;
		YbartL0->x[k++] = A_f[i + size1*j];
		nz_count++;
	    }
	}
	YbartL0->p[i+1] = nz_count;
    }
    QDLDL_Ltsolve_mat(size1,Lp,Li,Lx,size2,A_f);
    // Save Ybar_T
    nz_count = 0; k = 0;
    for (i=0; i<size1; i++){ // Loop over rows
	Ybar_T->p[i] = nz_count;
	for (j=0; j<size2; j++){ // For every column
	    if (fabs(A_f[i + j*size1]) > 10e-10){
		Ybar_T->i[k] = j;
		Ybar_T->x[k++] = A_f[i + j*size1];
		nz_count++;
	    }
	}
    }
    // ------------------------------------


    Ybar_T->p[size1] = nz_count;
    Ybar_T->n = size1;
    Ybar_T->m = size2;

    YbartL0->n = size1;
    YbartL0->m = size2;
    
    A_times_B_plus_I (size2, size1, 0, 1,
		      Ybar_T->x, Ybar_T->i, Ybar_T->p,
		      Lx, Li, Lp,
		      YbartL0->x, YbartL0->i, YbartL0->p,
		      size1*size2, 1.0);

    // Add the L0 and Ybar_TL0 values to the L matrix
    for (i=0; i<size1; i++){
	// Add L matrix
	for (j=Lp[i]; j<Lp[i+1]; j++){
	    if (fabs(Lx[j])>1e-10){
		L_csc->i[*count_L] = start_idx + Li[j];
		L_csc->x[(*count_L)++] = Lx[j];
	    }
	}
	// Add Aii*Ybar matrix under
	for (j=YbartL0->p[i]; j<YbartL0->p[i+1]; j++){
	    if (fabs(YbartL0->x[j])>1e-10){	    
		L_csc->i[*count_L] = start_idx + size1 + YbartL0->i[j];		
		L_csc->x[(*count_L)++] = YbartL0->x[j];
	    }
	}
	L_csc->p[start_idx  + i + 1] = *count_L;
    }
    for ( i = 0; i < size1; i++){    
	Dinv_arr[start_idx + i] = Dinv[i];
    }

    return 0;
}




static c_int pivot_final(const csc *X, //csc *L,
			 csc * L_csc,
			 c_float * Dinv_arr, 
			 c_int num_cols,
			 c_int * count_L,
			 c_int start_idx,
			 c_int * P,
			 c_int * perm_count,
			 c_int * A_count,
			 c_int nx){
    c_int factor_status;
    c_int i, j, k;
    c_int sum_Lnz;

    c_int Lp[X->n+1];
    c_int Li[num_cols*num_cols];
    c_float Lx[num_cols*num_cols];

    c_float Dinv[num_cols];
    c_float D[num_cols];

    
    c_int etree[num_cols];
    c_int Lnz[num_cols];

    c_int  iwork[num_cols*3];
    c_int  bwork[num_cols];
    c_float  fwork[num_cols];
    


    //L->n = X->n;
    //L->m = X->m;
    /*
      c_float *info;
      c_int amd_status;
      c_int  pP[num_cols]; // Permutation matrix
      c_int * Pinv;
      csc *Xnew;    
      #ifdef DLONG
      amd_status = amd_l_order(X->n, X->p, X->i, pP, (c_float *)OSQP_NULL, info);
      #else
      amd_status = amd_order(X->n, X->p, X->i, pP, (c_float *)OSQP_NULL, info);
      #endif
    
      Pinv = csc_pinv(pP, X->n);   
      Xnew = csc_symperm(X, Pinv, OSQP_NULL, 1);  // Xnew = PXP'

      for(i=0;i<num_cols;i++) {
      P[perm_count[0] + i] = A_count[0] +pP[i]; // + pP[i];
      }
      A_count[0] += num_cols;
      perm_count[0] += num_cols;
    */


    for(i=0;i<num_cols;i++) P[perm_count[0]++] = A_count[0]++;
    
    sum_Lnz = QDLDL_etree(X->n, X->p, X->i, iwork, Lnz, etree);
    factor_status = QDLDL_factor(X->n, X->p, X->i, X->x,
				 Lp, Li, Lx,
				 D, Dinv, Lnz,
				 etree, bwork, iwork, fwork);

    
    int nnz = 0, added=0, pointer;
    c_float Xx[num_cols];
    c_int Xi[num_cols];
    c_int Pi[num_cols];
    
    // P = [0, 3, 4, 5, 2, 6]
    // Meaning: first comes row number 0, then row number 3, then number 4
    /*
      for (i=0; i<nx; i++){
      nnz = 0;
      pointer = -1;
      for (j = L_csc->p[start_idx-nx+i]; j<L_csc->p[start_idx-nx+i+1]; j++){
      if(L_csc->i[j]>= start_idx){
      if(nnz==0) pointer = j;		
      Xi[nnz] = L_csc->i[j];
      Pi[nnz] = Pinv[L_csc->i[j]-start_idx]+start_idx;
      Xx[nnz] = L_csc->x[j];
      nnz++;
      }
      }
      if (nnz>0){
      added = 0;
      for (int k=0;k<num_cols;k++){
      if(added < nnz){
      for (j = 0; j<nnz; j++){
      if (start_idx+k==Pi[j]){
      // Which pointer should  use?
      L_csc->i[pointer] = Pi[j];
      L_csc->x[pointer] = Xx[j];
      added++;
      pointer++;
      }
      }
      }
      }
      }
	
      }
    */
    if (!factor_status){
	printf("ERROR IN FACTORIZATION: %i\n", factor_status);
	return factor_status;
    }
    // Save to L matrix
    for (i=0; i<num_cols; i++){
	for (j = Lp[i]; j<Lp[i+1]; j++){
	    L_csc->i[count_L[0]] = start_idx + Li[j];
	    L_csc->x[count_L[0]++] = Lx[j];
	}
	L_csc->p[start_idx + i+1] = count_L[0];
    }
    // Save to diagonal matrix
    for ( i = 0; i < num_cols; i++){    
	Dinv_arr[start_idx + i] = -Dinv[i];
    }

    return 0;
}


/**
 * Compute LDL factorization of matrix A
 * @param  P    Cost matrix to be factorized
 * @param  A    Constraint matrix to be factorized
 * @param  p    Private workspace
 * @param  nvar Number of QP variables
 * @return      exitstatus
 */
static c_int LDL_update_from_pivot(csc * Lmat, c_float * Dinv, c_int *Pmat, csc **X_even,
				   csc** Qmats, csc** Amats, 
				   c_float     sigma,
				   const c_float    *rho_inv_vec,
				   c_int nx,   c_int nu, c_int ny, c_int nt,
				   c_int deltaN, c_int Nnew, c_int N0, c_int Nmax){
    // N0 is the start value (0 for normal)
    // Nnew is the total value to reach
    // New size of L will depend on New and N0
    c_int status;
    c_int i;
    c_int Lmat_ptr = 0; // Where in Lmat to add elements
    c_int Lmat_col = 0; // What column in Lmat to add to
    c_int iter_start;
    c_int nx_ny = nx+ny;  
    c_int nx_nu = nx+nu;

    if (deltaN>0)  {            // Increasing horizon
        iter_start = Nnew-deltaN-1;
    } else if (deltaN<0)  {     // Reduced horizon
        iter_start = Nnew-1;
    } else return -1;
    // The column and pointer to continue from in L
    Lmat_col = Qmats[0]->n + (iter_start-N0)*(Qmats[1]->m + Amats[1]->m);
    Lmat_ptr = Lmat->p[Lmat_col];

    // Permutation matrix changed: Before A started at Nold*(nx+nu), now at Nnew*(nx+nu)
    // All A values gets new values. Q only gets new values if deltaN>0
    c_int perm_count = Qmats[0]->n;
    // Cannot know shift without assuming all Qi are same size
    //c_int A_count = Nnew*(Qmats[1]->m);
    c_int A_count = 0;
    for (int i=N0;i<=Nnew;i++)  A_count += Qmats[i]->m;
    c_int P_count = Qmats[0]->m + iter_start*(Qmats[1]->m);
    

    for(i=N0;i<Amats[0]->m;i++) Pmat[perm_count++] = A_count++;
    for (i=N0; i<iter_start; i++){
	perm_count += Qmats[i+1]->m;
	for(int j=0; j<Amats[i+1]->m; j++) Pmat[perm_count++] = A_count++;
    }
    // Temporary matrices    ===========================================================
    csc      *Aii = csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);    
    csc      *Ybar_T = csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);
    csc      *Xtemp = csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);    
    csc      *Aii_transpose =  csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);
    csc      *YbartL0 = csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);
    c_float  * Aij_f = (c_float*) c_calloc(nx_ny*nx_ny, sizeof(c_float));
    
    prea_int_vec_copy(X_even[iter_start]->p, Xtemp->p, X_even[iter_start]->n + 1);
    prea_int_vec_copy(X_even[iter_start]->i, Xtemp->i, X_even[iter_start]->p[X_even[iter_start]->n]);
    prea_vec_copy(X_even[iter_start]->x, Xtemp->x, X_even[iter_start]->p[X_even[iter_start]->n]);
    if (deltaN>0){
	for (c_int iter=iter_start+1;iter<Nnew;iter++){
	    if (pivot_odd(Xtemp, Ybar_T,  YbartL0,
			  Lmat, Dinv,
			  Xtemp->n, nx_nu, &Lmat_ptr, Lmat_col, nx, nu, ny)){
		return -1;
	    }
	    Lmat_col += Xtemp->n;
	    
	    // Compute next Xtemp: Qi - Ybar_T + I*sigma
	    copy_csc_plus_sigma(Qmats[iter], Xtemp, sigma);
	    A_minus_B(ny, ny+nx, Xtemp->n, 
		      Xtemp->p, Xtemp->i, Xtemp->x,
		      Ybar_T->p, Ybar_T->i, Ybar_T->x );

	    // Copy next RHS (Ai) to Aij_f
	    memset(Aij_f, 0, Amats[iter]->n*Amats[iter]->m*sizeof(c_float));
	    copy_csc_transpose(Amats[iter], Aij_f);
	    
	    if (pivot_even(Xtemp, Ybar_T,  YbartL0,
			   Aij_f, Lmat, Dinv,
			   Xtemp->n, nx_ny, &Lmat_ptr, Lmat_col, 
			   Pmat, &perm_count, &P_count)){
		printf("Error in Pivot even, iter = %i\n", iter);
		return -1;
	    }
	    
	    Lmat_col += Xtemp->n;
	    // Compute the next Xtemp:  -rhoinv - Ai*Ybar_ij
	    for(i=0;i<Amats[iter]->m;i++) Pmat[perm_count++] = A_count++; // Permutation: Ai	    
	    copy_csc_matrix(Amats[iter], Aii);
	    csc_to_csr(nx_ny, nx_nu,
		       Aii->p, Aii->i, Aii->x,
		       Aii_transpose->p, Aii_transpose->i, Aii_transpose->x);// Aii_transpose = Ai'
	    
	    Xtemp->n = nx_ny;
	    Xtemp->m = nx_ny;

	    status = A_times_B_plus_rho(nx_ny, 1,
					Ybar_T->x, Ybar_T->i, Ybar_T->p,
					Aii_transpose->x, Aii_transpose->i, Aii_transpose->p,
					Xtemp->x, Xtemp->i, Xtemp->p,
					nx_ny*nx_ny/2+nx_ny, &rho_inv_vec[(iter)*nx_ny]);
	    if (status) {
		printf("Could not complete X =  -rhoinv - Ai*Ybar_ij\n");
		return -1;
	    }
	    // Save all even X for later updates
	    if (X_even[iter] != OSQP_NULL){
		csc_spfree(X_even[iter]);
	    }
	    X_even[iter] = copy_csc_mat(Xtemp);
	}
    }
    
    // End with terminal cost and constraints
    pivot_odd(Xtemp, Ybar_T,  YbartL0,
	      Lmat, Dinv,
	      nx_ny, nx, &Lmat_ptr, Lmat_col, nx, nu, ny );
    Lmat_col += Xtemp->n;
    
    // Terminal cost  ------------------------------------------------
    copy_csc_plus_sigma(Qmats[Nnew], Xtemp, sigma);
    A_minus_B(ny, ny+nx, nx, 
	      Xtemp->p, Xtemp->i, Xtemp->x,
	      Ybar_T->p, Ybar_T->i, Ybar_T->x );

    memset(Aij_f, 0, nt*nx*sizeof(c_float));
    copy_csc_transpose(Amats[Nnew], Aij_f);
    
    pivot_even(Xtemp, Ybar_T,  YbartL0,
	       Aij_f, Lmat, Dinv,
               nx, nt, &Lmat_ptr, Lmat_col, 
	       Pmat, &perm_count, &P_count);
    Lmat_col += Xtemp->n;


    // Terminal constraint -------------------------------------------
    copy_csc_matrix(Amats[Nnew], Aii);
    csc_to_csr(nt, nx,
	       Aii->p, Aii->i, Aii->x,
	       Aii_transpose->p, Aii_transpose->i, Aii_transpose->x);
    
    status = A_times_B_plus_rho(nt, 1,
				Ybar_T->x, Ybar_T->i, Ybar_T->p,
				Aii_transpose->x, Aii_transpose->i, Aii_transpose->p,
				Xtemp->x, Xtemp->i, Xtemp->p,
				nt*nt/2+nt, &rho_inv_vec[(Nmax-1)*nx_ny]);
    if (status) {
	printf("Could not complete terminal X =  -rhoinv - Ai*Ybar_ij\n");
	return -1;
    }    
    Xtemp->n = nt;
    Xtemp->m = nt;

    pivot_final(Xtemp, Lmat, Dinv,
		nt, &Lmat_ptr, Lmat_col,
		Pmat, &perm_count, &A_count, nx );
    Lmat_col += Xtemp->n;

    Lmat->n = Lmat_col;
    Lmat->m = Lmat_col;

    csc_spfree(Aii);
    csc_spfree(Xtemp);
    csc_spfree(Ybar_T);
    
    csc_spfree(Aii_transpose);
    csc_spfree(YbartL0);
    free(Aij_f);

    return 0;
}

static c_int LDL_update_from_pivot_solver(qdldl_solver * s, csc **X_even, csc** Qmats, csc** Amats, 
					  c_float     sigma,
					  const c_float    *rho_inv_vec,
					  c_int nx,   c_int nu, c_int ny, c_int nt,
					  c_int deltaN, c_int Nnew, c_int Nmax){


    
    LDL_update_from_pivot(s->L, s->Dinv, s->P, X_even, 
			  Qmats, Amats, 
			  sigma, rho_inv_vec,
			  nx, nu,  ny,  nt,
			  deltaN, Nnew, 0, Nmax);
    
    update_problem_size(s, Nnew*(nx+nu), Nnew*(nx+ny)+nt);
    
}


/**
 * Compute LDL factorization of matrix A
 * @param  P    Cost matrix to be factorized
 * @param  A    Constraint matrix to be factorized
 * @param  p    Private workspace
 * @param  nvar Number of QP variables
 * @return      exitstatus (0 is good)
 */
static c_int LDL_factorize_recursive(csc * Lmat, c_float * Dinv, c_int *Pmat, csc **X_even,
				     csc** Qmats, csc** Amats, c_float sigma,
				     c_float  *rho_inv_vec,
				     c_int nx,   c_int nu, c_int ny, c_int nt,
				     c_int Nmax,  c_int N0, c_int Niter,
				     c_int * nnz_params){

    // N: Factorize to this N
    // Nmax: maximum allowed: used only in rho_invese to get last elements....
    // N0: start factorization from this element: 0 for recursive. 
    c_int i, j, k;
    c_int Lmat_ptr = 0; // Where in Lmat to add elements
    c_int Lmat_col = 0; // What column in Lmat to add to
    
    c_int perm_count=0, A_count=(Niter)*(nx+nu), P_count=0;
    
    c_int status;
    c_int nz_xeven = 0;

    c_int nx_ny = nx+ny;  
    c_int nx_nu = nx+nu;
    
    // Temporary matrices    ===========================================================
    // Be more careful with how many zeros are allocated?
    csc *Xtemp = csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);    
    csc * Ybar_T = csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);
    csc * Aii_transpose =  csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);
    csc * YbartL0 = csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);
    c_float * Aij_f = (c_float*) c_calloc(nx_ny*nx_ny, sizeof(c_float));
    
    // PIVOT POINT 0 =================================
    Lmat->p[0] = 0;
    //for(i=0;i<Q0->n;i++) Pmat[perm_count++] = P_count++; //  Q0 (even) nu
    copy_csc_transpose(Amats[0], Aij_f);
    // Assuming no non-zero elements on diagonal...
    copy_csc_plus_sigma(Qmats[0], Xtemp, sigma);

    pivot_even(Xtemp, Ybar_T,  YbartL0,
	       Aij_f, Lmat, Dinv,
               Xtemp->n, Amats[0]->m, &Lmat_ptr,  0, 
	       Pmat, &perm_count, &P_count);
    nnz_params[0] = Lmat_ptr;
    nnz_params[1] = Lmat_ptr;
    
    // Next matrix: X = -rhoinv - A0*Ybar_01  
    for(i=0;i<Amats[0]->m;i++) Pmat[perm_count++] = A_count++; // Permutation: A0

    csc_to_csr(Amats[0]->m, Amats[0]->n, Amats[0]->p, Amats[0]->i, Amats[0]->x,
	       Aii_transpose->p, Aii_transpose->i, Aii_transpose->x); 

    // Computes new X = Ybar_T*Aii_transpose + rhoinv*I
    Xtemp->n = Amats[0]->m;
    Xtemp->m = Amats[0]->m;
    status = A_times_B_plus_rho(Amats[0]->m, 1,
				Ybar_T->x, Ybar_T->i, Ybar_T->p,
				Aii_transpose->x, Aii_transpose->i, Aii_transpose->p,
				Xtemp->x, Xtemp->i, Xtemp->p,
				ceil(Xtemp->n*Xtemp->m/2+Xtemp->n)+10, &rho_inv_vec[(N0)*nx_ny]);
    
    if (status) {
	printf("Could not complete X =  -rhoinv - Ai*Ybar_ij\n");
	return -1;
    }



    // Save all even X for later updates
    X_even[N0] = copy_csc_mat(Xtemp);
    
    Lmat_col += Qmats[0]->n; // Start index in L matrix    
    for (c_int iter=1; iter<Niter; iter++){
	// Constraint pivot: X = -rhoinv*I - Ai*Ybar_ij
	pivot_odd(Xtemp, Ybar_T,  YbartL0, Lmat, Dinv,
		  Xtemp->n, nx_nu, &Lmat_ptr, Lmat_col, nx, nu, ny );
        Lmat_col += Xtemp->n;

	// Compute next Xtemp: Qi - Ybar_T
	//for(i=0;i<Qi->n;i++) Pmat[perm_count++] = P_count++; // Permutation: Qi (even) add Qi->n
	copy_csc_plus_sigma(Qmats[iter], Xtemp, sigma);
	// Subtract bottom nx rows of Ybar_T from top of X
	// Ybar_T = Ybar given in row major format - bottom nx rows = rightmost nx?
	A_minus_B(ny, ny+nx, Xtemp->n, 
		  Xtemp->p, Xtemp->i, Xtemp->x,
		  Ybar_T->p, Ybar_T->i, Ybar_T->x );

	memset(Aij_f, 0, Amats[iter]->n*Amats[iter]->m*sizeof(c_float));
	copy_csc_transpose(Amats[iter], Aij_f);
	        
	pivot_even(Xtemp, Ybar_T,  YbartL0,
		   Aij_f, Lmat, Dinv,
                   Xtemp->n, Amats[iter]->m, &Lmat_ptr, Lmat_col, 
		   Pmat, &perm_count, &P_count);
        Lmat_col += Xtemp->n;

	// Compute the next Xtemp:  -rhoinv - Ai*Ybar_ij
	for(i=0;i<Amats[iter]->m;i++) Pmat[perm_count++] = A_count++; // Permutation: Ai	
	csc_to_csr(Amats[iter]->m, Amats[iter]->n,
		   Amats[iter]->p, Amats[iter]->i, Amats[iter]->x,
		   Aii_transpose->p, Aii_transpose->i, Aii_transpose->x);// Aii_transpose = Ai'

	Xtemp->n = Amats[iter]->m;
	Xtemp->m = Amats[iter]->m;
	status = A_times_B_plus_rho(Amats[iter]->m, 1,
				    Ybar_T->x, Ybar_T->i, Ybar_T->p,
				    Aii_transpose->x, Aii_transpose->i, Aii_transpose->p,
				    Xtemp->x, Xtemp->i, Xtemp->p,
				    nx_ny*nx_ny/2+nx_ny, &rho_inv_vec[(iter)*nx_ny]);
	if (status) {
	    printf("Could not complete X =  -rhoinv - Ai*Ybar_ij\n");
	    return -1;
	}
        
        // Save all even X for later updates
	X_even[N0+iter] = copy_csc_mat(Xtemp);
	
    }
    nnz_params[1] = Lmat_ptr - nnz_params[0];

    // Second to last (pivot 2N, CONSTRAINT) ---------------------------------------------
    
    pivot_odd(Xtemp, Ybar_T,  YbartL0, Lmat, Dinv,
	      nx_ny, nx, &Lmat_ptr, Lmat_col, nx, nu, ny );
    Lmat_col += Xtemp->n;
    // Terminal cost constraint -------------------------------------
    //for(i=0;i<QN->n;i++) Pmat[perm_count++] = P_count++; // (even) add nx
    copy_csc_plus_sigma(Qmats[Niter], Xtemp, sigma);
    A_minus_B(ny, ny+nx, nx, Xtemp->p, Xtemp->i, Xtemp->x,
	      Ybar_T->p, Ybar_T->i, Ybar_T->x );

    memset(Aij_f, 0, nt*nx*sizeof(c_float));
    copy_csc_transpose(Amats[Niter], Aij_f);
    
    pivot_even(Xtemp, Ybar_T,  YbartL0,
	       Aij_f, Lmat, Dinv,
               nx, nt, &Lmat_ptr, Lmat_col,
	       Pmat, &perm_count, &P_count);
    Lmat_col += Xtemp->n;
    // Only include upper diagonal items...
    // For next X: need AN' and Ybar from previous solution


    //--------------------------------------------------------------------------------
    // Terminal constraint 
    // A = AN in R(nt * nx)
    //for(i=0;i<AN->m;i++) Pmat[perm_count++] = A_count++;    
    
    csc_to_csr(nt, nx,
               Amats[Niter]->p, Amats[Niter]->i, Amats[Niter]->x,
               Aii_transpose->p, Aii_transpose->i, Aii_transpose->x);
    
    A_times_B_plus_rho(nt, 1,
		       Ybar_T->x, Ybar_T->i, Ybar_T->p,
		       Aii_transpose->x, Aii_transpose->i, Aii_transpose->p,
		       Xtemp->x, Xtemp->i, Xtemp->p,
		       nt*nt/2+nt, &rho_inv_vec[(Nmax-1)*nx_ny]);
    Xtemp->n = nt;
    Xtemp->m = nt;

    pivot_final(Xtemp,  Lmat, Dinv, nt, &Lmat_ptr, Lmat_col, Pmat, &perm_count, &A_count,nx);
    nnz_params[0] = Lmat_ptr - nnz_params[1];
    nnz_params[1] = nnz_params[1]/(Niter);
    Lmat_col += Xtemp->n;
    // Fill the rest with zeros (not really necessary?)
    Lmat->p[Lmat_col+1] = Lmat_ptr;
    //for (c_int i = N*(nx+nu)+N*(nx+ny)+nt; i<Nmax*(nx+nu) + Nmax*(nx+ny)+nt;i++){
    //    Lmat->p[i] = Lmat_ptr;
    //}

    Lmat->n = Lmat_col;
    Lmat->m = Lmat_col;

    csc_spfree(Xtemp);
    csc_spfree(Ybar_T);
    csc_spfree(Aii_transpose);
    csc_spfree(YbartL0);
    free(Aij_f);

    return 0;

}

static c_int LDL_factorize_recursive_solver(qdldl_solver * s, csc **X_even,
					    csc ** Amats, csc **Qmats,  
					    c_float     sigma,
					    c_float    *rho_inv_vec,
					    c_int nx,   c_int nu, c_int ny, c_int nt,
					    c_int Nmax, c_int Niter,c_int * nnz_params){


    return LDL_factorize_recursive(s->L, s->Dinv, s->P, X_even, Qmats, Amats, 
				   sigma, rho_inv_vec, nx, nu,  ny,  nt, Nmax, 0, Niter, nnz_params);

}


void get_L_dimensions_from_solver(qdldl_solver * s, c_int *n, c_int *m, c_int *nnz){
    *n = s->L->n;
    *m = s->L->m;
    *nnz = s->L->p[s->L->n];
}

void get_L_dimensions(OSQPWorkspace* work, c_int *n, c_int *m, c_int *nnz){
    get_L_dimensions_from_solver(work->linsys_solver, n, m, nnz);
}


void compute_permutations(qdldl_solver * s, c_int n, c_int m, c_int N, 
			  const csc * Q0, const csc * Qi, const csc * QN, 
			  const csc * A0, const csc * Ai, const csc * Aij, const csc * AN){
    // PERMUTATION VECTOR P: A1 = A[P,P]  =============================================
    c_int i, perm_count=0, A_count=n, P_count=0;
    for(i=0;i<Q0->n;i++) {
        s->P[perm_count++] = P_count++;
    }
    for(i=0;i<A0->m;i++) {
	s->P[perm_count++] = A_count++;
    }
    
    for(c_int iter=0;iter<N-1;iter++){
        for(i=0;i<Qi->n;i++) s->P[perm_count++] = P_count++;
        for(i=0;i<Ai->m;i++) s->P[perm_count++] = A_count++;
    }
    for(i=0;i<QN->n;i++) s->P[perm_count++] = P_count++;
    for(i=0;i<AN->m;i++) s->P[perm_count++] = A_count++;
}


void compute_KKT_permutations(qdldl_solver * s, c_int n, c_int m,  c_int P_nnz_max, c_int A_nnz_max, 
			      c_int N, 
			      const csc * Q0, const csc * Qi, const csc * QN, 
			      const csc * A0, const csc * Ai, const csc * Aij, const csc * AN){
    c_int i, count=0, A_count=n, P_count=0;

    // Allocate vector of indices on the diagonal. Worst case it has m elements
    if (s->Pdiag_idx != OSQP_NULL) {
	(*s->Pdiag_idx) = c_malloc(n * sizeof(c_int));
	s->Pdiag_n     = 0; // Set 0 diagonal elements to start
    }
}

csc* form_KKT_permuted(const csc  *P,
		       const  csc *A,
		       c_int       format,
		       c_float     param1,
		       c_float    *param2,
		       c_int      *PtoKKT,
		       c_int      *AtoKKT,
		       c_int     **Pdiag_idx,
		       c_int      *Pdiag_n,
		       c_int      *param2toKKT) {
    c_int  nKKT, nnzKKTmax; // Size, number of nonzeros and max number of nonzeros
    // in KKT matrix
    csc   *KKT_trip, *KKT;  // KKT matrix in triplet format and CSC format
    c_int  ptr, i, j;       // Counters for elements (i,j) and index pointer
    c_int  zKKT = 0;        // Counter for total number of elements in P and in
    // KKT
    c_int *KKT_TtoC;        // Pointer to vector mapping from KKT in triplet form
    // to CSC

    // Get matrix dimensions
    nKKT = P->m + A->m;

    // Get maximum number of nonzero elements (only upper triangular part)
    nnzKKTmax = P->p[P->n] + // Number of elements in P
	P->m +       // Number of elements in param1 * I
	A->p[A->n] + // Number of nonzeros in A
	A->m;        // Number of elements in - diag(param2)

    // Preallocate KKT matrix in triplet format
    KKT_trip = csc_spalloc(nKKT, nKKT, nnzKKTmax, 1, 1);

    if (!KKT_trip) return OSQP_NULL;  // Failed to preallocate matrix

    // Allocate vector of indices on the diagonal. Worst case it has m elements
    if (Pdiag_idx != OSQP_NULL) {
	(*Pdiag_idx) = c_malloc(P->m * sizeof(c_int));
	*Pdiag_n     = 0; // Set 0 diagonal elements to start
    }

    // Allocate Triplet matrices
    // P + param1 I
    for (j = 0; j < P->n; j++) { // cycle over columns
	// No elements in column j => add diagonal element param1
	if (P->p[j] == P->p[j + 1]) {
	    KKT_trip->i[zKKT] = j;
	    KKT_trip->p[zKKT] = j;
	    KKT_trip->x[zKKT] = param1;
	    zKKT++;
	}

	for (ptr = P->p[j]; ptr < P->p[j + 1]; ptr++) { // cycle over rows
	    // Get current row
	    i = P->i[ptr];

	    // Add element of P
	    KKT_trip->i[zKKT] = i;
	    KKT_trip->p[zKKT] = j;
	    KKT_trip->x[zKKT] = P->x[ptr];

	    if (PtoKKT != OSQP_NULL) PtoKKT[ptr] = zKKT;  // Update index from P to
	    // KKTtrip

	    if (i == j) {                                 // P has a diagonal element,
		// add param1
		KKT_trip->x[zKKT] += param1;

		// If index vector pointer supplied -> Store the index
		if (Pdiag_idx != OSQP_NULL) {
		    (*Pdiag_idx)[*Pdiag_n] = ptr;
		    (*Pdiag_n)++;
		}
	    }
	    zKKT++;

	    // Add diagonal param1 in case
	    if ((i < j) &&                  // Diagonal element not reached
		(ptr + 1 == P->p[j + 1])) { // last element of column j
		// Add diagonal element param1
		KKT_trip->i[zKKT] = j;
		KKT_trip->p[zKKT] = j;
		KKT_trip->x[zKKT] = param1;
		zKKT++;
	    }
	}
    }

    if (Pdiag_idx != OSQP_NULL) {
	// Realloc Pdiag_idx so that it contains exactly *Pdiag_n diagonal elements
	(*Pdiag_idx) = c_realloc((*Pdiag_idx), (*Pdiag_n) * sizeof(c_int));
    }


    // A' at top right
    for (j = 0; j < A->n; j++) {                      // Cycle over columns of A
	for (ptr = A->p[j]; ptr < A->p[j + 1]; ptr++) {
	    KKT_trip->p[zKKT] = P->m + A->i[ptr];         // Assign column index from
	    // row index of A
	    KKT_trip->i[zKKT] = j;                        // Assign row index from
	    // column index of A
	    KKT_trip->x[zKKT] = A->x[ptr];                // Assign A value element

	    if (AtoKKT != OSQP_NULL) AtoKKT[ptr] = zKKT;  // Update index from A to
	    // KKTtrip
	    zKKT++;
	}
    }

    // - diag(param2) at bottom right
    for (j = 0; j < A->m; j++) {
	KKT_trip->i[zKKT] = j + P->n;
	KKT_trip->p[zKKT] = j + P->n;
	KKT_trip->x[zKKT] = -param2[j];

	if (param2toKKT != OSQP_NULL) param2toKKT[j] = zKKT;  // Update index from
	// param2 to KKTtrip
	zKKT++;
    }

    // Allocate number of nonzeros
    KKT_trip->nz = zKKT;

    // Convert triplet matrix to csc format
    if (!PtoKKT && !AtoKKT && !param2toKKT) {
	// If no index vectors passed, do not store KKT mapping from Trip to CSC/CSR
	if (format == 0) KKT = triplet_to_csc(KKT_trip, OSQP_NULL);
	else KKT = triplet_to_csr(KKT_trip, OSQP_NULL);
    }
    else {
	// Allocate vector of indices from triplet to csc
	KKT_TtoC = c_malloc((zKKT) * sizeof(c_int));

	if (!KKT_TtoC) {
	    // Error in allocating KKT_TtoC vector
	    csc_spfree(KKT_trip);
	    c_free(*Pdiag_idx);
	    return OSQP_NULL;
	}

	// Store KKT mapping from Trip to CSC/CSR
	if (format == 0)
	    KKT = triplet_to_csc(KKT_trip, KKT_TtoC);
	else
	    KKT = triplet_to_csr(KKT_trip, KKT_TtoC);

	// Update vectors of indices from P, A, param2 to KKT (now in CSC format)
	if (PtoKKT != OSQP_NULL) {
	    for (i = 0; i < P->p[P->n]; i++) {
		PtoKKT[i] = KKT_TtoC[PtoKKT[i]];
	    }
	}

	if (AtoKKT != OSQP_NULL) {
	    for (i = 0; i < A->p[A->n]; i++) {
		AtoKKT[i] = KKT_TtoC[AtoKKT[i]];
	    }
	}

	if (param2toKKT != OSQP_NULL) {
	    for (i = 0; i < A->m; i++) {
		param2toKKT[i] = KKT_TtoC[param2toKKT[i]];
	    }
	}

	// Free mapping
	c_free(KKT_TtoC);
    }

    // Clean matrix in triplet format and return result
    csc_spfree(KKT_trip);

    return KKT;
}



// Initialize LDL Factorization structure
// The result is p->L and p->D
c_int init_linsys_solver_qdldl_recursive(qdldl_solver ** sp, c_int n, c_int m, c_int P_nnz_max, c_int A_nnz_max, c_float sigma, const c_float * rho_vec, c_int polish){
    // Define Variables
    csc * KKT_temp;     // Temporary KKT pointer
    c_int i;            // Loop counter
    c_int n_plus_m;     // Define n_plus_m dimension

    // Allocate private structure to store KKT factorization
    qdldl_solver *s;
    s = c_calloc(1, sizeof(qdldl_solver));
    *sp = s;

    // Size of KKT
    s->n = n;
    s->m = m;
    n_plus_m = s->n + s->m;

    // Sigma parameter
    s->sigma = sigma;

    // Polishing flag
    s->polish = polish;

    // Link Functions
    s->solve = &solve_linsys_qdldl;


#ifndef EMBEDDED
    s->free = &free_linsys_solver_qdldl;
#endif

#if EMBEDDED != 1
    s->update_matrices = &update_linsys_solver_matrices_qdldl;
    s->update_rho_vec = &update_linsys_solver_rho_vec_qdldl;
    //s->update_LDL_factorization = &update_LDL_factorization_qdldl;
#endif

    // Assign type
    s->type = QDLDL_SOLVER;

    // Set number of threads to 1 (single threaded)
    s->nthreads = 1;

    // Sparse matrix L (lower triangular)
    // NB: We don not allocate L completely (CSC elements)
    //      L will be allocated during the factorization depending on the
    //      resulting number of elements.
    s->L = c_malloc(sizeof(csc));
    s->L->m = n_plus_m;
    s->L->n = n_plus_m;
    s->L->nz = -1;


    // Diagonal matrix stored as a vector D
    s->Dinv = (c_float *)c_malloc(sizeof(c_float) * n_plus_m);
    s->D    = (c_float *)c_malloc(sizeof(c_float) * n_plus_m);

    // Elimination tree workspace
    s->etree = (c_int *)c_malloc(n_plus_m * sizeof(c_int));
    s->Lnz   = (c_int *)c_malloc(n_plus_m * sizeof(c_int));
    // Permutation vector P: not defined yet......
    s->P     = (c_int *)c_malloc(n_plus_m * sizeof(c_int));
    
    // Working vector
    s->bp   = (c_float *)c_malloc(sizeof(c_float) * n_plus_m);

    // Solution vector
    s->sol  = (c_float *)c_malloc(sizeof(c_float) * n_plus_m);

    // Parameter vector
    s->rho_inv_vec = (c_float *)c_malloc(sizeof(c_float) * s->m);


    // Preallocate L matrix (Lx and Li are sparsity dependent)
    s->L->p = (c_int *)c_malloc((n_plus_m+1) * sizeof(c_int));
    s->L->p[0] = 0;

    // Lx and Li are sparsity dependent, so set them to
    // null initially so we don't try to free them prematurely
    s->L->i = (c_int *)c_malloc(sizeof(c_int)*((P_nnz_max+A_nnz_max)*2));     //OSQP_NULL;
    s->L->x = (c_float *)c_malloc(sizeof(c_float)*((P_nnz_max+A_nnz_max)*2)); //OSQP_NULL;

    // Preallocate workspace
    s->iwork = (c_int *)c_malloc(sizeof(c_int)*(30*n_plus_m));
    s->bwork = (c_int *)c_malloc(sizeof(c_int)*n_plus_m*50);
    s->fwork = (c_float *)c_malloc(sizeof(c_float)*n_plus_m*50);

    // Form and permute KKT matrix
    if (polish){ // Called from polish()
        // Use s->rho_inv_vec for storing param2 = vec(delta)
        for (i = 0; i < m; i++){
            s->rho_inv_vec[i] = sigma;
        }
    }
	
    else { // Called from ADMM algorithm	
        // Allocate vectors of indices
        s->PtoKKT = c_malloc((P_nnz_max) * sizeof(c_int));  // Number of nonzeros in cost
        s->AtoKKT = c_malloc((A_nnz_max) * sizeof(c_int));  // Number of nonzeros in constraint
        s->rhotoKKT = c_malloc((m) * sizeof(c_int));     // Number of constraints
        // Use p->rho_inv_vec for storing param2 = rho_inv_vec
        for (i = 0; i < m; i++){
            s->rho_inv_vec[i] = 1. / rho_vec[i];
	    
        }
    }
    
    if (polish){ // If KKT passed, assign it to KKT_temp
        // Polish, no need for KKT_temp
        csc_spfree(KKT_temp);
    }
    else { // If not embedded option 1 copy pointer to KKT_temp. Do not free it.
        s->KKT = KKT_temp;
    }


    // No error
    return 0;
}


c_int update_AP_matrices(OSQPDataRLDL * data_rldl, OSQPData * data, c_int Nold, c_int Nnew){

    c_int nz_P, nz_A;
    c_int P_row, P_col;
    c_int A_row, A_col;
    c_int nx = data_rldl->nx;
    c_int nu = data_rldl->nu;
    c_int ny = data_rldl->ny;
    c_int nt = data_rldl->nt;
    c_int nvar = Nnew*(nx+nu);
    c_int nconst = Nnew*(nx+ny)+nt;
    c_int Ntemp;

    csc * Q0 = data_rldl->Q0;
    csc * Qi = data_rldl->Qi;
    csc * QN = data_rldl->QN;
    csc * A0 = data_rldl->A0;
    csc * Ai = data_rldl->Ai;
    csc * Aij = data_rldl->Aij;
    csc * AN = data_rldl->AN;

    data->P->n = nvar;
    data->P->m = nvar;
    data->A->n = nvar;
    data->A->m = nconst;


    if (Nnew < Nold) {  // Return last few columns
        Ntemp = Nnew -1;
        P_row = nu + Ntemp*(nx+nu);
        P_col = nu + Ntemp*(nx+nu);
        A_row = Ntemp*(nx+ny);
        A_col = nu + Ntemp*(nx+nu);
        nz_P = data->P->p[P_col];
        nz_A = data->A->p[A_col];
    }else{ // N is increasing
        Ntemp = Nold -1;
        P_row = nu + Ntemp*(nx+nu);
        P_col = nu + Ntemp*(nx+nu);
        A_row = Ntemp*(nx+ny);
        A_col = nu + Ntemp*(nx+nu);
        nz_P = data->P->p[P_col];
        nz_A = data->A->p[A_col];
    }
    for (c_int iter=Ntemp;iter<Nnew-1;iter++){

        for (c_int i=0; i<Qi->n;i++){                  // For column in Qi
            // Move column i to column P_row + i
            for (c_int j=Qi->p[i];j<Qi->p[i+1];j++){   // Iterate over all rows in Qi
                data->P->i[nz_P] = P_row + Qi->i[j]; 
                data->P->x[nz_P] = Qi->x[j];
                nz_P++;
            }
            data->P->p[P_row+i+1] = nz_P; //  Pointer value for next is nz_P
        }
        P_row += Qi->m;

        for (c_int i=0; i<Aij->n;i++){
            data->A->p[A_col+i] = nz_A;
            for (c_int j=Aij->p[i];j<Aij->p[i+1];j++){
                data->A->i[nz_A] = A_row + Aij->i[j];
                data->A->x[nz_A++] = Aij->x[j];
            }
            for (c_int j=Ai->p[i];j<Ai->p[i+1];j++){
                data->A->i[nz_A] = A_row + Aij->m + Ai->i[j];
                data->A->x[nz_A++] = Ai->x[j];
            }
        }
        A_row += Aij->m;
        A_col += Aij->n;
    } // End iter
    data->P->n = nvar;
    // Add the terminal matrices
    for (c_int i=0; i<QN->n;i++){
        data->P->p[P_row+i] = nz_P;
        for (c_int j=QN->p[i];j<QN->p[i+1];j++){
            data->P->i[nz_P] = P_row + QN->i[j];
            data->P->x[nz_P++] = QN->x[j];
        }
    }
    P_row += QN->m;
    data->P->p[P_row] = nz_P;

    // Add the final AN->n columns
    for (c_int i=0; i<AN->n;i++){

        data->A->p[A_col+i] = nz_A;
        for (c_int j=Aij->p[i];j<Aij->p[i+1];j++){
            data->A->i[nz_A] = A_row + Aij->i[j];
            data->A->x[nz_A++] = Aij->x[j];
        }
        for (c_int j=AN->p[i];j<AN->p[i+1];j++){
            data->A->i[nz_A] = A_row + Aij->m + AN->i[j];
            data->A->x[nz_A++] = AN->x[j];
        }
    }

    A_row += AN->m;
    A_col += AN->n;

    data->A->p[A_col] = nz_A;


}


c_int setup_AP_matrices_initial_N(OSQPData * data, c_int Nmax, c_int N, c_int nx, c_int nu, c_int ny, c_int nt, 
				  csc * Q0, csc *Qi, csc *QN, 
				  csc * A0, csc *Ai, csc *Aij, csc *AN ){
    c_int nz=0;
    c_int nvar = Nmax*(nx+nu);
    c_int nconst = Nmax*(nx+ny)+nt;

    c_int P_nnz_max = Nmax*Qi->p[Qi->n] + Q0->p[Q0->n] + QN->p[QN->n]; // Maximum number of nonzeros in cost matrix
    c_int A_nnz_max = Nmax*(Ai->p[Ai->n] + Aij->p[Aij->n]) + A0->p[A0->n] + AN->p[AN->n]; // Maximum number of nonzeros in cost matrix

    data->P = csc_spalloc(data->n,data->n, P_nnz_max, 1, 0); // Constraint matrix    
    data->A = csc_spalloc(data->m, data->n, A_nnz_max, 1, 0); // Constraint matrix
    data->P->p[data->P->n] = P_nnz_max;
    data->A->p[data->A->n] = A_nnz_max;
    data->P->n = nvar;
    data->P->m = nvar;
    data->A->n = nvar;
    data->A->m = nconst;

    c_int nz_P = 0, nz_A = 0;
    c_int P_row=0;
    c_int A_row=0, A_col=0;
    for (c_int i=0; i<Q0->n;i++){ // For column 0->nu
        data->P->p[i] = nz_P;
        for (c_int j=Q0->p[i];j<Q0->p[i+1];j++){
            data->P->i[nz_P] = Q0->i[j];
            data->P->x[nz_P++] = Q0->x[j];
        }
        data->A->p[i] = nz_A;
        for (c_int j=A0->p[i];j<A0->p[i+1];j++){
            data->A->i[nz_A] = A0->i[j];
            data->A->x[nz_A++] = A0->x[j];
        }

    }
    P_row += Q0->m;
    A_col += A0->n;

    for (c_int iter=0;iter<N-2;iter++){
        for (c_int i=0; i<Qi->n;i++){
            data->P->p[P_row+i] = nz_P;
            for (c_int j=Qi->p[i];j<Qi->p[i+1];j++){
                data->P->i[nz_P] = P_row + Qi->i[j];
                data->P->x[nz_P++] = Qi->x[j];
            }
        }
        P_row += Qi->m;

        for (c_int i=0; i<Aij->n;i++){
            data->A->p[A_col+i] = nz_A;
            for (c_int j=Aij->p[i];j<Aij->p[i+1];j++){
                data->A->i[nz_A] = A_row + Aij->i[j];
                data->A->x[nz_A++] = Aij->x[j];
            }
            for (c_int j=Ai->p[i];j<Ai->p[i+1];j++){
                data->A->i[nz_A] = A_row + Aij->m + Ai->i[j];
                data->A->x[nz_A++] = Ai->x[j];
            }
        }
        A_row += Aij->m;
        A_col += Aij->n;

    } // End iter
    // Add last Q (for iter=N-2)
    for (c_int i=0; i<Qi->n;i++){
	data->P->p[P_row+i] = nz_P;
	for (c_int j=Qi->p[i];j<Qi->p[i+1];j++){
	    data->P->i[nz_P] = P_row + Qi->i[j];
	    data->P->x[nz_P++] = Qi->x[j];
	}
    }
    P_row += Qi->m;
    // Add last A (for iter=N-2, excluding Ai)
    for (c_int i=0; i<Aij->n;i++){
	data->A->p[A_col+i] = nz_A;
	for (c_int j=Aij->p[i];j<Aij->p[i+1];j++){
	    data->A->i[nz_A] = A_row + Aij->i[j];
	    data->A->x[nz_A++] = Aij->x[j];
	}
    }
    A_row += Aij->m;
    A_col += Aij->n;

    

    data->P->p[P_row] = nz_P;
    data->A->p[A_col] = nz_A;
    return 0;

}


c_int setup_AP_matrices(OSQPData * data, c_int Nmax, c_int N, c_int nx, c_int nu, c_int ny, c_int nt, 
                        csc * Q0, csc *Qi, csc *QN, 
                        csc * A0, csc *Ai, csc *Aij, csc *AN ){
    c_int nz=0;
    c_int nvar = Nmax*(nx+nu);
    c_int nconst = Nmax*(nx+ny)+nt;

    c_int P_nnz_max = Nmax*Qi->p[Qi->n] + Q0->p[Q0->n] + QN->p[QN->n]; // Maximum number of nonzeros in cost matrix
    c_int A_nnz_max = Nmax*(Ai->p[Ai->n] + Aij->p[Aij->n]) + A0->p[A0->n] + AN->p[AN->n]; // Maximum number of nonzeros in cost matrix
    //if (data->P) csc_spfree(data->P);
    //if (data->A) csc_spfree(data->A);
    
    data->P = csc_spalloc(data->n,data->n, P_nnz_max, 1, 0); // Constraint matrix    
    data->A = csc_spalloc(data->m, data->n, A_nnz_max, 1, 0); // Constraint matrix
    data->P->p[data->P->n] = P_nnz_max;
    data->A->p[data->A->n] = A_nnz_max;
    data->P->n = nvar;
    data->P->m = nvar;
    data->A->n = nvar;
    data->A->m = nconst;

    c_int nz_P = 0, nz_A = 0;
    c_int P_row=0;
    c_int A_row=0, A_col=0;

    for (c_int i=0; i<Q0->n;i++){ // For column 0->nu
        data->P->p[i] = nz_P;
        for (c_int j=Q0->p[i];j<Q0->p[i+1];j++){
            data->P->i[nz_P] = Q0->i[j];
            data->P->x[nz_P++] = Q0->x[j];
        }
        data->A->p[i] = nz_A;
        for (c_int j=A0->p[i];j<A0->p[i+1];j++){
            data->A->i[nz_A] = A0->i[j];
            data->A->x[nz_A++] = A0->x[j];
        }

    }
    P_row += Q0->m;
    A_col += A0->n;

    for (c_int iter=0;iter<N-1;iter++){
        for (c_int i=0; i<Qi->n;i++){
            data->P->p[P_row+i] = nz_P;
            for (c_int j=Qi->p[i];j<Qi->p[i+1];j++){
                data->P->i[nz_P] = P_row + Qi->i[j];
                data->P->x[nz_P++] = Qi->x[j];
            }
        }
        P_row += Qi->m;

        for (c_int i=0; i<Aij->n;i++){
            data->A->p[A_col+i] = nz_A;
            for (c_int j=Aij->p[i];j<Aij->p[i+1];j++){
                data->A->i[nz_A] = A_row + Aij->i[j];
                data->A->x[nz_A++] = Aij->x[j];
            }
            for (c_int j=Ai->p[i];j<Ai->p[i+1];j++){
                data->A->i[nz_A] = A_row + Aij->m + Ai->i[j];
                data->A->x[nz_A++] = Ai->x[j];
            }
        }
        A_row += Aij->m;
        A_col += Aij->n;

    } // End iter
    
    for (c_int i=0; i<QN->n;i++){
        data->P->p[P_row+i] = nz_P;
        for (c_int j=QN->p[i];j<QN->p[i+1];j++){
            data->P->i[nz_P] = P_row + QN->i[j];
            data->P->x[nz_P++] = QN->x[j];
        }
    }
    P_row += QN->m;
    data->P->p[P_row] = nz_P;

    // Add the final AN->n columns
    for (c_int i=0; i<AN->n;i++){
        data->A->p[A_col+i] = nz_A;
        for (c_int j=Aij->p[i];j<Aij->p[i+1];j++){
            data->A->i[nz_A] = A_row + Aij->i[j];
            data->A->x[nz_A++] = Aij->x[j];
        }
        for (c_int j=AN->p[i];j<AN->p[i+1];j++){
            data->A->i[nz_A] = A_row + Aij->m + AN->i[j];
            data->A->x[nz_A++] = AN->x[j];
        }
    }

    A_row += AN->m;
    A_col += AN->n;


    data->A->p[A_col] = nz_A;
    return 0;

}


c_int osqp_update_recursive(OSQPWorkspace* work, OSQPDataRLDL *data_rldl, c_int N) {
    work->data->n = N*(data_rldl->nx+data_rldl->nu);
    work->data->m = N*(data_rldl->nx+data_rldl->ny)+data_rldl->nt;

    c_int deltaN = N - data_rldl->N; 
    if((N>data_rldl->Nmax)){
        printf("ERROR: the new horizon is too large: %i > %i\n", N, data_rldl->Nmax);
        return -1;
    }
    if (deltaN == 0){
        printf("deltaN = 0, returning\n");
        return 0;
    }
    if(N<1){
        printf("ERROR: the new horizon is too small: %i < 1\n", N);
        return -1;
    }

    update_problem_size(work->linsys_solver,
			N*(data_rldl->nx+data_rldl->nu),
			N*(data_rldl->nx+data_rldl->ny)+data_rldl->nt);

    for (int i=data_rldl->N;i<N;i++) {
        data_rldl->Amats[i] = data_rldl->Ai;
        data_rldl->Qmats[i] = data_rldl->Qi;
    }

    data_rldl->Amats[N] = data_rldl->AN;
    data_rldl->Qmats[N] = data_rldl->QN;
    
    c_int result = LDL_update_from_pivot_solver(work->linsys_solver, data_rldl->X_even,
						data_rldl->Qmats,
						data_rldl->Amats, 
						work->settings->sigma, work->rho_inv_vec,
						data_rldl->nx,   data_rldl->nu, data_rldl->ny,
						data_rldl->nt,
						deltaN, N, data_rldl->Nmax);

    update_AP_matrices( data_rldl, work->data, data_rldl->N, N);
    data_rldl->N = N;

    return 0;

}

c_int osqp_setup_recursive(OSQPWorkspace** workp, OSQPDataRLDL *data_rldl, const OSQPSettings *settings, c_int Nmax, c_int N, c_int nx, c_int nu, c_int ny, c_int nt) {
    c_int exitflag;
    c_float epsilon = 1e-7;
    c_int nz=0,count = 0;

    OSQPWorkspace * work;
    OSQPData * data = &data_rldl->data;
    data_rldl->Nplus = 0;
    data_rldl->N = N;
    data_rldl->Nx = 0;    
    // Allocate empty workspace
    work = c_calloc(1, sizeof(OSQPWorkspace));
    if (!(work)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    *workp = work;
    // Copy problem data into workspace
    work->data = c_malloc(sizeof(OSQPData));
    if (!(work->data)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    work->data->n = data->n;
    work->data->m = data->m;


    c_int n = Nmax*(nx+nu);     // number of variables
    c_int m = Nmax*(nx+ny)+nt;  // number of constraints

    // Count the number of non-zeros
    csc* Q0 = data_rldl->Q0;
    csc* Qi = data_rldl->Qi;
    csc* QN = data_rldl->QN;

    
    csc* A0  = data_rldl->A0;
    csc* Ai  = data_rldl->Ai;
    csc* Aij = data_rldl->Aij;
    csc* AN  = data_rldl->AN;

    c_int P_nnz_max = Nmax*Qi->p[Qi->n] + Q0->p[Q0->n] + QN->p[QN->n]; // Maximum number of nonzeros in cost matrix
    c_int A_nnz_max = Nmax*(Ai->p[Ai->n] + Aij->p[Aij->n]) + A0->p[A0->n] + AN->p[AN->n]; // Maximum number of nonzeros in cost matrix
    c_int sum_Pnz = A0->p[A0->n] + AN->p[AN->n]  + Q0->p[Q0->n] + QN->p[QN->n]  + Nmax*(Ai->p[Ai->n] + Aij->p[Aij->n] + Qi->p[Qi->n]) ; // Total number of nonzeros

    // Compute the cost matrix P and the constraint matrix A
    // using the system matrices A0, Ai, Aij, AN, Q0, Qi, QN
    setup_AP_matrices(work->data, Nmax, N, nx, nu, ny, nt, 
		      data_rldl->Q0, data_rldl->Qi, data_rldl->QN, 
		      data_rldl->A0, data_rldl->Ai, data_rldl->Aij, data_rldl->AN );
    


    // Validate settings
    if (validate_settings(settings)) return osqp_error(OSQP_SETTINGS_VALIDATION_ERROR);
    // Copy settings
    work->settings = copy_settings(settings);
    if (!(work->settings)) return osqp_error(OSQP_MEM_ALLOC_ERROR);


    // Vectorized rho parameter
    work->rho_vec     = c_malloc(data->m * sizeof(c_float));
    work->rho_inv_vec = c_malloc(data->m * sizeof(c_float));
    if ( data->m && (!(work->rho_vec) || !(work->rho_inv_vec)) )
	return osqp_error(OSQP_MEM_ALLOC_ERROR);


    // Validate data: TODO: Add for recursive horizon
    //if (validate_data(data)) return osqp_error(OSQP_DATA_VALIDATION_ERROR);




    // Start and allocate directly timer
# ifdef PROFILING
    work->timer = c_malloc(sizeof(OSQPTimer));
    if (!(work->timer)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    osqp_tic(work->timer);
# endif /* ifdef PROFILING */

    // Cost
    //work->data->P = copy_csc_mat(data->P);
    work->data->q = vec_copy(data->q, data->n);

    if (!(work->data->P) || !(work->data->q)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    if (!(work->data->q)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    
    // Constraints
    //work->data->A = copy_csc_mat(data->A);
    if (!(work->data->A)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    work->data->l = vec_copy(data->l, data->m);
    work->data->u = vec_copy(data->u, data->m);

    if ( data->m && (!(work->data->l) || !(work->data->u)) )
	return osqp_error(OSQP_MEM_ALLOC_ERROR);
    

    work->data->P->n = N*(nx+nu);
    work->data->P->m = N*(nx+nu);
    work->data->A->n = N*(nx+nu);
    work->data->A->m = N*(nx+ny)+nt;


    // Type of constraints
    work->constr_type = c_calloc(data->m, sizeof(c_int));
    if (data->m && !(work->constr_type)) return osqp_error(OSQP_MEM_ALLOC_ERROR);

    // Allocate internal solver variables (ADMM steps)
    work->x        = c_calloc(data->n, sizeof(c_float));
    work->z        = c_calloc(data->m, sizeof(c_float));
    work->xz_tilde = c_calloc(data->n + data->m, sizeof(c_float));
    work->x_prev   = c_calloc(data->n, sizeof(c_float));
    work->z_prev   = c_calloc(data->m, sizeof(c_float));
    work->y        = c_calloc(data->m, sizeof(c_float));
    if (!(work->x) || !(work->xz_tilde) || !(work->x_prev))
	return osqp_error(OSQP_MEM_ALLOC_ERROR);
    if ( data->m && (!(work->z) || !(work->z_prev) || !(work->y)) )
	return osqp_error(OSQP_MEM_ALLOC_ERROR);

    // Initialize variables x, y, z to 0
    cold_start(work);

    // Primal and dual residuals variables
    work->Ax  = c_calloc(data->m, sizeof(c_float));
    work->Px  = c_calloc(data->n, sizeof(c_float));
    work->Aty = c_calloc(data->n, sizeof(c_float));

    // Primal infeasibility variables
    work->delta_y   = c_calloc(data->m, sizeof(c_float));
    work->Atdelta_y = c_calloc(data->n, sizeof(c_float));

    // Dual infeasibility variables
    work->delta_x  = c_calloc(data->n, sizeof(c_float));
    work->Pdelta_x = c_calloc(data->n, sizeof(c_float));
    work->Adelta_x = c_calloc(data->m, sizeof(c_float));

    if (!(work->Px) || !(work->Aty) || !(work->Atdelta_y) || !(work->delta_x) || !(work->Pdelta_x))
	return osqp_error(OSQP_MEM_ALLOC_ERROR);
    if ( data->m && (!(work->Ax) || !(work->delta_y) || !(work->Adelta_x)) )
	return osqp_error(OSQP_MEM_ALLOC_ERROR);



    
    // Perform scaling
    if (settings->scaling) {
	// Allocate scaling structure
	work->scaling = c_malloc(sizeof(OSQPScaling));
	if (!(work->scaling)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
	work->scaling->D    = c_malloc(data->n * sizeof(c_float));
	work->scaling->Dinv = c_malloc(data->n * sizeof(c_float));
	work->scaling->E    = c_malloc(data->m * sizeof(c_float));
	work->scaling->Einv = c_malloc(data->m * sizeof(c_float));
	if (!(work->scaling->D) || !(work->scaling->Dinv))
	    return osqp_error(OSQP_MEM_ALLOC_ERROR);
	if ( data->m && (!(work->scaling->E) || !(work->scaling->Einv)) )
	    return osqp_error(OSQP_MEM_ALLOC_ERROR);


	// Allocate workspace variables used in scaling
	work->D_temp   = c_malloc(data->n * sizeof(c_float));
	work->D_temp_A = c_malloc(data->n * sizeof(c_float));
	work->E_temp   = c_malloc(data->m * sizeof(c_float));
	// if (!(work->D_temp) || !(work->D_temp_A) || !(work->E_temp))
	//   return osqp_error(OSQP_MEM_ALLOC_ERROR);
	if (!(work->D_temp) || !(work->D_temp_A)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
	if (data->m && !(work->E_temp))           return osqp_error(OSQP_MEM_ALLOC_ERROR);

	// Scale data
	scale_data(work);
    } else {
	work->scaling  = OSQP_NULL;
	work->D_temp   = OSQP_NULL;
	work->D_temp_A = OSQP_NULL;
	work->E_temp   = OSQP_NULL;
    }

    // Set type of constraints
    set_rho_vec(work);

    // Load linear system solver
    if (load_linsys_solver(work->settings->linsys_solver)) return osqp_error(OSQP_LINSYS_SOLVER_LOAD_ERROR);

    data_rldl->Ytemp = OSQP_NULL;
    data_rldl->Uhat = OSQP_NULL;
    data_rldl->Vhat = OSQP_NULL;
    data_rldl->Yhat = OSQP_NULL;
    data_rldl->Lz = OSQP_NULL;
    data_rldl->Pz = OSQP_NULL;
    data_rldl->Dzinv = OSQP_NULL;
    data_rldl->Ly = OSQP_NULL;
    data_rldl->Dyinv = OSQP_NULL;
    data_rldl->P0 = OSQP_NULL;

    data_rldl->Qmats = malloc((Nmax+1)*sizeof(csc*));
    data_rldl->Amats = malloc((Nmax+1)*sizeof(csc*));
    
    data_rldl->Qmats[0] = Q0;
    data_rldl->Amats[0] = A0;
    for (int i=1;i<N;i++) {
        data_rldl->Qmats[i] = Qi;
        data_rldl->Amats[i] = Ai;
    }
    data_rldl->Qmats[N] = QN;
    data_rldl->Amats[N] = AN;

    // SET SIZES OF ALL MATRICES =======================================================

    data_rldl->X_even = c_malloc(Nmax*sizeof(csc*));
    for (int i = 0; i<Nmax;i++){
	data_rldl->X_even[i] = OSQP_NULL;
    }

    exitflag = init_linsys_solver_qdldl_recursive(&(work->linsys_solver), n, m, P_nnz_max*5, A_nnz_max, work->settings->sigma, work->rho_vec, work->pol);
    
    exitflag =  LDL_factorize_recursive_solver(work->linsys_solver, data_rldl->X_even,
					       data_rldl->Amats, data_rldl->Qmats,  
					       work->settings->sigma, work->rho_inv_vec,
					       nx,   nu, ny, nt,
					       Nmax, N, data_rldl->nnz);

    work->data->n = N*(data_rldl->nx+data_rldl->nu);
    work->data->m = N*(data_rldl->nx+data_rldl->ny)+data_rldl->nt;
    
    update_problem_size(work->linsys_solver,
			N*(data_rldl->nx+data_rldl->nu),
			N*(data_rldl->nx+data_rldl->ny)+data_rldl->nt);    
    
    if (exitflag) return osqp_error(exitflag);




    //compute_KKT_permutations(work->linsys_solver, N*(nx+nu), N*(nx+ny)+nt, P_nnz_max, A_nnz_max, 
    //                                     N, Q0, Qi, QN,  A0,  Ai, Aij, AN);

    // Initialize active constraints structure
    work->pol = c_malloc(sizeof(OSQPPolish));
    if (!(work->pol)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    work->pol->Alow_to_A = c_malloc(data->m * sizeof(c_int));
    work->pol->Aupp_to_A = c_malloc(data->m * sizeof(c_int));
    work->pol->A_to_Alow = c_malloc(data->m * sizeof(c_int));
    work->pol->A_to_Aupp = c_malloc(data->m * sizeof(c_int));
    work->pol->x         = c_malloc(data->n * sizeof(c_float));
    work->pol->z         = c_malloc(data->m * sizeof(c_float));
    work->pol->y         = c_malloc(data->m * sizeof(c_float));
    if (!(work->pol->x)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    if ( data->m && (!(work->pol->Alow_to_A) || !(work->pol->Aupp_to_A) ||
		     !(work->pol->A_to_Alow) || !(work->pol->A_to_Aupp) ||
		     !(work->pol->z) || !(work->pol->y)) )
	return osqp_error(OSQP_MEM_ALLOC_ERROR);

    // Allocate solution
    work->solution = c_calloc(1, sizeof(OSQPSolution));
    if (!(work->solution)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    work->solution->x = c_calloc(1, data->n * sizeof(c_float));
    work->solution->y = c_calloc(1, data->m * sizeof(c_float));
    if (!(work->solution->x))            return osqp_error(OSQP_MEM_ALLOC_ERROR);
    if (data->m && !(work->solution->y)) return osqp_error(OSQP_MEM_ALLOC_ERROR);

    // Allocate and initialize information
    work->info = c_calloc(1, sizeof(OSQPInfo));
    if (!(work->info)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    work->info->status_polish = 0;              // Polishing not performed
    update_status(work->info, OSQP_UNSOLVED);
# ifdef PROFILING
    work->info->solve_time  = 0.0;                   // Solve time to zero
    work->info->update_time = 0.0;                   // Update time to zero
    work->info->polish_time = 0.0;                   // Polish time to zero
    work->info->run_time    = 0.0;                   // Total run time to zero
    work->info->setup_time  = osqp_toc(work->timer); // Update timer information

    work->first_run         = 1;
    work->clear_update_time = 0;
    work->rho_update_from_solve = 0;
# endif /* ifdef PROFILING */
    work->info->rho_updates  = 0;                    // Rho updates set to 0
    work->info->rho_estimate = work->settings->rho;  // Best rho estimate
    // Print header

# ifdef PRINTING


    work->summary_printed = 0; // Initialize last summary  to not printed
# endif /* ifdef PRINTING */

    // If adaptive rho and automatic interval, but profiling disabled, we need to
    // set the interval to a default value
# ifndef PROFILING
    if (work->settings->adaptive_rho && !work->settings->adaptive_rho_interval) {
	if (work->settings->check_termination) {
	    // If check_termination is enabled, we set it to a multiple of the check
	    // termination interval
	    work->settings->adaptive_rho_interval = ADAPTIVE_RHO_MULTIPLE_TERMINATION *
		work->settings->check_termination;
	} else {
	    // If check_termination is disabled we set it to a predefined fix number
	    work->settings->adaptive_rho_interval = ADAPTIVE_RHO_FIXED;
	}
    }
# endif /* ifndef PROFILING */

    // Return exit flag
    return 0;
}


c_int osqp_allocate_combine_recursive(OSQPWorkspace** workp, OSQPDataRLDL *data_rldl,
				      const OSQPSettings *settings,
				      c_int Nmax, 
				      c_int nx, c_int nu, c_int ny, c_int nt) {
    OSQPWorkspace * work;
    
    // Allocate empty workspace
    work = c_calloc(1, sizeof(OSQPWorkspace));
    if (!(work)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    *workp = work;
    
    // Allocate data struct
    work->data = c_malloc(sizeof(OSQPData));
    if (!(work->data)) return osqp_error(OSQP_MEM_ALLOC_ERROR);

    c_int n_var = Nmax*(nx+nu);        // number of variables
    c_int n_constr = Nmax*(nx+ny)+nt;  // number of constraints
    c_int nx_ny = nx+ny;
    c_int nx_nu = nx+nu;


    // Validate settings
    if (validate_settings(settings)) return osqp_error(OSQP_SETTINGS_VALIDATION_ERROR);
    // Copy settings
    work->settings = copy_settings(settings);
    if (!(work->settings)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    
    // Maximum number of nonzeros in cost matrix
    c_int P_nnz_max = Nmax*data_rldl->Qi->p[data_rldl->Qi->n] +
	data_rldl->Q0->p[data_rldl->Q0->n] +
	data_rldl->QN->p[data_rldl->QN->n];
    // Maximum number of nonzeros in cost matrix
    c_int A_nnz_max = Nmax*(data_rldl->Ai->p[data_rldl->Ai->n] +
			    data_rldl->Aij->p[data_rldl->Aij->n]) +
	data_rldl->A0->p[data_rldl->A0->n] +
	data_rldl->AN->p[data_rldl->AN->n]; 

}


c_int osqp_setup_combine_recursive(OSQPWorkspace** workp, OSQPDataRLDL *data_rldl,
				   const OSQPSettings *settings,
				   c_int Nmax, c_int N, c_int Nplus,
				   c_int nx, c_int nu, c_int ny, c_int nt) {
    // Setup A, P matrices internally, and allow space to grow. 
    c_int exitflag;

    
    data_rldl->N = N+Nplus; // The current total horizon
    data_rldl->Nx = N;
    data_rldl->Nplus = Nplus;

    OSQPWorkspace * work;
    OSQPData * data = &data_rldl->data;

    // Allocate empty workspace
    work = c_calloc(1, sizeof(OSQPWorkspace));
    if (!(work)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    *workp = work;
    // Copy problem data into workspace
    work->data = c_malloc(sizeof(OSQPData));
    if (!(work->data)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    work->data->n = data->n;
    work->data->m = data->m;

    c_int n = Nmax*(nx+nu);     // number of variables
    c_int m = Nmax*(nx+ny)+nt;  // number of constraints
    c_int nx_ny = nx+ny;
    c_int nx_nu = nx+nu;
    // Count the number of non-zeros
    csc* Q0 = data_rldl->Q0;
    csc* Qi = data_rldl->Qi;
    csc* QN = data_rldl->QN;

    
    csc* A0  = data_rldl->A0;
    csc* Ai  = data_rldl->Ai;
    csc* Aij = data_rldl->Aij;
    csc* AN  = data_rldl->AN;

    c_int P_nnz_max = Nmax*Qi->p[Qi->n] + Q0->p[Q0->n] + QN->p[QN->n]; // Maximum number of nonzeros in cost matrix
    c_int A_nnz_max = Nmax*(Ai->p[Ai->n] + Aij->p[Aij->n]) + A0->p[A0->n] + AN->p[AN->n]; // Maximum number of nonzeros in cost matrix

    // Compute the cost matrix P and the constraint matrix A
    // using the system matrices A0, Ai, Aij, AN, Q0, Qi, QN
    setup_AP_matrices_initial_N(work->data, Nmax, N, nx, nu, ny, nt, 
				data_rldl->Q0, data_rldl->Qi, data_rldl->QN, 
				data_rldl->A0, data_rldl->Ai, data_rldl->Aij, data_rldl->AN );

    work->data->P->n =  (N-1)*(nx+nu)+nu;     // number of variables
    work->data->P->m =  (N-1)*(nx+nu)+nu;     // number of variables
    work->data->A->n =  (N-1)*(nx+nu)+nu;     // number of variables
    work->data->A->m =  (N-1)*(nx+ny);     // number of variables

    // Validate settings
    if (validate_settings(settings)) return osqp_error(OSQP_SETTINGS_VALIDATION_ERROR);
    // Copy settings
    work->settings = copy_settings(settings);
    if (!(work->settings)) return osqp_error(OSQP_MEM_ALLOC_ERROR);


    // Vectorized rho parameter
    work->rho_vec     = c_malloc(data->m * sizeof(c_float));
    work->rho_inv_vec = c_malloc(data->m * sizeof(c_float));
    if ( data->m && (!(work->rho_vec) || !(work->rho_inv_vec)) )
	return osqp_error(OSQP_MEM_ALLOC_ERROR);

    //data->P = copy_csc_mat(work->data->P);
    //data->A = copy_csc_mat(work->data->A);

    



    // Start and allocate directly timer
# ifdef PROFILING
    work->timer = c_malloc(sizeof(OSQPTimer));
    if (!(work->timer)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    osqp_tic(work->timer);
# endif /* ifdef PROFILING */

    // Cost    
    //work->data->P = copy_csc_mat(data->P);
    work->data->q = vec_copy(data->q, data->n);

    if (!(work->data->P) || !(work->data->q)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    if (!(work->data->q)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    
    // Constraints
    //work->data->A = copy_csc_mat(data->A);
    if (!(work->data->A)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    work->data->l = vec_copy(data->l, data->m);
    work->data->u = vec_copy(data->u, data->m);

    if ( data->m && (!(work->data->l) || !(work->data->u)) )
	return osqp_error(OSQP_MEM_ALLOC_ERROR);
    
    /*
      work->data->n = (N-1)*(nx+nu)+nu;
      work->data->m = (N-1)*(nx+ny);
      if (validate_data(work->data)) return osqp_error(OSQP_DATA_VALIDATION_ERROR);
      work->data->n = n;
      work->data->m = m;
    */
    // Type of constraints
    work->constr_type = c_calloc(data->m, sizeof(c_int));
    if (data->m && !(work->constr_type)) return osqp_error(OSQP_MEM_ALLOC_ERROR);

    // Allocate internal solver variables (ADMM steps)
    work->x        = c_calloc(data->n, sizeof(c_float));
    work->z        = c_calloc(data->m, sizeof(c_float));
    work->xz_tilde = c_calloc(data->n + data->m, sizeof(c_float));
    work->x_prev   = c_calloc(data->n, sizeof(c_float));
    work->z_prev   = c_calloc(data->m, sizeof(c_float));
    work->y        = c_calloc(data->m, sizeof(c_float));
    if (!(work->x) || !(work->xz_tilde) || !(work->x_prev))
	return osqp_error(OSQP_MEM_ALLOC_ERROR);
    if ( data->m && (!(work->z) || !(work->z_prev) || !(work->y)) )
	return osqp_error(OSQP_MEM_ALLOC_ERROR);

    // Initialize variables x, y, z to 0
    cold_start(work);

    // Primal and dual residuals variables
    work->Ax  = c_calloc(data->m, sizeof(c_float));
    work->Px  = c_calloc(data->n, sizeof(c_float));
    work->Aty = c_calloc(data->n, sizeof(c_float));

    // Primal infeasibility variables
    work->delta_y   = c_calloc(data->m, sizeof(c_float));
    work->Atdelta_y = c_calloc(data->n, sizeof(c_float));

    // Dual infeasibility variables
    work->delta_x  = c_calloc(data->n, sizeof(c_float));
    work->Pdelta_x = c_calloc(data->n, sizeof(c_float));
    work->Adelta_x = c_calloc(data->m, sizeof(c_float));

    if (!(work->Px) || !(work->Aty) || !(work->Atdelta_y) || !(work->delta_x) || !(work->Pdelta_x))
	return osqp_error(OSQP_MEM_ALLOC_ERROR);
    if ( data->m && (!(work->Ax) || !(work->delta_y) || !(work->Adelta_x)) )
	return osqp_error(OSQP_MEM_ALLOC_ERROR);

    

    
    // Perform scaling
    if (settings->scaling) {
	// Allocate scaling structure
	work->scaling = c_malloc(sizeof(OSQPScaling));
	if (!(work->scaling)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
	work->scaling->D    = c_malloc(data->n * sizeof(c_float));
	work->scaling->Dinv = c_malloc(data->n * sizeof(c_float));
	work->scaling->E    = c_malloc(data->m * sizeof(c_float));
	work->scaling->Einv = c_malloc(data->m * sizeof(c_float));
	if (!(work->scaling->D) || !(work->scaling->Dinv))
	    return osqp_error(OSQP_MEM_ALLOC_ERROR);
	if ( data->m && (!(work->scaling->E) || !(work->scaling->Einv)) )
	    return osqp_error(OSQP_MEM_ALLOC_ERROR);


	// Allocate workspace variables used in scaling
	work->D_temp   = c_malloc(data->n * sizeof(c_float));
	work->D_temp_A = c_malloc(data->n * sizeof(c_float));
	work->E_temp   = c_malloc(data->m * sizeof(c_float));
	// if (!(work->D_temp) || !(work->D_temp_A) || !(work->E_temp))
	//   return osqp_error(OSQP_MEM_ALLOC_ERROR);
	if (!(work->D_temp) || !(work->D_temp_A)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
	if (data->m && !(work->E_temp))           return osqp_error(OSQP_MEM_ALLOC_ERROR);

	// Scale data
	scale_data(work);
    } else {
	work->scaling  = OSQP_NULL;
	work->D_temp   = OSQP_NULL;
	work->D_temp_A = OSQP_NULL;
	work->E_temp   = OSQP_NULL;
    }

    // Set type of constraints
    set_rho_vec(work);   

    
    // Compute the P0 permutation vector
    data_rldl->P0 = (c_int *)c_malloc(sizeof(c_int)*(n+m));
    c_int P_count = 0;
    c_int cost_count = 0;
    c_int constr_count = data_rldl->N*(nx_nu);
    for (int j=0; j<nu;j++) data_rldl->P0[P_count++] = cost_count++;
    for (int i=0; i<data_rldl->N-1;i++){
	for (int j=0; j<nx_ny;j++) data_rldl->P0[P_count++] = constr_count++;
	for (int j=0; j<nx_nu;j++) data_rldl->P0[P_count++] = cost_count++;	
    }
    for (int j=0; j<nx_ny;j++) data_rldl->P0[P_count++] = constr_count++;
    for (int j=0; j<nx;j++)    data_rldl->P0[P_count++] = cost_count++;
    for (int j=0; j<nt;j++)    data_rldl->P0[P_count++] = constr_count++;
    


    // Load linear system solver
    if (load_linsys_solver(work->settings->linsys_solver)) return osqp_error(OSQP_LINSYS_SOLVER_LOAD_ERROR);

    
    // Initialize the linsys solver and compute the initial factorization P'XP = LDL'
    // Later compute Uhat and add below L
    exitflag = init_linsys_solver_qdldl_combined(&(work->linsys_solver),
						 work->data->P, work->data->A,
						 work->settings->sigma, work->rho_vec, work->pol, 
                                                 (nx+nu)*(N-1)+nu, (nx+ny)*(N-1),
						 (nx+nu)*Nmax, (nx+ny)*Nmax+nt,
						 P_nnz_max, A_nnz_max);
    compute_Px(work->linsys_solver, data_rldl);    

    // SET SIZES OF ALL MATRICES =======================================================
    data_rldl->X_even = c_malloc(Nmax*sizeof(csc*));
    for (int i = 0; i<Nmax;i++){
	data_rldl->X_even[i] = OSQP_NULL;
    }

    
    data_rldl->Qmats = malloc((Nmax+1)*sizeof(csc*));
    data_rldl->Amats = malloc((Nmax+1)*sizeof(csc*));
    
    data_rldl->Qmats[0] = Q0;
    data_rldl->Amats[0] = A0;
    for (int i=1;i<Nmax;i++) {
        data_rldl->Qmats[i] = Qi;
        data_rldl->Amats[i] = Ai;
    }
    data_rldl->Qmats[Nmax] = QN;
    data_rldl->Amats[Nmax] = AN;
       

    
    c_int increase_max_var = (nx+ny+nx+nu)*(Nmax-N)+nt;
    // Pz'*Z*Pz = Lz*Dz*Lz'
    data_rldl->Lz = csc_spalloc(increase_max_var,
				increase_max_var,
				(nx+ny+nx+nu)*(nx+ny+nx+nu)*(Nmax-N), 1, 0);
    
    data_rldl->Dzinv = (c_float *)c_malloc(sizeof(c_float)*increase_max_var);
    data_rldl->Pz = (c_int *)c_malloc(sizeof(c_int)*(increase_max_var));
    
    // Py'*Y*Py = Ly*Dy*Ly'
    
    data_rldl->Yhat = csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);

    
    data_rldl->Ytemp = (c_float *)c_calloc(nx_ny*nx_ny, sizeof(c_float));
    for (int i=0;i<nx_ny;i++){
	data_rldl->Ytemp[i+nx_ny*i] = -work->rho_inv_vec[(N-1)*nx_ny+i]; // TODO: should be N-1?
    }

    // Compute Uhat = U'P*Linv'*Dinv    
    data_rldl->Uhat = c_malloc(sizeof(csc));
    exitflag = compute_Uhat_solver(work->linsys_solver, Ai, data_rldl->Uhat, data_rldl->Ytemp);

    // Compute the Z factorization
    exitflag =  LDL_factorize_recursive(data_rldl->Lz, data_rldl->Dzinv, data_rldl->Pz, 
                                        data_rldl->X_even,
					&data_rldl->Qmats[data_rldl->Nmax - Nplus],
					&data_rldl->Amats[data_rldl->Nmax - Nplus],
					work->settings->sigma, work->rho_inv_vec,
					nx,   nu, ny, nt,
				        Nmax, data_rldl->Nx, Nplus, data_rldl->nnz);



    int ncol_z_max = (data_rldl->Nmax - data_rldl->Nx)*(nx_ny+nx_nu)+data_rldl->nt;
    
    // Compute the Vhat = V*Q*Minv'*Einv
    data_rldl->Vhat = csc_spalloc(nx_ny, ncol_z_max, nx_ny*ncol_z_max, 1, 0);
    compute_Vhat(data_rldl->Lz, data_rldl->Dzinv, data_rldl->Pz, Aij, data_rldl->Vhat,
		 data_rldl->Yhat, data_rldl->Ytemp);

    c_int factor_status;
    c_float D[nx_ny];
    c_int etree[nx_ny];
    c_int Lnz[nx_ny];
    c_int  iwork[nx_ny*3];
    c_int  bwork[nx_ny];
    c_float  fwork[nx_ny];
    c_int sum_Lnz;

    // LDL factorization
    sum_Lnz = QDLDL_etree(data_rldl->Yhat->n,
			  data_rldl->Yhat->p,
			  data_rldl->Yhat->i, iwork, Lnz, etree);

    if (sum_Lnz == -1){
	printf("Error in etree\n");
	return osqp_error(OSQP_DATA_VALIDATION_ERROR);
    }
    data_rldl->Ly = csc_spalloc(data_rldl->Yhat->n, data_rldl->Yhat->n, nx_ny*nx_ny, 1, 0);
    data_rldl->Dyinv = (c_float *)c_malloc(sizeof(c_float)*nx_ny);
    c_float * Dy = (c_float *)c_malloc(sizeof(c_float)*nx_ny);    

    factor_status = QDLDL_factor(data_rldl->Yhat->n,
				 data_rldl->Yhat->p,
				 data_rldl->Yhat->i,
				 data_rldl->Yhat->x,
				 data_rldl->Ly->p, data_rldl->Ly->i, data_rldl->Ly->x,
				 Dy, data_rldl->Dyinv, Lnz,
				 etree, bwork, iwork, fwork);



    if (factor_status == -1) {
	return osqp_error(OSQP_DATA_VALIDATION_ERROR);
    }
    

    // Fill in the rest of the L and D matrices:
    // L = [L, 0, 0; 0, M, 0; Uhat, Vhat, Ly], D = diag(D, E, Dx)  
    compute_L_matrix(work->linsys_solver, data_rldl);
    
    update_problem_size(work->linsys_solver,
			(data_rldl->N)*(data_rldl->nx+data_rldl->nu),
			(data_rldl->N)*(data_rldl->nx+data_rldl->ny)+data_rldl->nt);
    work->data->n = (data_rldl->N)*(data_rldl->nx+data_rldl->nu);
    work->data->m = (data_rldl->N)*(data_rldl->nx+data_rldl->ny)+data_rldl->nt;

    // TODO: sould be N
    update_AP_matrices( data_rldl, work->data, data_rldl->Nx-1, (data_rldl->N));

    if (exitflag) {
	return osqp_error(exitflag);
    }
    // Initialize active constraints structure
    work->pol = c_malloc(sizeof(OSQPPolish));
    if (!(work->pol)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    work->pol->Alow_to_A = c_malloc(data->m * sizeof(c_int));
    work->pol->Aupp_to_A = c_malloc(data->m * sizeof(c_int));
    work->pol->A_to_Alow = c_malloc(data->m * sizeof(c_int));
    work->pol->A_to_Aupp = c_malloc(data->m * sizeof(c_int));
    work->pol->x         = c_malloc(data->n * sizeof(c_float));
    work->pol->z         = c_malloc(data->m * sizeof(c_float));
    work->pol->y         = c_malloc(data->m * sizeof(c_float));
    if (!(work->pol->x)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    if ( data->m && (!(work->pol->Alow_to_A) || !(work->pol->Aupp_to_A) ||
		     !(work->pol->A_to_Alow) || !(work->pol->A_to_Aupp) ||
		     !(work->pol->z) || !(work->pol->y)) )
	return osqp_error(OSQP_MEM_ALLOC_ERROR);

    // Allocate solution
    work->solution = c_calloc(1, sizeof(OSQPSolution));
    if (!(work->solution)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    work->solution->x = c_calloc(1, data->n * sizeof(c_float));
    work->solution->y = c_calloc(1, data->m * sizeof(c_float));
    if (!(work->solution->x))            return osqp_error(OSQP_MEM_ALLOC_ERROR);
    if (data->m && !(work->solution->y)) return osqp_error(OSQP_MEM_ALLOC_ERROR);

    // Allocate and initialize information
    work->info = c_calloc(1, sizeof(OSQPInfo));
    if (!(work->info)) return osqp_error(OSQP_MEM_ALLOC_ERROR);
    work->info->status_polish = 0;              // Polishing not performed
    update_status(work->info, OSQP_UNSOLVED);
# ifdef PROFILING
    work->info->solve_time  = 0.0;                   // Solve time to zero
    work->info->update_time = 0.0;                   // Update time to zero
    work->info->polish_time = 0.0;                   // Polish time to zero
    work->info->run_time    = 0.0;                   // Total run time to zero
    work->info->setup_time  = osqp_toc(work->timer); // Update timer information

    work->first_run         = 1;
    work->clear_update_time = 0;
    work->rho_update_from_solve = 0;
# endif /* ifdef PROFILING */
    work->info->rho_updates  = 0;                    // Rho updates set to 0
    work->info->rho_estimate = work->settings->rho;  // Best rho estimate
    // Print header

# ifdef PRINTING


    work->summary_printed = 0; // Initialize last summary  to not printed
# endif /* ifdef PRINTING */

    // If adaptive rho and automatic interval, but profiling disabled, we need to
    // set the interval to a default value
# ifndef PROFILING
    if (work->settings->adaptive_rho && !work->settings->adaptive_rho_interval) {
	if (work->settings->check_termination) {
	    // If check_termination is enabled, we set it to a multiple of the check
	    // termination interval
	    work->settings->adaptive_rho_interval = ADAPTIVE_RHO_MULTIPLE_TERMINATION *
		work->settings->check_termination;
	} else {
	    // If check_termination is disabled we set it to a predefined fix number
	    work->settings->adaptive_rho_interval = ADAPTIVE_RHO_FIXED;
	}
    }
# endif /* ifndef PROFILING */
    if (Dy) free(Dy);

    // Return exit flag
    return 0;

}
c_int osqp_update_X_horizon(OSQPWorkspace* work, OSQPDataRLDL *data_rldl, c_int deltaN) {
    
}

c_int osqp_update_Z_horizon(OSQPWorkspace* work, OSQPDataRLDL *data_rldl, c_int deltaN) {
    c_int Nnew = data_rldl->N + deltaN;

    // Make sure that Nnew > N0
    if (Nnew < data_rldl->Nx){
	printf("ERROR: This N is too small for the recursive/combined MPC. Please setup the problem again. ");
	return -1;
    }

    // Update the total size of the problem
    c_int n_vars = Nnew*(data_rldl->nx + data_rldl->nu);
    c_int n_constr = Nnew*(data_rldl->nx + data_rldl->ny)+data_rldl->nt;

    if((Nnew > data_rldl->Nmax)){
        printf("ERROR: the new horizon is too large: %i > %i\n", Nnew, data_rldl->Nmax);
        return -1;
    }
    if (deltaN == 0)  return 0;

    // Update the solver size
    work->data->n = n_vars;
    work->data->m = n_constr;
    update_problem_size(work->linsys_solver, n_vars, n_constr);

    // Update the P0 permutation vector
    c_int nx_nu = data_rldl->nx+data_rldl->nu;
    c_int nx_ny = data_rldl->nx+data_rldl->ny;
    c_int P_count = data_rldl->nu;
    c_int cost_count = data_rldl->nu + (data_rldl->N-1)*nx_nu;
    c_int constr_count =  Nnew*(nx_nu);
    for (int i=0; i<data_rldl->N-1;i++){
	for (int j=0; j<nx_ny;j++) data_rldl->P0[P_count++] = constr_count++;
        P_count+=nx_nu;
    }
    for (int i=data_rldl->N-1; i<Nnew -1;i++){
	for (int j=0; j<nx_ny;j++) data_rldl->P0[P_count++] = constr_count++;
	for (int j=0; j<nx_nu;j++) data_rldl->P0[P_count++] = cost_count++;	
    }
    for (int j=0; j<nx_ny;j++)    data_rldl->P0[P_count++] = constr_count++;
    for (int j=0; j<data_rldl->nx;j++) data_rldl->P0[P_count++] = cost_count++;
    for (int j=0; j<data_rldl->nt;j++) data_rldl->P0[P_count++] = constr_count++;
    
    // Update the Z matrix with new horizon
    c_int result = LDL_update_from_pivot(data_rldl->Lz, data_rldl->Dzinv, data_rldl->Pz,
					 data_rldl->X_even,
					 &data_rldl->Qmats[data_rldl->Nmax-Nnew],
					 &data_rldl->Amats[data_rldl->Nmax-Nnew],
					 work->settings->sigma, work->rho_inv_vec,
					 data_rldl->nx, data_rldl->nu, data_rldl->ny,
					 data_rldl->nt,
					 deltaN, Nnew, data_rldl->Nx, data_rldl->Nmax);

    
    // Compute Vhat = V*Q*Minv'*Einv  and Yhat = V*Q*Minv'*Einv*Minv*Q'*V'
    compute_Vhat(data_rldl->Lz, data_rldl->Dzinv, data_rldl->Pz, data_rldl->Aij, data_rldl->Vhat,
		 data_rldl->Yhat, data_rldl->Ytemp);


    // LDL factorization of Yhat
    
    c_int   factor_status;
    c_int   etree[nx_ny];
    c_int   Lnz[nx_ny];
    c_int   iwork[nx_ny*3];
    c_int   bwork[nx_ny];
    c_float fwork[nx_ny];
    c_int   sum_Lnz;
    c_float Dy[nx_ny];
    sum_Lnz = QDLDL_etree(data_rldl->Yhat->n,
			  data_rldl->Yhat->p,
			  data_rldl->Yhat->i, iwork, Lnz, etree);    

    if (sum_Lnz == -1){
	printf("Error in etree\n");
	return osqp_error(OSQP_DATA_VALIDATION_ERROR);
    }
    factor_status = QDLDL_factor(data_rldl->Yhat->n,
				 data_rldl->Yhat->p,
				 data_rldl->Yhat->i,
				 data_rldl->Yhat->x,
				 data_rldl->Ly->p, data_rldl->Ly->i, data_rldl->Ly->x,
				 Dy, data_rldl->Dyinv, Lnz,
				 etree, bwork, iwork, fwork);

    if (factor_status == -1) {
	return osqp_error(OSQP_DATA_VALIDATION_ERROR);
    }

    // Fill in the rest of the L and D matrices:
    update_L_matrix(work->linsys_solver, data_rldl, deltaN);
    update_AP_matrices( data_rldl, work->data, data_rldl->Nx, Nnew);

    data_rldl->N = Nnew; 
    data_rldl->Nplus += deltaN; 
    return 0;
}


c_int osqp_solve_combine_recursive(OSQPWorkspace *work) {
    c_int exitflag = osqp_solve(work);
    return exitflag;
}




c_int osqp_solve_recursive(OSQPWorkspace *work) {
    c_int exitflag = osqp_solve(work);
    return exitflag;
}



