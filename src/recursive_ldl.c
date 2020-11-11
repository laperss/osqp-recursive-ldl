#include "recursive_ldl.h"
#include "osqp.h"
#include "auxil.h"
#include "util.h"
#include "scaling.h"
#include "error.h"
#include "qdldl_interface.h"

#include "cs_addon.h"



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

    /*
      printf("bounds = zeros(%i,%i); \n", work->data->m, 2);
      for (int i=0;i<work->data->m;i++){
      printf("bounds(%i,:) = [%9.5f,\t%9.5f];\n", i+1, work->data->m, work->data->l[i], work->data->u[i]);
      }
    */
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
    //printf("Update update_time\n");
    work->info->update_time += osqp_toc(work->timer);
#endif /* ifdef PROFILING */
    


    //printf("rhoinv1 = zeros(%2i,%i);\n", work->data->m, 1);
    //for(int i=0;i<work->data->m;i++){
    //    printf("rhoinv1(%2i,1) = %10.8f;\n", i+1, work->rho_inv_vec[i]);
    //}    
    //printf("Finished\n");
    return exitflag;
}


c_int update_problem_size(qdldl_solver * s, c_int n, c_int m){
    s->L->n = n+m;
    s->L->m = n+m;
    s->n = n;
    s->m = m;
}   






// This is for pivoting when the RHS is [0,-I; 0, 0]
// Then we can use that the final matrix has the same zero structure
static c_int pivot_odd( const csc *X, //csc *L,
			csc *Ybar_T,
			csc *YbartL0,
			c_int factor_status,
			csc *L_csc,
			c_float *Dinv_arr, 
			c_int size1, // number of rows in Aij == L->n == L->m == X->n == X->m
			c_int size2, // number of columns in Aij == Ybar_T->m == Ybar_TL0->m
			c_int *count_L,
			c_int start_idx,
			c_int nx, c_int nu, c_int ny){
    // nx_ny, nx, &count_L, nu+(nx_ny+nx_nu)*(Nnew-1), nx, nu, ny

    // Allocate memory for the L variable
    //csc *Ltemp = csc_spalloc(nx_ny, nx_ny,  nx_ny*nx_ny, 1, 0);    
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
    for (c_int i=0; i<nx;i++){
        Aij[ny + i + i*(size1)] = -1.0;
    }
    c_int i, j, k;
    c_int sum_Lnz;
    //L->n = X->n;
    //L->m = X->m;
    //if(X->n != size1)
//	printf("WARNING: L->n != size1\n");
    
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
    
    for(c_int i = 0; i < size1; i++) {
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
	    if (fabs(Aij[i + j*size1]) > 1e-10){
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
	    if (fabs(Aij[i + j*size1]) > 1e-10){
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
    amub_col_plus_I (nx, size1, ny, 1,
                     Ybar_T->x, Ybar_T->i, Ybar_T->p,
                     Lx, Li, Lp,
                     YbartL0->x, YbartL0->i, YbartL0->p,
                     200, -1.0);
    
    // Add Ltemp and Y0*Ltemp to large L matrix
    for (i=0; i<size1; i++){
	for (j=Lp[i]; j<Lp[i+1]; j++){
	    L_csc->i[*count_L] = start_idx + Li[j];
	    L_csc->x[(*count_L)++] = Lx[j];
	}
	for (j=YbartL0->p[i]; j<YbartL0->p[i+1]; j++){		
	    L_csc->i[*count_L] = start_idx + size1 + YbartL0->i[j];
	    L_csc->x[(*count_L)++] = YbartL0->x[j];
	}
	L_csc->p[start_idx+i+1] = *count_L;
    }
    for ( i = 0; i < size1; i++){    
	Dinv_arr[start_idx + i] = -Dinv[i]; // Negative
    }
    
    free(Aij);
}





static c_int pivot_even(const csc *X, //csc *L,
			csc *Ybar_T,
			csc *YbartL0,
			c_float *A_f,
			c_int factor_status,
			csc *L_csc,
			c_float *Dinv_arr, 
			c_int size1,
			c_int size2,
			c_int *count_L,
			int start_idx,
			c_int nx, c_int nu, c_int ny){
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
    
    amub_col_plus_I (size2, size1, 0, 1,
		     Ybar_T->x, Ybar_T->i, Ybar_T->p,
		     Lx, Li, Lp,
		     YbartL0->x, YbartL0->i, YbartL0->p,
		     200, 1.0);

    // Add the L0 and Ybar_TL0 values to the L matrix
    for (i=0; i<size1; i++){
	// Add L matrix
	for (j=Lp[i]; j<Lp[i+1]; j++){
	    L_csc->i[*count_L] = start_idx + Li[j];
	    L_csc->x[(*count_L)++] = Lx[j];
	}
	// Add Aii*Ybar matrix under
	for (j=YbartL0->p[i]; j<YbartL0->p[i+1]; j++){
	    L_csc->i[*count_L] = start_idx + size1 + YbartL0->i[j];		
	    L_csc->x[(*count_L)++] = YbartL0->x[j];
	}
	L_csc->p[start_idx  + i + 1] = *count_L;
    }
    for ( i = 0; i < size1; i++){    
	Dinv_arr[start_idx + i] = Dinv[i];
    }
}




static c_int pivot_final(const csc *X, //csc *L,
			 c_int factor_status,
			 csc * L_csc,
			 c_float * Dinv_arr, 
			 c_int num_cols,
			 c_int * count_L,
			 c_int start_idx){
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
    
    sum_Lnz = QDLDL_etree(X->n, X->p, X->i, iwork, Lnz, etree);
    factor_status = QDLDL_factor(X->n, X->p, X->i, X->x,
				 Lp, Li, Lx,
				 D, Dinv, Lnz,
				 etree, bwork, iwork, fwork);
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
}


/**
 * Compute LDL factorization of matrix A
 * @param  P    Cost matrix to be factorized
 * @param  A    Constraint matrix to be factorized
 * @param  p    Private workspace
 * @param  nvar Number of QP variables
 * @return      exitstatus
 */
static c_int LDL_update_from_pivot(qdldl_solver * s, csc *X_even, 
				   const csc *Q0, const csc *Qi, const csc *QN, 
				   const csc *A0,  const csc *Ai,  const csc *Aij,  const csc *AN,  
				   c_float     sigma,
				   const c_float    *rho_inv_vec,
				   c_int nx,   c_int nu, c_int ny, c_int nt,
				   c_int deltaN,c_int Nnew,c_int Nold,c_int Nmax){

    c_int factor_status;
    c_int i, j, k;
    c_int count_L;
    c_int Ntemp;
    c_int nx_ny = nx+ny;  
    c_int nx_nu = nx+nu;
    // Continue from this iteration
    if (deltaN>0)  {
        Ntemp = Nold-1;
        count_L = s->L->p[nu+(Ntemp)*(nx_ny+nx_nu)];
    } else if (deltaN<0)  { // Reduced horizon
        Ntemp = Nnew-1;
        count_L = s->L->p[nu+(Ntemp)*(nx_ny+nx_nu)];

    } else return -1;


    c_int nz_xeven = X_even->p[(Ntemp+1)*(nx_ny)-1]; // Continue from Ntemp?

    // Temporary matrices    ===========================================================
    csc      *Aii = csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);    
    csc      *Xtemp = csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);    
    csc      *Ybar_T = csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);
    csc      *Aii_transpose =  csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);
    csc      *YbartL0 = csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);
    c_float  * Aij_f = (c_float*) c_calloc(nx_ny*nx_ny, sizeof(c_float));

    // Continue from X value: -rhoinv - Ai*Ybar_i1 
    copy_from_offset(Xtemp, X_even, nx_ny, Ntemp*nx_ny-1); 

    if (deltaN>0){
	for (c_int iter=Ntemp;iter<Nnew-1;iter++){
	    pivot_odd(Xtemp, Ybar_T,  YbartL0,
		      factor_status, s->L, s->Dinv,
		      nx_ny, nx_nu, &count_L, nu+(nx_ny+nx_nu)*iter, nx, nu, ny);
	    
	    // Compute next Xtemp: Qi - Ybar_T + I*sigma
	    copy_csc_mat_cols_ident(Qi, Xtemp, 0, 0, nx_nu, nx_nu, sigma);
	    subtract_csc_matrix_cols(ny, ny+nx, Xtemp->n, 
				     Xtemp->p, Xtemp->i, Xtemp->x,
				     Ybar_T->p, Ybar_T->i, Ybar_T->x );

	    // Copy next RHS (Ai) to Aij_f
	    memset(Aij_f, 0, nx_ny*nx_ny*sizeof(c_float));
	    copy_csc_mat_transpose(Ai, Aij_f, 0, 0, nx_ny, nx_nu);
	    pivot_even(Xtemp, Ybar_T,  YbartL0,
		       Aij_f, factor_status, s->L, s->Dinv,
		       nx_nu, nx_ny, &count_L, nu+(nx_ny+nx_nu)*iter + nx_ny, nx, nu, ny );
	    
	    // Compute the next Xtemp:  -rhoinv - Ai*Ybar_ij
	    copy_csc_mat_cols(Ai, Aii, 0, 0, nx_ny, nx_nu);
	    csc_tocsr(nx_ny, nx_nu,
		      Aii->p, Aii->i, Aii->x,
		      Aii_transpose->p, Aii_transpose->i, Aii_transpose->x);// Aii_transpose = Ai'
	    Aii_transpose->m = Aii->n; 
	    Aii_transpose->n = Aii->m; 
	    
	    Xtemp->n = nx_ny;
	    Xtemp->m = nx_ny;
	    amub_col_plus_rho_upper_diag(nx_ny, nx_ny, 1,
				         Ybar_T->x, Ybar_T->i, Ybar_T->p,
				         Aii_transpose->x, Aii_transpose->i, Aii_transpose->p,
				         Xtemp->x, Xtemp->i, Xtemp->p,
				         200, &rho_inv_vec[(iter+1)*nx_ny]);
	    // Save all even X for later updates  
	    copy_mat_offset(X_even, Xtemp, nx_ny, (iter+1)*nx_ny, &nz_xeven);
	}
    }
    // End with terminal cost and constraints
    pivot_odd(Xtemp, Ybar_T,  YbartL0,
	      factor_status,
	      s->L, s->Dinv,
	      nx_ny, nx, &count_L, nu+(nx_ny+nx_nu)*(Nnew-1), nx, nu, ny );
    
    // Terminal cost constraint -------------------------------------
    // Q = QN in R(nx * nx)
    copy_csc_mat_cols_ident(QN, Xtemp, 0, 0, nx, nx, sigma);
    subtract_csc_matrix_cols(ny, ny+nx, nx, 
			     Xtemp->p, Xtemp->i, Xtemp->x,
			     Ybar_T->p, Ybar_T->i, Ybar_T->x );

    memset(Aij_f, 0, nt*nx*sizeof(c_float));
    copy_csc_mat_transpose(AN, Aij_f, 0, 0, nt, nx);
    
    pivot_even(Xtemp, Ybar_T,  YbartL0,
	       Aij_f, factor_status, s->L, s->Dinv,
               nx, nt, &count_L, nu+(nx_ny+nx_nu)*(Nnew-1) + nx_ny , nx, nu, ny);

    //--------------------------------------------------------------------------------
    // Terminal constraint 
    // A = AN in R(nt * nx)
    copy_csc_mat_cols(AN, Aii, 0, 0, nt, nx);
    csc_tocsr(nt, nx,
	      Aii->p, Aii->i, Aii->x,
	      Aii_transpose->p, Aii_transpose->i, Aii_transpose->x);
    Aii_transpose->m = Aii->n; 
    Aii_transpose->n = Aii->m; 
    
    amub_col_plus_rho_upper_diag(nt, nt, 1,
				 Ybar_T->x, Ybar_T->i, Ybar_T->p,
				 Aii_transpose->x, Aii_transpose->i, Aii_transpose->p,
				 Xtemp->x, Xtemp->i, Xtemp->p,
				 200, &rho_inv_vec[(Nmax-1)*nx_ny]);
    Xtemp->n = nt;
    Xtemp->m = nt;

    pivot_final(Xtemp, factor_status, s->L, s->Dinv,
		nt, &count_L, (nx_ny+nx_nu)*Nnew );

    // Fill the rest with zeros (not really necessary?)
    //for (c_int i = Nnew*(nx+nu)+Nnew*(nx+ny)+nt; i<Nmax*(nx+nu) + Nmax*(nx+ny)+nt;i++){
    //    s->L->p[i] = count_L;
    //}


    update_problem_size(s, Nnew*(nx+nu), Nnew*(nx+ny)+nt);

    s->L->n = Nnew*(nx+nu)+Nnew*(nx+ny)+nt;
    s->L->m = Nnew*(nx+nu)+Nnew*(nx+ny)+nt;
    //print_csc_matrix(s->L, "L");
    //printf("Dinv = zeros(%i,%i);\n",  Nnew*(nx+nu)+Nnew*(nx+ny)+nt,  Nnew*(nx+nu)+Nnew*(nx+ny)+nt);
    //for (i=0;i< Nnew*(nx+nu)+Nnew*(nx+ny)+nt;i++){
    //    printf("Dinv(%2i,%2i) = %9.5f;\n", i+1, i+1, s->Dinv[i]);
    //}



    csc_spfree(Aii);
    csc_spfree(Xtemp);
    csc_spfree(Ybar_T);
    csc_spfree(Aii_transpose);
    csc_spfree(YbartL0);
    free(Aij_f);

    return 0;
}




void get_L_dimensions_from_solver(qdldl_solver * s, c_int *n, c_int *m, c_int *nnz){
    *n = s->L->n;
    *m = s->L->m;
    *nnz = s->L->p[s->L->n];
    //printf("\tL(%i) in %i x%i\n", s->L->p[s->L->n], s->L->m, s->L->n);
    //print_csc_matrix(s->L,"L");
    //printf("Dval = zeros(%i, 1);\n", s->L->n);
    //for (int i=0;i<s->L->n;i++){
    //    printf("Dval(%i,1) = %12.8f;\n", i+1 , 1/s->Dinv[i]);
    //}
}

void get_L_dimensions(OSQPWorkspace* work, c_int *n, c_int *m, c_int *nnz){
    get_L_dimensions_from_solver(work->linsys_solver, n, m, nnz);
}

void get_permutation_matrix_from_solver(qdldl_solver * s){
    printf("Perm = zeros(%i, 1);\n", s->L->n);
    for(int i=0;i<s->L->n;i++){
        printf("Perm(%i,1) = %i;\n", i+1 , s->P[i]+1);
    }
}

void get_permutation_matrix(OSQPWorkspace* work){
    get_permutation_matrix_from_solver(work->linsys_solver);
}


/**
 * Compute LDL factorization of matrix A
 * @param  P    Cost matrix to be factorized
 * @param  A    Constraint matrix to be factorized
 * @param  p    Private workspace
 * @param  nvar Number of QP variables
 * @return      exitstatus (0 is good)
 */
static c_int LDL_factorize_recursive(qdldl_solver * s, csc *X_even, 
				     csc *Q0, csc *Qi, csc *QN, 
				     csc *A0,  csc *Ai,  csc *Aij,  csc *AN,  
				     c_float     sigma,
				     c_float    *rho_inv_vec,
				     c_int nx,   c_int nu, c_int ny, c_int nt,
				     c_int N,    c_int Nmax){

    c_int i, j, k;
    c_int count_L = 0;
    c_int factor_status;
    c_int nz_xeven = 0;

    c_int nx_ny = nx+ny;  
    c_int nx_nu = nx+nu;
    

    // Temporary matrices    ===========================================================
    // Be more careful with how many zeros are allocated?
    csc *Aii = csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);    
    csc *Xtemp = csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);    
    csc * Ybar_T = csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);
    csc * Aii_transpose =  csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);
    csc * YbartL0 = csc_spalloc(nx_ny, nx_ny, nx_ny*nx_ny, 1, 0);
    c_float * Aij_f = (c_float*) c_calloc(nx_ny*nx_ny, sizeof(c_float));
    
    // PIVOT POINT 0 =================================pivot_termin==================================

    copy_csc_mat_transpose(A0, Aij_f, 0, 0, nx_ny, nu);
    // Assuming no non-zero elements on diagonal...
    copy_csc_mat_cols_ident(Q0, Xtemp, 0, 0, nu, nu, sigma);
    //printf("* PIVOT POINT 0 -------------------\n"); 

    pivot_even(Xtemp, Ybar_T,  YbartL0,
	       Aij_f, factor_status, s->L, s->Dinv,
               nu, nx_ny, &count_L,  0 , nx, nu, ny);    

    // Next X value: X(2)  = -rhoinv - A0*Ybar_01  
    copy_csc_mat_cols(A0, Aii, 0, 0, nx_ny, nu);
    csc_tocsr(nx_ny, nu, Aii->p, Aii->i, Aii->x,
	      Aii_transpose->p, Aii_transpose->i, Aii_transpose->x); 
    Aii_transpose->n = Aii->m; // Number of rows in Aii
    Aii_transpose->m = Aii->n; // Number of rows in Aii
    Xtemp->n = nx_ny;
    Xtemp->m = nx_ny;    
    // Computes Ybar_T*Aii_transpose + rhoinv*I
    amub_col_plus_rho_upper_diag(nx_ny, nx_ny, 1,
				 Ybar_T->x, Ybar_T->i, Ybar_T->p,
				 Aii_transpose->x, Aii_transpose->i, Aii_transpose->p,
				 Xtemp->x, Xtemp->i, Xtemp->p,
				 200, rho_inv_vec);	

    // Save all even X for later updates
    copy_mat_offset(X_even, Xtemp, nx_ny,0, &nz_xeven);

    for (c_int iter=0;iter<N-1;iter++){
	// Constraint pivot: X = -rhoinv*I - Ai*Ybar_ij
	pivot_odd(Xtemp, Ybar_T,  YbartL0, factor_status,  s->L, s->Dinv,
		  nx_ny, nx_nu, &count_L, nu+(nx_ny+nx_nu)*iter, nx, nu, ny );
	    
	// Compute next Xtemp: Qi - Ybar_T	    
	copy_csc_mat_cols_ident(Qi, Xtemp, 0, 0, nx_nu, nx_nu, sigma);
	// Subtract bottom nx rows of Ybar_T from top of X
	// Ybar_T = Ybar given in row major format - bottom nx rows = rightmost nx?
	subtract_csc_matrix_cols(ny, ny+nx, Xtemp->n, 
				 Xtemp->p, Xtemp->i, Xtemp->x,
				 Ybar_T->p, Ybar_T->i, Ybar_T->x );

	memset(Aij_f, 0, nx_ny*nx_ny*sizeof(c_float));
	copy_csc_mat_transpose(Ai, Aij_f, 0, 0, nx_ny, nx_nu);
	    
	pivot_even(Xtemp, Ybar_T,  YbartL0,
		   Aij_f, factor_status, s->L, s->Dinv,
                   nx_nu, nx_ny, &count_L, nu+(nx_ny+nx_nu)*iter + nx_ny, nx, nu, ny );
	    
	// Compute the next Xtemp:  -rhoinv - Ai*Ybar_ij
	copy_csc_mat_cols(Ai, Aii, 0, 0, nx_ny, nx_nu); // Aii = Ai
	csc_tocsr(nx_ny, nx_nu,
		  Aii->p, Aii->i, Aii->x,
		  Aii_transpose->p, Aii_transpose->i, Aii_transpose->x);// Aii_transpose = Ai'
	Aii_transpose->m = Aii->n; 
	Aii_transpose->n = Aii->m; 
	    
	Xtemp->n = nx_ny;
	Xtemp->m = nx_ny;
	amub_col_plus_rho_upper_diag(nx_ny, nx_ny, 1,
				     Ybar_T->x, Ybar_T->i, Ybar_T->p,
				     Aii_transpose->x, Aii_transpose->i, Aii_transpose->p,
				     Xtemp->x, Xtemp->i, Xtemp->p,
				     200, &rho_inv_vec[(iter+1)*nx_ny]);

    
        // Save all even X for later updates    
        copy_mat_offset(X_even, Xtemp, nx_ny, (iter+1)*nx_ny, &nz_xeven);
    }

    // Second to last (pivot 2N, CONSTRAINT) ---------------------------------------------
    pivot_odd(Xtemp, Ybar_T,  YbartL0, factor_status, s->L, s->Dinv,
	      nx_ny, nx, &count_L,nu+(nx_ny+nx_nu)*(N-1), nx, nu, ny );

    // Terminal cost constraint -------------------------------------
    copy_csc_mat_cols_ident(QN, Xtemp, 0, 0, nx, nx, sigma);
    subtract_csc_matrix_cols(ny, ny+nx, nx, 
			     Xtemp->p, Xtemp->i, Xtemp->x,
			     Ybar_T->p, Ybar_T->i, Ybar_T->x );

    memset(Aij_f, 0, nt*nx*sizeof(c_float));
    copy_csc_mat_transpose(AN, Aij_f, 0, 0, nt, nx);
    
    pivot_even(Xtemp, Ybar_T,  YbartL0,
	       Aij_f, factor_status, s->L, s->Dinv,
               nx, nt, &count_L, nu+(nx_ny+nx_nu)*(N-1) + nx_ny, nx, nu, ny );

    
    // Only include upper diagonal items...
    // For next X: need AN' and Ybar from previous solution


    //--------------------------------------------------------------------------------
    // Terminal constraint 
    // A = AN in R(nt * nx)
    copy_csc_mat_cols(AN, Aii, 0, 0, nt, nx);
    csc_tocsr(nt, nx,
	      Aii->p, Aii->i, Aii->x,
	      Aii_transpose->p, Aii_transpose->i, Aii_transpose->x);
    Aii_transpose->m = Aii->n; 
    Aii_transpose->n = Aii->m; 
    
    amub_col_plus_rho_upper_diag(nt, nt, 1,
				 Ybar_T->x, Ybar_T->i, Ybar_T->p,
				 Aii_transpose->x, Aii_transpose->i, Aii_transpose->p,
				 Xtemp->x, Xtemp->i, Xtemp->p,
				 200, &rho_inv_vec[(Nmax-1)*nx_ny]);
    Xtemp->n = nt;
    Xtemp->m = nt;

    pivot_final(Xtemp,  factor_status, s->L, s->Dinv,
		nt, &count_L, (nx_ny+nx_nu)*N );



    // Fill the rest with zeros (not really necessary?)
    for (c_int i = N*(nx+nu)+N*(nx+ny)+nt; i<Nmax*(nx+nu) + Nmax*(nx+ny)+nt;i++){
        s->L->p[i] = count_L;
    }

    update_problem_size(s, N*(nx+nu), N*(nx+ny)+nt);
    

    s->L->n = N*(nx+nu)+N*(nx+ny)+nt;
    s->L->m = N*(nx+nu)+N*(nx+ny)+nt;

    csc_spfree(Aii);
    csc_spfree(Xtemp);
    csc_spfree(Ybar_T);
    csc_spfree(Aii_transpose);
    csc_spfree(YbartL0);
    free(Aij_f);

    return 0;

}

c_int compute_permutations(qdldl_solver * s, c_int n, c_int m, 
			   c_int N, 
			   const csc * Q0, const csc * Qi, const csc * QN, 
			   const csc * A0, const csc * Ai, const csc * Aij, const csc * AN){

    c_int i, count=0, A_count=n, P_count=0;
    // PERMUTATION VECTOR P: A1 = A[P,P]  =============================================
    for(i=0;i<Q0->n;i++) {
        s->P[count++] = P_count++;
    }
    for(i=0;i<A0->m;i++) s->P[count++] = A_count++;
    
    for(c_int iter=0;iter<N-1;iter++){
        for(i=0;i<Qi->n;i++) s->P[count++] = P_count++;
        for(i=0;i<Ai->m;i++) s->P[count++] = A_count++;
    }

    for(i=0;i<QN->n;i++) s->P[count++] = P_count++;
    for(i=0;i<AN->m;i++) s->P[count++] = A_count++;

}


c_int compute_KKT_permutations(qdldl_solver * s, c_int n, c_int m,  c_int P_nnz_max, c_int A_nnz_max, 
                               c_int N, 
                               const csc * Q0, const csc * Qi, const csc * QN, 
                               const csc * A0, const csc * Ai, const csc * Aij, const csc * AN){
    c_int i, count=0, A_count=n, P_count=0;

    // Allocate vector of indices on the diagonal. Worst case it has m elements
    if (s->Pdiag_idx != OSQP_NULL) {
	(*s->Pdiag_idx) = c_malloc(n * sizeof(c_int));
	s->Pdiag_n     = 0; // Set 0 diagonal elements to start
    }

    // PERMUTATION VECTOR P: A1 = A[P,P]  =============================================
    compute_permutations(s, n, m, 
                         N, Q0, Qi, QN,  A0,  Ai, Aij, AN);




    // PERMUTATION VECTORS ============================================================
    c_int  zKKT = 0;        // Counter for total number of elements in P and in
    // Initial cost Q0 ============================================================
    for (c_int j = 0; j < Q0->n; j++) { // cycle over columns in Q0
        // No elements in column j => add diagonal element param1
        if (Q0->p[j] == Q0->p[j + 1]) zKKT++;
        for (c_int ptr = Q0->p[j]; ptr < Q0->p[j + 1]; ptr++) { // cycle over rows
            if (s->PtoKKT != OSQP_NULL) s->PtoKKT[ptr] = zKKT;  
            if (Q0->i[ptr] == j) { // P has a diagonal element,
		// If index vector pointer supplied -> Store the index
		if (s->Pdiag_idx != OSQP_NULL) {
		    (s->Pdiag_idx)[s->Pdiag_n] = ptr;
		    (s->Pdiag_n)++;
		} 
            }
            zKKT++;
            // Add diagonal param1 in case last element of column j
            if ((Q0->i[ptr] < j) && (ptr + 1 == Q0->p[j + 1])) zKKT++;
        }
    }

    // Initial constraint A0 ============================================================
    for (c_int j = 0; j < A0->n; j++) {                      // Cycle over columns of A: should this be the rows? A'?
        for (c_int ptr = A0->p[j]; ptr < A0->p[j + 1]; ptr++) {
	    if (s->AtoKKT != OSQP_NULL) s->AtoKKT[ptr] = zKKT;  // Update index from A to KKTtrip
	    zKKT++;
        }

        if (s->rhotoKKT != OSQP_NULL) s->rhotoKKT[j] = zKKT;  // Update index from  param2 to KKTtrip
        zKKT++;
    }

    //  BEGIN ITERATIONS =================================================================
    for (c_int iter=0;iter<N-1;iter++){
	// cycle over columns in Aij/Qi
        for (c_int j = 0; j < Aij->n; j++) {
            // Aij  
            for (c_int ptr = Aij->p[j]; ptr < Aij->p[j + 1]; ptr++) {
                if (s->AtoKKT != OSQP_NULL) s->AtoKKT[ptr] = zKKT;  // Update index from A to KKTtrip
                zKKT++;
            }

            // Qi
            if (Qi->p[j] == Qi->p[j + 1]) zKKT++;
            for (c_int ptr = Qi->p[j]; ptr < Qi->p[j + 1]; ptr++) { // cycle over rows
                if (s->PtoKKT != OSQP_NULL) s->PtoKKT[ptr] = zKKT;  // Update index from P to
                zKKT++;
                // Add diagonal param1 in case last element of column j
                if ((Qi->i[ptr] < j) && (ptr + 1 == Qi->p[j + 1])) zKKT++;
            }
        }


	// cycle over columns in Ai/rhoinv
        // Ai
        for (c_int j = 0; j < Ai->n; j++) {                    
            // Ai: should this be the transpose?
            for (c_int ptr = Ai->p[j]; ptr < Ai->p[j + 1]; ptr++) {
                if (s->AtoKKT != OSQP_NULL) s->AtoKKT[ptr] = zKKT;  // Update index from A to KKTtrip
                zKKT++;
            }
            // rhoinv*I only adds one each column. 
            if (s->rhotoKKT != OSQP_NULL) s->rhotoKKT[j] = zKKT;  // Update index from  param2 to KKTtrip
            zKKT++;
        }




    }


    return 0;
}

c_int print_P_matrix(qdldl_solver * s){
    printf("* print_P_matrix()\n");
    for(c_int i=0;i<20;i++){
        //s->P[i] = 2*i;
        printf("P[%i] = %i\n", i, s->P[i]);

    }
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


    //print_csc_matrix(data->P,"Before");

    if (Nnew < Nold) {  // Return last few columns
        Ntemp = Nnew -2;
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

c_int setup_AP_matrices(OSQPData * data, c_int Nmax, c_int N, c_int nx, c_int nu, c_int ny, c_int nt, 
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
    if(!(N<=data_rldl->Nmax)){
        printf("ERROR: the new horizon is too large: %i\n", N);
        return -1;
    }
    if (deltaN == 0){
        printf("deltaN = 0, returning\n");
        return 0;
    }
    update_problem_size(work->linsys_solver, N*(data_rldl->nx+data_rldl->nu), N*(data_rldl->nx+data_rldl->ny)+data_rldl->nt);

    c_int result = LDL_update_from_pivot(work->linsys_solver, data_rldl->X_even, 
					 data_rldl->Q0, data_rldl->Qi, data_rldl->QN, 
					 data_rldl->A0,  data_rldl->Ai,  data_rldl->Aij,  data_rldl->AN,  
					 work->settings->sigma, work->rho_inv_vec,
					 data_rldl->nx,   data_rldl->nu, data_rldl->ny, data_rldl->nt,
					 deltaN, N, data_rldl->N, data_rldl->Nmax);

    update_AP_matrices( data_rldl, work->data, data_rldl->N, N);


    c_int exitflag =  compute_permutations(work->linsys_solver, N*(data_rldl->nx+data_rldl->nu), N*(data_rldl->nx+data_rldl->ny)+data_rldl->nt, 
					   N, data_rldl->Q0, data_rldl->Qi, data_rldl->QN,  data_rldl->A0,  data_rldl->Ai, data_rldl->Aij, data_rldl->AN);
    data_rldl->N = N;

}

c_int osqp_setup_recursive(OSQPWorkspace** workp, OSQPDataRLDL *data_rldl, const OSQPSettings *settings, c_int Nmax, c_int N, c_int nx, c_int nu, c_int ny, c_int nt) {
    printf("WARNING: Please note that in the current implementation, the Aij matrix is not used. Instead, it is assumed to have the shape [0,0;-I,0]. (FIX THIS)\n");
    c_int exitflag;
    c_float epsilon = 1e-7;
    c_int nz=0,count = 0;

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

    // SET SIZES OF ALL MATRICES =======================================================

    data_rldl->X_even = c_malloc(sizeof(csc));
    data_rldl->X_even->n = Nmax*(nx+ny);  // number of columns
    data_rldl->X_even->m = nx+ny;         // number of rows
    data_rldl->X_even->nz = -1;
    data_rldl->X_even->p = (c_int *)c_malloc((nx+ny)*Nmax * sizeof(c_int));
    data_rldl->X_even->i = (c_int *)c_malloc(sizeof(c_int) *(nx+ny)*(nx+ny)*Nmax);
    data_rldl->X_even->x = (c_float *)c_malloc(sizeof(c_float)*(nx+ny)*(nx+ny)*Nmax);
    data_rldl->X_even->p[0] = 0;

    exitflag = init_linsys_solver_qdldl_recursive(&(work->linsys_solver), n, m, P_nnz_max*5, A_nnz_max, work->settings->sigma, work->rho_vec, work->pol);
    
    exitflag =  LDL_factorize_recursive(work->linsys_solver, data_rldl->X_even,
					data_rldl->Q0, data_rldl->Qi, data_rldl->QN, 
					data_rldl->A0,  data_rldl->Ai,  data_rldl->Aij,  data_rldl->AN,  
					work->settings->sigma, work->rho_inv_vec,
					nx,   nu, ny, nt,
					N,    Nmax);
    if (exitflag) return osqp_error(exitflag);




    exitflag =  compute_KKT_permutations(work->linsys_solver, N*(nx+nu), N*(nx+ny)+nt, P_nnz_max, A_nnz_max, 
                                         N, Q0, Qi, QN,  A0,  Ai, Aij, AN);


    if (exitflag) return osqp_error(exitflag);

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





c_int osqp_solve_recursive(OSQPWorkspace *work, c_int N, c_int nx, c_int nu, c_int ny, c_int nt) {

    /*for (int i=0;i<  work->data->m; i++){
      printf("%3i\t%9.4f\t%9.4f\n", i,work->data->l[i], work->data->u[i] );
      }*/

    //print_csc_matrix(work->data->A, "A");
    //print_csc_matrix(work->data->P, "P");


    c_int exitflag = osqp_solve(work);
    /*for (int i=0;i<N;i++){
      printf("%f\n", work->solution->x[i]);
      }*/

    return exitflag;
}



