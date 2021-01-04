#include "cs_addon.h"

void copy_csc_mat_in_place(int An, const c_float *Ax, const c_int *Ai, const c_int*Ap,
                           c_float *Bx, c_int *Bi, c_int*Bp){
    prea_int_vec_copy(Ap, Bp, An + 1);
    prea_int_vec_copy(Ai, Bi, Ap[An]);
    prea_vec_copy(Ax, Bx, Ap[An]);
}

csc* copy_csc_mat_cols(const csc *A, csc *B) {
    if (!B) return OSQP_NULL;

    c_int i, column;
    c_int nz = 0;    
    for (column = 0; column < A->n; column++) {
        // For every row pointer
        B->p[column] = nz;
        for (i = A->p[column]; i < A->p[column + 1]; i++) {
            B->i[nz] = A->i[i];
            B->x[nz] = A->x[i];
            nz++;
        }
    }
    B->p[A->n] = nz;
    B->m = (A->m);
    B->n = (A->n);
    return B;
}


void copy_csc_mat_transpose(const csc *A, c_float * data) {
    c_int columns = A->n;
    c_int i, column, row;
    c_int nz = 0;
    for (column = 0; column < columns; column++) {
        for (i = A->p[column]; i < A->p[column + 1]; i++) {
            row = A->i[i];
            data[row*columns + (column)] = A->x[i];
        }
    }
}

void copy_mat_offset(csc *A, csc *B, c_int num_col, c_int column_offset, c_int * count){
    // Copy matrix B to matrix A, with an offset in column
    for (int i=0;i<num_col;i++){
        for(int p=B->p[i];p<B->p[i+1];p++){
            A->i[*count] = B->i[p];
            A->x[(*count)++] = B->x[p];
        }    
        A->p[i + column_offset] = *count;
    }

}


void copy_from_offset(const csc *A, csc *B, c_int num_col, c_int column_offset){
    c_int count = 0;
    A->p[0] = count;
    // Copy matrix B to matrix A, with an offset in column
    for (int i=column_offset;i<column_offset+num_col;i++){
        for(int p=B->p[i];p<B->p[i+1];p++){
            A->i[count] = B->i[p];
            A->x[count++] = B->x[p];
        }    
        A->p[i-column_offset+1] = count;
    }

}


void copy_csc_mat_cols_ident(const csc *A, csc *B, c_float offset) {
    c_int col, i;
    c_int nz = 0;
    
    for (col = 0; col < A->n; col++) {
        B->p[col] = nz;
        for (i = A->p[col]; i < A->p[col + 1]; i++) {
            if (A->i[i] == col) B->x[nz] = A->x[i] + offset;
            else                B->x[nz] = A->x[i];
            B->i[nz++] = A->i[i];
        } 
    }
    B->p[A->n] = nz;
    B->nz = nz;
    B->n = A->n;
    B->m = A->m;
}


void print_csc_matrix(csc *M, const char *name){
    c_int j, i, row_start, row_stop;
    c_int k = 0;

    c_print("%smat = zeros(%i,%i); \n", name, M->m, M->n);

    for (j = 0; j < M->n; j++) {
        row_start = M->p[j];
        row_stop  = M->p[j + 1];

        if (row_start == row_stop) continue;
        else {
            for (i = row_start; i < row_stop; i++) {
                c_print("%smat(%3u,%3u) = %23.19f;\n", name,
                        (int)M->i[i]+1, (int)j+1, M->x[k++]);
            }
        }
    }
}


void csr_to_csc(const c_int  n_row, const c_int  n_col, 
                const c_int * Ap, const c_int * Aj,  const c_float * Ax,
                c_int *  Bp,  c_int * Bi, c_float * Bx){
    c_int i, jj,  temp, dest;
    c_int row, col;
    c_int nnz = Ap[n_row];
    
    for (i=0; i<n_col;i++) Bp[i] = 0; 

    for (c_int  n = 0; n < nnz; n++){
        Bp[Aj[n]]++;
    }

    //cumsum the nnz per column to get Bp[]
    c_int cumsum = 0;
    for(col = 0; col < n_col; col++){     
        temp  = Bp[col];
        Bp[col] = cumsum;
        cumsum += temp;
    }
    Bp[n_col] = nnz; 

    for(row = 0; row < n_row; row++){
        for(jj = Ap[row]; jj < Ap[row+1]; jj++){
            col  = Aj[jj];
            dest = Bp[col];
            Bi[dest] = row;
            Bx[dest] = Ax[jj];
            Bp[col]++;
        }
    }  

    c_int last = 0;
    for(col = 0; col <= n_col; col++){
        temp  = Bp[col];
        Bp[col] = last;
        last    = temp;
    }

}   

void csc_to_csr(const c_int n_row, const c_int n_col, 
                const c_int * Ap, const c_int * Ai, const c_float * Ax,
                c_int * Bp, c_int * Bj, c_float * Bx){
    csr_to_csc(n_col, n_row, Ap, Ai, Ax, Bp, Bj, Bx);
}


int amub_col_plus_rho_upper_diag (int ncol, int values,
                                  c_float *Bx, c_int *Bi, c_int*Bp,
                                  c_float *Ax, c_int *Ai, c_int*Ap,
                                  c_float *Cx, c_int *Ci, c_int*Cp,
                                  int nzmax,
                                  c_float * rhoinv) {
    // Computes B*A + rho*I
    // The output matrix has dimensions ncol*ncol
    // Only include the upper diagonal values.    
    c_int ii, jj, jrow, jpos;
    c_int A_ptr, B_ptr;
    c_float coeff;
    c_int len = 0, diag_elem=0;   

    // Set work array = -1
    c_int  iw[ncol];
    for (int j=0; j < ncol; j++) iw[j] = -1;

    // Iterate over all columns in A
    for (ii=0; ii < ncol; ii++){
	Cp[ii] = len;	
	if (Ap[ii] ==  Ap[ii+1]){ // Only add the diagonal element
	    Ci[len] = ii;
	    Cx[len] = rhoinv[ii];
	    diag_elem++;	    
	    len++;
	}else{
	    for (A_ptr=Ap[ii]; A_ptr < Ap[ii+1]; A_ptr++){
		if (values){
		    coeff = Ax[A_ptr];
		}
		jj = Ai[A_ptr];

		for (B_ptr=Bp[jj]; B_ptr < Bp[jj+1]; B_ptr++){
                    // when jrow == ii, on diagonal, if jrow > ii, diagonal passed
		    jrow = Bi[B_ptr];
		    if (jrow <= ii){
			if (jrow > ii && diag_elem==ii){
			    Ci[len] = jrow;
			    Cx[len] = rhoinv[ii];
			    len++;
			    diag_elem++;
			}
			jpos = iw[jrow];
                        if (jpos == -1){ // Add new value
			    if (len > nzmax) return ii;
                            
			    Ci[len] = jrow;
			    iw[jrow]= len;
			    if (values)
				Cx[len]  = coeff*Bx[B_ptr];

			    if (jrow == ii){
				Cx[len] += rhoinv[ii];
				diag_elem++;
			    }
		    
			    len++;
			} else if (values){ // Add to old value
			    Cx[jpos] += coeff*Bx[B_ptr];
			}
		    }
		}
	    }
	}
	// Reset the work vector
	for (int k=Cp[ii]; k < len; k++) iw[Ci[k]] = -1;

    }
    Cp[ncol] = len;	
	
    return 0;
}


int amub_col (int nrow, int ncol, int values,
	      c_float *Bx, c_int *Bi, c_int*Bp,
	      c_float *Ax, c_int *Ai, c_int*Ap,
	      c_float *Cx, c_int *Ci, c_int*Cp,
	      int nzmax) {
    int ii, jj, jrow, jpos;
    c_int A_ptr, B_ptr;
    c_float coeff;
    int len = 0;
    Cp[0] = 0 ;

    // Set work array = -1
    c_int iw[ncol];
    for (int j=0; j < ncol; j++) iw[j] = -1;

    // Iterate over all columns in A
    for (ii=0; ii < ncol; ii++){
        for (A_ptr=Ap[ii]; A_ptr < Ap[ii+1]; A_ptr++){
            if (values){
                coeff = Ax[A_ptr];
            }
            jj   = Ai[A_ptr];
            for (B_ptr=Bp[jj]; B_ptr < Bp[jj+1]; B_ptr++){
                jrow = Bi[B_ptr];
                jpos = iw[jrow];
                if (jpos == -1){ // Add new value
                    if (len > nzmax){
                        return ii;
                    }
                    Ci[len] = jrow;
                    iw[jrow]= len;
                    if (values){
                        Cx[len]  = coeff*Bx[B_ptr];
                    }
                    len++;
                } else if (values){ // Add to old value
                    Cx[jpos] += coeff*Bx[B_ptr];
                }
            }
        }
        // Reset the work vector
        for (int k=Cp[ii]; k < len; k++) iw[Ci[k]] = -1;
        Cp[ii+1] = len;
    }
    return 0;
}




int amub_col_plus_I (c_int nrow, c_int ncol, c_int start_col,
		             c_int values,
                     c_float *Bx, c_int *Bi, c_int*Bp,
                     c_float *Ax, c_int *Ai, c_int*Ap,
                     c_float *Cx, c_int *Ci, c_int*Cp,
                     c_int nzmax, double scale) {
    // Solves B*(A+I), where A is lower diagonal without diagonal entries.
    //printf("amub_col_plus_I()\n");
    int ii, jj, jrow, jpos;
    c_int A_ptr, B_ptr;
    c_float coeff;
    int len = 0;
    Cp[0] = 0;

    // Set work array = -1
    c_int  iw[nrow];
    for (int j=0; j < nrow; j++) iw[j] = -1;

    // If A+I = I: Copy the B matrix
    if (Ap[ncol] == 0){
        copy_csc_mat_in_place(ncol, Bx, Bi, Bp, Cx, Ci, Cp);
        return 0;
    }
    
    // Iterate over all columns in A matrix
    for (ii=start_col; ii < ncol; ii++){
        // First add the diagonal element -----------------
        for (B_ptr=Bp[ii]; B_ptr < Bp[ii+1]; B_ptr++){
            jrow = Bi[B_ptr];
            jpos = iw[jrow];
            if (jpos == -1){ // Add new value
                if (len > nzmax) return ii;
                Ci[len] = jrow;
                iw[jrow]= len;
		//printf("I%i: Add new value %.10e\n", ii, scale*Bx[B_ptr]);
                if (values) Cx[len]  = scale*Bx[B_ptr];
                len++;
            } else if (values){ // Add to old value
                Cx[jpos] += scale*Bx[B_ptr];
            }
        } 
	    
        // Now add rest of the values --------------------
        for (A_ptr=Ap[ii]; A_ptr < Ap[ii+1]; A_ptr++){
            if (values){
                coeff = Ax[A_ptr];
            }
            jj   = Ai[A_ptr];
            for (B_ptr=Bp[jj]; B_ptr < Bp[jj+1]; B_ptr++){
                jrow = Bi[B_ptr];
                jpos = iw[jrow];
                if (jpos == -1){ // Add new value
                    if (len > nzmax) return ii;
                    Ci[len] = jrow;
                    iw[jrow]= len;
		    //printf("L%i: Add new value %.10e\n", ii, coeff*scale*Bx[B_ptr]);
                    if (values) Cx[len]  = coeff*scale*Bx[B_ptr];
                    len++;
                } else if (values){ // Add to old value
		    //printf("L%i: Add to old value %.10e\n", ii, coeff*scale*Bx[B_ptr]);
                    Cx[jpos] += coeff*scale*Bx[B_ptr];
                }
            }
        }
        // Reset the work vector
        for (int k=Cp[ii]; k < len; k++) iw[Ci[k]] = -1;
        Cp[ii+1] = len;
    }


    return 0;
}

int amub_row (int nrow, int ncol, int values,
              c_float *Ax, c_int *Ai, c_int*Ap,
              c_float *Bx, c_int *Bi, c_int*Bp,
              c_float *Cx, c_int *Ci, c_int*Cp,
              int nzmax) {
    int ii, jj, jcol, jpos;
    c_int A_ptr, B_ptr;
    c_float coeff;
    int len = 0;
    Cp[0] = 0 ;

    // Work array
    c_int  iw[nrow];
    for (int j=0; j< ncol; j++) iw[j] = 0;

    // Iterate over all rows in A    
    for (ii=0; ii < nrow; ii++){
	for (A_ptr=Ap[ii]; A_ptr < Ap[ii+1]; A_ptr++){
	    if (values){
		coeff = Ax[A_ptr];
	    }
	    jj   = Ai[A_ptr];
	    for (B_ptr=Bp[jj]; B_ptr < Bp[jj+1]; B_ptr++){
		jcol = Bi[B_ptr];
		jpos = iw[jcol];
		// Has this already been added?
		if (jpos == 0){
		    if (len > nzmax){
			return ii;
		    }
		    Ci[len] = jcol;
		    iw[jcol]= len;
		    if (values){
			Cx[len]  = coeff*Bx[B_ptr];
		    }
		    len++;
		} else if (values){
		    Cx[jpos] += coeff*Bx[B_ptr];
		}
	    } 
	} 
	for (int k=Cp[ii]; k < len; k++) iw[Ci[k]] = 0;
	Cp[ii+1] = len;
    } // end for ii
    return 0;
}


void subtract_csc_matrix_cols(int start_col, int stop_col, int n, 
                              c_int * Ap, c_int * Ai, c_float * Ax,
                              c_int * Bp, c_int * Bi, c_float * Bx ){
    // Subtract the columns between start_col and end_col of Bx to Ax.
    // Ax is upper triangular, Bx is not. 
    c_int j, i;
    c_int count=0;

    // Save new result in C matrix
    c_int   Cp[n+1];
    c_int   Ci[Ap[n]+Bp[n]]; 
    c_float Cx[Ap[n]+Bp[n]];

    // First: loop over initial columns
    int A_ptr, B_ptr, Amax, Bmax;
    int A_row, B_row;
    Cp[0] = 0;
        
    for (j = 0; j < stop_col-start_col; j++) {
        A_ptr = Ap[j];
        Amax  = Ap[j+1];
        B_ptr = Bp[start_col+j];
        Bmax  = Bp[start_col+j+1];
    
        while(A_ptr < Amax || (B_ptr < Bmax && Bi[B_ptr]<=j)){
            if (A_ptr < Amax)  A_row = Ai[A_ptr];
            else               A_row = n+1;
            
            if (B_ptr < Bmax && Bi[B_ptr]<=j)  B_row = Bi[B_ptr];
            else                               B_row = n+1;

            if (A_row  ==  B_row){
                Cx[count] = Ax[A_ptr]-Bx[B_ptr];
                Ci[count] = A_row;
                A_ptr++;
                B_ptr++;
                count++;
            }else if (A_row < B_row) {
                Ci[count] = A_row;
                Cx[count] = Ax[A_ptr];
                A_ptr++;
                count++;
            }else if (A_row > B_row){
                Ci[count] = B_row;
                Cx[count] = -Bx[B_ptr];
                B_ptr++;
                count++;
            }
        }
        Cp[j+1] = count;    
    }

    // Add the remainder of the columns 
    for (j=stop_col-start_col; j<n; j++) {    
        for (A_ptr=Ap[j]; A_ptr<Ap[j+1]; A_ptr++){
            Ci[count] = Ai[A_ptr];
            Cx[count++] = Ax[A_ptr];
        }
        Cp[j+1] = count;    
    }

    // Copy C matrix to A matrix
    for (i = 0; i <= n; i++)
        Ap[i] = Cp[i];
    for (i=0; i<count; i++){
        Ai[i] = Ci[i];
        Ax[i] = Cx[i];
    }
    
}



