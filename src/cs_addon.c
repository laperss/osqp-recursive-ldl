//#include "cs.h"
#include "cs_addon.h"



void array_to_csc(csc* C, c_float* matrix) {

  c_int nz = 0;
  c_float value;
  c_float epsilon = 1e-7;
  c_int rows = C->m;
  c_int count = 0;

  // For every column
  for(int i=0; i <C->n; i++){
    C->p[i] = count;
    for(int j=0;j<rows;j++){
      value = matrix[i*rows + j];
      if (fabs(value)>epsilon){
        C->i[count] = j;
        C->x[count++] = value;
      }
    }
  }
  C->p[C->n] = count;
  //assert(C->nz == count);

  return C;
}

void copy_csc_mat_in_place(int An, const c_float *Ax, const c_int *Ai, const c_int*Ap,
                           c_float *Bx, c_int *Bi, c_int*Bp){
  prea_int_vec_copy(Ap, Bp, An + 1);
  prea_int_vec_copy(Ai, Bi, Ap[An]);
  prea_vec_copy(Ax, Bx, Ap[An]);

}

csc* copy_csc_mat_cols(const csc *A, csc *B,
                       c_int row_start, c_int col_start,
                       c_int row_stop, c_int col_stop) {
    if (!B) return OSQP_NULL;

    c_int i, column;
    c_int nz = 0;    
    for (column = col_start; column < col_stop; column++) {
    // For every row pointer
    B->p[column-col_start] = nz;
    for (i = A->p[column]; i < A->p[column + 1]; i++) {
        if (A->i[i]>=row_start){
        B->i[nz] = A->i[i] - row_start;
        B->x[nz] = A->x[i];
        nz++;
        }
    }
    }
    B->p[col_stop-col_start] = nz;
    B->m = (row_stop-row_start);
    B->n = (col_stop-col_start);
    return B;
}




c_float *  copy_csc_submat(const csc *A,
                       c_int row_start, c_int col_start,
                       c_int row_stop, c_int col_stop) {
    // Save this into column format (FORTRAN)    
    //printf("WARNING: CHECK THAT THIS FUNCTION WORKS...\n");
    c_float * data = (c_float*) calloc((row_stop-row_start)*(col_stop-col_start), sizeof(c_float));

    c_int i, column, rows, columns;
    rows = row_stop-row_start;
    columns = col_stop-col_start;
    c_int nz = 0;
    // For every column
    for (column = col_start; column < col_stop; column++) {
    // For every row pointer
    for (i = A->p[column]; i < A->p[column + 1], A->i[i]<row_stop; i++) {
        data[A->i[i] + column*rows] = A->x[i];
    }
    }
    return data;

}

void copy_csc_mat_transpose(const csc *A, c_float * data,
                                 c_int row_start, c_int col_start,
                                 c_int row_stop, c_int col_stop) {
    // Save this into column format (FORTRAN)

    c_int i, column, rows, columns;
    c_int row;
    rows = row_stop-row_start;
    columns = col_stop-col_start;
    c_int nz = 0;
    // For every column
    for (column = col_start; column < col_stop; column++) {
    // For every row pointer
    for (i = A->p[column]; i < A->p[column + 1]; i++) {
        row = A->i[i];
        if (row>=row_start){
        data[(row-row_start)*columns + (column-col_start)] = A->x[i];
        }

        
    }
    }

    //printf("\n");
    return;
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


void copy_csc_mat_cols_ident(const csc *A, csc *B,
                             c_int row_start, c_int col_start,
                             c_int row_stop,   c_int col_stop, c_float offset) {

    c_int col, i;
    c_int k = 0;
  
    c_int nz = 0;
    for (col = col_start; col < col_stop; col++) {
        B->p[col-col_start] = nz;
        //printf("Col = %i,  p  = %i\n", col, nz);
        for (i = A->p[col]; i < A->p[col + 1]; i++) {
            if (A->i[i] == col) B->x[k] = A->x[i] + offset;
            else                B->x[k] = A->x[i];
            B->i[k++] = A->i[i] - row_start;
            //printf("nz++\n");
            nz++;
        } 
    }
    B->p[col_stop-col_start] = nz;
    B->nz = nz;
    B->n = col_stop-col_start;
    B->m = row_stop-row_start;
}


void print_csc_matrix(csc *M, const char *name){
  c_int j, i, row_start, row_stop;
  c_int k = 0;

  c_print("%smat = zeros(%i,%i); \n", name, M->m, M->n);
  // Loop over each column
  for (j = 0; j < M->n; j++) {
    row_start = M->p[j];
    row_stop  = M->p[j + 1];

    if (row_start == row_stop) continue;
    else {
      for (i = row_start; i < row_stop; i++) {
          c_print("%smat(%3u,%3u) = %12.9f;\n", name,   //\t\t[%i->%i]
                  (int)M->i[i]+1, (int)j+1, M->x[k++]); // , (int)row_start, (int)row_stop
      }
    }
  }
}


void print_csc_matrix_row_major(csc *M, const char *name)
{
    c_int j, i, column_start, column_stop;
    c_int k = 0;

    // Print name
    printf("print %s : (%i x %i)\n", name, M->n, M->m);

    for (j = 0; j < M->n; j++) {
    column_start = M->p[j];
    column_stop  = M->p[j + 1];
    

    if (column_start == column_stop) continue;
    else {
        for (i = column_start; i < column_stop; i++) {
        printf("\t[%3u,%3u] = %11.7f\t\t[%i => %i]\n",  (int)j, (int)M->i[i], M->x[k++], column_start, column_stop);
        }
    }
    }
    printf("\n");
}




void print_csc_matrix_cols(csc *M, const char *name, int start_col, int stop_col)
{
    c_int j, i, row_start, row_stop;
    c_int k = M->p[start_col];

    // Print name
    printf("Start column = %i,  stop column = %i\n", start_col, stop_col);
    printf("%smat =  zeros(%i , %i);\n", name, M->m, M->n);

    for (j = start_col; j < stop_col; j++) {
	row_start = M->p[j];
	row_stop  = M->p[j + 1];

	if (row_start == row_stop) continue;
	else {
	    for (i = row_start; i < row_stop; i++) {
		printf("%smat(%3u,%3u) = %12.8f;\n", name, (int)M->i[i]+1, (int)j-start_col+1, M->x[k++]);
	    }
	}
    }
    printf("\n");
}



void print_csc_matrix_dense(csc *M, const char *name, int cols, int rows)
{
  c_int j, i, row_start, row_stop;
  c_int k = 0;
  int nrow, ncol;

  nrow = M->m;
  ncol = M->n;
  printf("%s  =  (%i x %i)\n", name, nrow, ncol);
  
  c_float mat[nrow*cols];
  for (i=0; i<nrow*cols; i++) mat[i] = 0.0;

  // Loop over column
  printf("%s  =  (%i x %i [%i])\n", name, nrow, cols, ncol);
  for (j = 0; j < cols; j++) {
      row_start = M->p[j];
      row_stop  = M->p[j + 1];

      if (row_start == row_stop) continue;
      else {
	  for (i = row_start; i < row_stop; i++) {
	      
	      mat[M->i[k]*cols + j] = M->x[k++];
	  }
      }
  }

  for (i=0; i<rows; i++){
      for (j=0; j<cols; j++){
	  printf("%5.3f ", mat[i*rows+j]);
      }
      printf("\n");
  }
  
}



void csr_tocsc(const c_int  n_row, const c_int  n_col, 
           const c_int * Ap, const c_int * Aj,  const c_float * Ax,
           c_int *  Bp,  c_int * Bi, c_float * Bx){
    int i, jj, row, temp, dest;
    int nnz = Ap[n_row];
    for (i=0; i<n_col;i++) Bp[i] = 0; 

    for (c_int  n = 0; n < nnz; n++){
        Bp[Aj[n]]++;
    }

    //cumsum the nnz per column to get Bp[]
    for(c_int  col = 0, cumsum = 0; col < n_col; col++){     
        c_int  temp  = Bp[col];
        Bp[col] = cumsum;
        cumsum += temp;
    }
    Bp[n_col] = nnz; 

    for(c_int  row = 0; row < n_row; row++){
        for(c_int  jj = Ap[row]; jj < Ap[row+1]; jj++){
            c_int  col  = Aj[jj];
            c_int  dest = Bp[col];

            Bi[dest] = row;
            Bx[dest] = Ax[jj];
            Bp[col]++;
        }
    }  

    for(c_int  col = 0, last = 0; col <= n_col; col++){
        c_int  temp  = Bp[col];
        Bp[col] = last;
        last    = temp;
    }

}   

void csc_tocsr(const c_int n_row, const c_int n_col, 
               const c_int * Ap, const c_int * Ai, const c_float * Ax,
                     c_int * Bp, c_int * Bj, c_float * Bx){
    csr_tocsc(n_col, n_row, Ap, Ai, Ax, Bp, Bj, Bx);
}


int amub_col_plus_rho_upper_diag (int nrow, int ncol, int values,
				                  c_float *Bx, c_int *Bi, c_int*Bp,
				                  c_float *Ax, c_int *Ai, c_int*Ap,
				                  c_float *Cx, c_int *Ci, c_int*Cp,
				                  int nzmax,
				                  c_float * rhoinv) {
    // Computes B*A + rho*I
    // The output matrix has dimensions nrow*ncol
    // Only include the upper diagonal values.    
    int ii, jj, jrow, jpos;
    c_int A_ptr, B_ptr;
    c_float coeff;
    int len = 0, diag_elem=0;   

    // Work array
    c_int  iw[ncol];
    for (int j=0; j < ncol; j++) iw[j] = -1;

    // Iterate over all columns in A
    for (ii=0; ii < ncol; ii++){

	Cp[ii] = len;	
	if (Ap[ii] ==  Ap[ii+1]){
	    // Only add the diagonal element
	    Ci[len] = ii;
	    Cx[len] = rhoinv[ii];
	    //printf("0 ADD RHOINV[%i] = %f \n", ii, rhoinv[ii]);
	    diag_elem++;	    
	    len++;
	}else{
	    for (A_ptr=Ap[ii]; A_ptr < Ap[ii+1]; A_ptr++){
		if (values){
		    coeff = Ax[A_ptr];
		}
		jj   = Ai[A_ptr];
		//printf("\t* A row/B col = %i:   B_ptr = %i -> %i\n", jj, Bp[jj], Bp[jj+1]);

		for (B_ptr=Bp[jj]; B_ptr < Bp[jj+1]; B_ptr++){
                    // when jrow == ii, on diagonal, if jrow > ii, diagonal passed
		    jrow = Bi[B_ptr];
		    if (jrow <= ii){
			if (jrow > ii && diag_elem==ii){
			    Ci[len] = jrow;
			    Cx[len] = rhoinv[ii];
			    //printf("1 ADD RHOINV[%i] = %f\n", ii, rhoinv[ii]);			    
			    len++;
			    diag_elem++;
			    // must first add the diagonal element
			}
			jpos = iw[jrow];
			//printf("\t\t\t* jpos = iw[%i] = %i  =>  %f\n", jrow, iw[jrow], Bx[B_ptr]);
			// Has this already been added?
			if (jpos == -1){
			    if (len > nzmax){
				return ii;
			    }
			    Ci[len] = jrow;
			    iw[jrow]= len;
			    if (values){
				Cx[len]  = coeff*Bx[B_ptr];
			    }
			    if (jrow == ii){
				Cx[len] += rhoinv[ii];
				//printf("2 ADD RHOINV[%i,%i] = %f\n", jrow, ii, rhoinv[ii]);			    
				
				diag_elem++;
			    }
		    
			    len++;
			} else if (values){
			    Cx[jpos] += coeff*Bx[B_ptr];
			}
		    }
		}
	    }
	}
	// Reset the work vector
	for (int k=Cp[ii]; k < len; k++) {
	    iw[Ci[k]] = -1;
	}

    }
    Cp[ncol] = len;	
	
    return 0;
}



int amub_col_plus_rho (int nrow, int ncol, int values,
		       c_float *Bx, c_int *Bi, c_int*Bp,
		       c_float *Ax, c_int *Ai, c_int*Ap,
		       c_float *Cx, c_int *Ci, c_int*Cp,
		       int nzmax,
		       c_float * rhoinv) {
    // Computes B*A + rho*I
    // The output matrix has dimensions nrow*ncol
    int ii, jj, jrow, jpos;
    c_int A_ptr, B_ptr;
    c_float coeff;
    int len = 0, diag_elem=0;

    for (int i=0;i<80;i++){
	Ci[i] = 0;
	Cx[i] = 0.0;
    }

    

    // Work array
    c_int  iw[ncol];
    for (int j=0; j < ncol; j++) iw[j] = -1;

    // Iterate over all columns in A
    for (ii=0; ii < ncol; ii++){
	Cp[ii] = len;	
	if (Ap[ii] ==  Ap[ii+1]){
	    // Only add the diagonal element
	    Ci[len] = ii;
	    Cx[len] = rhoinv[ii];
	    diag_elem++;	    
	    len++;
	}else{
	    for (A_ptr=Ap[ii]; A_ptr < Ap[ii+1]; A_ptr++){
		if (values){
		    coeff = Ax[A_ptr];
		}
		jj   = Ai[A_ptr];

		for (B_ptr=Bp[jj]; B_ptr < Bp[jj+1]; B_ptr++){
		    jrow = Bi[B_ptr]; // when jrow == ii, on diagonal, if jrow > ii, diagonal passed
		    if (jrow > ii && diag_elem==ii){
			Ci[len] = jrow;
			Cx[len] = rhoinv[ii];
			len++;
			diag_elem++;
			// must first add the diagonal element
		    }
		    jpos = iw[jrow];
		    // Has this already been added?
		    if (jpos == -1){
			if (len > nzmax){
			    return ii;
			}
			Ci[len] = jrow;
			iw[jrow]= len;
			if (values){
			    Cx[len]  = coeff*Bx[B_ptr];
			}
			if (jrow == ii){
			    Cx[len] += rhoinv[ii];
			    diag_elem++;
			}
		    
			len++;
		    } else if (values){
			Cx[jpos] += coeff*Bx[B_ptr];
		    }
		}
	    }
	}
	// Reset the work vector
	for (int k=Cp[ii]; k < len; k++) {
	    iw[Ci[k]] = -1;
	}

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

    // Work array
    c_int  iw[ncol];
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
		        // Has this already been added?
		        if (jpos == -1){
		            if (len > nzmax){
			            return ii;
		            }
		            Ci[len] = jrow;
		            iw[jrow]= len;
		            if (values){
			            Cx[len]  = coeff*Bx[B_ptr];
		            }
		            len++;
		        } else if (values){
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




int amub_col_plus_I (c_int nrow, c_int ncol, c_int start_col, c_int values,
	                 c_float *Bx, c_int *Bi, c_int*Bp,
	                 c_float *Ax, c_int *Ai, c_int*Ap,
	                 c_float *Cx, c_int *Ci, c_int*Cp,
		             c_int nzmax, double scale) {
    // Solves B*(A+I), where A is lower diagonal without diagonal entries.
    //printf("amub_col_plus_I(): Solves B*(A+I), where A is lower diagonal without diagonal entries. \n");
    int ii, jj, jrow, jpos;
    c_int A_ptr, B_ptr;
    c_float coeff;
    int len = 0;
    Cp[0] = 0;
    for (int ii=0;ii<start_col;ii++) Cp[ii+1]=0; 

    // Work array: For each column: keeps track of if row has nonzero value
    c_int  iw[nrow];
    for (int j=0; j < nrow; j++) iw[j] = -1;

    // If A+I = I: Copy the B matrix
    if (Ap[ncol] == 0){
	    copy_csc_mat_in_place(ncol, Bx, Bi, Bp, Cx, Ci, Cp);
	    return 0;
    }
    
    // Iterate over all columns in A matrix
    for (ii=start_col; ii < ncol; ii++){
	    //printf("\nCOLUMN = %i:  %i->%i\n", ii, Bp[ii], Bp[ii+1]);
	    // First add the diagonal element -----------------


	    for (B_ptr=Bp[ii]; B_ptr < Bp[ii+1]; B_ptr++){
	        jrow = Bi[B_ptr];
	        jpos = iw[jrow];
	        // Has this already been added? This should always be the first?...
	        if (jpos == -1){
		        if (len > nzmax) return ii;
		        Ci[len] = jrow;
		        iw[jrow]= len;
		        if (values) Cx[len]  = scale*Bx[B_ptr];
                //printf("row: %i\tCx[%i] = %f\n", jrow, len,scale*Bx[B_ptr] );
		        len++;
            } else if (values){
                //printf("row: %i\tCx[%i] += %f\n", jrow, jpos,scale*Bx[B_ptr] );
	            Cx[jpos] += scale*Bx[B_ptr];
            }
	    } // End of for loop
	    
	    // Now add rest of the values --------------------
	    for (A_ptr=Ap[ii]; A_ptr < Ap[ii+1]; A_ptr++){
	        if (values){
		        coeff = Ax[A_ptr];
	        }
	        jj   = Ai[A_ptr];
	        for (B_ptr=Bp[jj]; B_ptr < Bp[jj+1]; B_ptr++){
		        jrow = Bi[B_ptr];
		        jpos = iw[jrow];
		        // Has this already been added?
		        if (jpos == -1){
		            if (len > nzmax){
			            return ii;
		            }
		            Ci[len] = jrow;
		            iw[jrow]= len;
		            if (values){
			            Cx[len]  = coeff*scale*Bx[B_ptr];
                        //printf("row: %i\tCx1[%i] = %f\n",jrow,  len,coeff*scale*Bx[B_ptr] );
		            }
		            len++;
		        } else if (values){
		            Cx[jpos] += coeff*scale*Bx[B_ptr];
                    //printf("row: %i\tCx2[%i] = %f\n",jrow,  jpos,coeff*scale*Bx[B_ptr] );
		        }
	        } // End for loop
	    } // End for loop
	    // Reset the work vector
	    for (int k=Cp[ii]; k < len; k++) iw[Ci[k]] = -1;
        Cp[ii+1] = len;
    } // End for loop
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


void csc_to_csr_with_diag(const c_int n_col,
                          const c_int n_row, 
                          const c_int * Ap, 
                          const c_int * Aj, 
                          const c_float * Ax,
                          c_int * Bp,
                          c_int * Bi,
                          c_float * Bx){
    const c_int nnz = Ap[n_col];
    int i, jj, row, column, cumsum, temp, last, dest; 
    int diag_elem = 0;
    for (i=0; i<n_row;i++) Bp[i] = 1; // Every row had the diagonal element

    // Loop over all indices and add to the corresponding row 

    for (int n = 0; n < nnz; n++){            
        Bp[Aj[n]]++;
    }

    //cumsum the nnz per row to get Bp[]
    for(row = 0, cumsum = 0; row < n_row; row++){     
        temp  = Bp[row];
        Bp[row] = cumsum;
        cumsum += temp;
    }
    Bp[n_row] = nnz; 

    for(column = 0; column < n_col; column++){
    // Add the diagonal element first
    dest = Bp[column]; // B pointer 
    Bi[dest] = column;
    Bx[dest] = 1.0;
    Bp[column]++;
        for(jj = Ap[column]; jj < Ap[column+1]; jj++){
            row  = Aj[jj];
            dest = Bp[row]; // B pointer 

            Bi[dest] = column;
            Bx[dest] = Ax[jj];

            Bp[row]++;
        }
    }  

    for(row = 0, last = 0; row <= n_row; row++){
        temp  = Bp[row];
        Bp[row] = last;
        last    = temp;
    }
}


void subtract_csc_matrix_cols(int start_col, int stop_col, int n, 
                              c_int * Ap, c_int * Ai, c_float * Ax,
                              c_int * Bp, c_int * Bi, c_float * Bx ){
    // Subtract the columns between start_col and end_col of Bx to Ax.
    // Ax is upper diagonal, Bx is not. 
    c_int j, i, k, nrows;
    c_int row, count=0;

    c_int Rowp[n];
    
    c_int Rowi[n];
    c_int Rowx[n];

    
    c_int Cp[n];
    c_int Ci[Ap[n]+Bp[n]]; // allocate memory
    c_float Cx[Ap[n]+Bp[n]]; // allocate memory

    // First: loop over initial columns
    int A_ptr, B_ptr, Amax, Bmax;
    int A_row, B_row;
    Cp[0] = 0;
    for (j = 0; j < stop_col-start_col; j++) {
    A_ptr = Ap[j];
    Amax = Ap[j+1];
    B_ptr = Bp[start_col+j];
    Bmax = Bp[start_col+j+1];
    
    while(A_ptr < Amax || (B_ptr < Bmax && Bi[B_ptr]<=j)){

        if (A_ptr < Amax) {
        A_row = Ai[A_ptr];
        }else{
        A_row = n+1;
        }
        if (B_ptr < Bmax && Bi[B_ptr]<=j){
        B_row = Bi[B_ptr] ;        
        }else{ 
        B_row = n+1;
        }


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

    // Loop over rest of the columns
    for (j=stop_col-start_col; j<n; j++) {    
    for (A_ptr=Ap[j]; A_ptr<Ap[j+1]; A_ptr++){
        Ci[count] = Ai[A_ptr];
        Cx[count++] = Ax[A_ptr];
    }
    Cp[j+1] = count;    
    }


    for (i = 0; i <= n; i++){
    Ap[i] = Cp[i];
    }
    for (i=0; i<count; i++){
    Ai[i] = Ci[i];
    Ax[i] = Cx[i];
    }
    
    
}


void apldia (int nrow, int job,
             c_float * Ax, c_int * Ai, c_int * Ap, c_int * diag,
             c_float * Bx, c_int * Bi, c_int *Bp) {
//----------------------------------------------------------------------
// Adds a diagonal matrix to a general sparse matrix:  B = A + Diag 
//----------------------------------------------------------------------
// nrow    = integer. The row dimension of A
// job   = integer. job indicator. Job=0 means get array b only
//         (i.e. assume that a has already been copied into array b,
//         or that algorithm is used in place. ) For all practical 
//         purposes enter job=0 for an in-place call and job=1 otherwise
// 
//         Note: in case there are missing diagonal elements in A, 
//         then the option job =0 will be ignored, since the algorithm 
//         must modify the data structure (i.e. Bi, Bp) in this 
//         situation.
// Ax,Ai, Ap   = Matrix A in compressed sparse row format.
// diag = diagonal matrix stored as a vector dig(1:n)
// Bx, Bi, Bp    = resulting matrix B in compressed sparse row sparse format.
// iw    = integer work array of length n. On return iw will
//         contain  the positions of the diagonal entries in the 
//         output matrix. (i.e., a(iw(k)), Ai(iw(k)), k=1,...n,
//         are the values/column indices of the diagonal elements 
//         of the output matrix. ). 
//
//----------------------------------------------------------------
    int test, nnz;
    c_int iw[nrow];
    c_int ii, k,  j, k1, k2, ko;
//     copy integer arrays into Bx's data structure if required

    if (job != 0){
    nnz = Ap[nrow];
    for(k=0; k< nnz; k++){
        Bi[k] = Ai[k];
        Bx[k]  = Ax[k] ;
    }
    for (k=0; k< nrow; k++){
        Bp[k] = Ap[k];
    }
    }
//  get positions of diagonal elements in data structure.
    diapos (nrow,Ai,Ap,iw);
//  count number of holes in diagonal and add diag(*) elements to
//  valid diagonal entries.
    int icount = 0;
    for(j=1; j<nrow; j++){
    if (iw[j] == 0){
        icount = icount+1;
    }else{
        Bx[iw[j]] = Ax[iw[j]] + diag[j] ;
    }
    }
//  if no diagonal elements to insert return
    if (icount == 0) return;

//  shift the nonzero elements if needed, to allow for created 
//  diagonal elements. 
    ko = Bp[nrow+1]+icount;

//  copy rows backward
    for(ii=nrow; ii>= 0; ii--){ 

//     go through  row ii
    k1 = Bp[ii];
    k2 = Bp[ii+1]-1 ;
    Bp[ii+1] = ko;
    test = (iw[ii] == 0) ;
    for (k = k2; k>=k1; k--){ 
            j = Bi[k];
        if (test && (j < ii)){
        test = 0; 
        ko--;
        Bx[ko] = diag[ii];
        Bi[ko] = ii;
        iw[ii] = ko;
        }
            ko--;
        Bx[ko] = Ax[k];
            Bi[ko] = j;
    }
//     diagonal element has not been added yet.
    if (test){
        ko--;
        Bx[ko] =  diag[ii] ;
        Bi[ko] = ii;
        iw[ii] = ko;
    }
    }
    Bp[0] = ko;
    return;

}

void diapos(c_int n, c_int *ja, c_int * ia, c_int *idiag){
//-----------------------------------------------------------------------
// returns the positions of the diagonal elements of a
// sparse matrix a, ja, ia, in the array idiag.
//-----------------------------------------------------------------------
// n    = integer. row dimension of the matrix a.
// a,ja, ia = matrix stored compressed sparse row format. a array skipped.
//
// idiag  = integer array of length n. The i-th entry of idiag 
//          points to the diagonal element a(i,i) in the arrays
//          a, ja. (i.e., a(idiag(i)) = element A(i,i) of matrix A)
//          if no diagonal element is found the entry is set to 0.
//----------------------------------------------------------------------//
    int k;
    for(int i=0; i< n; i++) idiag[i] = 0;

    for(int i=0; i<n; i++){
    for(k=ia[i]; k<ia[i+1]; k++){
        if (ja[k] == i) idiag[i] = k;
    }
    }
}



int aplb (int nrow,int ncol,int values,
          c_float *Ax, c_int *Ai, c_int*Ap,
          c_float *Bx, c_int *Bi, c_int*Bp,
          c_float *Cx, c_int *Ci, c_int*Cp,
          int nzmax,int * iw){
// performs the matrix sum  C = A+B. 
    int len = 0;
    int jpos, jcol;
    Cp[0] = 0; 
    for(int j=0; j < ncol; j++){
    iw[j] = 0;
    }

    for(int ii=0; ii <  nrow; ii++){
    
    for (int ka=Ap[ii]; ka <  Ap[ii+1]-1; ka++ ){
        jcol = Ai[ka];
            if (len > nzmax) return ii;
            Ci[len] = jcol;
        if (values) Cx[len]  = Ax[ka];
            iw[jcol]= len;
            len++;
    }
    for(int kb=Bp[ii]; kb <Bp[ii+1]-1; kb++){
          jcol = Bi[kb];
          jpos = iw[jcol];
          if (jpos == 0){
          if (len > nzmax) return ii;
          Ci[len] = jcol;
          if (values) Cx[len]  = Bx[kb];
          iw[jcol]= len;
          len++;
          }else if (values){
          Cx[jpos] = Cx[jpos] + Bx[kb];
          }
    }
    for(int k=Cp[ii]; k< len; k++){
        iw[Ci[k]] = 0;
    }
    Cp[ii+1] = len;
    }
    return 0;
}


csc* transpose(const csc *A){
    csc *B = csc_spalloc(A->n, A->m, A->p[A->n], 1, 0);
    return B;
}

