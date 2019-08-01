//============================================================================
// Name        : krylov_shur.cpp
// Author      : Guido Polles
// Version     :
// Copyright   :
// Description :
//============================================================================

#include <mpi.h>

#include <slepceps.h>

#include <cstdio>
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

#include "mylib/sparsematrix.hpp"

#if SLEPC_VERSION_MAJOR <= 3
#if SLEPC_VERSION_MINOR < 6 
#define EPS_ERROR_RELATIVE -1
#define MatCreateVecs(a, b, c) MatGetVecs(a, b, c)
#define EPSComputeError(a, b, c, d) EPSComputeResidualNorm(a, b, d)
#endif
#endif

typedef double real_t;

using namespace std;

static char help[] = "Usage: krylov_shur [slepc options] <matrix filename> <number of eigenvectors>\n";

int krylov_shur(mylib::SymmetricSparseMatrix<real_t>& input_matrix,
		int n_eigvec,
		vector< vector<real_t> >& eigvec,
		vector<real_t>& eigval)
{


  // Get problem size and row lenghts for faster allocation.
  PetscInt dim = input_matrix.size();
  PetscInt max_row_nnz = 0;
  PetscInt row_nnz[dim];
  for (int i = 0; i < dim; ++i){
    row_nnz[i] = input_matrix.getRow(i).size();
    max_row_nnz = max(row_nnz[i],max_row_nnz);
  }

  PetscErrorCode  ierr;
  int fargc=1;
  char** fargv = (char**) malloc(1*sizeof(char*));
  fargv[0] = (char*) malloc(2*sizeof(char));
  strcpy(fargv[0]," ");
  ierr = SlepcInitialize(&fargc,&fargv,(char*)0,help);CHKERRQ(ierr);

  // Set up the matrix
  Mat A;
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,dim,dim);CHKERRQ(ierr);
  ierr = MatSetType(A, MATSEQSBAIJ);CHKERRQ(ierr); // sequential symmetric block sparse matrix
  ierr = MatSeqAIJSetPreallocation(A, max_row_nnz, row_nnz);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
  ierr = MatSetOption(A,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_FALSE); // throws error if inserting on the wrong part

  // fill matrix
  PetscInt start,end;
  MatGetOwnershipRange(A,&start,&end);
  for(PetscInt row = start; row < end; row++) {
    int nnz=row_nnz[row];
    PetscInt    idxj[nnz];
    PetscScalar vals[nnz];
    for (int j = 0; j < nnz; ++j){
      idxj[j] = input_matrix.getRow(row).at(j).columnIndex;
      vals[j] = input_matrix.getRow(row).at(j).value;
    }
    MatSetValues(A,1,&row,nnz,idxj,vals,INSERT_VALUES);
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  // Solve the EP with SLEPc

  EPS eps;
  Vec eigen_vec_re, eigen_vec_im;
  PetscScalar eigen_val_re, eigen_val_im;
  PetscInt nconv,its;
  PetscInt nev = n_eigvec;

  cout << "  Call to Eigensolver:" << endl<< flush;

  ierr = MatCreateVecs(A,0,&eigen_vec_im);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,0,&eigen_vec_re);CHKERRQ(ierr);

  cout << "  > Preparing eigenproblem..." << endl<< flush;
  ierr = EPSCreate( PETSC_COMM_WORLD, &eps );CHKERRQ(ierr);
  ierr = EPSSetOperators( eps, A, 0 );CHKERRQ(ierr);
  ierr = EPSSetProblemType( eps, EPS_HEP );CHKERRQ(ierr);
  ierr = EPSSetDimensions( eps, nev, PETSC_DEFAULT, PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = EPSSetType( eps, EPSKRYLOVSCHUR );CHKERRQ(ierr);

  ierr = EPSSetWhichEigenpairs( eps, EPS_SMALLEST_MAGNITUDE);CHKERRQ(ierr);

  ierr = EPSSetFromOptions( eps );CHKERRQ(ierr);
  cout << "  > Solving eigenproblem..." << endl<< flush;
  #if SLEPC_VERSION_MAJOR <= 3
  #if SLEPC_VERSION_MINOR < 8 
    ierr = EPSMonitorSet(eps,EPSMonitorFirst,0,0);CHKERRQ(ierr);
  #endif
  #endif
  ierr = EPSSolve( eps );CHKERRQ(ierr);

  cout << "  > Done." << endl<< flush;

  // get results and print
  ierr = EPSGetConverged( eps, &nconv );CHKERRQ(ierr);
  ierr = EPSGetIterationNumber( eps, &its );CHKERRQ(ierr);

  cout << "  #### First 15 eigenvalues ####   "<< endl;
  cout << "  "  << nconv << "  converged eigenvectors in " << its << " iterations." << endl<< flush;
  cout << "\t" << setw(3) << "#" << setw(15) << "Re(l)" << setw(15) << "Im(l)"
       << setw(15) << "Residual" << endl<< flush;

  for (PetscInt j=0; j < nconv && j < 15 ; ++j) {
    PetscScalar res_norm;
    ierr = EPSGetEigenpair( eps, j,&eigen_val_re,&eigen_val_im, eigen_vec_re, eigen_vec_im);CHKERRQ(ierr);
    ierr = EPSComputeError( eps, j, EPS_ERROR_RELATIVE, &res_norm );CHKERRQ(ierr);
    cout << "\t" << setw(3) << j+1 << setw(15) << eigen_val_re
        << setw(15) << eigen_val_im << setw(15) << res_norm << endl<< flush;
  }

  eigvec.resize(nconv, vector<real_t>(dim) );
  eigval.resize(nconv,0.0);

  // save outputs
  PetscScalar* array;
  for (PetscInt j=0; j < nconv; ++j) {
      PetscReal res_norm;
      ierr = EPSGetEigenpair( eps, j,&eigen_val_re,&eigen_val_im,eigen_vec_re, eigen_vec_im);CHKERRQ(ierr);
      ierr = EPSComputeError( eps, j, EPS_ERROR_RELATIVE, &res_norm );CHKERRQ(ierr);

      // output: eigenvalues
      VecGetArray(eigen_vec_re,&array);
      eigval[j] = eigen_val_re;
      for (PetscInt k=0; k<dim; ++k) {
          eigvec[j][k] = array[k];
      }
  }

  // destroy workspace
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&eigen_vec_im);CHKERRQ(ierr);
  ierr = VecDestroy(&eigen_vec_re);CHKERRQ(ierr);
  ierr = SlepcFinalize();

  free(fargv[0]);
  free(fargv);

  return nconv;
}
