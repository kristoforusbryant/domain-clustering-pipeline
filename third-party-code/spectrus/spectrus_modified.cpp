// LAUNCH WITH:
// ./main.out <trajectory file> 
// OR:
// ./main.out <pdb file> <num normal modes> <num decimation>:
//
// INPUT (launch options)
// - a trajectory file (or a set of conformations), in .xyz form
// - alternatively, a pdb file for the ENM. In this case, the number of modes and the decimation number must be
//   provided as well, and the pdb file _has_to_be_prepared_already_decimated_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <exception>
#include <map>
#include <math.h>
#include <numeric>
#include <ctime>
#include <queue> 

#include "mylib/mylib.h"

using namespace std;
using namespace mylib;

#include "krylov_shur.hpp"		/* Routine for solving eig. problem for sparse matrices, using SLEPc */
#include "my_malloc.hpp"  		/* Routines for memory allocation */
#include "k_medoids.hpp"		/* Routine which performs k-medoids on a set of points */

#define VERSION_MAJOR 0
#define VERSION_MINOR 8

typedef std::map<std::string,std::string> dict;

dict read_parameters(const char* fname);
vector<int> top_k_index(vector<double> row, int k); 

void compute_Q_dist_matrix(int N_ATOMS, int q_max, double** x, double **mtrx_sigma, double *sigma_max_p);
void find_sign_groups(int N_ATOMS, int Q, double** x, int* groups);
void print_matrix(vector<vector<double> > matrix, char filename[256], int n_atoms);

template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T> >& v) {
    std::size_t total_size = 0;
    for (size_t i = 0; i < v.size(); ++i)
      total_size += v[i].size();
    std::vector<T> result;
    result.reserve(total_size);
    for (size_t i = 0; i < v.size(); ++i)
      result.insert(result.end(), v[i].begin(), v[i].end());
    return result;
}



/// MAIN ///

int main(int argc, char* argv[])
{

  printf("MODIFIED FROM SPECTRUS v. %d.%d\n TO TAKE DIRECT COUPLING MATRIX",VERSION_MAJOR,VERSION_MINOR);
  if(argc < 3 || argc > 3){
    cerr << "  Usage: ./main.out direct_coupling_matrix" << endl;
    return 1;
  }

  char* fname = argv[1]; // name of direct coupling matrix 
  string n_atoms_str = argv[2]; 
  int n_atoms = atoi(n_atoms_str); // in this modified context, this is number of residues  

  if(argc==3) {
    cout << "  Distance fluctuation matrix computed from the file:" << endl;
    cout << "  " << argv[1] << endl << endl;
  }

  int max_dom, min_dom, kmed_iter, keep_k;
  double temp;
  string line;
  FILE *fp;
  char test_filename[256], command[256];

  // read parameters from file
  dict params = read_parameters("INPUT_PARAMS.DAT");
  try{
    // n_atoms = atoi(params.at("N_ATOMS")); // means n of residues in the modified context
    max_dom = atoi(params.at("N_DOM_MAX")); 
    min_dom = atoi(params.at("N_DOM_MIN")); 
    kmed_iter = atoi(params.at("KMED_ITER")); 
    keep_k = atoi(params.at("KEEP_K"));
  } catch (exception& e){
    cerr << "  Some options missing in parameter file." << endl;
    cerr << e.what() << endl;
    return 1;
  } 
  

  // read direct coupling matrix
  // TODO: can use some error handling to ensure the size of the matrix
  vector< vector<double> > sim_matrix_raw(n_atoms, vector<double>(n_atoms));
  vector< vector<double> > sim_matrix(n_atoms, vector<double>(n_atoms, 0.));

  // setting up raw sim matrix
  ifstream fin(fname); 
  if (fin.is_open())
  {
  for (int i=0; i < (n_atoms -1); ++i){
     for (int j=(i + 1); j < n_atoms; ++j){
        getline(fin, line, '\n');
	temp = stod(line); 
 
	if (temp < 0.) { 
	  sim_matrix_raw[i][j] = 0.;
	  sim_matrix_raw[j][i] = 0.;
	}
	else{
	  sim_matrix_raw[i][j] = temp;
	  sim_matrix_raw[j][i] = temp; 
	}
     }
  }
  }
  fin.close();

  // sprintf(test_filename, "test/results/%s_parsed_mtx.txt", fname);
  // print_matrix(sim_matrix_raw, test_filename, n_atoms);  
 
  vector<int> top_k(keep_k);
  vector<double> temp_vec(n_atoms, 0.);
  for (int i=0; i < n_atoms; ++i){
     // reset temporary array to 0
     for (int i=0; i < n_atoms; ++i){temp_vec[i] = 0.;}
     // keep the top keep_k interactions
     top_k = top_k_index(sim_matrix_raw[i], keep_k);
     //for (int j=0; j < keep_k; ++j) printf(" %d,", top_k[j]); 
     //printf("\n");
     
     // populate temporary vector with top keep_k interactions 
     for (int j=0; j < keep_k; ++j){temp_vec[top_k[j]] = sim_matrix_raw[i][top_k[j]];}
     // populate temporary vector with consecutive elements 
     if(i < (n_atoms-1)){
	if (sim_matrix_raw[i][i+1] < 0.000001) {temp_vec[i+1] = 0.000001;}
	else{ temp_vec[i+1] = sim_matrix_raw[i][i+1]; }
     }
     if(i > 0){
	if (sim_matrix_raw[i][i-1] < 0.000001) {temp_vec[i-1] = 0.000001;}
	else{ temp_vec[i-1] = sim_matrix_raw[i][i-1]; }
     }
     //for (int j=0; j < n_atoms; ++j) printf(" %f,", temp_vec[j]);
     //printf(" \n"); 
     
     // assign temporary vector value to similarity matrix row 
     for (int j=0; j < n_atoms; ++j) { 
	if (temp_vec[j] >= 0.00000001) {
	  sim_matrix[i][j] = temp_vec[j];
	  sim_matrix[j][i] = temp_vec[j]; 
        }
     }
  }
 
  // sprintf(test_filename, "test/results/%s_sim_mtx.txt", fname);
  // print_matrix(sim_matrix, test_filename, n_atoms);

  // compute the sum of each row (degree matrix)
  vector<double> diag(n_atoms,0.0);
  for (int i=0; i<n_atoms; ++i){
    for (int j=0; j<n_atoms; ++j){
      diag[i] += sim_matrix[i][j];
    }
  }
  
  // check diagonal consistency
  for (int i=0; i<n_atoms; ++i){
    if(diag[i] <= 0.) {
      fprintf(stderr, "# ERROR! Null row detected in the sparse matrix. \n");
      fprintf(stderr, "# ERROR! Row %d has only zeroes. Stopped. \n", i);
      exit(1);
    }
  }
  
  // build the SYMMETRIC LAPLACIAN L_sym = I - D^{-1/2} S D^{-1/2}
  mylib::SymmetricSparseMatrix<double> laplacian_matrix(n_atoms);
  for(int i=0; i<n_atoms; ++i) {
    laplacian_matrix.insert( i, i, 1. );
    for(int j=0; j < n_atoms; ++j){
      if (j < i && sim_matrix[i][j] > 0.000001) laplacian_matrix.insert(i,j, -sim_matrix[i][j] / sqrt(diag[i]*diag[j])); 
    }
  }
 
  
  // create "/results" directory
    // system("if [ ! -d results ]; then mkdir results; fi");
  
  // spectral analysis of the Laplacian of the similarity matrix.
  std::vector< std::vector<double> > eigvecs;
  std::vector<double> eigvals;
  int n_requested_eigv = max_dom + 2; // why do we request 2 more eigenvalues than the maximum number of domain? 
  if (n_requested_eigv > n_atoms) n_requested_eigv = n_atoms;
  int n_converged_eigv = krylov_shur(laplacian_matrix, n_requested_eigv, eigvecs, eigvals);

  if (n_converged_eigv < n_requested_eigv) {
    cerr << "# ERROR: diagonalization routine was not able "
            "to converge all the requested eigenvalues" << endl;
    return 1;
  }
  
  // create results directory
  sprintf(command, "rm -rf results_%s; mkdir results_%s", fname, fname); 
  system(command);

  // check: print the Laplacian eigenspace 
  sprintf(command, "mkdir results_%s/eigenspace", fname); 
  system(command);
  for(int l=0; l<max_dom; ++l){
    char stringa[100];
    sprintf(stringa, "results_%s/eigenspace/eigenvector_%d.dat", fname, l);
    fp = fopen(stringa, "w");
    for(int j=0; j<n_atoms; ++j)
      fprintf(fp,"%d %le \n", j, eigvecs[l][j]);
    fclose(fp);
  }
  sprintf(test_filename, "results_%s/eigenspace/eigenvalues.dat", fname); 
  fp = fopen(test_filename, "w");
  for(int l=0; l<n_atoms; ++l) fprintf(fp,"%d %le \n", l, eigvals[l]);
  fclose(fp);

  // check the number of null eigenvalues
  int nNullEigval = 0;
  while (eigvals[nNullEigval] < 1e-15) ++nNullEigval;

  // if nNullEigval>1, the graph is disconnected and spectral clustering
  // is not well defined anymore (how to deal with degenerate eigenvectors?).
  // In that case, just output the clusterization into the disconnected parts.
  if( nNullEigval > 1 ) {
    fprintf(stderr,"# WARNING: %d disconnected parts have been found:\n", nNullEigval);
    fprintf(stderr,"           The subdivision into these disconnected parts will be returned instead. \n");
    fprintf(stderr,"           You can try again the domain decomposition by retaining only one single part. \n");
    min_dom = nNullEigval;
    max_dom = min_dom;
  }

  /************ cycle on Q (number of domains) *************/
  
  double **points_on_qsphere = d2t(n_atoms,max_dom);
  int *groups = i1t(n_atoms);
  int *label_kmed = i1t(n_atoms);
  double q_score_median, q_score_mean;
  char filename[256];

  for(int q=min_dom; q<=max_dom; ++q) {
    
    // normalization of the points on the unitary Q-dimensional sphere
    int n_zero_norm = 0;
    for(int i=0; i<n_atoms; ++i) {
      double norm = 0.;

      // CARE: now for eigenvectors the first index is the eigenvector
      // ordinal number, the second is the atom number
      for(int j = 0; j < q; ++j) {
        norm += pow(eigvecs[j][i],2);
      }
      norm = sqrt(norm);
      if(norm == 0.) { 
	++n_zero_norm;
	continue;
      }
      for(int j = 0; j < q; ++j) {
        points_on_qsphere[i][j] = eigvecs[j][i]/norm;
      }
    }

    // if we have zero norm points throws a warning and jump to next Q
    if(n_zero_norm) {
      fprintf(stderr, "# WARNING: Some eigenvectors of the kernel space are missing! \n");
      fprintf(stderr, "# WARNING: K-medoids for Q = %d will be skipped. \n",q);
      fprintf(stderr, "# WARNING: (number of atoms with all components equal to zero = %d.) \n", n_zero_norm);
      continue;
    }
    

    // print the coordinates of the points on the Q-hypersphere
    fp = fopen("k_medoids_coordinates.dat","w");
    for(int i = 0; i < n_atoms; ++i) {
      for(int j = 0; j < q; ++j) {
	fprintf(fp,"%le ",points_on_qsphere[i][j]);
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
          
    // find atoms with the same sequence of signs in the eigvectors
    find_sign_groups(n_atoms, q, points_on_qsphere, groups);

    time_t now = time(0);
    tm* localtm = localtime(&now);
    cout << "  kmedoids started on: " << asctime(localtm) << endl;
    k_medoids(n_atoms, q, q, 1, q, kmed_iter);
    now = time(0);
    localtm = localtime(&now);
    cout << "  kmedoids finished  on: " << asctime(localtm) << endl << endl;
    
    // import the cluster profile (labels)
    sprintf(filename, "k_medoids_output/output_k_%d/cluster_profile.dat",q);
    fp = fopen(filename,"r");
    for(int i = 0; i < n_atoms; ++i) {
      fscanf(fp, "%*d \t%d \n", &label_kmed[i]);
      if(label_kmed[i]<0) label_kmed[i] = -1;
    }
    fclose(fp);

    // import the quality score
    double prefactor, prefactor2;
    fp = fopen("k_medoids_output/quality_score.dat","r");
    fscanf(fp,"%*d %le %le %le %le\n", &q_score_median, &q_score_mean, &prefactor, &prefactor2);
    fclose(fp);
 
    // move the directory into results
    sprintf(command, "mv k_medoids_output results_%s/k_medoids_output_%d", fname, q);
    system(command);
    system("rm -f k_medoids_coordinates.dat k_medoids_sign_groups.dat");
    
    /*** PRINT THE RESULTS ***/
    
    sprintf(filename, "results_%s/quality_score.dat", fname);
    fp = fopen(filename,"a");
    fprintf(fp, "%4d %le %le %le %le\n", q, q_score_median, q_score_mean, prefactor, prefactor2);
    fclose(fp);

    // prints the final clusterization (in .pdb format)
    sprintf(filename, "results_%s/final_labels_kmed-%d.dat", fname, q);
    FILE *fp_labels = fopen(filename,"w");
    for(int j = 0; j < n_atoms; ++j) ++label_kmed[j];
    for(int j = 0; j < n_atoms; ++j) {
      int i = j + 1;
      fprintf(fp_labels, "%d %d \n", i,label_kmed[j]);
    }
    fclose(fp_labels);
    
  }

  free_d2t(points_on_qsphere);
  free_i1t(label_kmed);
  free_i1t(groups);

  return 0;
}


/******************************/
/***** SUBROUTINES ************/
/******************************/


int read_last_frame(const char* fname, int N_ATOMS, 
                    std::vector<mylib::Vector3d>& pos_atom, int* n_frames_p){  // very dumb way to do it, but whatever

  FILE *fp = fopen(fname,"r");
  int n_frames,exit_flag=0;

  for(n_frames=0; ; ++n_frames) {
    // read one conformation (or frame) from the file
    for(int i=0; i<N_ATOMS; ++i) {
      if(fscanf(fp,"%lf %lf %lf",&pos_atom[i][0],&pos_atom[i][1],&pos_atom[i][2] )==EOF) {
	// check for a correct EOF
	if(i==0) {
	  exit_flag = 1;
	  break;
	}
	else {
	  fprintf(stderr, "# ERROR in reading conformation file!\n");
	  exit(1);
	}
      }
    }
    if(exit_flag) break;
  }
  fclose(fp);
  *n_frames_p = n_frames;

  return 0;
}


/******************************/
/******************************/


dict read_parameters(const char* fname){

  dict params;

  ifstream fin(fname);
  if (!fin.good()){
    cerr<< "# Cannot read parameter file" << endl;
    exit(1);
  }
  string key,val;
  while(!(fin >> key >> val).fail()) params[key]=val;

  return params;
}


vector<int> top_k_index(vector<double> row, int k){
  vector<int> top_k(k); 
  priority_queue < pair<double, int> > q; 
  for (int i = 0; i < row.size(); ++i){
    q.push(pair<double,int>(row[i],i));
  } 
  for (int i=0; i < k; ++i){
    top_k[i] = q.top().second;   
    q.pop();
  }
  return top_k;
}

void print_matrix(vector <vector<double> >  matrix, char filename[256], int n_atoms) { 
  FILE *fp = fopen(filename,"w");
  for(int i=0; i<n_atoms; ++i) {
    for(int j=0; j < n_atoms; ++j) fprintf(fp,"%f,", matrix[i][j]);
    fprintf(fp, "\n");  
  }
  fclose(fp);
}


/******************************/


void compute_Q_dist_matrix(int N_ATOMS, int q_max, 
			   double** x, double **mtrx_sigma, double *sigma_max_p) {

  double sigma_max,tmp;

  printf("  DISTANCE MATRIX FROM TOP %d EIGENVECTORS \n", q_max);
  printf("  ---------------------------------------- \n");

  sigma_max = 0.;
  for(int i=0; i<N_ATOMS; ++i) {
    for(int j=i+1; j<N_ATOMS; ++j) {
      tmp = 0.;
      // distance on a unit hypersphere
      for(int q=0; q<q_max; ++q) tmp+=x[i][q]*x[j][q];
      if(tmp>1.) tmp=0.;
      else if(tmp<-1.) tmp=3.1415926535;
      else tmp = acos(tmp);
      // OR euclidean distance
//       for(int q=0; q<q_max; ++q) tmp+=pow(x[i][q]-x[j][q],2);
//       tmp = sqrt(tmp);
      mtrx_sigma[i][j] = tmp;
      if(tmp>sigma_max) sigma_max=tmp;
    }
  }  

  // outputs
  *sigma_max_p = sigma_max;

  printf("  Computed. \n\n");
}


/******************************/


void find_sign_groups(int N_ATOMS, int Q, double** first_eigvecs, int* groups) {

  for(int i=0; i<N_ATOMS; ++i) groups[i]=-1;
  for(int i=0; i<N_ATOMS; ++i) {
    if(groups[i]==-1) groups[i]=i;
    else continue;
    for(int j=i+1,in_the_same_group; j<N_ATOMS; ++j) {
      if(groups[j]!=-1) continue;
      in_the_same_group = 1;
      for(int q=1; q<Q; ++q) {
        if((first_eigvecs[i][q]*first_eigvecs[j][q]) < 0.) {
          in_the_same_group = 0;
          break;
        }
      }
      if(in_the_same_group) groups[j]=i;
    }
  }
  FILE* fp = fopen("k_medoids_sign_groups.dat","w");
  for(int i=0; i<N_ATOMS; ++i) fprintf(fp,"%d\n",groups[i]);
  fclose(fp);

}

