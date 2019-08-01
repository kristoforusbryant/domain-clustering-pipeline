//============================================================================
// Author       : Luca Ponzoni, Cristian Micheletti
// Version      : 1.0
// Description  : k-medoids algorithm
// Input        : k_medoids_coordinates.dat, 
// Date		: Jan 2015
// References	: Hastie et al., "The elements of statistical learning", 
//		  pg. 468, Springer
//============================================================================

/* The input file 'k_medoids_coordinates.dat' contains the points' coordinates
(which have to be normalized vectors!),
instead of the distance matrix (in order to save memory on the disk). 
The distance matrix (i.e. the matrix of the distances on the 
hypersphere) is then computed inside here. */

/* The initial cluster profile for each iteration is not randomly chosen, 
but it's chosen by taking into account the "sign groups" imported from file 
(so that, at least at the beginning, two cluster centres don't belong
to the same group). 
Moreover, after each iteration, some "merge and split" steps are taken, as far as
the dissimilarity keep decreasing.*/


#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include "prefactor.hpp"
#include "prefactor2.hpp"

#include "my_malloc.hpp"
#include "mylib/matrix.hpp"

#define NMAXelements (1000)

using namespace std;

void read_data(mylib::SymmetricMatrix<float>& dist_matr, char labels[][20], int n_elements);
void read_coordinates(mylib::SymmetricMatrix<float>& dist_matr, char labels[][20], int n_elements, int Q);
void output_results(mylib::SymmetricMatrix<float>& dist_matr,int *cluster_centres, int *cluster_profile, 
                    int n_elements, int k, char labels[][20], double dissim);
void dump_regrouped_distance_matrix(mylib::SymmetricMatrix<float>& dist_matr,int *cluster_profile,  
                                    int n_elements, int k);
void copy_cluster_centres_a_to_b(int *cluster_centres_a, int *cluster_centres_b, int k);
void assign_elements_to_centres_and_rerank(mylib::SymmetricMatrix<float>& dist_matr,int *cluster_centres, 
                                           int *cluster_profile, int n_elements, int k);
void assign_new_centres(mylib::SymmetricMatrix<float>& dist_matr,int *cluster_centres, int *cluster_profile, 
                        int n_elements, int k);
int unchanged_assignment(int *cluster_centres, int *new_cluster_centres, int k);
double within_cluster_dissimilarity(mylib::SymmetricMatrix<float>& dist_matr,int *cluster_centres, 
                                    int * cluster_profile, int n_elements, int k);
double inter_cluster_dissimilarity(mylib::SymmetricMatrix<float>& dist_matr,int *cluster_centres, 
                                   int * cluster_profile, int n_elements, int k, double *C_inter, double *C_ratio, double *C_ratio2, 
                                   double *C_ratio3, double *C_ratio4);
double minimize_inter_cluster_dissimilarity(mylib::SymmetricMatrix<float>& dist_matr,int *cluster_centres, 
                                            int * cluster_profile, int n_elements, int k, double *C_inter, double *C_ratio, double *C_ratio2, 
                                            double *C_ratio3, double *C_ratio4, double *threshold);
void random_within_cluster_dissimilarities(mylib::SymmetricMatrix<float>& dist_matr,int *cluster_centres, 
                                           int n_elements, int k, double *mean_C, double *delta_C);
void random_improvement_of_cluster_centres(mylib::SymmetricMatrix<float>& dist_matr,int *cluster_centres, 
                                           int n_elements, int k, int n_trials);
void random_initial_cluster_centres(int *cluster_centres, int n_elements,int k, int *sign_groups);
int are_cluster_centres_in_different_sign_groups(int *cluster_centres, int k, int *sign_groups);
void MergeNSplit(int n_elements, int k, mylib::SymmetricMatrix<float>& dist_matr, int *cluster_profile, 
                 int *cluster_centres, int *sign_groups);
double compute_median_2(std::vector<double>& vec);
bool does_file_exist(const char *fileName);

long int idum= (-1);
double ran2 (long *idum);


int k_medoids(int n_elements, int kmin, int kmax, int kstep, int n_coord, int loopmax) {

// performs k-medoids into 'kmin' to 'kmax' clusters, with a pace 'kstep'. 
// 'loopmax' = # of iterations (i.e. choices of the initial medoids).
// 'n_coord' = # of point coordinates as in 'k_medoids_coordinates.dat'


  int *cluster_profile;
  int *new_cluster_centres, *n_cluster_elements, *best_cluster_centres;
  int i, k, iteration, major_loop;
  //int best_loop = -1;
  int do_MergeNSplit;
  double C, mean_C, delta_C, Cbest, C_old;
  double C_inter, C_ratio, C_ratio2, C_ratio3, C_ratio4;
  double C_ratio4_mean, C_ratio4_var, threshold;
  int loopmax_cont, any_bad_centre;
  char string[100];
  FILE *fp;
  
  
  /* check for input files and directories */
  struct stat st;
  if(stat("k_medoids_output", &st) == 0) {
    fprintf(stderr, "# ERROR: directory 'k_medoids_output' already exists. Exiting.\n\n");
    exit(1);
  }
  else system("mkdir k_medoids_output");
  if(!does_file_exist("k_medoids_coordinates.dat")) {
    fprintf(stderr, "# ERROR: file 'k_medoids_coordinates.dat' not found.\n\n");
    exit(1);
  } 
  if(!does_file_exist("k_medoids_sign_groups.dat")) {
    fprintf(stderr, "# WARNING: file 'k_medoids_sign_groups.dat' with initial groups not found.\n\n");
  }


  /* allocations and definitions */
  mylib::SymmetricMatrix<float> distance_matrix(n_elements);
  cluster_profile = i1t(n_elements);
  new_cluster_centres = i1t(n_elements);
  n_cluster_elements = i1t(n_elements);
  best_cluster_centres = i1t(n_elements);
  int cluster_centres[n_elements];
  char labels[n_elements][20];
  
  
  /* read from file the point coordinates, then compute the distance matrix */
  read_coordinates(distance_matrix, labels, n_elements, n_coord);
  
  /* in case you want to read from file the distance matrix directly, 
  do something like this: */
//   read_data(distance_matrix, labels, n_elements);
  /* NB: ONLY SYMMETRIC MATRICES! */
  /* input file: 'k_medoids_distances.dat' */
  
  /* read sign groups from file, if available */
  int* sign_groups = i1t(n_elements);
  fp = fopen("k_medoids_sign_groups.dat","r");
  if(fp == NULL) {
    fprintf(stderr,"# WARNING: file 'k_medoids_sign_groups.dat' with initial groups not found.\n");
    sign_groups[0] = -1;
  } 
  else {
    for(i=0; i<n_elements; ++i)
      fscanf(fp,"%d",&sign_groups[i]);
    fclose(fp);
  }

  /* output file with various order parameters */
  fp = fopen("k_medoids_output/order_parameters.dat","w");
  fprintf(fp,"   k \tC \t\t|C-<C>|/DeltaC \t<C> \t\tDeltaC \t\tC_inter \tC_inter/C \t");
  fprintf(fp,"C_ratio \tC_ratio2 \tC_ratio3 \tC/C_inter \tC_ratio4 \tloopmax_cont \tC_ratio5 \n");
  fclose(fp);
  
  
  for(k=kmin; k <= kmax; k+=kstep){

    sprintf(string, "mkdir k_medoids_output/output_k_%d", k);
    system(string);

    /* initial choice */
    for(i=0; i<k; ++i) {
      cluster_centres[i] = i;
    }
    copy_cluster_centres_a_to_b(cluster_centres,best_cluster_centres,k);
    
    Cbest = -1.;
    C_old = Cbest;
    do_MergeNSplit = 0;
    
    C_ratio4_mean = 0.;
    C_ratio4_var = 0.;
    loopmax_cont = 0;

    /* Collect statistics on random scores */
//     printf("K: %3d Random reference for within_cluster_dissimilarity %lf +- %lf\n",k,mean_C,delta_C);
    mean_C = 999.;
    delta_C = -666.;



    /********** Random optimization of the WITHIN cluster dissimilarities **********/
    
    for(major_loop=0; major_loop<loopmax; ++major_loop){

      if(do_MergeNSplit == 1) {
       MergeNSplit(n_elements,k,distance_matrix,cluster_profile,cluster_centres,sign_groups);
       assign_new_centres(distance_matrix,cluster_centres,cluster_profile,n_elements,k);
     }
     else if(sign_groups[0] != -1) {
        /* pick up k random medoids, in different sign groups */
      random_initial_cluster_centres(cluster_centres,n_elements,k,sign_groups);
        /* fill the cluster profile */
      assign_elements_to_centres_and_rerank(distance_matrix,cluster_centres,cluster_profile,n_elements,k);
    }
    else {
	/* (no info about sign groups) */
	/* do a random exploration of possible clusters and collect statistics about mean
	and spread of the order parameter: within_cluster_dissimilarities 
	The latter corresponds to the sum of distances of each member from the nearest centre.
	After calling the routine, the cluster_centres have been randomised */
     random_within_cluster_dissimilarities(distance_matrix,cluster_centres,n_elements,k,&mean_C,&delta_C);
	/* Now we optimise the centres by substituting centres at random and keeping the
	results if it improves the order parameter */
//	  random_improvement_of_cluster_centres(distance_matrix,cluster_centres,n_elements,k,((n_elements*k)/100));
   }

   any_bad_centre = 0;
   for(int z=0; z<k; ++z) {
     if(cluster_profile[cluster_centres[z]]!=z) {
       any_bad_centre=1; 
       break;
     }
   }
   if(any_bad_centre == 0) {
     ++loopmax_cont;
     C = inter_cluster_dissimilarity(distance_matrix,cluster_centres,cluster_profile,n_elements,k,
                                     &C_inter,&C_ratio,&C_ratio2,&C_ratio3,&C_ratio4);
     C_ratio4_mean += C/C_inter;
     C_ratio4_var += (C/C_inter)*(C/C_inter);
   }


      /* Now we iteratively assign centres and cluster members until convergence */

      /********** Iterative optimization of the WITHIN cluster dissimilarities **********/

   iteration=0; copy_cluster_centres_a_to_b(cluster_centres,new_cluster_centres,k);
   do {
    copy_cluster_centres_a_to_b(new_cluster_centres,cluster_centres,k);      
    assign_elements_to_centres_and_rerank(distance_matrix,cluster_centres,cluster_profile,n_elements,k);
    assign_new_centres(distance_matrix,new_cluster_centres,cluster_profile,n_elements,k);

    C = within_cluster_dissimilarity(distance_matrix,new_cluster_centres,cluster_profile,n_elements,k);
    ++iteration;
  } while(unchanged_assignment(cluster_centres,new_cluster_centres,k)==0);


  if(C<(Cbest-1.0e-8) || Cbest<0.) {
    printf("  K = %d. Loop: %d. Best Cluster dissimilarity: %le \n",k,major_loop,C);
    Cbest = C;
    copy_cluster_centres_a_to_b(new_cluster_centres,best_cluster_centres,k);
    //best_loop = major_loop;
  }

      /* Now we have certainly ended up in a minimum of the order parameter. Let's repeat the
         cycle so to explore other possible minima, by either MergeNSplit or random initialization */

  if(C<C_old || do_MergeNSplit==0) {
   do_MergeNSplit = 1;
   C_old = C;
 } else {
   do_MergeNSplit = 0;
 }

}


copy_cluster_centres_a_to_b(best_cluster_centres,cluster_centres,k);      
assign_elements_to_centres_and_rerank(distance_matrix,cluster_centres,cluster_profile,n_elements,k);
C = inter_cluster_dissimilarity(distance_matrix,cluster_centres,cluster_profile,n_elements,k,
                                &C_inter,&C_ratio,&C_ratio2,&C_ratio3,&C_ratio4);


    /* eliminate the ambiguous points */
    /* (to turn off this part, set threshold > 1.) */
threshold = 1.1;
do {
  C = minimize_inter_cluster_dissimilarity(distance_matrix,cluster_centres,cluster_profile,n_elements,k,
                                           &C_inter,&C_ratio,&C_ratio2,&C_ratio3,&C_ratio4,&threshold);
} while(threshold > 0.);

output_results(distance_matrix,cluster_centres,cluster_profile,n_elements,k, labels, C);
//    dump_regrouped_distance_matrix(distance_matrix,cluster_profile,  n_elements, k);

C_ratio4_mean = C_ratio4_mean/loopmax_cont;
C_ratio4_var = C_ratio4_var/loopmax_cont - C_ratio4_mean*C_ratio4_mean;

fp = fopen("k_medoids_output/order_parameters.dat","a");
fprintf(fp,"%4d \t%le \t%le \t%le \t%le \t%le \t%le \t%le \t%le \t%le \t%le \t%le \t%8d \t%le\n",
        k, C, fabs(C-mean_C)/delta_C,
        mean_C, delta_C,
        C_inter, C_inter/C, C_ratio, C_ratio2, C_ratio3, 
        C/C_inter, C_ratio4,
        loopmax_cont,(C_ratio4_mean-C/C_inter)/C_ratio4_var);
fclose(fp);

}

printf("\n");

free_i1t(cluster_profile);
free_i1t(new_cluster_centres);
free_i1t(n_cluster_elements);
free_i1t(best_cluster_centres);
free_i1t(sign_groups);

return 0;
}


/******************************/
/******************************/

void  read_data(mylib::SymmetricMatrix<float>& distance_matrix, char labels[][20], int n){

  FILE *fp;
  int i, j;
  double temp;

  fp = fopen("k_medoids_distances.dat","r");

  for(i=0; i < n; i++){
//     fscanf(fp,"%s",&(labels[i]));
    // instead of reading the labels from file, print the sequence of numbers from 0:
    sprintf(labels[i],"%d",i);
  }
  
  /* SYMMETRIC MATRICES ONLY! */
  while (fscanf(fp,"%d %d %le",&i,&j,&temp)!=EOF) {
    distance_matrix(i,j)=temp;
  }
  fclose(fp);
  for(i=0; i<n; ++i) distance_matrix(i,i) = 0.;

}

/******************************/

void  read_coordinates(mylib::SymmetricMatrix<float>& distance_matrix, char labels[][20], int n, int Q){

//   FILE *fp;
  int i, j;
  double temp;

  
  for(i=0; i < n; i++){
//     fscanf(fp,"%s",&(labels[i]));
    // instead of reading the labels from file, print the sequence of numbers from 0:
    sprintf(labels[i],"%d",i);
  }

  // define a vector for coordinates
  std::vector< vector<double> > coordinates(n);
  for (i=0; i<n; ++i){
    coordinates[i].resize(Q);
  }

  // read coordinates from file
  ifstream fin("k_medoids_coordinates.dat");
  for(i=0; i<n; ++i){
    for(int q=0; q<Q; ++q)
      fin >> coordinates[i][q];
  }
  fin.close();

  // build the distance matrix
  printf("  DISTANCE MATRIX FROM TOP %d EIGENVECTORS \n", Q);
  printf("  ---------------------------------------- \n");
  for(i=0; i<n; ++i) {
    for(j=i+1; j<n; ++j) {
      temp = 0.;
      // distance on a unit hypersphere
      for(int q=0; q<Q; ++q) temp+=coordinates[i][q]*coordinates[j][q];
        if(temp>1.) temp=0.;
      else if(temp<-1.) temp=3.1415926535;
      else temp = acos(temp);
      distance_matrix(i,j)=temp;    
    }
  }
  printf("  Computed. \n\n");
  for(i=0; i<n; ++i) distance_matrix(i,i)=0.;

}

/******************************/

void  output_results(mylib::SymmetricMatrix<float>& distance_matrix,int *cluster_centres, int *cluster_profile, 
                     int n_elements, int k, char labels[][20], double dissim){

  int i,j, *populations; 
  FILE *fp, *fp2;
  char name[100];
  
  populations=i1t(k);
  for(i=0; i < k; i++){
    populations[i]=0;
  }


//   sprintf(name,"k_medoids_output/output_k_%d/cluster_members.dat",k);
//   fp2 = fopen(name,"w");
  
  for(i=0; i < k; i++){
    for(j=0 ; j < n_elements; j++){
      if (cluster_profile[j]==i){
        populations[i]++;
//         fprintf(fp2,"%d \t%d \n",i,j);
      }
    }
  }

//   fclose(fp2);


  sprintf(name,"k_medoids_output/output_k_%d/summary.dat",k);
  fp = fopen(name,"w");

  fprintf(fp,"Within_cluster_dissimilarity:  %le\n\n",dissim);
  fprintf(fp,"cluster: \tpopulation: \tmedoid id: \n");
  
  for(i=0; i < k; i++){
    fprintf(fp,"%d \t\t%d \t\t%d\n",i, populations[i], cluster_centres[i]);  
  }

  fclose(fp);


  sprintf(name,"k_medoids_output/output_k_%d/medoids_distances.dat",k);
  fp = fopen(name,"w");

  for(i=0; i < k; i++){
    for(j=0; j < k; j++){
      fprintf(fp,"%d \t%d \t%le \n",i,j,distance_matrix(cluster_centres[i],cluster_centres[j]));
    }
  }

  fclose(fp);
  
  
  sprintf(name,"k_medoids_output/output_k_%d/cluster_profile.dat",k);
  fp2 = fopen(name,"w");
  
  for(i=0 ; i < n_elements; ++i){
    fprintf(fp2,"%d \t%d \n",i,cluster_profile[i]);
  }  

  fclose(fp2);
  
  
  sprintf(name,"k_medoids_output/output_k_%d/medoids.dat",k);
  fp = fopen(name,"w");

  for(i=0; i < k; i++){
    fprintf(fp,"%d \t%d\n", cluster_profile[cluster_centres[i]], cluster_centres[i]);  
  }

  fclose(fp);

  
  free_i1t(populations);
  
}

/******************************/

void dump_regrouped_distance_matrix(mylib::SymmetricMatrix<float>& distance_matrix, 
                                    int *cluster_profile,  int n_elements, int k){

  int i,j, l, *mapping;
  FILE *fp;
  char name[100];

  sprintf(name,"k_medoids_output/output_k_%d/regrouped_distance_matrix.gnuplot",k);
  fp = fopen(name,"w");

  mapping=i1t(n_elements);
  l=0;
  for(i=0; i < k; i++){
    for(j=0 ; j < n_elements; j++){
      if (cluster_profile[j]==i) {
        mapping[l]=j;
        l++;
      }

    }
  }

  for(i=0; i < n_elements; i++){
    for(j=0 ; j < n_elements; j++){
      fprintf(fp, "%3d %3d %lf\n",i,j,distance_matrix(mapping[i],mapping[j]));
    }
    fprintf(fp,"\n");
  }
  fclose(fp);


  free_i1t(mapping);
}

/******************************/

void copy_cluster_centres_a_to_b(int *cluster_centres_a, int *cluster_centres_b, int k){
  int i;
  for(i=0; i < k; i++){
    cluster_centres_b[i]=cluster_centres_a[i];
  }
}


/******************************/

void assign_elements_to_centres_and_rerank(mylib::SymmetricMatrix<float>& distance_matrix, 
                                           int *cluster_centres, int *cluster_profile, int n_elements, int k){

  double d;
  int i,j, *population, *new_rank, *assigned, n_assigned, max, best = -1;
  int *temp_centres;

  /* for each member find the closest representative and attach it to its cluster */

  population=i1t(k);
  new_rank=i1t(n_elements);
  assigned=i1t(n_elements);
  temp_centres=i1t(n_elements);

  for(j=0; j < k; j++){
    population[j]=0;
    assigned[j]=0;
  }


  for(i=0; i <  n_elements; i++){
    d=100000000000.0;
    for(j=0; j < k; j++){
      if (d > distance_matrix(i,cluster_centres[j])){
        d = distance_matrix(i,cluster_centres[j]);
        cluster_profile[i]=j;
      }
    }
    population[ cluster_profile[i]]++;
  }


  /* reorder centres for decreasing population */

  n_assigned=0;
  do{
    max =-1;
    for(i=0; i < k; i++){
      if (assigned[i]==1) continue;
      if (population[i]>max){
        max=population[i];
        best=i;
      }
    }
    temp_centres[n_assigned]= cluster_centres[best];
    new_rank[best]=n_assigned;
    assigned[best]=1;
    n_assigned++;
  } while (n_assigned < k);
  
  for(i=0; i < k; i++){
    cluster_centres[i]=temp_centres[i];
  }
  
  /* and relabel the cluster belonging accordingly */

  for(i=0; i <  n_elements; i++){
    cluster_profile[i]=new_rank[cluster_profile[i]];
  }
  free_i1t(population);
  free_i1t(new_rank);
  free_i1t(assigned);
  free_i1t(temp_centres);
}

/******************************/

void assign_new_centres(mylib::SymmetricMatrix<float>& distance_matrix, 
                        int *cluster_centres, int *cluster_profile, int n_elements, int k){


  double d, temp;
  int i, j, l, best=-1, *cluster_members, n_cluster_members;
  
  cluster_members=i1t(n_elements);
  for(j=0; j < k; j++){

    /* identify the elements of the jth cluster */
    n_cluster_members=0;
    for(i=0; i < n_elements; i++){
      if (cluster_profile[i]==j)  {
        cluster_members[n_cluster_members]=i;
        n_cluster_members++;
      }
    }

    temp=0.0;
    for(l=0; l < n_cluster_members; l++){
      temp+= distance_matrix(cluster_centres[j],cluster_members[l]);
    }


    /* find the one that has minimal total distance from the other members */
    d=100000000000.0;
    for(i=0; i < n_cluster_members; i++){
      temp=0.0;
      for(l=0; l < n_cluster_members; l++){
        temp+= distance_matrix(cluster_members[i],cluster_members[l]);
      }
      
      if (temp < d){
        d=temp;
        best= cluster_members[i];
      }
    }
    
    /* and make it the new cluster representative */
    cluster_centres[j]=best;

  }

 //cluster_members=i1t(NMAXelements);
  free_i1t(cluster_members);

}

/******************************/

int unchanged_assignment(int *cluster_centres, int *new_cluster_centres, int k){

  int j;

  for(j=0; j < k; j++){
    if (cluster_centres[j]!= new_cluster_centres[j]) return(0);
  }
  
  return(1);

}


/******************************/

double within_cluster_dissimilarity(mylib::SymmetricMatrix<float>& distance_matrix, 
                                    int *cluster_centres, int * cluster_profile, int n_elements, int k){

  double temp;
  int i, representative;


  temp=0.0;
  
  for(i=0 ; i < n_elements; i++){
    representative=cluster_centres[cluster_profile[i]];
    temp += distance_matrix(i,representative);
  }

  return(temp);
}

/******************************/

double inter_cluster_dissimilarity(mylib::SymmetricMatrix<float>& distance_matrix, int *cluster_centres, 
                                   int * cluster_profile, int n_elements, int k, double *C_inter_p, double *C_ratio_p, 
                                   double *C_ratio2_p, double *C_ratio3_p, double *C_ratio4_p){

  double temp,d, min_dist,dist,C_inter, C_ratio,C_ratio2,C_ratio3,C_ratio4;
  int i,j,representative,label;

  temp = 0.0;
  C_inter = 0.;
  C_ratio = 0.;
  C_ratio2 = 0.;
  C_ratio3 = 0.;
  C_ratio4 = 0.;
  
  for(i=0 ; i < n_elements; ++i) {
    label = cluster_profile[i];
    representative = cluster_centres[label];
    d = distance_matrix(i,representative);
    temp += d;
    
    min_dist = -1.;
    for(j=0; j<k; ++j) {
      if(j==label) continue;
      representative = cluster_centres[j];
      dist = distance_matrix(i,representative);
      if(dist<min_dist || min_dist<0.) min_dist=dist;
    }
    C_inter += min_dist;
    
    if(d>0.) C_ratio += min_dist/d;
    C_ratio2 += (min_dist-d)/(min_dist+d);
    C_ratio3 += (min_dist-d)/(min_dist);
    C_ratio4 += d/min_dist;
  }

  *C_inter_p = C_inter;
  *C_ratio_p = C_ratio/n_elements;
  *C_ratio2_p = C_ratio2/n_elements;
  *C_ratio3_p = C_ratio3/n_elements;
  *C_ratio4_p = C_ratio4/n_elements;
  
  return(temp);
}


/******************************/


double compute_median_2(std::vector<double>& vec){

  std::vector<double> v = vec;
  std::sort(v.begin(),v.end());

  if( v.size()%2==0 ) return (v[v.size()/2] + v[v.size()/2-1]) /2.0;
  else return v[v.size()/2];
}


/******************************/

double minimize_inter_cluster_dissimilarity(mylib::SymmetricMatrix<float>& distance_matrix, 
                                            int *cluster_centres, int * cluster_profile, int n_elements, int k, double *C_inter_p, 
                                            double *C_ratio_p, double *C_ratio2_p, double *C_ratio3_p, double *C_ratio4_p,
                                            double *threshold_p){

  double temp,d, min_dist,dist;
  double C_inter,C_ratio,C_ratio2,C_ratio3,C_ratio4;
  int i,j,representative,label;
  double y,max_ratio;
  double prefactor, prefactor2;
  int i_max_ratio=-1, n_outliers, n_elements_new;
  FILE *fp;

  
  temp = 0.0;
  C_inter = 0.;
  C_ratio = 0.;
  C_ratio2 = 0.;
  C_ratio3 = 0.;
  C_ratio4 = 0.;
  max_ratio = 0.;
  n_outliers = 0;
  
  // vector of ratios (dist_to_1st_medoid/dist_to_2nd_medoid)
  vector<double> ratios(n_elements);
  int l=0; 
  
  for(i=0 ; i < n_elements; ++i) {
    label = cluster_profile[i];
    if(label < 0) {
      ++n_outliers;
      continue;
    }
    representative = cluster_centres[label];
    d = distance_matrix(i,representative);
    temp += d;
    
    min_dist = -1.;
    for(j=0; j<k; ++j) {
      if(j==label) continue;
      representative = cluster_centres[j];
      dist = distance_matrix(i,representative);
      if(dist<min_dist || min_dist<0.) min_dist=dist;
    }
    C_inter += min_dist;
    if(min_dist==0.) fprintf(stderr, "# WARNING: incorrect 2nd-medoid distance found.\n");
    
    if(d>0.) C_ratio += min_dist/d;
    C_ratio2 += (min_dist-d)/(min_dist+d);
    C_ratio3 += (min_dist-d)/(min_dist);
    y = d/min_dist;
    C_ratio4 += y;
    if(y > max_ratio) {
      max_ratio = y;
      i_max_ratio = i;
    }
    
    // store the ratios (dist_to_1st_medoid/dist_to_2nd_medoid)
    ratios[l] = d/min_dist; 
    ++l;
  }
  
  n_elements_new = n_elements - n_outliers;
  
  // resize the ratios vector
  ratios.resize(n_elements_new);

  if(max_ratio > *threshold_p) {
    cluster_profile[i_max_ratio] = -(cluster_profile[i_max_ratio] + 1);
  }
  else {
    *threshold_p=-1.;

    // compute the median of the ratios
    double median;
    median = compute_median_2(ratios);

    // compute the mean of the ratios
    double mean = 0.;
    for(i=0; i<n_elements_new; ++i) 
      mean += ratios[i];
    mean = mean/n_elements_new;
    
    // import the prefactor
    if (k < 2) {
      fprintf(stderr, "# ERROR! Invalid number of domains \n");
      exit(1);
    }
    else if (k > n_elements/20) {
      fprintf(stderr, "# WARNING! Number of domains is too large compared to the number of elements. \n");
      fprintf(stderr, "#          The prefactor may be not exact. \n");
    }
    else if (k > 400) {
      fprintf(stderr, "# WARNING! An approximate prefactor value will be used. \n");
    }

    //prefactor
    std::map<int,double> prefactor_table = get_prefactors();
    std::map<int,double> prefactor2_table = get_prefactors2();
    prefactor = (k<=400)? prefactor_table[k] : prefactor_table[400];
    prefactor2 = (k<=400)? prefactor2_table[k] : prefactor2_table[400];

    // print the quality score (from either the median and the mean of ratios)
    // by using the prefactor imported from file
    fp = fopen("k_medoids_output/quality_score.dat","a");
    if(fp==NULL) {fprintf(stderr,"Cannot open %s\n","k_medoids_output/quality_score.dat"); exit(1);}
    fprintf(fp, "%4d %le %le %le %le\n", k, prefactor/median, prefactor2/mean, prefactor, prefactor2);
    fclose(fp);
    
//    /* print the set of all ratios (dist_1/dist_2) */
//    char name[100];
//    sprintf(name,"k_medoids_output/output_k_%d/ratios.dat",k);
//    FILE *fp = fopen(name,"w");
//    for(i=0 ; i < n_elements; ++i) {
//      label = cluster_profile[i];
//      if(label < 0) continue;
//      representative = cluster_centres[label];
//      d = distance_matrix(i,representative);
//      
//      min_dist = -1.;
//      for(j=0; j<k; ++j) {
//	if(j==label) continue;
//	representative = cluster_centres[j];
//	dist = distance_matrix(i,representative);
//	if(dist<min_dist || min_dist<0.) min_dist=dist;
//      }
//      fprintf(fp, "%d %le %le %le \n", i, d/min_dist, d, min_dist);
//    }
//    fclose(fp);
  }

  *C_inter_p = C_inter;
  *C_ratio_p = C_ratio/n_elements_new;
  *C_ratio2_p = C_ratio2/n_elements_new;
  *C_ratio3_p = C_ratio3/n_elements_new;
  *C_ratio4_p = C_ratio4/n_elements_new;

  vector<double>().swap(ratios);

  return(temp);
}

/******************************/

void  random_within_cluster_dissimilarities(mylib::SymmetricMatrix<float>& distance_matrix, 
                                            int *cluster_centres, int n_elements, int k, double *mean_C, double *delta_C){


  int i,j,l,m,ok, loop;
  int *new_cluster_centres, *cluster_profile;
  double C, sum, sum2, ave, delta;
  

  new_cluster_centres=i1t(n_elements);
  cluster_profile=i1t(n_elements);
  sum=sum2=0.0;
//for(loop=0; loop < (50*n_elements); loop++){
//for(loop=0; loop < (10*k); loop++){
  for(loop=0; loop < 10; loop++){
    copy_cluster_centres_a_to_b(cluster_centres,new_cluster_centres,k);
    

    for(m=0; m < k; m++){
      do{
        i = (int)floor(ran2(&idum)*k); /* pick an existing centre */
        j = (int)floor(ran2(&idum)*n_elements); /* and reassign it */
        ok=1;
        for(l=0;l< k;l++) {
          if (l==i) continue;
          if(j==new_cluster_centres[l]) ok=0;
        }
      }while(ok==0);
      new_cluster_centres[i]=j; /* we have reassigned it */
    }

    /* compute the within_cluster_dissimilarity_score */
    assign_elements_to_centres_and_rerank(distance_matrix,new_cluster_centres,cluster_profile,n_elements,k);
    C= within_cluster_dissimilarity(distance_matrix,new_cluster_centres,cluster_profile,n_elements,k);
    sum+=C;
    sum2+=C*C;
  }
  copy_cluster_centres_a_to_b(new_cluster_centres,cluster_centres,k);

  sum=sum/loop;
  sum2=sum2/loop;

  ave=sum;
  delta= sqrt(sum2-sum*sum);
  
  (*mean_C)=ave;
  (*delta_C)=delta;
  free_i1t(new_cluster_centres);
  free_i1t(cluster_profile);
}

/******************************/

void  random_improvement_of_cluster_centres(mylib::SymmetricMatrix<float>& distance_matrix, 
                                            int *cluster_centres, int n_elements, int k, int n_trials){

  int i,j,l,m,ok, loop;
  int *new_cluster_centres, *cluster_profile;
  double C, Cbest;
  
  new_cluster_centres=i1t(n_elements);
  cluster_profile=i1t(n_elements);
  assign_elements_to_centres_and_rerank(distance_matrix,cluster_centres,cluster_profile,n_elements,k);
  Cbest= within_cluster_dissimilarity(distance_matrix,cluster_centres,cluster_profile,n_elements,k);

  for(loop=0; loop < n_trials; loop++){
    copy_cluster_centres_a_to_b(cluster_centres,new_cluster_centres,k);
    
    m = (int) 1+ floor(ran2(&idum)*(k-1));
    for(m=0; m < k; m++){
      do{
        i = (int)floor(ran2(&idum)*k); /* pick an existing centre */
        j = (int)floor(ran2(&idum)*n_elements); /* and reassign it */
        ok=1;
        for(l=0;l< k;l++) {
          if (l==i) continue;
          if(j==new_cluster_centres[l]) ok=0;
        }
      }while(ok==0);
      new_cluster_centres[i]=j; /* we have reassigned it */
    }
    
    /* compute the within_cluster_dissimilarity_score */
    assign_elements_to_centres_and_rerank(distance_matrix,new_cluster_centres,cluster_profile,n_elements,k);
    C= within_cluster_dissimilarity(distance_matrix,new_cluster_centres,cluster_profile,n_elements,k);
    
    /* if we have improved over the previous one we keep the new centres */

    if (C < Cbest){
      Cbest=C;
      copy_cluster_centres_a_to_b(new_cluster_centres,cluster_centres,k);
    }
  }

  free_i1t(new_cluster_centres);
  free_i1t(cluster_profile);

}


/******************************/


void random_initial_cluster_centres(int *cluster_centres, int n_elements, int k, int *sign_groups) {

  int flag;
  int c_q,c_qq;

  for(int q=0; q<k; ++q){
    do {
      flag = 0;
      /* assign a tentative random element to cluster_centres[q] */
      c_q = (int)floor(ran2(&idum)*n_elements);
      /* CHECK: */
      for(int qq=0; qq<q; ++qq) {
        c_qq = cluster_centres[qq];
        /* check whether it was already picked */
        if(c_q == c_qq) flag=1;
        /* or whether its sign group was already picked */
        if(sign_groups[c_q] == sign_groups[c_qq]) flag=1;
      } 
    } while(flag==1);
    /* medoid accepted */
    cluster_centres[q] = c_q;
  }

}

/******************************/

int are_cluster_centres_in_different_sign_groups(int *cluster_centres, int k, int *sign_groups) {

  int i,j,ii,check=1;
  
  for(i=0; i<k; ++i) {
    ii = sign_groups[cluster_centres[i]];
    for(j=i+1; j<k; ++j) {
      if(sign_groups[cluster_centres[j]]==ii) {
       check=0;
       break;
     }
     if(check==0) break;
   }
 }

 return check;
}

/******************************/

bool does_file_exist(const char *fileName) {

  std::ifstream infile(fileName);
  return infile.good();
}

/******************************/

void MergeNSplit(int n_elements, int k, mylib::SymmetricMatrix<float>& distance_matrix, 
                 int *cluster_profile, int *cluster_centres, int *sign_groups) {

  int i,j,n, *cluster_0_elements;
  int split_centre_1, split_centre_2, smallest_cluster_centre, nearest_cluster;
  double min_distance;
  
  cluster_0_elements = i1t(n_elements);
  n = 0;
  for(i=0; i<n_elements; ++i) {
    if(cluster_profile[i] == 0) {
      // store which elements are in the cluster 0
      cluster_0_elements[n] = i;
      ++n;
    }
    else {
      // shift the cluster labels by 1
      ++cluster_profile[i];
    }
  }
  
  // split the cluster 0 (which is the biggest in size) in two:
  /* find randomly the first centre of the split cluster */
  i = (int) (ran2(&idum)*n);
  split_centre_1 = cluster_0_elements[i];
  /* find the second centre of the split cluster */
  do {
    i = (int)(ran2(&idum)*n);
    split_centre_2 = cluster_0_elements[i];
  } while(split_centre_2 == split_centre_1);
  
  // assign each element in the former cluster 0 to the same cluster as the nearest one 
  // between split_centre_1 and split_centre_2
  for(i=0; i<n; ++i) {
    j = cluster_0_elements[i];
    if(distance_matrix(j,split_centre_1) < distance_matrix(j,split_centre_2)) {
      cluster_profile[j] = 0;
    }
    else cluster_profile[j]=1;
  }
  
  // in order to avoid empty clusters, each cluster must have
  // at least one element, i.e. itself
  cluster_profile[split_centre_1] = 0;
  cluster_profile[split_centre_2] = 1;
  
  // update the cluster centres
  smallest_cluster_centre = cluster_centres[k-1];
  for(i=k-1; i>1; --i) {
    cluster_centres[i] = cluster_centres[i-1];
  }
  cluster_centres[1] = split_centre_2;
  cluster_centres[0] = split_centre_1;

  
  // merge the smallest cluster to the nearest one
  /* find the nearest cluster */
  min_distance = distance_matrix(smallest_cluster_centre,cluster_centres[k-1]);
  nearest_cluster = k - 1;
  for(i=k-2; i>=0; --i) {
    if(distance_matrix(smallest_cluster_centre,cluster_centres[i]) < min_distance)
      nearest_cluster = i;
  }
  /* merge the two clusters */
  j = cluster_profile[cluster_centres[nearest_cluster]];
  for(i=0; i<n_elements; ++i) {
    if(cluster_profile[i] == k)
      cluster_profile[i] = j;
  }

  free_i1t(cluster_0_elements);
}

/******************************/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


double ran2 (long *idum)
{
  int j;
  long k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  double temp;


  if (*idum <= 0)
  {
    if (-(*idum) < 1)
     *idum = 1;
   else
     *idum = -(*idum);
   idum2 = (*idum);
   for (j = NTAB + 7; j >= 0; j--)
   {
     k = (*idum) / IQ1;
     *idum = IA1 * (*idum - k * IQ1) - k * IR1;
     if (*idum < 0)
       *idum += IM1;
     if (j < NTAB)
       iv[j] = *idum;
   }
   iy = iv[0];
 }
 k = (*idum) / IQ1;
 *idum = IA1 * (*idum - k * IQ1) - k * IR1;
 if (*idum < 0)
  *idum += IM1;
k = idum2 / IQ2;
idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
if (idum2 < 0)
  idum2 += IM2;
j = iy / NDIV;
iy = iv[j] - idum2;
iv[j] = *idum;
if (iy < 1)
  iy += IMM1;
if ((temp = AM * iy) > RNMX)
  return RNMX;
else
  return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX



/******************************/
