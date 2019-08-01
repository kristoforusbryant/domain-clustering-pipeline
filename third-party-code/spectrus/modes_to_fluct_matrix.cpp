//============================================================================
// Name        : ModesToFluctMatrix.cpp
// Author      : Luca Ponzoni, Guido Polles
// Version     : 1.0
// Description : builds a correlation matrix, then a distance fluctuation
//		 matrix from normal modes (saved in a folder named "eigenspace/").
// Input       : a folder with modes in the parent directory ("../eigenspace/")
//============================================================================

#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <iomanip>

#include "mylib/mylib.h"

using namespace std;
using namespace mylib;



int ModesToFluctMatrix(char* filename, int num_atoms, int decimation, int num_modes, mylib::ContactMap& cm,
	 std::vector< std::vector<double> >& fluct) {

  // make up for null modes
  num_modes += 6;

  size_t num_dec_atoms = num_atoms/decimation;

  NormalModes normal_modes;
  normal_modes.readFromFile( string("../eigenspace/").append(filename),
          num_modes,
          num_atoms);

  // put inverse eigenvalues in lambda
  vector<double> lambda(num_modes);
  for (size_t mode = 6; mode < num_modes; ++mode){
    lambda[mode]=1.0/normal_modes.getEigenvalue(mode);
  }

  // put the values of eigenvectors relative to the atoms
  // which are not decimated into dx
  vector< vector<double> > dx(num_modes);
  for (size_t mode = 6; mode < num_modes; ++mode){
    vector<Vector3d> eigenvec = normal_modes.getEigenvector(mode);
    dx[mode]=vector<double>((num_atoms/decimation)*3, 0.0);
    for (size_t i = 0; i < num_atoms/decimation; ++i){
      dx[mode][i*3]   = eigenvec[i*decimation].X;
      dx[mode][i*3+1] = eigenvec[i*decimation].Y;
      dx[mode][i*3+2] = eigenvec[i*decimation].Z;
    }
  }

//   // opens output file
//   fstream fout("correlation_matrix.bindat",ios::out | ios::binary);
//   size_t num_elements = 0;
//   size_t line_size = num_dec_atoms*3;
//   fout.write(reinterpret_cast<char*>(&line_size),sizeof(size_t));
// 
//   // writes the correlation matrix
//   for (size_t i = 0; i < num_dec_atoms*3; ++i){
//     for (size_t j = 0; j < num_dec_atoms*3; ++j){
//       if(cm.isNeighbor(i,j)) {
// 	float el = 0.0;
// 	for (size_t mode = 6; mode < num_modes; ++mode){
// 	  el += dx[mode][i]*dx[mode][j]*lambda[mode];
// 	}
// 	// write i,j element
// 	fout.write(reinterpret_cast<char*>(&el),sizeof(float));
//       }
//       ++num_elements;
//     }
//   }
//   fout.close();


  // compute distance fluctuations between atoms in the contact map
  double entry;
  int j,ii,jj;
  for (int i=0; i<num_dec_atoms; ++i){
    for (int z=0; z<cm[i].size(); ++z){
      j = cm[i][z];
      entry = 0.;
      for(int mu=0; mu<3; ++mu) {
	ii = 3*i + mu;
	jj = 3*j + mu;
	// given the correlation matrix M, the formula would be:
	// entry += M[ii][ii] + M[jj][jj] - 2.*M[jj][ii];
	// Here, we compute the correlation matrix M "on-the-fly"
	for (size_t mode = 6; mode < num_modes; ++mode){
	  entry += dx[mode][ii]*dx[mode][ii]*lambda[mode];
	  entry += dx[mode][jj]*dx[mode][jj]*lambda[mode];
	  entry -= 2.*dx[mode][jj]*dx[mode][ii]*lambda[mode];
	}
      }
      fluct[i][z] = sqrt(entry);
    }
  }


  return 0;
}
