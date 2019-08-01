#ifndef __MODES_TO_FLUCT_MATRIX_H__
#define __MODES_TO_FLUCT_MATRIX_H__
// builds a correlation matrix, then a distance fluctuation
// matrix from normal modes (saved in a folder named "eigenspace/").

#include "mylib/contactmap.hpp"
#include <vector>

int ModesToFluctMatrix(char* filename, int num_atoms, int decimation, int num_modes, mylib::ContactMap& cm,
                       std::vector< std::vector<double> >& fluct);

#endif
