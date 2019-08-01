
#ifndef __KRYLOV_SHUR_HPP__
#define __KRYLOV_SHUR_HPP__

#include "mylib/sparsematrix.hpp"

int krylov_shur(mylib::SymmetricSparseMatrix<double>& input_matrix,
		int n_eigvec,
		std::vector< std::vector<double> >& eigvec,
		std::vector<double>& eigval);

#endif
