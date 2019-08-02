/*
 * NormalModes.h
 *
 *  Created on: Jul 5, 2012
 *      Author: gpolles
 */

#ifndef __NORMALMODES_HPP__
#define __NORMALMODES_HPP__

#include <vector>
#include "Vector3d.hpp"

namespace mylib{

class NormalModes
{
public:
  NormalModes() :  _numModes(0) {}
  virtual ~NormalModes() {}

  void    readFromFile(const std::string& basename, size_t numModes, size_t dim){
    char fname[basename.length()+100] ;
    sprintf(fname,"%s%s", basename.c_str(), "_eigenvalues.dat");
    std::ifstream eigvalFile(fname);
    if (!eigvalFile.is_open()){
      std::cerr <<  "Fatal error. Could not open " << fname << std::endl;
      exit(1);
    }

    try{
      _eigenvalues.reserve(numModes);
    }catch(std::bad_alloc const&){
      std::cerr <<  "NormalMode eigenvalues memory allocation failed!" << std::endl << std::flush;
      exit(1);
    }

    while(eigvalFile.good()){
      if (_eigenvalues.size() == numModes) break;
      int index;
      double val;
      eigvalFile >> index >> val;
      _eigenvalues.push_back(val);
    }

    if(_eigenvalues.size() < numModes){
      std::cerr <<  "Error. Could read only " <<  _eigenvalues.size() << " eigenvalues."<< std::endl;
      exit(1);
    }

    try{
      _eigenvectors.resize(numModes);
      for (size_t i = 0; i < numModes; ++i) {

        sprintf(fname, "%s%s%lu%s",basename.c_str(),"_eigenvector_",i,".dat");
        std::ifstream eigvecFile(fname);
        if (!eigvecFile.is_open()){
          std::cerr <<  "Fatal error. Could not open " << fname << std::endl<< std::flush;
          exit(1);
        }

        _eigenvectors[i].reserve(dim);
        for (size_t j = 0; j < dim; ++j) {
          int index;
          Vector3d p3d;
          eigvecFile >> index >> p3d;
          _eigenvectors[i].push_back(p3d);

        }

        if(eigvecFile.fail()){
          std::cerr <<  "Fatal error reading " << fname << std::endl << std::flush;
          exit(1);
        }

      }

    }catch(std::bad_alloc const&){
      std::cerr <<  "NormalMode eigenvector memory allocation failed!" << std::endl << std::flush;
      exit(1);
    }


    _numModes = numModes;

  }

  double  getEigenvalue(size_t i) const{  return _eigenvalues[i]; }
  const std::vector<Vector3d>& getEigenvector(size_t i) const {return _eigenvectors[i];}

  size_t  getNumModes() const { return _numModes; }

  void clear(){  _eigenvalues.resize(0); _eigenvectors.resize(0); _numModes = 0; }

private:
  std::vector<double> _eigenvalues;
  std::vector<std::vector<Vector3d> > _eigenvectors;
  size_t _numModes;

};



}//namespace mylib

#endif /* __NORMALMODES_HPP__ */
