/*    sparsematrix.hpp

Library for sparse squared matrices.
Note that insertion requires sorting of elements
and every access requires a binary search, this is NOT
fast.

*/


#ifndef __SPARSEMATRIX_HPP__
#define __SPARSEMATRIX_HPP__

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>

#define MATRIX_VERSION 1

namespace mylib{

enum MATRIX_TYPE{
  MT_NONE,
  MT_GENERAL,
  MT_SPARSE,
  MT_SYMMETRIC,
  MT_SPARSESYMMETRIC
};



struct SparseMatrixHeader{
  int size;
  int version;
  int precision;
  MATRIX_TYPE type;
};

template <class real_t>
struct CscSparseMatrix{
  int dimension; /* dimension of the matrix */
  int numNonZero; /* number of nonzero values */
  std::vector<int> icols; /* index of the columns */
  std::vector<int> prows; /* index at which a new row start */
  std::vector<real_t> values; /* array of values */
};

template <class real_t>
struct SparseMatrixElement{
  size_t columnIndex;
  real_t value;
};

template <class real_t>
inline bool SparseElementsSortingFunction(const SparseMatrixElement<real_t>& a, const SparseMatrixElement<real_t>& b) {
      return (a.columnIndex<b.columnIndex);
}

template <class real_t>
int sparseMatrixBinFind(size_t a, const std::vector<SparseMatrixElement<real_t> >& v){
  int imin=0; int imax=v.size()-1;
  while(imax>=imin){
    int imid = (imin+imax)/2;
    if (a==v[imid].columnIndex){
      return imid;
    }else if (v[imid].columnIndex < a){
      imin = imid+1;
    }else{
      imax = imid-1;
    }
  }
  return -1;
}

template <class real_t>
inline int sparseFind(size_t a, const std::vector<SparseMatrixElement<real_t> >& v){
  for (size_t i = 0; i < v.size(); ++i) {
    if(a==v[i].columnIndex) return i;
  }
  return -1;
}

template <class real_t = double>
class SparseMatrix
{
public:
  SparseMatrix() : _size(0), _numNonZero(0){}
  SparseMatrix(size_t size) : _size(size), _numNonZero(0){
    _rows.resize(size);
  }
  virtual ~SparseMatrix() {}

  real_t& insert(size_t i, size_t j){
    typename std::vector<SparseMatrixElement<real_t> >::iterator it;
    for ( it=_rows[i].begin(); it!=_rows[i].end(); ++it){
      if (j==(*it).columnIndex){
        std::cerr << "SparseMatrix: Insert() called to existing object ("<<i<<","<<j<<")"<< std::endl;
        exit(1);
      }
      if ((*it).columnIndex>j) break;
    }
    it=_rows[i].insert(it,SparseMatrixElement<real_t>());
    (*it).columnIndex = j;
    (*it).value = 0;
    ++_numNonZero;
    return (*it).value;
  }

  void insert(size_t i, size_t j, real_t v){
    insert(i,j) = v;
  }

  real_t& operator() (size_t i, size_t j){
    int pos=sparseMatrixBinFind(j,_rows[i]);
    if(pos==-1){
      std::cerr << "SparseMatrix: operator() called to non-existing object ("<<i<<","<<j<<")"<< std::endl;
      exit(1);
    }
    return _rows[i][pos].value;
  }
  
  const real_t& operator() (size_t i, size_t j) const{
    int pos=sparseMatrixBinFind(j,_rows[i]);
    if(pos==-1){
      std::cerr << "SparseMatrix: operator() called to non-existing object ("<<i<<","<<j<<")"<< std::endl;
      exit(1);
    }
    return _rows[i][pos].value;
  }

  virtual real_t& at(size_t i, size_t j){ return operator()(i,j);  }

  void getCsc(struct CscSparseMatrix<real_t>& csc){

    csc.dimension = static_cast<int>(_size);
    csc.numNonZero = static_cast<int>(_numNonZero);
    csc.values.reserve(_numNonZero);
    csc.icols.reserve(_numNonZero);
    csc.prows.reserve(_size+1);

    for (int i = 0; i < static_cast<int>(_size); ++i) {
      sort(_rows[i].begin(),_rows[i].end(),SparseElementsSortingFunction);
    }
    int num_ptr = 0;
    for (int i = 0; i < static_cast<int>(_size); ++i) {
      csc.prows.push_back(num_ptr) ;
      for (int j = 0; j < static_cast<int>(_rows[i].size()); ++j) {
        csc.values.push_back(_rows[i][j].value);
        csc.icols.push_back(_rows[i][j].columnIndex);
        ++num_ptr;
      }
    }
    csc.prows.push_back(_numNonZero);
  }

  virtual bool isDefined(size_t i, size_t j) const{
    int pos=sparseMatrixBinFind(j,_rows[i]);
    if(pos==-1){
      return false;
    }
    return true;
  }

  virtual bool isDefined(size_t i, size_t j, real_t* value) {
    int pos=sparseMatrixBinFind(j,_rows[i]);
    if(pos==-1){
      value = NULL;
      return false;
    }
    *value = _rows[i][pos].value;
    return true;
  }
  
  std::vector<SparseMatrixElement<real_t> >& getRow(size_t i){
    return _rows[i];
  }

  const std::vector<SparseMatrixElement<real_t> >& getRow(size_t i) const{
    return _rows[i];
  }

  size_t getNumNonZero() const {
    return _numNonZero;
  }

  size_t size() const {
    return _size;
  }

  void clear(){
    _size =0;
    _numNonZero = 0;
    for(size_t i = 0; i< _rows.size(); ++i){
      std::vector<SparseMatrixElement<real_t> >().swap(_rows[i]);
    }
  }

  void remove(size_t i, size_t j){
    int pos=sparseMatrixBinFind(j,_rows[i]);
    if (pos==-1){
      std::cerr << "SparseMatrix: remove() called to non-existing object (" << i << "," << j << ")" << std::endl;
      exit(1);
    }
    _rows[i].erase(_rows[i].begin() + pos);
    --_numNonZero;
  }

  void reserve(size_t n){  // reserve space for n elements on each row
    for (size_t i = 0; i < n; ++i) {
      try{
        _rows[i].reserve(n);
      }catch(std::bad_alloc const&){
        std::cerr << "Sparse matrix reserve() failed." << std::endl << std::flush;
        exit(1);
      }
    }
  }

  virtual void readFromFile(const char* fname){

    std::ifstream in(fname, std::ios::binary);
    if (!in.is_open()) {
      std::cerr << "SparseMatrix::readFromFile(): Unable to open file: " <<fname <<std::endl;
      return;
    }
    int header_size;
    in.read((char*)&header_size, sizeof(int));

    in.seekg (0, in.beg);
    SparseMatrixHeader header;
    in.read((char*)&header, header_size);


    if(header.type != this->type()) {
      std::cerr << "SparseMatrix::readFromFile(): Uncompatible file format.\n";
      return;
    }

    if(sizeof(real_t)!=header.precision){
      std::cerr << "SparseMatrix::readFromFile(): wrong precision.\n";
      return;
    }

    in.read((char*)&_size,sizeof(_size));
    in.read((char*)&_numNonZero,sizeof(_numNonZero));
    if(!in.good()) {
      std::cerr << "SparseMatrix::readFromFile(): Reading error 1.\n";
      return;
    }
    try{
        _rows.resize(_size);
        for (size_t i = 0; i < _size; ++i)
        {
          size_t num_ptr = 0;
          in.read((char*)&num_ptr,sizeof(num_ptr));
          _rows[i].resize(num_ptr);
        }
        for (size_t i = 0; i < _size; ++i){
          in.read((char*)&_rows[i][0],sizeof(SparseMatrixElement<real_t>)*_rows[i].size());
        }
      }catch (std::bad_alloc const&) {
        return;
      }
  }

  virtual void saveToFile(const char* fname){
    std::ofstream out(fname, std::ios::binary);

    SparseMatrixHeader header;
    header.type = this->type();
    header.version = MATRIX_VERSION;
    header.precision = sizeof(real_t);
    header.size = sizeof(header);

    out.write((char*)&header, sizeof(header));

    for (int i = 0; i < static_cast<int>(_size); ++i) {
      sort(_rows[i].begin(),_rows[i].end(),SparseElementsSortingFunction<real_t>);
    }

    out.write((char*)&_size,sizeof(_size));
    out.write((char*)&_numNonZero,sizeof(_numNonZero));
    size_t num_ptr = 0;
    for (size_t i = 0; i < _size; ++i){
      num_ptr= _rows[i].size();
      out.write((char*)&num_ptr,sizeof(num_ptr));
    }
    for (size_t i = 0; i < _size; ++i){
      out.write((char*)&_rows[i][0],sizeof(SparseMatrixElement<real_t>)*_rows[i].size());
    }
    out.close();
  }

  void updateNumNonZero(){
    _numNonZero = 0;
    for (size_t i = 0; i < _rows.size(); ++i) {
      _numNonZero+=_rows[i].size();
    }
  }

  virtual MATRIX_TYPE type() {return MT_SPARSE;}

public:
  std::vector<std::vector<SparseMatrixElement<real_t> > > _rows;
  size_t _size;
  size_t _numNonZero;
};



/**************************
 * Symmetric Sparse Matrix
 **************************/
template <class real_t = double>
class SymmetricSparseMatrix : public SparseMatrix<real_t>
{
public:
  SymmetricSparseMatrix() : SparseMatrix<real_t>() {}
  virtual ~SymmetricSparseMatrix() {}
  SymmetricSparseMatrix(size_t dim) : SparseMatrix<real_t>(dim){}
  SymmetricSparseMatrix(SparseMatrix<real_t>& m) : SparseMatrix<real_t>(m.size()){
    SparseMatrix<real_t>::_numNonZero = 0;
    for (size_t i = 0; i < m.size(); ++i) {
      for (size_t j = 0; j < m.getRow(i).size(); ++j) {
        if (m.getRow(i)[j].columnIndex >= i){
          SparseMatrix<real_t>::_rows[i].push_back(m.getRow(i)[j]) ;
          ++SparseMatrix<real_t>::_numNonZero;
        }
      }
    }
  }


  real_t& insert(size_t i, size_t j){
    if (j<i) {
      return SparseMatrix<real_t>::insert(j,i);
    }else{
      return SparseMatrix<real_t>::insert(i,j);
    }
  }

  void insert(size_t i, size_t j, real_t v){
    if (j<i) {
      SparseMatrix<real_t>::insert(j,i) = v;
    }else{
      SparseMatrix<real_t>::insert(i,j) = v;
    }
  }


  inline real_t& operator() (size_t i, size_t j){
    if (j<i) {
      return SparseMatrix<real_t>::operator ()(j,i);
    }else{
      return SparseMatrix<real_t>::operator ()(i,j);
    }
  }

  inline const real_t& operator() (size_t i, size_t j) const{
    if (j<i) {
      return SparseMatrix<real_t>::operator ()(j,i);
    }else{
      return SparseMatrix<real_t>::operator ()(i,j);
    }
  }

  virtual real_t& at(size_t i, size_t j) {return operator()(i,j);};

  inline bool isDefined(size_t i, size_t j) const{
    if (j<i) {
      return SparseMatrix<real_t>::isDefined(j,i);
    }else{
      return SparseMatrix<real_t>::isDefined(i,j);
    }
  }

  virtual MATRIX_TYPE type(){ return MT_SPARSESYMMETRIC; }

};

template <class T>
inline std::ostream& operator<<(std::ostream& s,const SparseMatrix<T>& m){
  for(size_t i=0; i<m.size(); ++i){
    for(size_t j=0; j<m.size(); ++j){
      if(m.isDefined(i,j)){
        s << std::setw(10) << std::setprecision(3) << std::scientific << m(i,j) << " ";
      }else{
        s << std::setw(10) << std::setprecision(3) << std::scientific << 0.000 << " ";
      }
    }
    s << std::endl;
  }
  return s;
}

template <class T>
inline std::ostream& operator<<(std::ostream& s,const SymmetricSparseMatrix<T>& m){
  for(size_t i=0; i<m.size(); ++i){
    for(size_t j=0; j<m.size(); ++j){
      if(m.isDefined(i,j)){
        s << std::setw(10) << std::setprecision(3) << std::scientific << m(i,j) << " ";
      }else{
        s << std::setw(10) << std::setprecision(3) << std::scientific << 0.000 << " ";
      }
    }
    s << std::endl;
  }
  return s;
}


} // namespace mylib
#endif // __SPARSEMATRIX_HPP__
