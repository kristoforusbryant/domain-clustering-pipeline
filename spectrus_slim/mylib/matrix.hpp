/*
 * matrix.hpp
 *
 *  Created on: Jul 9, 2014
 *      Author: gpolles
 */

#ifndef MATRIX_HPP_
#define MATRIX_HPP_
#include <vector>
#include <iomanip>
#include "Vector3d.hpp"

namespace mylib{

template <class T>
class SymmetricMatrix{
public:
  SymmetricMatrix(): N(0) {}
  SymmetricMatrix(size_t n): N(n) {
    try{
      data.resize(N);
      for (size_t i = 0; i < N; ++i){
        data[i].resize(i+1);
      }
    }catch(std::bad_alloc const&){
      std::cerr << "SymmetricMatrix resize() failed." << std::endl << std::flush;
      exit(1);
    }
  }

  T& operator()(size_t i,size_t j){
    return (j<=i) ? data[i][j]:data[j][i];
  }

  void resize(size_t n){
    try{
      data.resize(n);
      for (size_t i = 0; i < n; ++i){
        data[i].resize(i+1);
      }
    }catch(std::bad_alloc const&){
      std::cerr << "SymmetricMatrix resize() failed." << std::endl << std::flush;
      exit(1);
    }
    N=n;
  }

  void read(const char* fname){
    std::ifstream fin(fname, std::ios::binary);
    size_t n;
    fin.read((char*) &n,sizeof(size_t));
    resize(n);
    for (size_t i = 0; i < n; ++i){
      fin.read((char*) &(data[i][0]),(i+1)*sizeof(T));
    }
  }

  void read(const std::string fname) {read(fname.c_str());}

  void write(const char* fname) const{
    std::ofstream fout(fname, std::ios::binary);
    fout.write((char*) &N,sizeof(size_t));
    for (size_t i = 0; i < N; ++i){
      fout.write((char*) &(data[i][0]),(i+1)*sizeof(T));
    }
  }

  void write(const std::string fname) const {write(fname.c_str());}
  size_t size() {return N;}

private:
  std::vector<std::vector<T> > data;
  size_t N;
};

template <class T>
class Matrix {
public:
  Matrix(): N(0),M(0) {}
  Matrix(size_t n,size_t m): N(n),M(m) {
    try{
      data.resize(n);
      for (size_t i = 0; i < n; ++i){
        data[i].resize(m);
      }
    }catch(std::bad_alloc const&){
      std::cerr << "Matrix resize() failed." << std::endl << std::flush;
      exit(1);
    }
  }

  void reshape(size_t n, size_t m){
    try{
      data.resize(n);
      for (size_t i = 0; i < n; ++i){
        data[i].resize(m);
      }
    }catch(std::bad_alloc const&){
      std::cerr << "Matrix reshape() failed." << std::endl << std::flush;
      exit(1);
    }
    N=n;
    M=m;
  }

  std::vector<T> operator*(const std::vector<T>& v) const{
    if (M!=v.size()){
      std::cerr << "Matrix and vector have incompatible shape." << std::endl << std::flush;
      exit(1);
    }
    std::vector<T> w(N);
    for (size_t i = 0; i < N; ++i){
      for (size_t j = 0; j < M; ++j){
        w[i] += data[i][j]*v[j];
      }
    }
    return w;
  }

  std::vector<T>& operator[](size_t i){
    return data[i];
  }

  const std::vector<T>& operator[](size_t i) const{
    return data[i];
  }

  T& operator()(size_t i,size_t j){
    return data[i][j];
  }

  const T& operator()(size_t i,size_t j) const{
    return data[i][j];
  }

  
  size_t num_rows() const{
    return N;
  }

  size_t num_cols() const{
    return M;
  }

  //TODO: to_pointer
public:
  size_t N;
  size_t M;
private:
  std::vector<std::vector<T> > data;
};

template <class T>
inline std::ostream& operator<<(std::ostream& s,const Matrix<T>& m){
  for(size_t i=0; i<m.num_rows(); ++i){
    for(size_t j=0; j<m.num_cols(); ++j){
      s << std::setw(10) << std::setprecision(3) << std::scientific << m[i][j] << " ";
    }
    s << std::endl;
  }
  return s;
}


/////////////////
/// MATRIX 3d ///
/////////////////

class Matrix3d{
public:
  Matrix3d(){
    data.resize(3);
  }

  Vector3d operator*(const Vector3d& v) const{
    Vector3d w;
    for (size_t i = 0; i < 3; ++i){
      for (size_t j = 0; j < 3; ++j){
        w[i] += data[i][j]*v[j];
      }
    }
    return w;
  }

  Matrix3d operator+(const Matrix3d& B) const{
    Matrix3d C;
    for (size_t i = 0; i < 3; ++i){
      C[i] = data[i]+B[i];
    }
    return C;
  }
  Matrix3d operator-(const Matrix3d& B) const{
    Matrix3d C;
    for (size_t i = 0; i < 3; ++i){
      C[i] = data[i]-B[i];
    }
    return C;
  }

  Matrix3d operator-() const{
    Matrix3d C;
    for (size_t i = 0; i < 3; ++i){
      C[i] = -data[i];
    }
    return C;
  }

  Matrix3d& operator+=(const Matrix3d& B) {
    for (size_t i = 0; i < 3; ++i) data[i]+=B[i];
    return *this;
  }
  Matrix3d& operator-=(const Matrix3d& B) {
    for (size_t i = 0; i < 3; ++i) data[i]-=B[i];
    return *this;
  }

  Matrix3d operator*(double s) const{
    Matrix3d C;
    for (size_t i = 0; i < 3; ++i){
      C[i] = data[i]*s;
    }
    return C;
  }
  Matrix3d operator/(double s) const{
    Matrix3d C;
    for (size_t i = 0; i < 3; ++i){
      C[i] = data[i]/s;
    }
    return C;
  }

  Matrix3d& operator*=(double s) {
    for (size_t i = 0; i < 3; ++i) data[i]*=s;
    return *this;
  }
  Matrix3d& operator/=(double s) {
    for (size_t i = 0; i < 3; ++i) data[i]/=s;
    return *this;
  }

  Vector3d& operator[](size_t i){
    return data[i];
  }

  const Vector3d& operator[](size_t i) const{
    return data[i];
  }

  double det() const{
    return data[0][0]*(data[1][1]*data[2][2]-data[2][1]*data[1][2])
          -data[0][1]*(data[1][0]*data[2][2]-data[1][2]*data[2][0])
          +data[0][2]*(data[1][0]*data[2][1]-data[1][1]*data[2][0]);
  }

  Matrix3d inv() const{
    Matrix3d inv;
    inv[0][0] = (data[1][1] * data[2][2] - data[2][1] * data[1][2]) ;
    inv[0][1] = (data[0][2] * data[2][1] - data[0][1] * data[2][2]) ;
    inv[0][2] = (data[0][1] * data[1][2] - data[0][2] * data[1][1]) ;
    inv[1][0] = (data[1][2] * data[2][0] - data[1][0] * data[2][2]) ;
    inv[1][1] = (data[0][0] * data[2][2] - data[0][2] * data[2][0]) ;
    inv[1][2] = (data[1][0] * data[0][2] - data[0][0] * data[1][2]) ;
    inv[2][0] = (data[1][0] * data[2][1] - data[2][0] * data[1][1]) ;
    inv[2][1] = (data[2][0] * data[0][1] - data[0][0] * data[2][1]) ;
    inv[2][2] = (data[0][0] * data[1][1] - data[1][0] * data[0][1]) ;
    double invdet = data[0][0]*inv[0][0] + data[0][1]*inv[1][0] + data[0][2]*inv[2][0];
    return inv/invdet;
  }

  Matrix3d& trans(){
    std::swap(data[0][1],data[1][0]);
    std::swap(data[0][2],data[2][0]);
    std::swap(data[2][1],data[1][2]);
    return *this;
  }

 private:
  std::vector<Vector3d > data;
 
};



inline Matrix3d operator*(const Matrix3d& A, const Matrix3d& B){
  Matrix3d C;
  for (size_t i = 0; i < 3; ++i){
    for (size_t j = 0; j < 3; ++j){
      for (size_t q = 0; q < 3; ++q){
        C[i][j] += A[i][q]*B[q][j];
      }
    }
  }
  return C;
}

inline Matrix3d operator*(double s, const Matrix3d& A){return A*s;}

template <class T>
Matrix<T> operator*(const Matrix<T>& A, const Matrix<T>& B){
  if (A.M!=B.N ){
    std::cerr << "operator*(Matrix&, Matrix&): matrices have incompatible shape." << std::endl << std::flush;
    exit(1);
  }

  Matrix<T> C(A.N,B.M);
  for (size_t i = 0; i < A.N; ++i){
    for (size_t j = 0; j < B.M; ++j){
      for (size_t q = 0; q < A.M; ++q){
        C[i][j] += A[i][q]*B[q][j];
      }
    }
  }
  return C;
}

inline Matrix3d x_rot_matrix(double angle){
  Matrix3d R;
  double cost = cos(angle);
  double sint = sin(angle);
  R[0][0] = 1;
  R[1][1] = R[2][2] = cost;
  R[2][1] = sint;
  R[1][2] = -sint;
  return R;
}

inline Matrix3d y_rot_matrix(double angle){
  Matrix3d R;
  double cost = cos(angle);
  double sint = sin(angle);
  R[1][1] = 1;
  R[0][0] = R[2][2] = cost;
  R[2][0] = sint;
  R[0][2] = -sint;
  return R;
}

inline Matrix3d z_rot_matrix(double angle){
  Matrix3d R;
  double cost = cos(angle);
  double sint = sin(angle);
  R[2][2] = 1;
  R[0][0] = R[1][1] = cost;
  R[1][0] = sint;
  R[0][1] = -sint;
  return R;
}

inline std::ostream& operator<<(std::ostream& s,const Matrix3d& m){
  for(int i=0; i<3; ++i){
      s << m[i] << std::endl;
  }
  return s;
}


class SymmetricMatrix3f{
public:
  SymmetricMatrix3f() : data(std::vector<float>(6,0.0f)) {}

  float& operator()(size_t i,size_t j){
    if(j<i) return (data[j=0? i: i+j+1]);
    return (data[i=0? j: i+j+1]);
  }
  const float& operator()(size_t i,size_t j)const{
    if(j<i) return (data[j=0? i: i+j+1]);
    return (data[i=0? j: i+j+1]);
  }

  SymmetricMatrix3f operator +(const SymmetricMatrix3f& B) const{
    SymmetricMatrix3f C;
    for (size_t i=0; i<6;++i) C.data[i] = data[i]+B.data[i];
    return C;
  }
  SymmetricMatrix3f operator -(const SymmetricMatrix3f& B) const {
    SymmetricMatrix3f C;
    for (size_t i=0; i<6;++i) C.data[i] = data[i]+B.data[i];
    return C;
  }

  SymmetricMatrix3f operator *(float c) const {
    SymmetricMatrix3f C;
    for (size_t i=0; i<6;++i) C.data[i] = data[i]*c;
    return C;
  }
  SymmetricMatrix3f operator /(float c) const {
    SymmetricMatrix3f C;
    for (size_t i=0; i<6;++i) C.data[i] = data[i]/c;
    return C;
  }

  SymmetricMatrix3f& operator *=(float c) {
    for (size_t i=0; i<6;++i) data[i]*=c;
    return *this;
  }
  SymmetricMatrix3f& operator /=(float c) {
    for (size_t i=0; i<6;++i) data[i]/=c;
    return *this;
  }



  std::vector<float> data;
};

// operations between vectors
template <class T>
inline T dotProduct(const std::vector<T>& a, const std::vector<T>& b){
  T dp = 0;
  size_t dim = a.size();
  if(dim!=b.size()) {
    std::cerr << "ERROR: dotProduct(const std::vector<T>& a, const std::vector<T>& b): incompatible size" << std::endl;
    return 0;
  }
  for (size_t i = 0; i<dim; ++i){
    dp+=a[i]*b[i];
  }
  return dp;
}

template <class T>
inline T vec_normSQ(const std::vector<T>& a){
  T n2 = 0;
  size_t dim = a.size();
  for (size_t i = 0; i<dim; ++i){
    n2+=a[i]*a[i];
  }
  return n2;
}

template <class T>
inline T vec_norm(const std::vector<T>& a){
  return sqrt(vec_normSQ(a));
}

template <class T> 
inline std::vector<T> operator*(const std::vector<T>& v, T a){
  std::vector<T> w = v;
  for (size_t i =0; i< v.size(); ++i) w[i]*=a;
  return w;
}

template <class T> 
inline std::vector<T> operator*(T a, const std::vector<T>& v){
  return v*a;
}


template <class T> 
inline std::vector<T> operator+(const std::vector<T>& v, const std::vector<T>& w){
  if(v.size()!=w.size()) {
    std::cerr << "ERROR: operator + (const std::vector<T>& a, const std::vector<T>& b): incompatible size" << std::endl;
    exit(1);
  }
  std::vector<T> u = v;
  for (size_t i =0; i< v.size(); ++i) u[i] = v[i]+w[i];
  return u;
}

template <class T> 
inline std::vector<T> operator-(const std::vector<T>& v, const std::vector<T>& w){
  if(v.size()!=w.size()) {
    std::cerr << "ERROR: operator - (const std::vector<T>& a, const std::vector<T>& b): incompatible size" << std::endl;
    exit(1);
  }
  std::vector<T> u(v.size());
  for (size_t i =0; i< v.size(); ++i) u[i] = v[i]-w[i];
  return u;
}

template <class T> 
inline std::vector<T> operator-(const std::vector<T>& v){
  std::vector<T> u = v;
  for (size_t i =0; i< v.size(); ++i) u[i] = -u[i];
  return u;
}

template <class T> 
inline std::ostream& operator<<(std::ostream& s,const std::vector<T>& v){
  for(size_t i=0; i<v.size(); ++i){
      s << v[i] << " ";
  }
  return s;
}

} // namespace mylib

#endif /* MATRIX_HPP_ */
