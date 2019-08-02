/*
  * Vector3d.h
 *
 *  Created on: Jul 2, 2012
 *      Author: gpolles
 */

#ifndef POINT3D_H_
#define POINT3D_H_

#include <iostream>
#include <vector>
#include <cmath>

namespace mylib{

inline double SQ(double i){ return (i)*(i) ;}


class Vector3d
{
  friend std::ostream& operator <<(std::ostream &out, const Vector3d& p);
  friend std::istream& operator >>(std::istream &in, Vector3d& p);


public:
  double X;
  double Y;
  double Z;
  public:
  Vector3d(){X=0;Y=0;Z=0;}
  Vector3d(double ax, double ay, double az) {X=ax;Y=ay;Z=az;}
  Vector3d(const Vector3d& a) {X = a.X ; Y = a.Y ; Z = a.Z;}


  // Comparison operators
  bool operator ==(const Vector3d& a) { if(a.X == X && a.Y==Y && a.Z==Z) return true; return false;}
  bool operator !=(const Vector3d& a) { if(a.X == X && a.Y==Y && a.Z==Z) return false; return true;}

  // vector-vector operators
  inline Vector3d operator-(const Vector3d& a) const{
      return Vector3d(X-a.X,Y-a.Y,Z-a.Z);
  }
  inline Vector3d operator-() const{
        return Vector3d(-X,-Y,-Z);
  }
  inline Vector3d operator+(const Vector3d& a) const{
      return Vector3d(X+a.X,Y+a.Y,Z+a.Z);
  }
  inline Vector3d operator*(const Vector3d& a) const{
        return Vector3d(X*a.X,Y*a.Y,Z*a.Z);
  }
  inline Vector3d operator/(const Vector3d& a) const{
        return Vector3d(X/a.X,Y/a.Y,Z/a.Z);
  }

  // Self-vector operators
  inline Vector3d& operator-=(const Vector3d& a) {
      X-=a.X; Y-=a.Y; Z-=a.Z;
      return *this;
  }
  inline Vector3d& operator+=(const Vector3d& a) {
        X+=a.X; Y+=a.Y; Z+=a.Z;
        return *this;
  }
  inline Vector3d& operator/=(const Vector3d& a) {
      X/=a.X; Y/=a.Y; Z/=a.Z;
      return *this;
  }
  inline Vector3d& operator*=(const Vector3d& a) {
        X*=a.X; Y*=a.Y; Z*=a.Z;
        return *this;
  }


  // Self-scalar operators
  inline Vector3d& operator/=(const double a) {
      X/=a; Y/=a; Z/=a;
      return *this;
  }
  inline Vector3d& operator*=(const double a) {
    X*=a; Y*=a; Z*=a;
    return *this;
  }
  inline Vector3d& operator-=(const double a) {
      X-=a; Y-=a; Z-=a;
      return *this;
  }
  inline Vector3d& operator+=(const double a) {
    X+=a; Y+=a; Z+=a;
    return *this;
  }

  // Vector-scalar operators
  inline Vector3d operator*(double a) const{
        return Vector3d(X*a,Y*a,Z*a);
  }
  inline Vector3d operator/(double a) const{
        return Vector3d(X/a,Y/a,Z/a);
  }
  inline Vector3d operator+(double a) const{
        return Vector3d(X+a,Y+a,Z+a);
  }
  inline Vector3d operator-(double a) const{
        return Vector3d(X-a,Y-a,Z-a);
  }




  double& operator[](size_t i) {
    return *(&X+i);
  }
  const double& operator[](size_t i) const{
    return *(&X+i);
  }

  inline double dotProduct(const Vector3d& p) const {return X*p.X + Y*p.Y + Z*p.Z;}
  inline Vector3d crossProduct(const Vector3d& p) const {
    return Vector3d( Y*p.Z - Z*p.Y,
                    Z*p.X - X*p.Z,
                    X*p.Y - Y*p.X);
  }


  inline double distance(const Vector3d& p) const {return sqrt(SQ(p.X-X)+SQ(p.Y-Y)+SQ(p.Z-Z)); }
  inline double distanceSQ(const Vector3d& p) const {return SQ(p.X-X)+SQ(p.Y-Y)+SQ(p.Z-Z); }

  inline double normSQ() const {return X*X + Y*Y + Z*Z;}
  inline double norm() const {return sqrt(X*X + Y*Y + Z*Z);}

  inline Vector3d& normalize(){
    if (X==0 && Y==0 && Z==0) return *this;
    double n = 1.0/sqrt(X*X + Y*Y + Z*Z);
    X *= n;Y *= n; Z *= n;
    return *this;
  }

  inline double angle(const Vector3d& p) const {
    double x = dotProduct(p)/norm()/p.norm();
    return acos( x );
  }
  inline double cosAngle(const Vector3d& p) const {
    return dotProduct(p)/norm()/p.norm();
  }

  Vector3d& rotate(const Vector3d& axis, const double angle){
    Vector3d normaxis=axis/axis.norm();
    Vector3d par = normaxis*(dotProduct(normaxis));
    Vector3d k = normaxis.crossProduct(*this);
    *this = par + (*this-par)*cos(angle) + k*sin(angle);
    return *this;
  }

  Vector3d& projectOn(const Vector3d& u){
    *this = u*(dotProduct(u)/u.normSQ());
    return *this;

  }

  Vector3d& projectOrthogonally(const Vector3d& u){ // project on orthogonal plane
    *this -= u*(dotProduct(u)/u.normSQ());
    return *this;
  }
};

class box_t: public std::vector<Vector3d>{
  public:
    box_t() {resize(2);}
    Vector3d lengths() const{
      return at(1)-at(0);
    }
};

inline double norm(const Vector3d& v) {return sqrt(v.X*v.X + v.Y*v.Y + v.Z*v.Z);}
inline double normSQ(const Vector3d& v) {return (v.X*v.X + v.Y*v.Y + v.Z*v.Z);}

inline double distance(const Vector3d& p,const Vector3d& q)  {return sqrt(SQ(p.X-q.X)+SQ(p.Y-q.Y)+SQ(p.Z-q.Z)); }
inline double distanceSQ(const Vector3d& p,const Vector3d& q)  {return SQ(p.X-q.X)+SQ(p.Y-q.Y)+SQ(p.Z-q.Z); }

inline std::ostream& operator <<(std::ostream &out, const Vector3d& p){
      out << " " << p.X << " " <<  p.Y << " " << p.Z ;
      return out;
}

inline std::istream& operator >>(std::istream &in, Vector3d& p){
        in >> p.X >>  p.Y >> p.Z;
        return in;
}

inline Vector3d operator *(const double x, const Vector3d& v) {
  return v*x;
}

inline double dotProduct(const Vector3d& p,const Vector3d& q) {
  return q.X*p.X + q.Y*p.Y + q.Z*p.Z;
}

inline Vector3d crossProduct(const Vector3d& v,const Vector3d& p) {
  return v.crossProduct(p);
}

inline Vector3d rotate(const Vector3d& v, const Vector3d& axis, const double angle){
  Vector3d normaxis=axis/norm(axis);
  Vector3d par = dotProduct(v,normaxis)*normaxis;
  Vector3d k = crossProduct(normaxis,v);
  return par + cos(angle)*(v-par) + sin(angle)*k ;
}

inline Vector3d normalize(const Vector3d& v){
  return v/v.norm();
}

inline Vector3d distanceVectorPbc(const Vector3d& p,const Vector3d& q,const Vector3d& box_l){
  Vector3d d = (q-p)/box_l;
  d.X = ( d.X - round(d.X) )*box_l.X;
  d.Y = ( d.Y - round(d.Y) )*box_l.Y;
  d.Z = ( d.Z - round(d.Z) )*box_l.Z;
  return d;
}
inline double distanceSqPbc(const Vector3d& p,const Vector3d& q,const Vector3d& box_l){
  return distanceVectorPbc(p,q,box_l).normSQ();
}
inline double distancePbc(const Vector3d& p,const Vector3d& q,const Vector3d& box_l){
  return sqrt(distanceVectorPbc(p,q,box_l).normSQ());
}
inline Vector3d PbcWrap(const Vector3d& p, const box_t& box){
  Vector3d box_l = box.lengths();
  Vector3d d = p/box_l;
  d.X = ( d.X - round(d.X) )*box_l.X;
  d.Y = ( d.Y - round(d.Y) )*box_l.Y;
  d.Z = ( d.Z - round(d.Z) )*box_l.Z;
  return box[0]+d;
}

inline double cosAngle(const Vector3d& p,const Vector3d& q) {
    return dotProduct(p,q)/norm(p)/norm(q);
}
inline double angle(const Vector3d& p,const Vector3d& q) {
    return acos(dotProduct(p,q)/norm(p)/norm(q));
}

inline Vector3d x_rotation(const Vector3d& v, double angle){
  double sint=sin(angle);
  double cost=cos(angle);
  return Vector3d(v.X,v.Y*cost+v.Z*sint,v.Y*(-sint)+v.Z*cost);
}

inline Vector3d y_rotation(const Vector3d& v, double angle){
  double sint=sin(angle);
  double cost=cos(angle);
  return Vector3d(v.X*cost+v.Z*(-sint),v.Y,v.X*sint+v.Z*cost);
}

inline Vector3d z_rotation(const Vector3d& v, double angle){
  double sint=sin(angle);
  double cost=cos(angle);
  return Vector3d(v.X*cost+v.Y*sint,v.X*(-sint)+v.Y*cost,v.Z);
}

inline Vector3d parallelProjection(const Vector3d& v, const Vector3d& u){
  return u*(dotProduct(v,u)/normSQ(u));
}


inline Vector3d orthogonalProjection(const Vector3d& v, const Vector3d& u){
  return v-u*(dotProduct(v,u)/normSQ(u));
}

namespace{
const Vector3d x_axis = Vector3d(1,0,0);
const Vector3d y_axis = Vector3d(0,1,0);
const Vector3d z_axis = Vector3d(0,0,1);
}
} // namespace mylib
#endif /* POINT3D_H_ */
