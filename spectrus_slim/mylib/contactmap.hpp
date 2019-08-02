#ifndef CONTACTMAP_HPP
#define CONTACTMAP_HPP

#include <vector>
#include <algorithm>

#include "Vector3d.hpp"
#include "verletlist.hpp"
#include "containerutils.hpp"

namespace mylib{

class ContactMap{
 public:
  ContactMap(std::vector<Vector3d>& crd, double cutoff){
    _cm.resize(crd.size());
    box_t box = get_enclosing_box(crd);
    box[0]-=0.01*cutoff;
    box[1]+=0.01*cutoff;
    VerletList3D vl(crd,box,cutoff);
    for (size_t i=0; i<crd.size(); ++i){
      _cm[i] = vl.getNeighbors(crd, i, cutoff);
      std::sort(_cm[i].begin(),_cm[i].end());
    }
  }

  std::vector<size_t>& getNeighbors(size_t i){
    return _cm[i];
  }

  std::vector<size_t>& operator[](size_t i){
    return _cm[i];
  }

  const std::vector<size_t>& getNeighbors(size_t i) const {
    return _cm[i];
  }

  const std::vector<size_t>& operator[](size_t i) const{
    return _cm[i];
  }

  bool isNeighbor(size_t i,size_t j) const{
    if(i==j) return true;
    return is_in(i,_cm[j]);
  }

  size_t size() const{ return _cm.size();}

  void clear() {std::vector<std::vector<size_t> >().swap(_cm);}
 public:
  std::vector<std::vector<size_t> > _cm;
};

class ContactMapPBC{
 public:
  ContactMapPBC(std::vector<Vector3d>& crd, const box_t& box, double cutoff){
    _cm.resize(crd.size());
    if(norm(box.lengths())==0)
      std::cerr << "WARNING: ContactMapPBC box of 0 size" << std::endl;
    VerletList3D vl(crd,box,cutoff);
    for (size_t i=0; i<crd.size(); ++i){
      _cm[i] = vl.getNeighborsPBC(crd, i, cutoff);
    }
  }

  std::vector<size_t>& getNeighbors(size_t i){
    return _cm[i];
  }

  std::vector<size_t>& operator[](size_t i){
    return _cm[i];
  }

  const std::vector<size_t>& getNeighbors(size_t i) const{
    return _cm[i];
  }

  const std::vector<size_t>& operator[](size_t i) const{
    return _cm[i];
  }

  size_t size() const{ return _cm.size();}
 public:
  std::vector<std::vector<size_t> > _cm;
};


} // namespace mylib
#endif // CONTACTMAP_HPP
