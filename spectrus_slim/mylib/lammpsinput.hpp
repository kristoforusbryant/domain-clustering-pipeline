#ifndef __LAMMPSINPUT_HPP__
#define __LAMMPSINPUT_HPP__

#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include "Vector3d.hpp"
#include "stringutils.hpp"

namespace mylib{

class LammpsInput{
public:
  LammpsInput() : num_atoms(0), num_atom_types(0), num_mol(0) {}
  LammpsInput(const char* fname) : num_atoms(0), num_atom_types(0) {read(fname);}
  void read(const char* fname){
    std::ifstream fin(fname);
    if (!fin.good()){
      std::cerr << "Unable to read file " << fname << std::endl << std::flush;
      exit(2);
    }
    char line[2048];


    while(fin.getline(line,2047).good()){
      if (is_empty_line(line)) continue;
      if (strstr(line,"atoms")!=NULL){
        std::sscanf(line,"%zu",&num_atoms);
        coords.resize(num_atoms);
        mol_id.resize(num_atoms);
        atom_type.resize(num_atoms);
        continue;
      }
      if (strstr(line,"atom types")!=NULL){
        std::sscanf(line,"%zu",&num_atom_types);
        continue;
      }

      if (strstr(line,"xlo xhi")!=NULL){
        std::sscanf(line,"%lf %lf",&(box[0].X),&(box[1].X));
        continue;
      }

      if (strstr(line,"ylo yhi")!=NULL){
        std::sscanf(line,"%lf %lf",&(box[0].Y),&(box[1].Y));
        continue;
      }

      if (strstr(line,"zlo zhi")!=NULL){
        std::sscanf(line,"%lf %lf",&(box[0].Z),&(box[1].Z));
        continue;
      }

      // read coordinates
      if (strstr(line,"Atoms")!=NULL) {
        for (size_t i = 0; i<num_atoms; ++i){
          fin.getline(line,2047);
          while (is_empty_line(line)) {fin.getline(line,2047);}
          double x,y,z;
          int id,molid,atype;
          std::sscanf(line,"%d %d %d %lf %lf %lf",&id,&molid,&atype,&x,&y,&z);
          --id;

          mol_id[id]=molid;
          atom_type[id]=atype;
          coords[id].X = x;
          coords[id].Y = y;
          coords[id].Z = z;
          if (atoms_of_mol.size() < molid) atoms_of_mol.resize(molid);
          atoms_of_mol[molid-1].push_back(id);
        }
        num_mol = *std::max_element(mol_id.begin(),mol_id.end());
        continue;
      }
    }
  }
public:
  size_t num_atoms;
  size_t num_atom_types;
  size_t num_mol;

  std::vector<Vector3d> coords;
  std::vector<size_t> mol_id;
  std::vector<size_t> atom_type;
  std::vector<std::vector<size_t> > atoms_of_mol;

  box_t box;

};

} // namespace mylib
#endif
