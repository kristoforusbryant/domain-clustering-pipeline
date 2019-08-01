#ifndef __LAMMPSTRAJ_HPP__
#define __LAMMPSTRAJ_HPP__

#include <vector>
#include <fstream>
#include <cstring>

#include <sstream>

#include "Vector3d.hpp"
#include "stringutils.hpp"

namespace mylib{

class LammpsFrame{
public:
  LammpsFrame(): num_atoms(0), timestep(0) {}

  Vector3d box_sizes() {return box[1]-box[0];}
  std::vector<Vector3d> coords;
  std::vector<Vector3d> image_idx;
  size_t num_atoms;
  size_t timestep;
  box_t box;

};

class LammpsTraj: public std::vector<LammpsFrame> {
public:
  LammpsTraj() {}
  LammpsTraj(const char* fname) {read(fname);}
public:
  void read(const char* fname){
    std::ifstream fin(fname);
    if (!fin.good()){
      std::cerr << "Unable to read file " << fname << std::endl << std::flush;
      exit(1);
    }
    char line[2048];

    while(fin.getline(line,2047).good()){
      while (is_empty_line(line)) {fin.getline(line,2047);}
      if (strstr(line,"ITEM: TIMESTEP")){
        push_back(LammpsFrame());
        fin.getline(line,2047);
        size_t ts = 0;
        std::sscanf(line,"%zu",&ts);
        back().timestep = ts;
        continue;
      }
      if (strstr(line,"ITEM: NUMBER OF ATOMS")){
        fin.getline(line,2047);
        size_t na = 0;
        std::sscanf(line,"%zu",&na);
        back().num_atoms = na;
        back().coords.resize(na);
        back().image_idx.resize(na);
        continue;
      }
      if (strstr(line,"ITEM: BOX BOUNDS")){
        for (int i=0; i<3; ++i){
          fin.getline(line,2047);
          std::sscanf(line,"%lf %lf", &(back().box[0][i]),&(back().box[1][i]));
        }
        continue;
      }

      if (strstr(line,"ITEM: ATOMS")){
        int posx,posy,posz =-1;
        int posix,posiy,posiz =-1;
        int posid= -1;
        std::stringstream ss(line+11);
        std::string s;
        int idx = 0;
        while(!(ss>>s).fail()){
          if(s.compare("id")==0) posid = idx;
          else if(s.compare("x")==0 || s.compare("ux")==0) posx = idx;
          else if(s.compare("y")==0 || s.compare("uy")==0) posy = idx;
          else if(s.compare("z")==0 || s.compare("uz")==0) posz = idx;
          else if(s.compare("ix")==0) posix = idx;
          else if(s.compare("iy")==0) posiy = idx;
          else if(s.compare("iz")==0) posiz = idx;
          ++idx;
        }

        if (posid==-1){
          std::cerr << "LammpsTraj::read(\""<<fname<<"\") error: Trajectory with no ids." << std::endl << std::flush;
          return;
        }

        if (posx==-1 || posy==-1 || posz==-1){
          std::cerr << "LammpsTraj::read(\""<<fname<<"\") error: Trajectory with no coordinates I can read." << std::endl << std::flush;
          return;
        }


        bool images_flag = (posix!=-1 && posiy!=-1 && posiz!=-1);


        for (size_t i=0; i<back().num_atoms; ++i){
          fin.getline(line,2047);
          std::stringstream ssl(line);
          double val[idx];
          for ( size_t j = 0; j<idx; ++j) ssl>>val[j];
          size_t id = size_t(val[posid])-1;
          back().coords[id].X = val[posx];
          back().coords[id].Y = val[posy];
          back().coords[id].Z = val[posz];
          if (images_flag){
            back().image_idx[id][0] = int(val[posix]);
            back().image_idx[id][1] = int(val[posiy]);
            back().image_idx[id][2] = int(val[posiz]);
          }
        }
        continue;
      }
    }
  }

};



} // namespace mylib

inline std::istream& operator>>(std::istream& fin, mylib::LammpsFrame& f){
  char line[2048];

  while(fin.getline(line,2047).good()){

    while (mylib::is_empty_line(line)) {
      fin.getline(line,2047);
    }

    if (strstr(line,"ITEM: TIMESTEP")){
      fin.getline(line,2047);
      size_t ts = 0;
      std::sscanf(line,"%zu",&ts);
      f.timestep = ts;
      continue;
    }

    if (strstr(line,"ITEM: NUMBER OF ATOMS")){
      fin.getline(line,2047);
      size_t na = 0;
      std::sscanf(line,"%zu",&na);
      f.num_atoms = na;
      f.coords.resize(na);
      f.image_idx.resize(na);
      continue;
    }

    if (strstr(line,"ITEM: BOX BOUNDS")){
      for (int i=0; i<3; ++i){
        fin.getline(line,2047);
        std::sscanf(line,"%lf %lf", &(f.box[0][i]),&(f.box[1][i]));
      }
      continue;
    }

    if (strstr(line,"ITEM: ATOMS")){
      int posx,posy,posz =-1;
      int posix,posiy,posiz =-1;
      int posid= -1;
      std::stringstream ss(line+11);
      std::string s;
      int idx = 0;
      while(!(ss>>s).fail()){
        if(s.compare("id")==0) posid = idx;
        else if(s.compare("x")==0 || s.compare("ux")==0) posx = idx;
        else if(s.compare("y")==0 || s.compare("uy")==0) posy = idx;
        else if(s.compare("z")==0 || s.compare("uz")==0) posz = idx;
        else if(s.compare("ix")==0) posix = idx;
        else if(s.compare("iy")==0) posiy = idx;
        else if(s.compare("iz")==0) posiz = idx;
        ++idx;
      }

      if (posid==-1){
        std::cerr << "operator >>(istream, LammpsTraj) error: Trajectory with no ids." << std::endl << std::flush;
        return fin;
      }

      if (posx==-1 || posy==-1 || posz==-1){
        std::cerr << "operator >>(istream, LammpsTraj) error: Trajectory with no coordinates I can read." << std::endl << std::flush;
        return fin;
      }


      bool images_flag = (posix!=-1 && posiy!=-1 && posiz!=-1);


      for (size_t i=0; i<f.num_atoms; ++i){
        fin.getline(line,2047);
        std::stringstream ssl(line);
        double val[idx];
        for ( size_t j = 0; j<idx; ++j) ssl>>val[j];
        size_t id = size_t(val[posid])-1;
        f.coords[id].X = val[posx];
        f.coords[id].Y = val[posy];
        f.coords[id].Z = val[posz];
        if (images_flag){
          f.image_idx[id][0] = int(val[posix]);
          f.image_idx[id][1] = int(val[posiy]);
          f.image_idx[id][2] = int(val[posiz]);
        }
      }
      return fin;
    }
  }
  return fin;
}

#endif
