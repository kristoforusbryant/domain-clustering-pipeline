/********* datafile.phh **************
Class for generic ascii datafile I/O

Author: Guido Polles
Date: 30/05/2014
*************************************/

#ifndef __DATAFILE_HPP__
#define __DATAFILE_HPP__

#define MAX_LINE_LEN 10000

#include <iostream>
#include <vector>
#include <cstring>
#include <cstdio>
#include <string>
#include <fstream>

#include "stringutils.hpp"


namespace mylib{




class Datafile{
 public:
  Datafile() : _num_col(0) {}
  Datafile(const char* fname){ read( fname ); }

  void read(const char* fname){
    std::ifstream fin(fname);
    char line[MAX_LINE_LEN];
    std::vector<std::string> svals;
    std::vector<double> lvals;
    strcpy(line,"");

    // read first line and sets once and for all the number of fields
    while( strlen(line) == 0 || line[0] == '#') fin.getline(line,MAX_LINE_LEN-1);
    svals = strsplit(line," \n\t");
    _num_col = svals.size();

    _data.push_back(std::vector<double>(_num_col));
    lvals.resize(svals.size());
    double val;
    for (size_t i=0 ; i<svals.size(); ++i){
      std::sscanf(svals[i].c_str(),"%lf",&val);
      _data[0][i] = val;
    }

    // read other lines
    while (fin.getline(line,MAX_LINE_LEN-1).good()){
      if ( strlen(line) == 0 || line[0] == '#') continue;
      svals = strsplit(line," \n\t");
      if (_num_col!=svals.size()){
        std::cerr<< "Datafile: read() error: num_fields!=svals.size() <" << line << ">" << std::endl;
        return;
      }
      _data.push_back(std::vector<double>(_num_col));
      for (size_t i=0 ; i<_num_col; ++i){
        if(std::sscanf(svals[i].c_str(),"%lf",&val)!=1){
          std::cerr<< "Datafile: read() error: sscanf() failed <" << svals[i].c_str() << ">" << std::endl;
          return;
        }
        _data.back()[i] = val;
      }
    }
  }

  std::vector<double>& operator[] (size_t i){ return _data[i]; }
  const std::vector<double>& operator[] (size_t i) const{ return _data[i]; }

  void setData(std::vector< std::vector<double> >& d){ _data = d; }
  std::vector< std::vector<double> >& getData(){ return _data; }
  void clearData() { _data.resize(0); _num_col = 0;}

  void write(const char* fname){
    std::ofstream fout(fname);
    if(!fout.good()){
      std::cerr<< "Datafile: write() error: unable to open file <" << fname << ">" << std::endl;
      return;
    }
    for (size_t i=0; i<_data.size(); ++i){
      for (size_t j=0; j<_num_col; ++j){
        fout << _data[i][j] << " ";
      }
      fout << std::endl;
    }
  }



  size_t size(){ return _data.size(); }
  size_t ncol(){ return _num_col; }
 private:
  std::vector< std::vector<double> > _data;
  size_t _num_col;


  friend std::ostream& operator <<(std::ostream &out, const Datafile& p);
};


inline std::ostream& operator <<(std::ostream &out, const Datafile& p){
  for (size_t i=0; i<p._data.size(); ++i){
    for (size_t j=0; j<p._num_col; ++j){
      out << p._data[i][j] << " ";
    }
    out << std::endl;
  }
  return out;
}

}// namespace mylib

#endif
