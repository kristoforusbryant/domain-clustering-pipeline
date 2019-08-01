/*
 * PdbIoWrapper.h
 *
 *  Created on: Oct 30, 2012
 *      Author: gpolles
 */

#ifndef PDBIOWRAPPER_H_
#define PDBIOWRAPPER_H_

#include <string>
#include <cstring>
#include <fstream>
#include "Vector3d.hpp"
#include "stringutils.hpp"

namespace mylib {

class PdbAtomInfo{
public:
  friend std::ostream& operator <<(std::ostream &out, const PdbAtomInfo& p);
  friend std::istream& operator >>(std::istream &in, PdbAtomInfo& p);

  PdbAtomInfo() : serial(0), name("CA"), altloc(' '), res_name("GLY"),
                  chain_id(' '), res_id(0), iCode(' '), x(0.0), y(0.0), z(0.0),
                  occupancy(0.0), beta(0.0) {
  }
  std::string toAscii();
  void setCoords(const Vector3d& v){
    x = v.X; y = v.Y ; z = v.Z;
  }

public:
  int serial;
  std::string name;
  char altloc;
  std::string res_name;
  char chain_id;
  int res_id;
  char iCode;
  double x,y,z;
  double occupancy;
  double beta;
  std::string element;
  std::string charge;
};

/*
int serial;
  char name[5];
  char altloc;
  char res_name[5];
  char chain_id;
  int res_id;
  char iCode;
  double x,y,z;
  double occupancy;
  double beta;
  char element[3];
  char charge[3];
  */

// Parsing and writing functions

inline void generate_pdb_atom_line(char* line, const PdbAtomInfo& info) {
  sprintf(line,"ATOM  ");
  sprintf(&line[6],"%5d ",info.serial%100000);
  std::string nstr(" ");
  if(info.name.length() < 4 )
    nstr+= info.name;
  else
    nstr = info.name;
  sprintf(&line[12],"%-4s",nstr.substr(0,4).c_str());
  sprintf(&line[16],"%c",info.altloc);
  sprintf(&line[17],"%-4s",info.res_name.substr(0,4).c_str());
  sprintf(&line[21],"%c",info.chain_id);
  sprintf(&line[22],"%4d",info.res_id%10000);
  sprintf(&line[26],"%c   ",info.iCode);
  sprintf(&line[30],"%8.3lf%8.3lf%8.3lf",info.x,info.y,info.z);
  sprintf(&line[54],"%6.2lf%6.2lf",info.occupancy,info.beta);
  sprintf(&line[54],"          %2s%2s",info.element.substr(0,2).c_str(),info.charge.substr(0,2).c_str());


}

inline std::string generate_pdb_atom_line(const PdbAtomInfo& info) {
  char line[100];
  generate_pdb_atom_line(line,info);
  return std::string(line);
}

inline int parse_pdb_atom_line(char* line, PdbAtomInfo* info) {

  char tmp_buffer[10];
  *info = PdbAtomInfo();

  int num_elements_read = 0;

  int len = strlen(line);
  if (len<6) return 0;

  //  7 - 11        Integer       serial       Atom serial number.
  if (len<11) return num_elements_read;
  char* p = line+6;
  strncpy(tmp_buffer,p,5);
  tmp_buffer[5]='\0';
  if (sscanf(tmp_buffer,"%d",&(info->serial))){
    ++num_elements_read;
  }

  //13 - 16        Atom          name         Atom name.
  if (len<16) return num_elements_read;
  p = line + 12;
  strncpy(tmp_buffer,p,4);
  tmp_buffer[4]='\0';
  info->name = tmp_buffer;
  trim(info->name);
  if ((info->name).length()){
    ++num_elements_read;
  }

  // 17             Character     altLoc       Alternate location indicator.
  if (len<17) return num_elements_read;
  info->altloc = *(line + 16);

  // 18 - 21        Residue name  res_name      Residue name.
  if (len<21) return num_elements_read;
  p = line + 17;
  strncpy(tmp_buffer,p,4);
  tmp_buffer[4]='\0';
  info->res_name = tmp_buffer;
  trim(info->res_name);
  if ((info->res_name).length()){
    ++num_elements_read;
  }

  // 22             Character     chain_id      Chain identifier.
  if (len<22) return num_elements_read;
  info->chain_id = *(line + 21);

  // 23 - 26        Integer       res_id       Residue sequence number.
  if (len<26) return num_elements_read;
  p = line+22;
  strncpy(tmp_buffer,p,4);
  tmp_buffer[4]='\0';
  if (sscanf(tmp_buffer,"%d",&(info->res_id))){
    ++num_elements_read;
  }

  //27             AChar         iCode        Code for insertion of residues.
  if (len<27) return num_elements_read;
  info->iCode = *(line + 26);

  // 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.

  if (len<38) return num_elements_read;
  p = line+30;
  strncpy(tmp_buffer,p,8);
  tmp_buffer[8]='\0';
  if (sscanf(tmp_buffer,"%lf",&(info->x))){
    ++num_elements_read;
  }

  // 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
  if (len<46) return num_elements_read;
  p = line+38;
  strncpy(tmp_buffer,p,8);
  tmp_buffer[8]='\0';
  if (sscanf(tmp_buffer,"%lf",&(info->y))){
    ++num_elements_read;
  }

  // 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
  if (len<54) return num_elements_read;
  p = line+46;
  strncpy(tmp_buffer,p,8);
  tmp_buffer[8]='\0';
  if (sscanf(tmp_buffer,"%lf",&(info->z))){
    ++num_elements_read;
  }


  // 55 - 60        Real(6.2)    occupancy     Occupancy.
  if (len<60) return num_elements_read;
  p = line+54;
  strncpy(tmp_buffer,p,6);
  tmp_buffer[6]='\0';
  if (sscanf(tmp_buffer,"%lf",&(info->occupancy))){
    ++num_elements_read;
  }

  // 61 - 66        Real(6.2)    beta    Temperature factor.
  if (len<66) return num_elements_read;
  p = line+60;
  strncpy(tmp_buffer,p,6);
  tmp_buffer[6]='\0';
  if (sscanf(tmp_buffer,"%lf",&(info->beta))){
    ++num_elements_read;
  }

  //77 - 78        LString(2)   element       Element symbol, right-justified.
  if (len<54) return num_elements_read;
  p = line+76;
  strncpy(tmp_buffer,p,2);
  tmp_buffer[2]='\0';
  info->element = tmp_buffer;
  trim(info->element);
  if ((info->element).length()){
    ++num_elements_read;
  }

  // 79 - 80        LString(2)   charge        Charge on the atom.
  if (len<80) return num_elements_read;
  p = line+78;
  strncpy(tmp_buffer,p,2);
  tmp_buffer[2]='\0';
  info->charge = tmp_buffer;
  trim(info->charge);
  if ((info->charge).length()){
    ++num_elements_read;
  }

  return num_elements_read;
}

inline std::string PdbAtomInfo::toAscii(){
  return generate_pdb_atom_line(*this);
}

class PdbIoWrapper {
public:
  PdbIoWrapper() {strcpy(_line,"");}
  virtual ~PdbIoWrapper() {}

public:
  void openAsInput(const char* fname);
  void openAsOutput(const char* fname);
  void close();

  bool nextAtom(PdbAtomInfo* info);
  bool writeAtom(const PdbAtomInfo& info);


private:
  std::fstream _stream;
  char _line[2000];

};

class PdbFile {
 public:
  PdbFile(){}
  PdbFile(const char* fname){ read(fname);}

  void read(const char* fname){
    clear();
    PdbIoWrapper pdb;
    pdb.openAsInput(fname);
    PdbAtomInfo info;
    while(pdb.nextAtom(&info)) {
      coords.push_back( Vector3d(info.x,info.y,info.z) );
      serial.push_back(info.serial);
      atom_name.push_back(info.name);
      altloc.push_back(info.altloc);
      chain_id.push_back(info.chain_id);
      iCode.push_back(info.iCode);
      res_name.push_back(info.res_name);
      res_seq.push_back(info.res_id);
      occupancy.push_back(info.occupancy);
      beta.push_back(info.beta);
      element.push_back(info.element);
      charge.push_back(info.charge);
    }
    pdb.close();
  }

  void clear(){
    coords.resize(0);
    serial.resize(0);
    atom_name.resize(0);
    altloc.resize(0);
    chain_id.resize(0);
    iCode.resize(0);
    res_name.resize(0);
    res_seq.resize(0);
    occupancy.resize(0);
    beta.resize(0);
    element.resize(0);
    charge.resize(0);
  }

  void close(){ clear(); }

  size_t size() const { return coords.size(); }

  PdbAtomInfo operator[] (size_t i){
    PdbAtomInfo info;
    info.serial = serial[i];
    info.x = coords[i].X;
    info.y = coords[i].Y;
    info.z = coords[i].Z;
    info.altloc = altloc[i];
    info.chain_id = chain_id[i];
    info.iCode = iCode[i];
    info.res_id = res_seq[i];
    info.occupancy = occupancy[i];
    info.beta = beta[i];
    info.name = atom_name[i];
    info.res_name = res_name[i];
    info.element = element[i];
    info.charge = charge[i];
    return info;
  }

 public:
  std::vector<Vector3d> coords;
  std::vector<int> serial;
  std::vector<std::string> atom_name;
  std::vector<char> altloc;
  std::vector<char> chain_id;
  std::vector<char> iCode;
  std::vector<std::string> res_name;
  std::vector<int> res_seq;
  std::vector<float> occupancy;
  std::vector<float> beta;
  std::vector<std::string> element;
  std::vector<std::string> charge;


};



inline void PdbIoWrapper::openAsInput(const char* fname) {
  _stream.open(fname,std::fstream::in);
  if(!_stream.is_open()){
    std::cerr << "Error: could not open input file: " << fname << std::endl << std::flush;
  }
}

inline void PdbIoWrapper::openAsOutput(const char* fname) {
  _stream.open(fname,std::fstream::out);
  if(!_stream.is_open()){
    std::cerr << "Error: could not open output file: " << fname << std::endl << std::flush;
  }
}

inline void PdbIoWrapper::close() {
  _stream.close();
}



inline bool PdbIoWrapper::nextAtom(PdbAtomInfo* info) {
  while(_stream.getline(_line,256)){
    if(!strncmp (_line, "ATOM", 4)){
      parse_pdb_atom_line(_line,info);
      return true;
    }
  }
  return false;
}

inline bool PdbIoWrapper::writeAtom(const PdbAtomInfo& info) {
  generate_pdb_atom_line(_line,info);
  _stream << _line<< std::endl;
  return _stream.fail();
}

inline std::ostream& operator <<(std::ostream &out, const PdbAtomInfo& p){
  out << generate_pdb_atom_line(p);
  return out;
}
inline std::istream& operator >>(std::istream &in, PdbAtomInfo& p){
  char line[256];
  in.getline(line,255);
  parse_pdb_atom_line(line,&p);
  return in;
}

}//namespace mylib
#endif /* PDBIOWRAPPER_H_ */
