/*
 * PdbIoWrapper.h
 *
 *  Created on: Oct 30, 2012
 *      Author: gpolles
 */

#ifndef PDBIOWRAPPER_H_
#define PDBIOWRAPPER_H_

#include <cstring>
#include <fstream>
#include <Vector3d.hpp>

class PdbAtomInfo{
public:
  PdbAtomInfo() : serial(0), altloc(' '),
                  chainId(' '), resSeq(0), iCode(' '), x(0.0), y(0.0), z(0.0),
                  occupancy(0.0), tempFactor(0.0) {
    strcpy(name,"");
    strcpy(resName,"");
    strcpy(element,"");
    strcpy(charge,"");
  }
public:
  int serial;
  char name[5];
  char altloc;
  char resName[5];
  char chainId;
  int resSeq;
  char iCode;
  double x,y,z;
  double occupancy;
  double tempFactor;
  char element[3];
  char charge[3];
};

class PdbIoWrapper
{
public:
  PdbIoWrapper() {}
  virtual ~PdbIoWrapper() {}

public:
  void openAsInput(const char* fname);
  void openAsOutput(const char* fname);
  void close();

  int parse_pdb_atom_line(char* line, PdbAtomInfo* info);
  void generate_pdb_atom_line(char* line, const PdbAtomInfo& info);

  bool nextAtom(PdbAtomInfo* info);
  bool writeAtom(const PdbAtomInfo& info);


private:
  std::fstream _stream;
  char _line[2000];
  char tmp_buffer[10];

};

class PdbFile
{
 public:
  PdbFile(){}
  PdbFile(const char* fname){ read(fname);}
  
  void read(const char* fname){
    clear();
    PdbIoWrapper pdb;
    pdb.openAsInput(fname);
    PdbAtomInfo info;
    while(pdb.nextAtom(&info)) {
      _coords.push_back( Vector3d(info.x,info.y,info.z) );
      _atom_names.push_back(info.name);
      _chain_ids.push_back(info.chainId);
    }
    pdb.close();
  }
  
  void clear(){
    _coords.resize(0);
    _atom_names.resize(0);
    _chain_ids.resize(0);
  }
  
  void close(){ clear(); }
  
  std::vector<Vector3d> getCoords() const { return _coords;}
  std::vector<std::string> getAtomNames() const { return _atom_names; }
  std::vector<char> getChainIDs() const { return _chain_ids; }
  
 private:
  std::vector<Vector3d> _coords;
  std::vector<std::string> _atom_names;
  std::vector<char> _chain_ids;
  
};


int PdbIoWrapper::parse_pdb_atom_line(char* line, PdbAtomInfo* info) {

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
  if (sscanf(tmp_buffer,"%s",info->name)){
    ++num_elements_read;
  }

  // 17             Character     altLoc       Alternate location indicator.
  if (len<17) return num_elements_read;
  info->altloc = *(line + 16);

  // 18 - 21        Residue name  resName      Residue name.
  if (len<21) return num_elements_read;
  p = line + 17;
  strncpy(tmp_buffer,p,4);
  tmp_buffer[4]='\0';
  if (sscanf(tmp_buffer,"%s",info->resName)){
    ++num_elements_read;
  }

  // 22             Character     chainID      Chain identifier.
  if (len<22) return num_elements_read;
  info->chainId = *(line + 21);

  // 23 - 26        Integer       resSeq       Residue sequence number.
  if (len<26) return num_elements_read;
  p = line+22;
  strncpy(tmp_buffer,p,4);
  tmp_buffer[4]='\0';
  if (sscanf(tmp_buffer,"%d",&(info->resSeq))){
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

  // 61 - 66        Real(6.2)    tempFactor    Temperature factor.
  if (len<66) return num_elements_read;
  p = line+60;
  strncpy(tmp_buffer,p,6);
  tmp_buffer[6]='\0';
  if (sscanf(tmp_buffer,"%lf",&(info->tempFactor))){
    ++num_elements_read;
  }

  //77 - 78        LString(2)   element       Element symbol, right-justified.
  if (len<54) return num_elements_read;
  p = line+76;
  strncpy(tmp_buffer,p,2);
  tmp_buffer[2]='\0';
  if (sscanf(tmp_buffer,"%s",info->element)){
    ++num_elements_read;
  }

  // 79 - 80        LString(2)   charge        Charge on the atom.
  if (len<80) return num_elements_read;
  p = line+78;
  strncpy(tmp_buffer,p,2);
  tmp_buffer[2]='\0';
  if (sscanf(tmp_buffer,"%s",info->charge)){
    ++num_elements_read;
  }

  return num_elements_read;
}

void PdbIoWrapper::openAsInput(const char* fname) {
  _stream.open(fname,std::fstream::in);
  if(!_stream.is_open()){
    std::cerr << "Error: could not open input file: " << fname << std::endl << std::flush;
  }
}

void PdbIoWrapper::openAsOutput(const char* fname) {
  _stream.open(fname,std::fstream::out);
  if(!_stream.is_open()){
    std::cerr << "Error: could not open output file: " << fname << std::endl << std::flush;
  }
}

void PdbIoWrapper::close() {
  _stream.close();
}

void PdbIoWrapper::generate_pdb_atom_line(char* line, const PdbAtomInfo& info) {
  sprintf(line,"ATOM  ");
  sprintf(&line[6],"%5d ",info.serial%100000);
  sprintf(&line[12],"%4s",info.name);
  sprintf(&line[16],"%c",info.altloc);
  sprintf(&line[17],"%4s",info.resName);
  sprintf(&line[21],"%c",info.chainId);
  sprintf(&line[22],"%4d",info.resSeq%10000);
  sprintf(&line[26],"%c   ",info.iCode);
  sprintf(&line[30],"%8.3lf%8.3lf%8.3lf",info.x,info.y,info.z);
  sprintf(&line[54],"%6.2lf%6.2lf",info.occupancy,info.tempFactor);
  sprintf(&line[54],"          %2s%2s",info.element,info.charge);


}

bool PdbIoWrapper::nextAtom(PdbAtomInfo* info) {
  while(_stream.getline(_line,256)){
    if(!strncmp (_line, "ATOM", 4)){
      parse_pdb_atom_line(_line,info);
      return true;
    }
  }
  return false;
}

bool PdbIoWrapper::writeAtom(const PdbAtomInfo& info) {
  generate_pdb_atom_line(_line,info);
  _stream << _line<< std::endl;
  return _stream.fail();
}

#endif /* PDBIOWRAPPER_H_ */
