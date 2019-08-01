//============================================================================
// Name        : ElasticNetworkDynamics.cpp
// Author      : Guido Polles
// Version     :
// Copyright   : 
// Description : CARE! It is an ANM simulation, not an ENM
//============================================================================

// Compile with: icpc -std=c++11 -openmp ElasticNetworkDynamics.cpp

#include <iostream>
#include <random>
#include <omp.h>
#include <cmath>
#include <map>
#include <fstream>
#include <sstream>
#include <chrono>


#include "mylib/Vector3d.hpp"
#include "mylib/contactmap.hpp"


typedef  double real_t;
typedef mylib::Vector3d vec3;


using namespace std;

int omp_n_proc = 1; // default number of processors;


// normal random number generator
template <class T>
class NormalRNG
{
public:
  NormalRNG(T mean, T sigma, int seed = time(0))
      : distribution(normal_distribution<T>(mean,sigma)),
        generator(std::mt19937(seed))
    {}
  T operator() () { return distribution(generator);}
private:
  normal_distribution<T> distribution;
  mt19937 generator;
  double padding[2]; //avoid false sharing - probably not useful or even wrong.
};

vector<NormalRNG<double>* > omp_rng;


typedef std::map<std::string,std::string> dict;
dict read_parameters(const char* fname){
  dict params;
  ifstream fin(fname);
  if (!fin.good()){
    cerr << "Warning: cannot read parameter file " << fname << endl;
    cerr << "Falling back to default values" << endl;
    return params;
  }
  string key,val;
  while(!(fin >> key >> val).fail()) params[key]=val;
  return params;
}


void compute_forces(vector<vec3>& f,
                    const vector<vec3>& x,
                    const mylib::ContactMap& cm,
                    const vector<int>& n_neigh,
                    const vector<vector<vec3> >& d0)
{
  int n_atoms = x.size();
  //static int static_chunk = n_atoms/omp_n_proc +1;

  for (int i = 0; i < n_atoms; ++i){
      f[i].X=0;f[i].Y=0;f[i].Z=0;
  }

  #pragma omp parallel for num_threads(omp_n_proc) schedule(static)
  for (int i = 0; i < n_atoms; ++i){
    for (size_t j = 0; j < n_neigh[i] ; ++j){
      size_t k = cm[i][j];
      vec3 d = d0[i][j] - (x[k]-x[i]);
      f[i] -= d;
    }
  }
}

void update_positions(vector<vec3>& x, const vector<vec3>& f, double dt){

  int n_atoms = x.size();
  const double sqdt=sqrt(2*dt);

  #pragma omp parallel for num_threads(omp_n_proc) schedule(static)
  for (int i = 0; i < n_atoms; ++i){
    int iproc = omp_get_thread_num();
    NormalRNG<double>& rng = *(omp_rng[iproc]);
    vec3 gauNoise(rng(),rng(),rng());
    x[i] += dt*f[i]+sqdt*gauNoise;
  }
}

//int stoi(string s){std::istringstream iss(s);
//                   int i; iss >> i; return i;}
//double stof(string s){std::istringstream iss(s);
//                   double i; iss >> i; return i;}

int main() {
  
  // read parameters from file or set to default
  dict param = read_parameters("params.input");
  size_t n_steps    = param.count("steps")  ? stoi(param["steps"])  : 1000;
  double cutoff     = param.count("cutoff") ? stof(param["cutoff"]) : 10.0; // interaction cutoff
  double dt         = param.count("dt")     ? stof(param["dt"])     : 0.01;
  size_t write_steps= param.count("write")  ? stoi(param["write"])  :  100;
  string ifname     = param.count("infile") ? param["infile"]       : "in.xyz";
  string ofname     = param.count("outfile")? param["outfile"]      : "out.xyz";
  size_t skip       = param.count("skip")   ? stoi(param["skip"])   : 1;
  omp_n_proc        = param.count("nproc")  ? stoi(param["nproc"])  : 1;


  // TODO: write parameters
  cout << "steps: " << n_steps << endl;
  cout << "cutoff: " << cutoff << endl;
  cout << "dt: " << dt << endl;
  cout << "write: " << write_steps << endl;
  cout << "infile: " << ifname << endl;
  cout << "outfile: " << ofname << endl;
  cout << "skip: " << skip << endl;
  cout << "nproc: " << omp_n_proc << endl;

  // prepare random generators with different seeds for each processor
  omp_rng.resize(omp_n_proc);
  for (size_t i = 0; i < omp_n_proc; ++i){
    omp_rng[i] = new NormalRNG<double>(0.0,1.0,time(0)^i) ;
  }

  // read starting configuration in xyz format - the vmd one
  char buffer[256];
  ifstream inf(ifname.c_str());
  inf.getline(buffer,255);
  size_t n_atoms = atoi(buffer);
  inf.getline(buffer,255);
  vector<vec3> x(n_atoms);
  for (size_t i = 0; i < n_atoms; ++i){
    string aname; // not gonna use it
    inf.getline(buffer,255);
    std::istringstream iss(buffer);
    iss >> aname >> x[i];
  }

  // generate interactions map
  mylib::ContactMap cm(x,cutoff);
  cout << "finished with cm" << endl;

  // save starting distances
  vector<int> n_neigh(n_atoms,0);
  vector<vector<vec3> > d0(n_atoms);
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < n_atoms; ++i){
    ++n_neigh[i];
    for (size_t j = 0; j < cm[i].size(); ++j){
      size_t k = cm[i][j];
      d0[i].push_back(x[k]-x[i]);
    }
  }

  cout << "distances calculated" << endl;

  // write first frame
  size_t wrt_atoms = (n_atoms-1)/skip +1;
  ofstream outf(ofname.c_str());
  outf << n_atoms << endl << endl;
  for (size_t i = 0; i < n_atoms; i+=skip){
    outf << "CA " << x[i] << endl;
  }


  // execution time stuff
  std::chrono::time_point<std::chrono::system_clock> start,end;
  std::chrono::duration<double> elapsed_seconds;
  start = std::chrono::system_clock::now();

  // main loop
  vector<vec3> f = vector<vec3>(n_atoms,vec3(0,0,0)); // forces
  for (size_t t = 1; t <= n_steps; ++t){
    compute_forces(f,x,cm,n_neigh,d0);
    update_positions(x,f,dt);
    if(t%write_steps == 0){
      outf << n_atoms << endl << endl;
      for (size_t i = 0; i < n_atoms; i+=skip){
        outf << "CA " << x[i] << endl;
      }
      cout << "step: " << t << endl;
    }
  }
  outf.close();

  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  cout << ">>>> Simulation time (s):" << elapsed_seconds.count() << endl;

  // clean up random generators
  for (size_t i = 0; i < omp_n_proc; ++i){
    delete omp_rng[i];
  }

  return 0;
}
