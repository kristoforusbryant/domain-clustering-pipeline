#include <iostream>
#include <algorithm>
#include "mylib/mylib.h"

using namespace std;

// sort class
class sort_class{
  public:
  sort_class(vector<mylib::Vector3d>& coords, int from) : crd(&coords), i(from) {}
  bool operator() (int j, int k){ return mylib::distance( crd->at(i),crd->at(j) ) < mylib::distance( crd->at(i),crd->at(k) );}
  vector <mylib::Vector3d>* crd;
  int i;
};

  
int smart_decimation(vector<mylib::Vector3d>& coords, int decimation, int cutoff, vector<bool>& keep) {
  // return value = -1 >> error: impossible to decimate, probably because 
  //                      i)  cutoff is to stringent,
  //                      ii) decimation >= # beads.
  // otherwise         >> return value = # retained beads
  
  
  
  keep.resize(coords.size(), false);
  vector<bool> visited(coords.size(), false);
    
  // generate contact map for given cutoff
  mylib::ContactMap cm(coords,cutoff);
  
  // check we have enough contacts
  for (int i = 0; i < cm.size(); ++i) {
    if (cm[i].size() < decimation) {
      // cerr << "Error: Number of neighs for atom " << i << " (" << cm[i].size() << ") is less than decimation factor." << endl;
      return -1;
    }
  }
 
  // sort contact map
  for (int i = 0; i < cm.size(); ++i) 
    std::sort(cm[i].begin(), cm[i].end(), sort_class(coords, i));
  

  int i = 0; /* current bead */
  int keep_count = 0;
  while (i != -1) {
    // NB: the first bead (with index 0) is always retained
    keep[i] = true;
    keep_count++;
    visited[i] = true;
    
    // the closest n=decimation beads are not retained
    for (int j = 0; j < decimation; ++j) {
      visited[ cm[i][j] ] = true; 
    }

    // find next atom to retain as the next closest to i which was not visited
    i = -1;
    for (int j = decimation; j < cm[i].size(); ++j) {
      if (!visited[ cm[i][j] ]) {
        i = cm[i][j];
        break;
      }
    }

    // if we're out of neighbors, get the first non visited site 
    if (i == -1) {
      for (int j = 0; j < coords.size(); ++j) {
        if (!visited[j]) {
          i = j;
          break;
        }
      }
    }
    
  }

  
  return keep_count;
}
