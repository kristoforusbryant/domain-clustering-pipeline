#ifndef __VERLETLIST_HPP__
#define __VERLETLIST_HPP__

#include <cmath>
#include <vector>
#include <mylib/Vector3d.hpp>

// TODO: use a mask instead of computing the cells distance every time.

namespace mylib{
const size_t MAXCELLS=1000000;
inline void FALLBACK_MAXCELL(size_t& a, size_t& b, size_t& c) {a=100;b=100;c=100;}


inline std::vector<size_t> cutoff_to_ncells(double cutoff, const box_t& box);



inline int iabs(int a) {return a<0? -a :a;}

class VerletList3D{
public:
  VerletList3D(): num_points(0) {dim[0] = 0; dim[1] = 0; dim[2] = 0;}
  VerletList3D(const std::vector<Vector3d>& crd, const box_t& pbox, double cutoff): num_points(0) {
    std::vector<size_t> ncells = cutoff_to_ncells(cutoff,  pbox);
    generateList(crd,pbox,ncells[0],ncells[1],ncells[2]);
  }
  VerletList3D(const std::vector<Vector3d>& crd, const box_t& pbox, const std::vector<size_t>& ncells): num_points(0) {
    generateList(crd,pbox,ncells[0],ncells[1],ncells[2]);
  }
  VerletList3D(const std::vector<Vector3d>& crd, const box_t& pbox, size_t n_cell_x, size_t n_cell_y, size_t n_cell_z): num_points(0) {
    generateList(crd,pbox,n_cell_x,n_cell_y,n_cell_z);
  }
  void generateList(const std::vector<Vector3d>& crd, const box_t& pbox, size_t n_cell_x, size_t n_cell_y, size_t n_cell_z){
    box = pbox;
    size_t tot_cell = n_cell_x*n_cell_y*n_cell_z;
    if (tot_cell>MAXCELLS){
      std::cerr << "VerletList3D: generateList(): exceded max number of cells MAXCELLS=("<<MAXCELLS<<"). Falling back to 100,100,100.\n" << std::flush;
      FALLBACK_MAXCELL(n_cell_x,n_cell_y,n_cell_z);
    }

    box_sizes = (box[1]-box[0]);
    intv = Vector3d(box_sizes.X/n_cell_x,box_sizes.Y/n_cell_y,box_sizes.Z/n_cell_z);

    // ALLOCATION of CELLS

    cells.resize(n_cell_x);
    for (size_t i = 0; i<n_cell_x; ++i){
      cells[i].resize(n_cell_y);
      for (size_t j = 0; j<n_cell_y; ++j){
        cells[i][j].resize(n_cell_z);
        for (size_t k = 0; k<n_cell_z; ++k){
          cells[i][j][k].resize(0);
        }
      }
    }
    dim[0] = n_cell_x; dim[1] = n_cell_y; dim[2] = n_cell_z;


    for (size_t i = 0; i<crd.size(); ++i){
      if (crd[i].X < box[0].X || crd[i].X > box[1].X ||
          crd[i].Y < box[0].Y || crd[i].Y > box[1].Y ||
          crd[i].Z < box[0].Z || crd[i].Z > box[1].Z ){
        std::cerr << "VerletList3D: generateList(): coordinate " << i << "out of box\n" << std::flush;
        continue;
      }
      size_t cx = int( floor( (crd[i].X-box[0].X)/intv.X ) ) % n_cell_x;
      size_t cy = int( floor( (crd[i].Y-box[0].Y)/intv.Y ) ) % n_cell_y;
      size_t cz = int( floor( (crd[i].Z-box[0].Z)/intv.Z ) ) % n_cell_z;
      cells[cx][cy][cz].push_back(i);
      ++num_points;
    }
  }

  void preallocate(size_t n_cell_x, size_t n_cell_y, size_t n_cell_z, size_t n_items ){
    size_t tot_cell = n_cell_x*n_cell_y*n_cell_z;
    if (tot_cell>MAXCELLS){
      std::cerr << "VerletList3D: generateList(): exceded max number of cells MAXCELLS=("<<MAXCELLS<<"). Falling back to 100,100,100.\n" << std::flush;
      FALLBACK_MAXCELL(n_cell_x,n_cell_y,n_cell_z);
    }
    cells.resize(n_cell_x);
    for (size_t i = 0; i<n_cell_x; ++i){
      cells[i].resize(n_cell_y);
      for (size_t j = 0; j<n_cell_y; ++j){
        cells[i][j].resize(n_cell_z);
        for (size_t k = 0; k<n_cell_y; ++k){
          cells[i][j][k].reserve(n_items);
        }
      }
    }
    dim[0] = n_cell_x; dim[1] = n_cell_y; dim[2] = n_cell_z;
  }

  double getAveragePointsPerCell(){
    double totcell = dim[0]*dim[1]*dim[2];
    if (totcell==0){
      std::cerr << "VerletList3D: getAveragePointsPerCell(): null dimension\n" << std::flush;
      return 0;
    }
    return double(num_points)/totcell;
  }

  // TODO :  bool isNeighbor() ? {  }

  std::vector<size_t> getNeighborsPBC(const std::vector<Vector3d>& crd, size_t i, double cutoff){
    size_t nx = ceil(cutoff/intv.X);
    size_t ny = ceil(cutoff/intv.Y);
    size_t nz = ceil(cutoff/intv.Z);

    int cx = int( floor( (crd[i].X-box[0].X)/intv.X ) ) % dim[0];
    int cy = int( floor( (crd[i].Y-box[0].Y)/intv.Y ) ) % dim[1];
    int cz = int( floor( (crd[i].Z-box[0].Z)/intv.Z ) ) % dim[2];

    std::vector<size_t> neigh;
    double cutoff2 = cutoff*cutoff;
    Vector3d intv2 = intv*intv;

    for (int ii = cx - int(nx) ; ii <= cx + int(nx); ++ii){
      size_t ix = ii>=0 ? ii % dim[0] : ii + dim[0]*( (-ii/dim[0]) +1 );
      for (int jj = cy - int(ny) ; jj <= cy + int(ny); ++jj){
        size_t iy = jj>=0 ? jj % dim[1] : jj + dim[1]*( (-jj/dim[1]) +1 );
        for (int kk = cz - int(nz) ; kk <= cz + int(nz); ++kk){
          size_t iz = kk>=0 ? kk % dim[2] : kk + dim[2]*( (-kk/dim[2]) +1 );
          int dx = iabs(cx-ii); int dy = iabs(cy-jj); int dz = iabs(cz-kk);
          double cdist=0.0;
          if(dx > 1) cdist+=SQ(dx-1)*intv2.X;
          if(dy > 1) cdist+=SQ(dy-1)*intv2.Y;
          if(dz > 1) cdist+=SQ(dz-1)*intv2.Z;
          if ( cdist > cutoff2 ) continue;
          // now check all potential neighbors against our guy
          std::vector<size_t>& ccell = cells[ix][iy][iz];
          for (size_t j = 0; j<ccell.size(); ++j){
            if (ccell[j]==i) continue;
            if (distanceSqPbc(crd[i],crd[ccell[j]],box_sizes) < cutoff2) neigh.push_back(ccell[j]);
          }
        }
      }
    }


    return neigh;
  }

  std::vector<size_t> getNeighbors(const std::vector<Vector3d>& crd, size_t i, double cutoff){
    size_t nx = ceil(cutoff/intv.X);
    size_t ny = ceil(cutoff/intv.Y);
    size_t nz = ceil(cutoff/intv.Z);

    int cx = int( floor( (crd[i].X-box[0].X)/intv.X ) ) ;
    int cy = int( floor( (crd[i].Y-box[0].Y)/intv.Y ) ) ;
    int cz = int( floor( (crd[i].Z-box[0].Z)/intv.Z ) ) ;

    // create mask for cells
    bool mask[nx+1][ny+1][nz+1];
    for (size_t dx = 0 ; dx <= nx; ++dx){
      for (size_t dy = 0 ; dy <= ny; ++dy){
        for (size_t dz = 0 ; dz <= nz; ++dz){
          double cdist=0.0;
          if(dx > 1) cdist+=SQ(dx-1)*intv2.X;
          if(dy > 1) cdist+=SQ(dy-1)*intv2.Y;
          if(dz > 1) cdist+=SQ(dz-1)*intv2.Z;
          if ( cdist > cutoff2 ) mask[dx][dy][dz] = false;
          else mask[dx][dy][dz] = true;
        }
      }
    }

    std::vector<size_t> neigh;
    double cutoff2 = cutoff*cutoff;
    Vector3d intv2 = intv*intv;


    for (int ii = cx - int(nx) ; ii <= cx + int(nx); ++ii){
      if (ii<0 || ii>=dim[0]) continue;
      size_t ix = ii;
      for (int jj = cy - int(ny) ; jj <= cy + int(ny); ++jj){
        if (jj<0 || jj>=dim[1]) continue;
        size_t iy = jj;
        for (int kk = cz - int(nz) ; kk <= cz + int(nz); ++kk){
          if (kk<0 || kk>=dim[2]) continue;
          size_t iz = kk;
          int dx = iabs(cx-ix); int dy = iabs(cy-iy); int dz = iabs(cz-iz);
          /*double cdist=0.0;
          if(dx > 1) cdist+=SQ(dx-1)*intv2.X;
          if(dy > 1) cdist+=SQ(dy-1)*intv2.Y;
          if(dz > 1) cdist+=SQ(dz-1)*intv2.Z;
          if ( cdist > cutoff2 ) continue;
          */
          
          if (mask[dx][dy][dz]==false) continue;
          
          // now check all potential neighbors against our guy
          std::vector<size_t>& ccell = cells[ix][iy][iz];
          for (size_t j = 0; j<ccell.size(); ++j){
            if (ccell[j]==i) continue;
            if (distanceSQ(crd[i],crd[ccell[j]]) < cutoff2) neigh.push_back(ccell[j]);
          }
        }
      }
    }
    return neigh;
  }

  std::vector<std::vector<size_t> > getNeighborsMap(const std::vector<Vector3d>& crd, double cutoff){
    std::vector<std::vector<size_t> > neighs(crd.size());
    size_t nx = ceil(cutoff/intv.X);
    size_t ny = ceil(cutoff/intv.Y);
    size_t nz = ceil(cutoff/intv.Z);
    double cutoff2 = cutoff*cutoff;
    Vector3d intv2 = intv*intv;

    // create mask for cells
    bool mask[nx+1][ny+1][nz+1];
    for (size_t dx = 0 ; dx <= nx; ++dx){
      for (size_t dy = 0 ; dy <= ny; ++dy){
        for (size_t dz = 0 ; dz <= nz; ++dz){
          double cdist=0.0;
          if(dx > 1) cdist+=SQ(dx-1)*intv2.X;
          if(dy > 1) cdist+=SQ(dy-1)*intv2.Y;
          if(dz > 1) cdist+=SQ(dz-1)*intv2.Z;
          if ( cdist > cutoff2 ) mask[dx][dy][dz] = false;
          else mask[dx][dy][dz] = true;
        }
      }
    }

    for (size_t i = 0; i < crd.size(); ++i){
      
      int cx = int( floor( (crd[i].X-box[0].X)/intv.X ) ) ;
      int cy = int( floor( (crd[i].Y-box[0].Y)/intv.Y ) ) ;
      int cz = int( floor( (crd[i].Z-box[0].Z)/intv.Z ) ) ;

      for (int ii = cx - int(nx) ; ii <= cx + int(nx); ++ii){
        if (ii<0 || ii>=dim[0]) continue;
        size_t ix = ii;
        for (int jj = cy - int(ny) ; jj <= cy + int(ny); ++jj){
          if (jj<0 || jj>=dim[1]) continue;
          size_t iy = jj;
          for (int kk = cz - int(nz) ; kk <= cz + int(nz); ++kk){
            if (kk<0 || kk>=dim[2]) continue;
            size_t iz = kk;
            int dx = iabs(cx-ix); int dy = iabs(cy-iy); int dz = iabs(cz-iz);
            if (mask[dx][dy][dz]==false) continue;
            // now check all potential neighbors against our guy
            std::vector<size_t>& ccell = cells[ix][iy][iz];
            for (size_t j = 0; j<ccell.size(); ++j){
              if (ccell[j]>=i) continue;
              if (distanceSQ(crd[i],crd[ccell[j]]) < cutoff2) {
                neighs[i].push_back(ccell[j]);
                neighs[ccell[j]].push_back(i);
              }
            }
          }
        }
      }
    }
    return neighs;
  }

public:
  std::vector<std::vector<std::vector<std::vector<size_t> > > > cells;
  Vector3d clen;
  box_t box;
  size_t num_points;
  size_t dim[3];
  Vector3d box_sizes;
  Vector3d intv;
};

inline box_t get_enclosing_box(std::vector<Vector3d>& crd){
  box_t b;
  if(crd.size()==0) return b;
  b[0] = crd[0];
  b[1] = crd[0];
  for (size_t i = 0; i < crd.size(); ++i){
    if(b[1].X < crd[i].X) b[1].X = crd[i].X;
    if(b[1].Y < crd[i].Y) b[1].Y = crd[i].Y;
    if(b[1].Z < crd[i].Z) b[1].Z = crd[i].Z;
    if(b[0].X > crd[i].X) b[0].X = crd[i].X;
    if(b[0].Y > crd[i].Y) b[0].Y = crd[i].Y;
    if(b[0].Z > crd[i].Z) b[0].Z = crd[i].Z;
  }
  return b;
}

inline std::vector<size_t> cutoff_to_ncells(double cutoff, const box_t& box){
  double x_size = box[1].X-box[0].X;
  double y_size = box[1].Y-box[0].Y;
  double z_size = box[1].Z-box[0].Z;
  std::vector<size_t> nc(3);
  nc[0] = floor(x_size/cutoff*2);
  nc[1] = floor(y_size/cutoff*2);
  nc[2] = floor(z_size/cutoff*2);
  return nc;
}

} // namespace mylib
#endif
