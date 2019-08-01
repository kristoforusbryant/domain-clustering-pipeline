
ifndef MPICXX
MPICXX=mpic++
endif

SLEPC_DIR=/gpfs/ysm/project/kbo6/conda_envs/spectrus
# Modified in 2019/07/17
PETSC_DIR=/gpfs/ysm/project/kbo6/conda_envs/spectrus
# Modified in 2019/07/17

#LFLAGS = -L ${PETSC_DIR}/${PETSC_ARCH}/lib -L ${SLEPC_DIR}/${SLEPC_ARCH}/lib
LFLAGS = -L ${SLEPC_DIR}/${SLEPC_ARCH}/lib -L ${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${SLEPC_DIR}/${SLEPC_ARCH}/lib,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib

INCLUDES = -I ${PETSC_DIR}/${PETSC_ARCH}/include -I ${PETSC_DIR}/include -I ./include -I ${SLEPC_DIR}/${SLEPC_ARCH}/include -I ${SLEPC_DIR}/include

CXXFLAGS=-c -Wall -O3 -Wno-sign-compare -Wno-unused-result
LDFLAGS= -lslepc -lpetsc 
SOURCES=spectrus_modified.cpp my_malloc.cpp krylov_shur.cpp k_medoids.cpp # modes_to_fluct_matrix.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=spectrus

ifndef SLEPC_DIR
$(error SLEPC_DIR is not set)
endif
ifndef PETSC_DIR
$(error PETSC_DIR is not set)
endif

all: $(SOURCES) $(EXECUTABLE)
    
$(EXECUTABLE): $(OBJECTS) 
	$(MPICXX) $(LFLAGS) $(LDFLAGS) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(MPICXX) $(INCLUDES) $(CXXFLAGS) $< -o $@

clean:
	rm -f *.o
	rm -f spectrus:
