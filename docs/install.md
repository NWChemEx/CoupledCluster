
The installation instructions for this repository are the same as those for the [TAMM Library](https://github.com/NWChemEx-Project/TAMM).

Installation
=============
- [Software Requirements](https://tamm.readthedocs.io/en/latest/prerequisites.html)

- [Build Instructions](https://tamm.readthedocs.io/en/latest/install.html)

Build instructions for a quick start
=====================================
## Step 1
```
git clone https://github.com/NWChemEx-Project/TAMM.git
cd TAMM && mkdir build && cd build
```
- ### A detailed list of the cmake build options available are listed [here](https://tamm.readthedocs.io/en/latest/install.html)
```
CC=gcc CXX=g++ FC=gfortran cmake -DCMAKE_INSTALL_PREFIX=<exachem-install-path> -DBUILD_LIBINT=ON -DUSE_GAUXC=ON ..
make -j4 install
```

## Step 2
```
git clone https://github.com/ExaChem/exachem.git
cd exachem && mkdir build && cd build
CC=gcc CXX=g++ FC=gfortran cmake -DCMAKE_INSTALL_PREFIX=<exachem-install-path> -DBUILD_LIBINT=ON -DUSE_GAUXC=ON ..
make -j4
```

## `NOTE:` The cmake configure line in Steps 1 and 2 should be the same.


Building via Spack
------------------
```
spack repo add ./spack
spack install coupledcluster [+cuda]
```

Running the code
=====================
- SCF  
`export CC_EXE=$REPO_INSTALL_PATH/bin/HartreeFock`  

- CCSD  
`export CC_EXE=$REPO_INSTALL_PATH/bin/CD_CCSD`  

- CCSD(T)   
`export CC_EXE=$REPO_INSTALL_PATH/bin/CCSD_T`

### General run:
```
export OMP_NUM_THREADS=1
export INPUT_FILE=$REPO_ROOT_PATH/inputs/ozone.json

mpirun -n 2 $CC_EXE $INPUT_FILE
```

### On Summit:
```
export PAMI_IBV_ENABLE_DCT=1
export PAMI_ENABLE_STRIPING=1
export PAMI_IBV_ADAPTER_AFFINITY=1
export PAMI_IBV_DEVICE_NAME="mlx5_0:1,mlx5_3:1"
export PAMI_IBV_DEVICE_NAME_1="mlx5_3:1,mlx5_0:1"

export GA_PROGRESS_RANKS_DISTRIBUTION_PACKED=1
export GA_NUM_PROGRESS_RANKS_PER_NODE=6

export INPUT_FILE=$REPO_ROOT_PATH/inputs/ubiquitin_dgrtl.json

jsrun -a12 -c12 -g6 -r1 -dpacked $CC_EXE $INPUT_FILE
```
