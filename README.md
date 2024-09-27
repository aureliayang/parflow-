# ParFlow-CoLM
Coupling ParFlow and the latest Common Land Model

### Citing Parflow

If you want the DOI for a specific release see:
[Zendo](https://zenodo.org/search?page=1&size=20&q=parflow&version)

A generic DOI that always links to the most current release :
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4816884.svg)](https://doi.org/10.5281/zenodo.4816884)

If you use ParFlow in a publication and wish to cite a paper reference
please use the following that describe model physics:

* Ashby S.F. and R.D. Falgout, Nuclear Science and Engineering 124:145-159, 1996
* Jones, J.E. and C.S. Woodward, Advances in Water Resources 24:763-774, 2001
* Kollet, S.J. and R.M. Maxwell, Advances in Water Resources 29:945-958, 2006
* Maxwell, R.M. Advances in Water Resources 53:109-117, 2013

If you use ParFlow coupled to CLM in a publication, please also cite
two additional papers that describe the coupled model physics:

* Maxwell, R.M. and N.L. Miller, Journal of Hydrometeorology 6(3):233-247, 2005
* Kollet, S.J. and R.M. Maxwell, Water Resources Research 44:W02402, 2008

## Quick Start on Unix/Linux

Important note for users that have built with Autoconf, the CMake
configure process is one step by default.  Most builds of of ParFlow
are on MPP architectures or workstations where the login node and
compute nodes are same architecture the default build process builds
both the ParFlow executable and tools with the same compilers and
libraries in one step.  This will hopefully make building easier for
the majority of users.  It is still possible to build the two
components separately; see instruction below for building pftools and
pfsimulator separately.

CMake supports builds for several operating systems and IDE tools
(like Visual Studio on Windows and XCode on MacOS).  The ParFlow team
has not tested building on platforms other than Linux; there will
likely be some issues on other platforms.  The ParFlow team welcomes
bug reports and patches if you attempt other builds.

### Step 1: Setup

Decide where to install ParFlow and associated libraries.

Set the environment variable `PARFLOW_DIR` to the chosen location:

For bash:

```shell
   export PARFLOW_DIR=/home/snoopy/parflow
```   

For csh and tcsh:

```shell
   setenv PARFLOW_DIR /home/snoopy/parflow
```

### Step 2: Extract the Source

Extract the source files from the compressed tar file.

Obtain the release from the ParFlow GitHub web site:

https://github.com/parflow/parflow/releases

and extract the release.  Here we assume you are building in new
subdirectory in your home directory:

```shell
   mkdir ~/parflow 
   cd ~/parflow 
   tar -xvf ../parflow.tar.gz
```

Note the ParFlow tar file will be have a different name based on the
version number.

If you are not using GNU tar or have a very old version GNU tar you
will need to uncompress the file first:

```shell
   mkdir ~/parflow 
   cd ~/parflow 
   gunzip ../parflow.tar.gz
   tar -xvf ../parflow.tar
```

### Step 3: Running CMake to configure ParFlow

CMake is a utility that sets up makefiles for building ParFlow.  CMake
allows setting of compiler to use and other options.  First create a
directory for the build.  It is generally recommend to build outside
of the source directory to make it keep things clean.  For example,
restarting a failed build with a separate build directory simply
involves removing the build directory.

#### Building with the ccmake GUI

You can control build options for ParFlow using the ccmake GUI.

```shell
   mkdir build
   cd build
   ccmake ../parflow 
```
At a minimum, you will want to set the CMAKE_INSTALL_PREFIX value to the same thing
as PARFLOW_DIR was set to above.  Other variables should be set as desired.

After setting a variable 'c' will configure `ParFlow.  When you are
completely done setting configuration options, use 'g' to generate the
configuration and exit ccmake.

If you are new to CMake, the creators of CMake provide some additional ccmake usage notes here:

https://cmake.org/runningcmake/

#### Building with the cmake command line

CMake may also be configured from the command line using the cmake
command. Instructions to build with different accelerator backends are found from the following documents: [CUDA, KOKKOS](README-GPU.md), [OpenMP](README-OPENMP.md). The default will configure a sequential version of ParFlow
using MPI libraries.  CLM is being enabled.

```shell
   mkdir build
   cd build
   cmake ../parflow \
   	 -DCMAKE_INSTALL_PREFIX=${PARFLOW_DIR} \
   	 -DPARFLOW_HAVE_CLM=ON
```

If TCL is not installed in the standard locations (/usr or /usr/local)
you need to specify the path to the tclsh location:

```shell
	-DTCL_TCLSH=${PARFLOW_TCL_DIR}/bin/tclsh8.6
```

Building a parallel version of ParFlow requires the communications
layer to use must be set.  The most common option will be MPI.  Here
is a minimal example of an MPI build with CLM:

```shell
   mkdir build
   cd build
   cmake ../parflow \
      	 -DCMAKE_INSTALL_PREFIX=${PARFLOW_DIR} \
   	 -DPARFLOW_HAVE_CLM=ON \
	 -DPARFLOW_AMPS_LAYER=mpi1
```

Here is a more complex example where location of various external
packages are being specified and some features are being enabled:

```shell
   mkdir build
   cd build
   cmake ../parflow \
        -DPARFLOW_AMPS_LAYER=mpi1 \
	-DHYPRE_ROOT=${PARFLOW_HYPRE_DIR} \
	-DHDF5_ROOT=${PARFLOW_HDF5_DIR} \
	-DSILO_ROOT=${PARFLOW_SILO_DIR} \
	-DCMAKE_BUILD_TYPE=Debug \
	-DPARFLOW_ENABLE_TIMING=TRUE \
	-DPARFLOW_HAVE_CLM=ON \
	-DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}
```

### Step 4: Building and installing

Once CMake has configured and created a set of Makefiles; building is
easy:

```shell
   cd build
   make 
   make install
```

