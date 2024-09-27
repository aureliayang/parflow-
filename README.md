# ParFlow-CoLM
Coupling ParFlow and the latest Common Land Model

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

### Step 4: Building and installing

Once CMake has configured and created a set of Makefiles; building is
easy:

```shell
   cd build
   make 
   make install
```

