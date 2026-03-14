
Usage:

1. compile with make [-e debug=no];
2. generate the mesh with easymesh as:
	     easymesh D
3. run with ./main;
4. view the obtained numerical solution with OpenDX;



  

INSTALLATION OF AFEPACK 

System environment：ubuntu22.04

Table of Contents


1. Environment:
2. Some issue:
3. Flowchart
.. 1. deal.II installation
.. 2. AFEPack Installation:


1 Environment:

sudo apt-get install cmake

sudo apt-get install g++ gcc gfortran fort77

sudo apt-get install libboost-all-dev

sudo apt-get install libtbb-dev

sudo apt-get install automake

2 Some issue(openmpi):


Here we directly install the followings from the library:
sudo apt-get install build-essential binutils libopenmpi-dev openmpi-doc openmpi-bin 
sudo apt-get install liblapack-dev libblas-dev libboost-all-dev libpcre3-dev

Check if the installation was successful. After completing the installation, you can run the following commands to verify that OpenMPI has been installed and is working properly. For example, Run a simple example:
mpicc hello.c -o hello
mpirun -np 2 ./hello

 Generally, if openmpi is installed successfully from the apt source, the link would be established by default. Enter mpicc --version(mpicxx --version) to check the link. Otherwise establish the link
mpic++ -> /usr/bin/mpic++.openmpi*,   mpicxx -> /usr/bin/mpic++.openmpi*


3 Flowchart

3.1 deal.II installation

It is recommended to install deal.II in /usr/local/deal.II.

  ----------
  • cd /usr/lib/x86_64-linux-gnu/cmake/Boost-1.74.0
  Add the following line of command before line 240 of the BoostConfig.cmake file:
  if(POLICY CMP0057) cmake_policy(SET CMP0057 NEW) endif()
   • If this line is not added, errors related to Boost may occur.


  • Unpack dealii.8.1.0.tar.gz
  • Follow README.md
  • Using following command
  ┌────
	mkdir build
	cd build
	cmake -DCMAKE_INSTALL_PREFIX=/usr/local/deal.II -DDEAL_II_WITH_MPI=ON -DDEAL_II_WITH_UMFPACK=OFF -DDEAL_II_WITH_TRILINOS=OFF -	DDEAL_II_WITH_ZLIB=OFF -DDEAL_II_COMPONENT_EXAMPLES=OFF -DDEAL_II_WITH_PETSC=OFF -DDEAL_II_COMPONENT_MESH_CONVERT=OFF -	DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_FC_COMPILER="mpif90" ..
	make -j
        sudo make install -j
  └────

   • In CMakeList.txt, Put following sentences after CMAKE_ MINIMUM_
    REQUIRED …
  ┌────
  │ CMAKE_POLICY(SET CMP0026 OLD)
  │ CMAKE_POLICY(SET CMP0037 OLD)
  │ CMAKE_POLICY(SET CMP0042 NEW)
  └────

  • During the compiling, we have following issues

    • In /lac/sparsity_ pattern.h:1505: error 'upper-bound' is not a
      member of 'std'

      Reason: need the header #include <algorithm>

    • In /lac/vector.templates.h:806: error: 'numerical_limits' is not a
      member of 'std'

      Reason: need the header #include <limits>

    • In /source/base/parameter_ handler.cc: 1278: error: cannot convert
      'boost::optional<std::__ cxx11::basic_ string<char> >' to 'bool'
      in return

      How: In 1278 line: change to "return bool (p.get_
      optinal<std::string>("value"));"

    • /include/openmpi-gcc12/mpi.h: 342: error: expected
      primary-expression before 'static_ assert'

      How: actual error place is: deal.II/source/base/mpi.cc: 253: MPI_
      Type_ struct

      Take a look at:
      <https://www.open-mpi.org/faq/?category=mpi-removed#mpi-1-mpi-type-struct>

      In Line 253: change the function to "MPI_ Type_ create_ struct".



3.2 AFEPack Installation:


  • Unpack the package

  • Implement the following:
    
  ┌───
  │ $ sudo ln -sf /usr/lib/x86_64-linux-gnu/openmpi/include/* /usr/local/include
  │ $ CC="mpicc" CXX="mpicxx" CPPFLAGS="-I/usr/local/include -I/usr/local/deal.II/include -I/usr/local/deal.II/include/deal.II" LDFLAGS="-L/usr/local/lib -L/usr/local/deal.II/lib"  ./configure --prefix=/usr/local/AFEPack
  │ $ make 
  │ $ make install
  └────
 • Environment variables
┌───
│ $ vim ~/.bashrc
      Add the followings at the end of this file.
      export AFEPACK_PATH="/usr/local/AFEPack/include/AFEPack"
     export AFEPACK_TEMPLATE_PATH="$AFEPACK_PATH/template/triangle:$AFEPACK_PATH/template/rectangle:$AFEPACK_PATH/template/interval:$AFEPACK_PATH/template/twin_triangle:$AFEPACK_PATH/template/tetrahedron:$AFEPACK_PATH/template/twin_tetrahedron:$AFEPACK_PATH/template/four_tetrahedron"

     export LD_LIBRARY_PATH="/usr/local/deal.II/lib:/usr/local/AFEPack/lib"

  $ source ~/.bashrc
└────

• Copy the easymesh file from this folder to /usr/bin
