# SKiES - Solver of Kinetic Equation for Solids

SKiES software v. 1.0.0.

SKiES is a parallel C++/Fortran code for first principles transport properties calculations. It uses the LOVA method for solution of kinetic Boltzmann equation proposed by P. Allen.
Current version relies on preliminary calculation with Quantum Espresso and EPW codes.
EPW in turn relies on Wannier90 calculations to obtain all the quantities required next in maximally localized Wannier function representation. SKiES uses quantities in this representation and routines for Wannier interpolation as provided in EPW as an external library to obtain electron and phonon energies, matrix elements of electron-phonon interaction and electron velocities on coarse grids in Bloch representation. With these quantities obtained SKiES calculates transport Eliasberg spectral functions and with them finally electrical resistivity, thermal conductivity and Seebeck coefficient as functions of temperature. See the deatils in **examples** folder.

## Build and Install Quantum Espresso - EPW

QE must be built and installed with the use of cmake/ccmake because SKiES uses the same build system and exploits qeConfig.cmake and MbdConfig.cmake files produced by cmake after the installation of QE. It is higly recommended to first configure the project with ccmake as follows.

Create **build** folder and launch ccmake:
```bash
cd /path/to/src/qunatum-espresso (A file with the main CMakeLists.txt QE file, e.g /home/user/soft/qe)
mkdir build
cd build
ccmake ..
```
Configure project settings pressing [c].
```
 CMAKE_BUILD_TYPE                Release                                 # choose the build type (pressing Enter)
 CMAKE_INSTALL_PREFIX            /home/user/soft/qe/install              # choose the install folder (an example
                                                                         # given here)
 QE_ENABLE_MPI                   OFF                                     # it is recommended to disable MPI in
                                                                         # the pathced version of QE to avoid
                                                                         # some conflicts with the TBB
                                                                         # library used in SKiES (no MPI routines
                                                                         # are actually used in EPW Wannier
                                                                         # interpolation)
 QE_FFTW_VENDOR                  Internal                                # We highly recommend set this to 
                                                                         # "Internal" for the subsequent correct
                                                                         # integration of FFTW library into SKiES # via CMake configuration files

```

## Build and Install SKiES

SKiES is a CMake-built project.
It is higly recommended to first configure the project with ccmake as follows.

Create **build** folder and launch ccmake:
```bash
mkdir build
cd build
ccmake /path/to/src/skies
```
Configure project settings pressing [c].
```
 CMAKE_BUILD_TYPE                Release                                        # choose the build type
 CMAKE_INSTALL_PREFIX            /home/user/soft/skies                          # choose the install folder
                                                                                # (an example given)
 qe_DIR                          /home/user/soft/qe/install/lib/cmake/qe        # path to the qeConfig.cmake
                                                                                # installed by cmake in QE
                                                                                # (an example given)
 Mbd_DIR                         /home/user/soft/qe/install/lib/cmake/mbd       # path to the MbdConfig.cmake
                                                                                # installed by cmake in QE
                                                                                # (an example given)
 SKIES_ENABLE_DOC                ON                                             # create documentation
 SKIES_ENABLE_EXAMPLES           ON                                             # install some ready-to-launch
                                                                                # examples
 SKIES_ENABLE_TBB                ON                                             # install TBB parallel version of
                                                                                # the code
 SKIES_ENABLE_TESTS              ON                                             # build and install tests
```

**IMPORTANT NOTE**. The Fortran compiler used to build QE project must be the same as used to build SKiES. Also it must be compatible with C++ compiler used. If there are some issues at configuration step, try pressing [t] and manually changing  MPI_CXX_COMPILER and MPI_Fortran_COMPILER.

After configuration step press [c] again until there is [g] option is available. Press [g] to generate the project. 

Right after lauch building with make:
```bash
make -j4 install
```

Alternative way (not using ccmake UI) is to do configuration step is only using cmake with additional -D flags in command line, e.g.
```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -Dqe_DIR=/home/user/soft/qe/install/lib/cmake/qe -S /path/to/src/skies
```
