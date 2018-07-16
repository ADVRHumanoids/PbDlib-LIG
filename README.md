# PbDLib C++ Library

Copyright (C) 2015, Davide De Tommaso, Leonel Rozo, Tohid Alizadeh, Milad Malekzadeh, João Silvério, Sylvain Calinon, Danilo Bruno, Martijn Zeestraten, Ioannis Havoutis and Daniel Berio. Pbdlib has been jointly developed at the Italian Institute of Technology, Genoa, Italy and at the Idiap Research Institute, Martigny, Switzerland, under a GPL Version 3 license.

See also https://gitlab.idiap.ch/rli/pbdlib-matlab for a Matlab/GNU-Octave version of the library containing additional experimental functionalities.

### Dependencies

PbDLib requires Armadillo 5.4 for linear algebra operations, which runs faster if Lapack and OpenBLAS are also installed. PbDLib can be compiled with CMake.

Armadillo 5.4 can be downloaded from [here](http://sourceforge.net/projects/arma/files/armadillo-5.400.2.tar.gz/download).

To install the dependencies run:
```
sudo apt-get install cmake liblapack3 liblapack-dev libopenblas-dev
```
 
### Installation instructions

```
cd pbdlib-sandbox
mkdir build
cd build
cmake ..
make
sudo make install
```

A GUI can be used after installation of the library: 
https://gitlab.idiap.ch/rli/pbdlib_gui

### Test (assuming a build folder was created for cmake install)

```
cd examples
./test_gmm
```

### Using PbDLib in your program

In order to use PbDLib in your program, its CMakeLists.txt must be edited to add the include directories of PbDLib and Armadillo, as well as to link it with both libraries.

#### 1) Add PbDLib and Armadillo include directories

The following commands add the include directories:
```
find_package(Armadillo 5.4 REQUIRED)

include_directories(include 
                    ${ARMADILLO_INCLUDE_DIRS}
                    )
```

The installation of PbDLib copies its header files to /usr/local/include, so the command include_directories(include) adds this directory to the search path.

Header files can then be included in source code using:

```
#include "pbdlib/gmm.h"     // includes header file of the GMM class
```

#### 2) Linking your program with PbD and Armadillo

In order to link your program with PbDLib and Armadillo, add the following commands to CMakeLists.txt:
```
find_library(ARMADILLO_LIBRARIES armadillo)

target_link_libraries(yourProgram pbd
                                  ${ARMADILLO_LIBRARIES}
                                )
```
# PbDlib-LIG
# PbDlib-LIG
# PbDlib-LIG
