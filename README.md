## Introduction :

This project contains the C++ reference implementation of the Cascading Metric Tree (CMT).

The CMT data structure is presented in the paper "The Cascading Metric Tree",
by Jeffrey K. Uhlmann and Miguel R. Zuniga (submitted 2021). The paper is under
review; please email the authors for preprints or further information.

The CMT C++ class file is CMTree.h in directory CMT/src/search. Much of this class
is a straightforward implementation of the pseudocode presented in the paper.
The project also contains an implementation of a standard metric tree (see file
SPMTree.h) and various routines used to test and compare the trees and generate
the data published in the paper.

The directory CMT/src/app contains the top-level (function main() ) source code of
applications that create or read datasets and then use the metric trees to perform
various metric queries on the datasets. Features of the query search results that
are needed to measure and study search performance are usually then recorded to a
file. There are also some applications used to test query result correctness, by
comparing the results produced by one tree against another or against brute force.

## Getting the source code :

```
git clone --recursive https://github.com/ngs333/CMT.git
cd CMT
git submodule update --init --recursive
```

(If you are contributing to this project,  instead of cloning from
 https://github.com/ngs333/CMT.git  as above, you should clone from your own fork)

## Compiling:

The project uses CMake. It is recommended that you compile in a directory purposely
set aside for compiling. Say this directory is to be called "build". Consider
performing these operations :

```
cd CMT/src
mkdir build
cd build
cmake ..
make
```

Note that the executables will be automatically placed in directory CMT/bin.

## Using the VS Code IDE: 

Technically no particular IDE is required to change, compile debug and add code,
but development has largely been with Microsoft VS Code and that is what we recommend.
Open with Code (with "File -> Open Folder") the directory CMT and the file explorer
should show the direcotry structure.

The hidden directory CMT/src/.vscode has some visual code files with current project
preferences.

Also we recommend that these Extensions be installed:
a) C/C++ for Visual Studio Code (by Microsoft)
b) C/C++ Extension Pack (by Microsoft)
c) C/C++ Themes (by Microsoft)
d) CMake Tools (by Microsoft)
e) CMake  (by TWXS)


## Compiler Information:

This software has been compiled and tested with g++ 9.3 on Linux using C++ standard c++17.

## Contributing to the project:

Currently we are only accepting software contributions via pull requests. Please
fork your own copies of the project and make official pull requests for contributing.